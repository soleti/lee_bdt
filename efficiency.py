#!/usr/local/bin/python3

import math
import ROOT
import sys
from array import array
from glob import glob

from settings import printProgressBar, distance, binning, labels
from settings import is_fiducial, is_active, variables, spectators
from settings import pre_cuts, load_bdt, apply_cuts, fill_kin_branches
from settings import BINS, BINS, BDT_TYPES, N_UNI, ELECTRON_THRESHOLD, PROTON_THRESHOLD
from settings import PROTON_MASS, ELECTRON_MASS

from proton_energy import length2energy

ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptFit(0)
ROOT.gStyle.SetPalette(87)
ROOT.gStyle.SetNumberContours(99)


class Efficiency:
    def __init__(self, tot, selected={}, var="", pot=0, weights=[]):
        self.tot = tot
        self.pot = pot
        self.selected = selected
        self.efficiency = {}
        self.efficiency_err = {}
        self.variable = var
        self.weights = weights
        for s in selected:
            eff_value = selected[s].Integral()/tot.Integral()
            self.efficiency[s] = eff_value
            self.efficiency_err[s] = 1./tot.Integral() * math.sqrt(selected[s].Integral()*(1-eff_value))

    def tefficiency(self, num, uni=0):
        teff = ROOT.TEfficiency(self.selected[num], self.tot)
        return teff


def check_reco_fidvol(chain_nue):
    track_fidvol = True

    for i_tr in range(chain_nue.n_tracks):
        track_start = [
            chain_nue.track_start_x[i_tr],
            chain_nue.track_start_y[i_tr],
            chain_nue.track_start_z[i_tr]
        ]
        track_end = [
            chain_nue.track_end_x[i_tr],
            chain_nue.track_end_y[i_tr],
            chain_nue.track_end_z[i_tr]
        ]

        track_fidvol = track_fidvol and is_fiducial(
            track_start) and is_fiducial(track_end)
        if not track_fidvol:
            break

    shower_fidvol = True

    for i_sh in range(chain_nue.n_showers):
        shower_start = [
            chain_nue.shower_start_x[i_sh],
            chain_nue.shower_start_y[i_sh],
            chain_nue.shower_start_z[i_sh]
        ]

        shower_fidvol = shower_fidvol and is_fiducial(shower_start)
        if not shower_fidvol:
            break

    neutrino_vertex = [chain_nue.vx, chain_nue.vy, chain_nue.vz]
    neutrino_fidvol = is_fiducial(neutrino_vertex)

    return track_fidvol and shower_fidvol and neutrino_fidvol

def check_cc0pinp(chain_nue):

    protons = 0
    electrons = 0
    photons = 0
    pions = 0
    electron_energy = 0
    proton_energy = 0

    for i, energy in enumerate(chain_nue.nu_daughters_E):
        if chain_nue.nu_daughters_pdg[i] == 2212:
            proton_energy += energy - 0.938

            p_vertex = [chain_nue.nu_daughters_start_v[i][0],
                        chain_nue.nu_daughters_start_v[i][1],
                        chain_nue.nu_daughters_start_v[i][2]]

            p_end = [chain_nue.nu_daughters_end_v[i][0],
                     chain_nue.nu_daughters_end_v[i][1],
                     chain_nue.nu_daughters_end_v[i][2]]

            if energy - 0.938 > PROTON_THRESHOLD:
                protons += 1
            if not is_fiducial(p_vertex) or not is_fiducial(p_end):
                proton = 0
                break

        if chain_nue.nu_daughters_pdg[i] == 11:
            electron_energy += energy
            e_vertex = [chain_nue.nu_daughters_end_v[i][0],
                        chain_nue.nu_daughters_end_v[i][1],
                        chain_nue.nu_daughters_end_v[i][2]]

            if energy - 0.51e-3 > ELECTRON_THRESHOLD and is_fiducial(e_vertex):
                electrons += 1

        if chain_nue.nu_daughters_pdg[i] == 22:
            # if energy > 0.035:
            photons += 1

        if chain_nue.nu_daughters_pdg[i] == 111 or chain_nue.nu_daughters_pdg[i] == 211:
            # if energy > 0.06:
            pions += 1

    return electrons == 1 and pions == 0 and protons >= 1 and photons == 0 and chain_nue.ccnc == 0


def draw_eff(efficiencies, var, name, descriptions, colors):
    c_eff = ROOT.TCanvas("c_eff_%s_%s" % (var, name))
    c_eff.SetTopMargin(0.18)

    leg = ROOT.TLegend(0.097, 0.844, 0.897, 0.958)
    leg.SetNColumns(2)
    plots = []

    for i, e in enumerate(efficiencies):
        teff = efficiencies[e].tefficiency(name)
        teff.SetLineColor(colors[i])
        teff.SetMarkerStyle(0)
        teff.SetLineWidth(2)
        value = efficiencies[e].efficiency[name]*100
        err = efficiencies[e].efficiency_err[name]*100

        plots.append(teff)

        if i == 0:
            teff.Draw("AP")
        else:
            teff.Draw("P SAME")
        c_eff.Draw()

        leg.AddEntry(teff,
                     "%s: (%.1f #pm %.1f) %%" % (descriptions[i], value, err),
                     "le")
        teff.GetPaintedGraph().SetMinimum(0)
        teff.GetPaintedGraph().SetMaximum(1)

    c_eff.Draw()

    leg.Draw()
    ROOT.gPad.Update()
    return c_eff, plots, leg

def efficiency(files_path, name, effs=["nu_E"]):
    nue_cosmic = glob(files_path+"/newfilter/*.root")
    chain_nue = ROOT.TChain("robertoana/pandoratree")
    chain_filter = ROOT.TChain("nueFilter/filtertree")
    chain_pot = ROOT.TChain("nueFilter/pot")

    # print(nue_cosmic)
    for f in nue_cosmic:
        chain_filter.Add(f)
        chain_nue.Add(f)
        chain_pot.Add(f)

    scale = 10
    entries = int(chain_filter.GetEntries() / scale)
    pot_entries = int(chain_pot.GetEntries() / scale)

    total_pot = 0
    for i in range(pot_entries):
        chain_pot.GetEntry(i)
        total_pot += chain_pot.pot

    print("Sample %s" % name)
    print("Entries", pot_entries)
    print("Total pot", total_pot)

    h_tot = {}
    h_selected = {}
    h_selected_numu = {}
    h_selected_precuts = {}
    h_selected_cuts = {}
    h_selected_bdt = {}
    h_selected_contained = {}

    for variable in effs:
        h_tot[variable] = ROOT.TH1F("h_tot_%s_%s" % (variable, name),
                                    ";%s;Efficiency" % labels[variable].split(";")[1],
                                    binning[variable][0],
                                    binning[variable][1],
                                    binning[variable][2])
        h_selected[variable] = ROOT.TH1F("h_selected_%s_%s" % (variable, name),
                                         ";%s;Efficiency" % labels[variable].split(";")[1],
                                         binning[variable][0],
                                         binning[variable][1],
                                         binning[variable][2])
        h_selected_numu[variable] = ROOT.TH1F("h_selected_numu_%s_%s" % (variable, name),
                                              ";%s;Efficiency" % labels[variable].split(";")[1],
                                              binning[variable][0],
                                              binning[variable][1],
                                              binning[variable][2])
        h_selected_precuts[variable] = ROOT.TH1F("h_selected_precuts_%s_%s" % (variable, name),
                                                 ";#%s;Efficiency" % labels[variable].split(";")[1],
                                                 binning[variable][0],
                                                 binning[variable][1],
                                                 binning[variable][2])
        h_selected_contained[variable] = ROOT.TH1F("h_selected_contained_%s_%s" % (variable, name),
                                                   ";#%s;Efficiency" % labels[variable].split(";")[1],
                                                   binning[variable][0],
                                                   binning[variable][1],
                                                   binning[variable][2])
        h_selected_cuts[variable] = ROOT.TH1F("h_selected_cuts_%s_%s" % (variable, name),
                                              ";%s;Efficiency" % labels[variable].split(";")[1],
                                              binning[variable][0],
                                              binning[variable][1],
                                              binning[variable][2])
        h_selected_bdt[variable] = ROOT.TH1F("h_selected_bdt_%s_%s" % (variable, name),
                                             ";%s;Efficiency" % labels[variable].split(";")[1],
                                             binning[variable][0],
                                             binning[variable][1],
                                             binning[variable][2])

    var_dict = dict(variables + spectators)
    ROOT.TMVA.Tools.Instance()
    reader = ROOT.TMVA.Reader(":".join([
        "!V",
        "Silent",
        "Color"]))
    load_bdt(reader)

    vars = {}

    i_nue = -1
    for i_evt in range(entries):
        printProgressBar(i_evt, entries, prefix="Progress:",
                        suffix="Complete", length=20)
        chain_filter.GetEntry(i_evt)

        eNp = check_cc0pinp(chain_filter)

        true_neutrino_vertex = [chain_filter.true_vx_sce,
                                chain_filter.true_vy_sce,
                                chain_filter.true_vz_sce]
        true_neutrino_vertex_nosce = [chain_filter.true_vx,
                                      chain_filter.true_vy,
                                      chain_filter.true_vz]

        if chain_filter.passed:
            i_nue += 1

        if not (is_fiducial(true_neutrino_vertex) and eNp):
            continue

        for i, pdg in enumerate(chain_filter.nu_daughters_pdg):
            if pdg == 11:
                px = chain_filter.nu_daughters_end_v[i][0] - chain_filter.nu_daughters_start_v[i][0]
                py = chain_filter.nu_daughters_end_v[i][1] - chain_filter.nu_daughters_start_v[i][1]
                pz = chain_filter.nu_daughters_end_v[i][2] - chain_filter.nu_daughters_start_v[i][2]
                theta = math.degrees(math.acos(pz/math.sqrt(px**2+py**2+pz**2)))
                phi = math.degrees(math.atan(py/px))
                break

        vars["nu_E"] = chain_filter.nu_energy
        vars["theta"] = theta
        vars["phi"] = phi

        flux_weights = [1]#[1]*len(chain_nue.flux_weights[0])
        # for fl in chain_nue.flux_weights:
        #     flux_weights = [a*b for a, b in zip(flux_weights, fl)]

        for v in effs:
            weight = chain_filter.bnbweight
            h_tot[v].Fill(vars[v], weight)

        if not chain_filter.passed:
            continue

        chain_nue.GetEntry(i_nue)
        contaminated = chain_nue.cosmic_fraction < 0.5 and chain_nue.category == 7
        selected = (chain_nue.category == 2 or contaminated) and chain_nue.passed == 1
        if not selected:
            continue

        # If there are no tracks we require at least two showers
        showers_2_tracks_0 = True
        if chain_nue.n_tracks == 0 and chain_nue.n_showers == 1:
            showers_2_tracks_0 = False
        if not showers_2_tracks_0:
            continue

        for v in effs:
            weight = chain_nue.bnbweight
            h_selected[v].Fill(vars[v], weight)

        fill_kin_branches(chain_nue, 1, var_dict, "nue", False)

        # if not int(var_dict["numu_score"][0]) == 0:
        #     continue

        for v in effs:
            for u in range(N_UNI):
                weight = chain_nue.bnbweight
                h_selected_numu[v].Fill(vars[v], weight)

        if not check_reco_fidvol(chain_nue):
            continue

        for v in effs:
            for u in range(N_UNI):
                # * chain_nue.genie_weights[0][u] * chain_nue.bnbweight
                weight = chain_nue.bnbweight
                h_selected_contained[v].Fill(vars[v], weight)

        if not pre_cuts(var_dict):
            continue

        bdt_values = {}
        for bdt_name in BDT_TYPES:
            bdt_values[bdt_name] = reader.EvaluateMVA("BDT%s" % bdt_name)

        for v in effs:
            weight = chain_nue.bnbweight
            h_selected_precuts[v].Fill(vars[v], weight)
            if apply_cuts(bdt_values, var_dict, bdt=False, manual=True):
                h_selected_cuts[v].Fill(vars[v], weight)
            if apply_cuts(bdt_values, var_dict, bdt=True, manual=False):
                h_selected_bdt[v].Fill(vars[v], weight)

    eff_objects = {}

    for v in effs:
        selected = {"topology": h_selected[v],
                    "precuts": h_selected_precuts[v],
                    "cuts": h_selected_cuts[v],
                    "bdt": h_selected_bdt[v],
                    "contained": h_selected_contained[v]}
        eff_objects[v] = (Efficiency(h_tot[v], selected, variable, total_pot))


    print()
    print("Integral", h_tot["nu_E"].Integral())
    return eff_objects
