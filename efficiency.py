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


def stacked(eff, h_stack, colors):
    h_cats = []
    leg = ROOT.TLegend(0.1, 0.79, 0.9, 0.91)

    for cat in eff.selected:
        h_cat = eff.selected[cat].Clone()
        e = eff.efficiency[cat]
        h_cat.Divide(eff.tot)
        if colors:
            h_cat.SetFillColor(ROOT.TColor.GetColor(colors[cat]))
        h_cat.SetLineColor(ROOT.TColor.GetColor("#555555"))
        h_cat.SetLineWidth(1)
        if cat == "passed":
            h_cat.SetLineWidth(3)
            h_cat.SetLineColor(ROOT.TColor.GetColor(1))
            h_cat.SetFillColor(ROOT.TColor.GetColor(0))
        h_cats.append(h_cat)
        leg.AddEntry(h_cat, "%s: %.1f%%" % (h_cat.GetName(), e*100), "f")
        h_stack.Add(h_cat)

    return h_cats, leg


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


def lepton_angles(chain_filter):
    theta, phi = -999, -999

    for i, pdg in enumerate(chain_filter.nu_daughters_pdg):
        if pdg == 11:
            px = chain_filter.nu_daughters_p[i][0]
            py = chain_filter.nu_daughters_p[i][1]
            pz = chain_filter.nu_daughters_p[i][2]
            theta = math.degrees(math.acos(pz/math.sqrt(px**2+py**2+pz**2)))
            phi = math.degrees(math.atan2(py,px))
            break

    return theta, phi


def efficiency(files_path, eff_variables=[], scale=1):
    nue_cosmic = glob(files_path+"/*.root")
    chain_nue = ROOT.TChain("robertoana/pandoratree")
    chain_filter = ROOT.TChain("nueFilter/filtertree")
    chain_pot = ROOT.TChain("nueFilter/pot")

    # print(nue_cosmic)
    for f in nue_cosmic:
        chain_filter.Add(f)
        chain_nue.Add(f)
        chain_pot.Add(f)

    entries = int(chain_filter.GetEntries() / scale)
    pot_entries = int(chain_pot.GetEntries() / scale)

    total_pot = 0
    for i in range(pot_entries):
        chain_pot.GetEntry(i)
        total_pot += chain_pot.pot

    print("Entries", pot_entries)
    print("%.2e POT" % total_pot)

    categories = ["passed", "quality cuts", "CC #nu_{#mu} selected", "not contained",
                  "cosmic selected", "1 shower",  "no showers",  "no flash", "no data products"]

    h_tot = {}

    for v in eff_variables:
        h_tot[v] = ROOT.TH1F("h_tot_%s" % v,
                             labels[v],
                             binning[v][0],
                             binning[v][1],
                             binning[v][2])

    h_selected = {}
    for v in eff_variables:
        h_sel = {}
        for c in categories:
            h_sel[c] = ROOT.TH1F("%s" % c,
                                 labels[v],
                                 binning[v][0],
                                 binning[v][1],
                                 binning[v][2])
        h_selected[v] = h_sel

    i_nue = -1
    eff_vars = {}
    var_dict = dict(variables + spectators)
    ROOT.TMVA.Tools.Instance()
    reader = ROOT.TMVA.Reader(":".join([
        "!V",
        "Silent",
        "Color"]))
    load_bdt(reader)

    for i_evt in range(entries):
        printProgressBar(i_evt, entries, prefix="Progress:",
                         suffix="Complete", length=20)
        chain_filter.GetEntry(i_evt)

        eNp = check_cc0pinp(chain_filter)

        true_neutrino_vertex = [chain_filter.true_vx_sce,
                                chain_filter.true_vy_sce,
                                chain_filter.true_vz_sce]

        if chain_filter.passed:
            i_nue += 1

        if not (is_fiducial(true_neutrino_vertex) and eNp):
            continue

        eff_vars["energy"] = chain_filter.nu_energy
        eff_vars["theta"], eff_vars["phi"] = lepton_angles(chain_filter)

        flux_weights = [1]  # [1]*len(chain_nue.flux_weights[0])
        # for fl in chain_nue.flux_weights:
        #     flux_weights = [a*b for a, b in zip(flux_weights, fl)]
        weight = chain_filter.bnbweight
        for v in eff_variables:
            h_tot[v].Fill(eff_vars[v], weight)

            if not chain_filter.passed:
                if chain_filter.selection_result == 2:
                    h_selected[v]["no showers"].Fill(eff_vars[v], weight)
                if chain_filter.selection_result == 4:
                    h_selected[v]["no flash"].Fill(eff_vars[v], weight)
                if chain_filter.selection_result == 5:
                    h_selected[v]["no data products"].Fill(eff_vars[v], weight)
                continue

            chain_nue.GetEntry(i_nue)
            contaminated = chain_nue.cosmic_fraction < 0.5 and chain_nue.category == 7
            selected = (chain_nue.category == 2 or contaminated) and chain_nue.passed == 1
            if not selected:
                h_selected[v]["cosmic selected"].Fill(eff_vars[v], weight)
                continue

            # If there are no tracks we require at least two showers
            showers_2_tracks_0 = True
            if chain_nue.n_tracks == 0 and chain_nue.n_showers == 1:
                showers_2_tracks_0 = False

            if not showers_2_tracks_0:
                h_selected[v]["1 shower"].Fill(eff_vars[v], weight)
                continue

            if chain_nue.numu_passed == 1:
                h_selected[v]["CC #nu_{#mu} selected"].Fill(eff_vars[v], weight)
                continue

            if not check_reco_fidvol(chain_nue):
                h_selected[v]["not contained"].Fill(eff_vars[v], weight)
                continue

            fill_kin_branches(chain_nue, 1, var_dict, "nue", False)

            if not pre_cuts(var_dict):
                h_selected[v]["quality cuts"].Fill(eff_vars[v], weight)
                continue

            h_selected[v]["passed"].Fill(eff_vars[v], weight)

            # bdt_values = {}
            # for bdt_name in BDT_TYPES:
            #     bdt_values[bdt_name] = reader.EvaluateMVA("BDT%s" % bdt_name)

        #     if apply_cuts(bdt_values, var_dict, bdt=False, manual=True):
        #         h_selected["cuts"].Fill(eff_vars[v], weight)
        #     if apply_cuts(bdt_values, var_dict, bdt=True, manual=False):
        #         h_selected["bdt"].Fill(eff_vars[v], weight)

    eff_objects = {}
    for v in eff_variables:
        eff_objects[v] = Efficiency(h_tot[v], h_selected[v], v, total_pot)
    return eff_objects

