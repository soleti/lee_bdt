#!/usr/local/bin/python3

import math
import statistics
from glob import glob

import ROOT

from settings import printProgressBar, binning, labels
from settings import is_fiducial, variables, spectators
from settings import pre_cuts, load_bdt, apply_cuts, fill_kin_branches
from settings import N_UNI, ELECTRON_THRESHOLD, PROTON_THRESHOLD
from settings import PROTON_MASS, ELECTRON_MASS


ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptFit(0)
ROOT.gStyle.SetPalette(87)
ROOT.gStyle.SetNumberContours(99)


class Efficiency:
    def __init__(self, tot, selected={}, var="", pot=0, total_sys_err=0, sys_err=[]):
        self.tot = tot
        self.pot = pot
        self.selected = selected
        self.efficiency = {}
        self.efficiency_err = {}
        self.variable = var
        self.sys_err = sys_err
        self.total_sys_err = total_sys_err
        for sel in selected:
            eff_value = selected[sel].Integral()/tot.Integral()
            self.efficiency[sel] = eff_value
            self.efficiency_err[sel] = 1./tot.Integral() * math.sqrt(selected[sel].Integral()*(1-eff_value))

    def tefficiency(self, num):
        teff = ROOT.TEfficiency(self.selected[num], self.tot)
        return teff

    def tefficiency_err(self, num):
        teff = self.tefficiency(num)
        teff_sys = self.selected[num].Clone()
        teff_sys.Divide(self.tot)

        for i_bin in range(1, self.selected[num].GetNbinsX()+1):
            stat_error = teff.GetEfficiencyErrorLow(i_bin) + teff.GetEfficiencyErrorUp(i_bin)
            teff_sys.SetBinError(i_bin, math.sqrt(stat_error**2 + self.sys_err[i_bin - 1]**2))

        return teff_sys


    def draw(self, systematics=True):
        c_eff_sys = ROOT.TCanvas("c_eff_sys_%s" % self.tot.GetName())
        eff_stat = self.tefficiency("passed")
        eff_sys = self.tefficiency_err("passed")
        eff_stat.SetLineWidth(2)
        eff_stat.Draw("ap")
        eff_sys.SetLineColor(1)

        if systematics:
            eff_sys.Draw("e same")
            eff_stat.Draw("p same")

        l_eff = ROOT.TLegend(0.097, 0.854, 0.9, 0.9)

        if systematics:
            l_eff.AddEntry(eff_stat,
                        "Selection efficiency: (%.1f #pm %.1f (stat) #pm %.1f (sys)) %%"
                        % (self.efficiency["passed"] * 100, self.efficiency_err["passed"] * 100, self.total_sys_err * 100),
                        "le")
        else:
            l_eff.AddEntry(eff_stat,
                           "Selection efficiency: (%.1f #pm %.1f) %%"
                           % (self.efficiency["passed"] * 100, self.efficiency_err["passed"] * 100),
                           "le")
        l_eff.Draw()
        c_eff_sys.SetTopMargin(0.18)

        return c_eff_sys, eff_stat, eff_sys, l_eff


def draw_true(histo):
    tot = histo.Clone()
    tot.SetLineColor(ROOT.TColor.GetColor("#62b570"))
    tot.SetLineWidth(3)
    leg = ROOT.TLegend(0.097, 0.854, 0.9, 0.9)
    leg.SetNColumns(2)
    leg.AddEntry(tot,
                 "Generated #nu_{e} CC0#pi-Np: %.1f events" % tot.Integral(),
                 "f")
    c = ROOT.TCanvas("c_tot")
    tot.Draw("hist")
    c.SetTopMargin(0.16)

    return c, tot, leg


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
            h_cat.SetLineColor(1)
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
            proton_energy += energy - PROTON_MASS

            p_vertex = [chain_nue.nu_daughters_start_v[i][0],
                        chain_nue.nu_daughters_start_v[i][1],
                        chain_nue.nu_daughters_start_v[i][2]]

            p_end = [chain_nue.nu_daughters_end_v[i][0],
                     chain_nue.nu_daughters_end_v[i][1],
                     chain_nue.nu_daughters_end_v[i][2]]

            if energy - PROTON_MASS > PROTON_THRESHOLD:
                protons += 1
            if not is_fiducial(p_vertex) or not is_fiducial(p_end):
                protons = 0
                break

        if chain_nue.nu_daughters_pdg[i] == 11:
            electron_energy += energy
            e_vertex = [chain_nue.nu_daughters_end_v[i][0],
                        chain_nue.nu_daughters_end_v[i][1],
                        chain_nue.nu_daughters_end_v[i][2]]

            if energy - ELECTRON_MASS > ELECTRON_THRESHOLD and is_fiducial(e_vertex):
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
            phi = math.degrees(math.atan2(py, px))
            break

    return theta, phi


def eff_total_sys_error(h_tot_sys, h_passed_sys):
    ratios = []
    for u in range(N_UNI):
        ratio = h_passed_sys[u].Integral()/h_tot_sys[u].Integral()
        ratios.append(ratio)

    return statistics.stdev(ratios)


def eff_sys_error(h_tot_sys, h_passed_sys):
    sys_error = []

    for i_bin in range(1, h_tot_sys[0].GetNbinsX()+1):
        err = []
        for u in range(N_UNI):
            if h_tot_sys[u].GetBinContent(i_bin):
                err.append(h_passed_sys[u].GetBinContent(i_bin)/h_tot_sys[u].GetBinContent(i_bin))

        if len(err) > 1:
            sys_error.append(statistics.stdev(err))
        else:
            sys_error.append(0)

    return sys_error


def efficiency(files_path, eff_variables=[], systematics=False, scale=1):
    nue_cosmic = glob(files_path+"/*.root")
    chain_nue = ROOT.TChain("robertoana/pandoratree")
    chain_filter = ROOT.TChain("nueFilter/filtertree")
    chain_pot = ROOT.TChain("nueFilter/pot")

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
    h_tot_sys = {}

    for v in eff_variables:
        h_tot[v] = ROOT.TH1F("h_tot_%s" % v,
                             labels[v],
                             binning[v][0],
                             binning[v][1],
                             binning[v][2])
        h_tot_uni = []
        for u in range(N_UNI):
            h_tot_uni.append(ROOT.TH1F("h_tot_%s_%i" % (v, u),
                                       labels[v],
                                       binning[v][0],
                                       binning[v][1],
                                       binning[v][2]))
        h_tot_sys[v] = h_tot_uni

    h_selected = {}
    h_passed_sys = {}
    for v in eff_variables:
        h_sel = {}
        for c in categories:
            h_sel[c] = ROOT.TH1F("%s" % c,
                                 labels[v],
                                 binning[v][0],
                                 binning[v][1],
                                 binning[v][2])
        h_selected[v] = h_sel

        h_passed_uni = []
        for u in range(N_UNI):
            h_passed_uni.append(ROOT.TH1F("h_tot_%s_%i" % (v, u),
                                          labels[v],
                                          binning[v][0],
                                          binning[v][1],
                                          binning[v][2]))
        h_passed_sys[v] = h_passed_uni

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

        if systematics:
            flux_weights = [1]*len(chain_filter.flux_weights[0])
            for fl in chain_filter.flux_weights:
                flux_weights = [a*b for a, b in zip(flux_weights, fl)]
            genie_weights = chain_filter.genie_weights[0]
        weight = chain_filter.bnbweight

        for v in eff_variables:
            h_tot[v].Fill(eff_vars[v], weight)

            for u in range(N_UNI):
                h_tot_sys[v][u].Fill(eff_vars[v], weight * flux_weights[u] * genie_weights[u])

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
            for u in range(N_UNI):
                h_passed_sys[v][u].Fill(eff_vars[v], weight * flux_weights[u] * genie_weights[u])

            # bdt_values = {}
            # for bdt_name in BDT_TYPES:
            #     bdt_values[bdt_name] = reader.EvaluateMVA("BDT%s" % bdt_name)

        #     if apply_cuts(bdt_values, var_dict, bdt=False, manual=True):
        #         h_selected["cuts"].Fill(eff_vars[v], weight)
        #     if apply_cuts(bdt_values, var_dict, bdt=True, manual=False):
        #         h_selected["bdt"].Fill(eff_vars[v], weight)

    eff_objects = {}
    for v in eff_variables:
        sys_err = eff_sys_error(h_tot_sys[v], h_passed_sys[v])
        total_sys_err = eff_total_sys_error(h_tot_sys[v], h_passed_sys[v])
        eff_objects[v] = Efficiency(h_tot[v], h_selected[v], v, total_pot, total_sys_err, sys_err)

    return eff_objects
