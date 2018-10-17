#!/usr/local/bin/python3

import math
import ROOT

from glob import glob
from bdt_common import x_start, x_end, y_start, y_end, z_start, z_end
from bdt_common import printProgressBar, bins, bins2, distance, is_fiducial, is_active
from bdt_common import N_UNI, variables, spectators, pre_cuts, load_bdt, bdt_types, apply_cuts
from fill import fill_kin_branches
from array import array
from proton_energy import length2energy

ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptFit(0)
ROOT.gStyle.SetPalette(87)
ROOT.gStyle.SetNumberContours(99)

# glob("mc_nue_ubxsec/*.root")
nue_cosmic = glob("data_files/mc_nue_sys/*.root")
chain_nue = ROOT.TChain("robertoana/pandoratree")
chain_filter = ROOT.TChain("nueFilter/filtertree")
chain_pot = ROOT.TChain("nueFilter/pot")

for f in nue_cosmic:
    chain_nue.Add(f)
    chain_filter.Add(f)
    chain_pot.Add(f)

total_pot = 0
for i in range(chain_pot.GetEntries()):
    chain_pot.GetEntry(i)
    total_pot += chain_pot.pot
print(total_pot)

entries = int(chain_nue.GetEntries() / 1)
print("entries", entries)

PROTON_THRESHOLD = 0.020
ELECTRON_THRESHOLD = 0.020

h_tot = []
h_selected = []
h_selected_numu = []
h_selected_cuts = []
h_selected_precuts = []
h_selected_bdt = []

for u in range(N_UNI):
    h_tot.append(ROOT.TH1F("h_tot_%i" % u,
                           ";#nu_{e} energy [GeV];Efficiency",
                           len(bins) - 1,
                           bins2))
    h_selected.append(ROOT.TH1F("h_selected_%i" % u,
                                ";#nu_{e} energy [GeV];Efficiency",
                                len(bins) - 1,
                                bins2))
    h_selected_numu.append(ROOT.TH1F("h_selected_numu_%i" % u,
                                     ";#nu_{e} energy [GeV];Efficiency",
                                     len(bins) - 1,
                                     bins2))
    h_selected_precuts.append(ROOT.TH1F("h_selected_precuts_%i" % u,
                                        ";#nu_{e} energy [GeV];Efficiency",
                                        len(bins) - 1,
                                        bins2))
    h_selected_cuts.append(ROOT.TH1F("h_selected_cuts_%i" % u,
                                     ";#nu_{e} energy [GeV];Efficiency",
                                     len(bins) - 1,
                                     bins2))
    h_selected_bdt.append(ROOT.TH1F("h_selected_bdt_%i" % u,
                                    ";#nu_{e} energy [GeV];Efficiency",
                                    len(bins) - 1,
                                    bins2))
var_dict = dict(variables + spectators)
ROOT.TMVA.Tools.Instance()
reader = ROOT.TMVA.Reader(":".join([
    "!V",
    "!Silent",
    "Color"]))
load_bdt(reader)

for i_evt in range(entries):
    printProgressBar(i_evt, entries, prefix="Progress:",
                     suffix="Complete", length=20)
    chain_nue.GetEntry(i_evt)

    protons = 0
    electrons = 0
    photons = 0
    pions = 0
    electron_energy = 0
    proton_energy = 0

    for i, energy in enumerate(chain_nue.nu_daughters_E):
        if chain_nue.nu_daughters_pdg[i] == 2212:
            proton_energy += energy - 0.938

            p_vertex = [chain_nue.nu_daughters_vx[i],
                        chain_nue.nu_daughters_vy[i],
                        chain_nue.nu_daughters_vz[i]]

            p_end = [chain_nue.nu_daughters_endx[i],
                     chain_nue.nu_daughters_endy[i],
                     chain_nue.nu_daughters_endz[i]]

            if energy - 0.938 > PROTON_THRESHOLD:
                protons += 1
            if not is_fiducial(p_vertex) or not is_fiducial(p_end):
                proton = 0
                break

        if chain_nue.nu_daughters_pdg[i] == 11:
            electron_energy += energy
            e_vertex = [chain_nue.nu_daughters_vx[i],
                        chain_nue.nu_daughters_vy[i],
                        chain_nue.nu_daughters_vz[i]]
            e_end = [chain_nue.nu_daughters_endx[i],
                     chain_nue.nu_daughters_endy[i],
                     chain_nue.nu_daughters_endz[i]]

            if energy - 0.51e-3 > ELECTRON_THRESHOLD and is_fiducial(e_vertex) and is_active(e_end):
                electrons += 1

        if chain_nue.nu_daughters_pdg[i] == 22:
            # if energy > 0.035:
            photons += 1

        if chain_nue.nu_daughters_pdg[i] == 111 or chain_nue.nu_daughters_pdg[i] == 211:
            # if energy > 0.06:
            pions += 1

    eNp = electrons == 1 and pions == 0 and protons >= 1

    true_neutrino_vertex = [chain_nue.true_vx_sce,
                            chain_nue.true_vy_sce,
                            chain_nue.true_vz_sce]
    true_neutrino_vertex_nosce = [chain_nue.true_vx,
                                  chain_nue.true_vy,
                                  chain_nue.true_vz]

    # if not eNp: print(electrons, photons, protons, pions)

    if is_fiducial(true_neutrino_vertex) and eNp and chain_nue.ccnc == 0:
        chain_filter.GetEntry(i_evt)

        flux_weights = [1]*len(chain_nue.flux_weights[0])
        for fl in chain_nue.flux_weights:
            flux_weights = [a*b for a, b in zip(flux_weights, fl)]


        for u in range(N_UNI):
            weight = chain_nue.genie_weights[0][u] * chain_nue.bnbweight

            if weight > 100:
                weight = chain_nue.bnbweight
            if u == 0:
                weight = chain_nue.bnbweight

            h_tot[u].Fill(chain_nue.nu_E, weight)

        contaminated = chain_nue.cosmic_fraction > 0.5 and chain_nue.category == 7
        selected = chain_nue.category == 2 or contaminated
        if not selected:
            continue


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

        if not track_fidvol:
            continue

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

        if not shower_fidvol:
            continue

        neutrino_vertex = [chain_nue.vx, chain_nue.vy, chain_nue.vz]

        if not is_fiducial(neutrino_vertex):
            continue

        # If there are no tracks we require at least two showers
        showers_2_tracks_0 = True
        if chain_nue.n_tracks == 0 and chain_nue.n_showers == 1:
            showers_2_tracks_0 = False
        if not showers_2_tracks_0:
            continue

        for u in range(N_UNI):
            weight = chain_nue.genie_weights[0][u] * chain_nue.bnbweight
            if weight > 100:
                weight = chain_nue.bnbweight
            if u == 0:
                weight = chain_nue.bnbweight
            h_selected[u].Fill(chain_nue.nu_E, weight)

        fill_kin_branches(chain_nue, 1, var_dict, "nue", True)

        if not int(var_dict["numu_score"][0]) == 0:
            continue

        for u in range(N_UNI):
            weight = chain_nue.genie_weights[0][u] * chain_nue.bnbweight
            if weight > 100:
                weight = chain_nue.bnbweight
            if u == 0:
                weight = chain_nue.bnbweight
            h_selected_numu[u].Fill(chain_nue.nu_E, weight)

        if not pre_cuts(var_dict):
            continue

        bdt_values = {}
        for bdt_name in bdt_types:
            bdt_values[bdt_name] = reader.EvaluateMVA("BDT%s" % bdt_name)

        for u in range(N_UNI):
            weight = chain_nue.genie_weights[0][u] * chain_nue.bnbweight
            if weight > 100:
                weight = chain_nue.bnbweight
            if u == 0:
                weight = chain_nue.bnbweight

            h_selected_precuts[u].Fill(chain_nue.nu_E, weight)

            if apply_cuts(bdt_values, var_dict, bdt=False, manual=True):
                h_selected_cuts[u].Fill(chain_nue.nu_E, weight)
            if apply_cuts(bdt_values, var_dict, bdt=True, manual=False):
                h_selected_bdt[u].Fill(chain_nue.nu_E, weight)

h_eff_2d = ROOT.TH2F("h_eff_2d",
                     ";#nu_{e} energy [GeV];Efficiency",
                     len(bins) - 1,
                     bins2,
                     N_UNI, 0, 1)
sys_err = [0] * (h_selected[0].GetNbinsX() + 1)
eff_cv = ROOT.TEfficiency(h_selected[0], h_tot[0])
for u in range(1, N_UNI):
    if not h_selected[u].Integral() > 0:
        continue

    eff_u = ROOT.TEfficiency(h_selected[u], h_tot[u])

    for i_bin in range(1, h_selected[u].GetNbinsX() + 1):
        value = eff_u.GetEfficiency(i_bin)
        center = h_selected[u].GetBinCenter(i_bin)

        diff = (eff_cv.GetEfficiency(i_bin) -
                eff_u.GetEfficiency(i_bin))

        # if n == "reco_energy":
        #     bin_width = h_sys[n][u].GetBinWidth(bin)
        #     diff = (h_cv[n].GetBinContent(bin) / (bin_width / 0.05) -
        #             h_sys[n][u].GetBinContent(bin) / (bin_width / 0.05))
        #     value /= bin_width / 0.05

        h_eff_2d.Fill(center, value)
        sys_err[i_bin] += diff**2


if __name__ == "__main__":
    c_e = ROOT.TCanvas("c_e")
    # h_eff_2d.Draw("colz")
    h_eff_2d.GetZaxis().SetRangeUser(0, 50)
    eff_e = ROOT.TEfficiency(h_selected[0], h_tot[0])
    eff_e.SetMarkerStyle(0)
    eff_e.SetLineColor(ROOT.kRed + 1)
    eff_e.SetLineWidth(2)
    eff_e.Draw("AP")
    # eff_e.GetPaintedGraph().GetYaxis().SetRangeUser(0.001, 1)

    eff_e_numu = ROOT.TEfficiency(h_selected_numu[0], h_tot[0])
    eff_e_numu.SetMarkerStyle(0)
    eff_e_numu.SetLineColor(ROOT.kBlue + 1)
    eff_e_numu.SetLineWidth(2)
    eff_e_numu.Draw("P SAME")

    eff_e_precuts = ROOT.TEfficiency(h_selected_precuts[0], h_tot[0])
    eff_e_precuts.SetMarkerStyle(0)
    eff_e_precuts.SetLineColor(ROOT.kGreen + 1)
    eff_e_precuts.SetLineWidth(2)
    eff_e_precuts.Draw("P SAME")

    eff_e_cuts = ROOT.TEfficiency(h_selected_cuts[0], h_tot[0])
    eff_e_cuts.SetMarkerStyle(0)
    eff_e_cuts.SetLineColor(ROOT.kOrange + 1)
    eff_e_cuts.SetLineWidth(2)
    eff_e_cuts.Draw("P SAME")

    eff_e_bdt = ROOT.TEfficiency(h_selected_bdt[0], h_tot[0])
    eff_e_bdt.SetMarkerStyle(0)
    eff_e_bdt.SetLineColor(ROOT.kMagenta + 1)
    eff_e_bdt.SetLineWidth(2)
    eff_e_bdt.Draw("P SAME")

    leg = ROOT.TLegend(0.4,0.64,0.83,0.84)
    leg.AddEntry(eff_e, "Topology and flash requirements: %.1f %%" % (h_selected[0].Integral()/h_tot[0].Integral()*100), "le")
    leg.AddEntry(eff_e_numu, "#nu_{#mu} rejection: %.1f %%" % (h_selected_numu[0].Integral()/h_tot[0].Integral()*100), "le")
    leg.AddEntry(eff_e_precuts, "Quality precuts: %.1f %%" % (h_selected_precuts[0].Integral()/h_tot[0].Integral()*100), "le")
    leg.AddEntry(eff_e_cuts, "Rectangular cuts: %.1f %%" % (h_selected_cuts[0].Integral()/h_tot[0].Integral()*100), "le")
    leg.AddEntry(eff_e_bdt, "BDT: %.1f %%" % (h_selected_bdt[0].Integral()/h_tot[0].Integral()*100), "le")

    leg.Draw()
    c_e.Update()

    input()
