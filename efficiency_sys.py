#!/usr/local/bin/python3

import math
import ROOT

from glob import glob
from bdt_common import x_start, x_end, y_start, y_end, z_start, z_end, printProgressBar, bins, bins2, distance, is_fiducial, is_active, N_UNI, variables, spectators
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
h_selected_precuts = []

for u in range(N_UNI):
    h_tot.append(ROOT.TH1F("h_tot_%i" % u, ";#nu_{e} energy [GeV];Efficiency", len(bins) - 1, bins2))
    h_selected.append(ROOT.TH1F("h_selected_%i" % u, ";#nu_{e} energy [GeV];Efficiency", len(bins) - 1, bins2))
    h_selected_numu.append(ROOT.TH1F("h_selected_numu_%i" % u, ";#nu_{e} energy [GeV];Efficiency", len(bins) - 1, bins2))
    h_selected_precuts.append(ROOT.TH1F("h_selected_precuts_%i" % u, ";#nu_{e} energy [GeV];Efficiency", len(bins) - 1, bins2))


variables = dict(variables + spectators)

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
            weight = flux_weights[u] * chain_nue.genie_weights[0][u] * chain_nue.bnbweight

            if weight > 100:
                weight = chain_nue.bnbweight
            if u == 0:
                weight = chain_nue.bnbweight

            h_tot[u].Fill(chain_nue.nu_E, weight)

        if not chain_nue.passed:
            continue

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
            weight = flux_weights[u] * chain_nue.genie_weights[0][u] * chain_nue.bnbweight
            if weight > 100:
                weight = chain_nue.bnbweight
            if u == 0:
                weight = chain_nue.bnbweight
            h_selected[u].Fill(chain_nue.nu_E, weight)

        fill_kin_branches(chain_nue, 1, variables, "nue")
        hits = variables["track_hits"][0] > 5 and variables["shower_hits"][0] > 5 and variables["total_hits_y"][0] > 0 and variables["total_hits_u"][0] > 0 and variables["total_hits_v"][0] > 0
        sh_id = int(variables["shower_id"][0])
        shower_track_energy = variables["total_shower_energy"][0] > 0.01 and variables["total_track_energy_length"][0] > 0 and variables["shower_energy"][sh_id] > 0.01

        for u in range(N_UNI):
            weight = flux_weights[u] * chain_nue.genie_weights[0][u] * chain_nue.bnbweight
            if weight > 100:
                weight = chain_nue.bnbweight
            if u == 0:
                weight = chain_nue.bnbweight

            if not shower_track_energy:
                continue
            if not hits:
                continue

            h_selected_precuts[u].Fill(chain_nue.nu_E, weight)

            if chain_nue.numu_passed == 1:
                continue

            h_selected_numu[u].Fill(chain_nue.nu_E, weight)


h_eff_2d = ROOT.TH2F("h_eff_2d", ";#nu_{e} energy [GeV];Efficiency", len(bins) - 1, bins2, N_UNI, 0, 1)
sys_err = [0] * (h_selected[0].GetNbinsX() + 1)
eff_cv = ROOT.TEfficiency(h_selected[0], h_tot[0])
for u in range(1, N_UNI):
    if not h_selected[u].Integral() > 0:
        continue

    eff_u = ROOT.TEfficiency(h_selected[u], h_tot[u])

    for bin in range(1, h_selected[u].GetNbinsX() + 1):
        value = eff_u.GetEfficiency(bin)
        center = h_selected[u].GetBinCenter(bin)

        diff = (eff_cv.GetEfficiency(bin) -
                eff_u.GetEfficiency(bin))

        # if n == "reco_energy":
        #     bin_width = h_sys[n][u].GetBinWidth(bin)
        #     diff = (h_cv[n].GetBinContent(bin) / (bin_width / 0.05) -
        #             h_sys[n][u].GetBinContent(bin) / (bin_width / 0.05))
        #     value /= bin_width / 0.05

        h_eff_2d.Fill(center, value)
        sys_err[bin] += diff**2


if __name__ == "__main__":
    c_e = ROOT.TCanvas("c_e")
    h_eff_2d.Draw("colz")
    h_eff_2d.GetZaxis().SetRangeUser(0, 50)
    eff_e = ROOT.TEfficiency(h_selected[0], h_tot[0])
    eff_e.SetMarkerStyle(0)
    eff_e.SetLineColor(ROOT.kRed + 1)
    eff_e.SetLineWidth(2)
    eff_e.Draw("P SAME")

    eff_e_numu = ROOT.TEfficiency(h_selected_numu[0], h_tot[0])
    eff_e_numu.SetMarkerStyle(21)
    eff_e_numu.SetMarkerColor(ROOT.kBlue + 1)
    # eff_e_numu.Draw("P SAME")

    eff_e_precuts = ROOT.TEfficiency(h_selected_precuts[0], h_tot[0])
    eff_e_precuts.SetMarkerStyle(22)
    eff_e_precuts.SetMarkerColor(ROOT.kGreen + 1)
    # eff_e_precuts.Draw("P SAME")
    c_e.Update()

    input()