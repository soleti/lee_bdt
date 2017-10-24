#!/usr/local/bin/python3

import math
import ROOT
from glob import glob
from bdt_common import x_start, x_end, y_start, y_end, z_start, z_end

ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetPalette(87)
ROOT.gStyle.SetNumberContours(99)


def is_fiducial(point):
    ok_y = point[1] > y_start + 20 and point[1] < y_end - 20
    ok_x = point[0] > x_start + 10 and point[0] < x_end - 10
    ok_z = point[2] > z_start + 10 and point[2] < z_end - 50
    return ok_y and ok_x and ok_z


def is_active(point):
    ok_y = point[1] > y_start and point[1] < y_end
    ok_x = point[0] > x_start and point[0] < x_end
    ok_z = point[2] > z_start and point[2] < z_end
    return ok_y and ok_x and ok_z


nue_cosmic = glob("mc_nue_old/*/Pandora*.root")
#nue_cosmic = glob("softmerge.root")
chain_nue = ROOT.TChain("robertoana/pandoratree")

for f in nue_cosmic:
    chain_nue.Add(f)

e_energy = ROOT.TEfficiency("e_energy",
                            ";#nu_{e} energy [GeV];Fraction", 20, 0, 2)
ep_energy = ROOT.TEfficiency("e_energy",
                             ";#nu_{e} energy [GeV];Efficiency #times Purity",
                             20, 0, 2)

e_proton = ROOT.TEfficiency("e_proton",
                            ";p kinetic energy [GeV];#epsilon #times P_{reco}",
                            20, 0, 0.5)

h_dedx_electron = ROOT.TH1F("h_dedx_electron",
                            "Electrons;dE/dx [MeV/cm];Area normalized",
                            50, 0, 5)
h_dedx_photon = ROOT.TH1F("h_dedx_photon",
                          "Photons;dE/dx [MeV/cm];Area normalized",
                          50, 0, 5)

p_energy = ROOT.TEfficiency("p_energy",
                            ";#nu_{e} energy [GeV];Purity",
                            20, 0, 2)
p_dist_energy = ROOT.TEfficiency("p_dist_energy",
                                 ";#nu_{e} energy [GeV];Purity",
                                 20, 0, 2)
ep_dist_energy = ROOT.TEfficiency("p_dist_energy",
                                  ";#nu_{e} energy [GeV];Efficiency #times \
                                  Purity",
                                  20, 0, 2)

l_e_proton = ROOT.TH2F("l_e_proton",
                       ";Reco. track length [cm];True p kinetic energy [GeV]",
                       100, 0, 50, 100, 0, 0.5)

h_dist = ROOT.TH1F("h_dist", ";Distance [cm];N. Entries / 0.2 cm", 50, 0, 10)
h_dist_nosce = ROOT.TH1F("h_dist_nosce",
                         ";Distance [cm];N. Entries / 0.2 cm",
                         50, 0, 10)

h_x_diff = ROOT.TH1F("h_x_diff",
                     ";#Delta x [cm]; N. Entries / 0.2 cm", 50, -5, 5)
h_y_diff = ROOT.TH1F("h_y_diff",
                     ";#Delta y [cm]; N. Entries / 0.2 cm", 50, -5, 5)
h_z_diff = ROOT.TH1F("h_z_diff",
                     ";#Delta z [cm]; N. Entries / 0.2 cm", 50, -5, 5)

h_x = ROOT.TH1F("h_x",
                     ";#Delta x [cm]; N. Entries / 0.2 cm", 50, x_start-40, x_end+40)
h_y = ROOT.TH1F("h_y",
                     ";#Delta y [cm]; N. Entries / 0.2 cm", 50, y_start-40, y_end+40)
h_z = ROOT.TH1F("h_z",
                     ";#Delta z [cm]; N. Entries / 0.2 cm", 50, z_start-40, z_end+40)


fiducial = 0
eNp_events = 0

passed = 0
not_passed = 0
flash_passed = 0
reco_ok = 0
vertex_ok = 0
noflash = 0
entries = chain_nue.GetEntries()
cc_events = 0

for i in range(entries):
    chain_nue.GetEntry(i)

    protons = 0
    electrons = 0
    photons = 0
    pions = 0
    electron_energy = 0
    proton_energy = 0
    for i, energy in enumerate(chain_nue.nu_daughters_E):
        if chain_nue.nu_daughters_pdg[i] == 2212:
            proton_energy += energy
            if energy - 0.938 > 0.000005:
                protons += 1

        if chain_nue.nu_daughters_pdg[i] == 11:
            electron_energy += energy
            if energy > 0.000003:
                electrons += 1

        if chain_nue.nu_daughters_pdg[i] == 22:
            # if energy > 0.035:
            photons += 1

        if chain_nue.nu_daughters_pdg[i] == 111:
            # if energy > 0.06:
            pions += 1

    eNp = electrons == 1 and photons == 0 and pions == 0 and protons > 0

    true_neutrino_vertex = [chain_nue.true_vx_sce,
                            chain_nue.true_vy_sce,
                            chain_nue.true_vz_sce]

    if chain_nue.category != 4:
        cc_events += 1

    if not eNp: print(electrons, photons, protons, pions)

    if eNp and chain_nue.nu_E > 0.1:
        eNp_events += 1
        h_x.Fill(true_neutrino_vertex[0])
        h_y.Fill(true_neutrino_vertex[1])
        h_z.Fill(true_neutrino_vertex[2])

        if is_fiducial(true_neutrino_vertex):
            fiducial += 1

            p = False
            p_track = False
            p_shower = False

            if protons == chain_nue.nu_matched_tracks:
                p_track = True

            if electrons == chain_nue.nu_matched_showers:
                p_shower = True

            if p_track and p_shower:
                p = True

            proton_energy = sum([chain_nue.nu_daughters_E[i]
                                 for i, pdg in
                                 enumerate(chain_nue.nu_daughters_pdg)
                                 if pdg == 2212])

            electron_energy = max([chain_nue.nu_daughters_E[i]
                                   for i, pdg in
                                   enumerate(chain_nue.nu_daughters_pdg)
                                   if abs(pdg) == 11])

            neutrino_vertex = [chain_nue.vx, chain_nue.vy, chain_nue.vz]

            true_neutrino_vertex = [chain_nue.true_vx_sce,
                                    chain_nue.true_vy_sce,
                                    chain_nue.true_vz_sce]
            true_neutrino_vertex_nosce = [chain_nue.true_vx,
                                          chain_nue.true_vy,
                                          chain_nue.true_vz]

            dist = 0

            if chain_nue.passed:
                dist = math.sqrt(sum([(t - r) ** 2
                                      for t, r in zip(neutrino_vertex,
                                                      true_neutrino_vertex)]))

                passed += 1
                if p:
                    reco_ok += 1
                if chain_nue.distance < 5:
                    vertex_ok += 1

                h_dist.Fill(dist)
                h_dist_nosce.Fill(math.sqrt(
                    sum([(t - r)**2 for t, r in zip(
                        neutrino_vertex, true_neutrino_vertex_nosce)])
                ))

                h_x_diff.Fill(neutrino_vertex[0] - true_neutrino_vertex[0])
                h_y_diff.Fill(neutrino_vertex[1] - true_neutrino_vertex[1])
                h_z_diff.Fill(neutrino_vertex[2] - true_neutrino_vertex[2])

                p_energy.Fill(p_track and p_shower, chain_nue.nu_E)
                p_dist_energy.Fill(chain_nue.distance < 5, chain_nue.nu_E)

            else:
                not_passed += 1
                if chain_nue.flash_passed != 1:
                    noflash += 1

            if protons == 1:
                e_proton.Fill(chain_nue.passed and p_track and p_shower,
                              proton_energy - 0.938)

            if chain_nue.passed and protons == 1 and dist < 2:
                l_e_proton.Fill(chain_nue.track_len[0], proton_energy - 0.938)

            ep_energy.Fill(chain_nue.passed and p_track and p_shower,
                           chain_nue.nu_E)
            ep_dist_energy.Fill(chain_nue.passed and chain_nue.distance < 5,
                                chain_nue.nu_E)
            e_energy.Fill(chain_nue.passed, chain_nue.nu_E)


print("Entries", entries)
print("CC", cc_events)
print("1eNp", eNp_events)
print("1eNp + Is fiducial", fiducial)

print("Passed", passed)
print("Not passed", not_passed)

eff = passed / fiducial
eff_err = math.sqrt((eff * (1 - eff)) / fiducial)

p_reco = reco_ok / passed
p_reco_err = math.sqrt((p_reco * (1 - p_reco)) / passed)
p_vertex = vertex_ok / passed
p_vertex_err = math.sqrt((p_vertex * (1 - p_vertex)) / passed)

ep_reco = reco_ok / fiducial
ep_reco_err = math.sqrt((ep_reco * (1 - ep_reco)) / fiducial)
ep_vertex = vertex_ok / fiducial
ep_vertex_err = math.sqrt((ep_vertex * (1 - ep_vertex)) / fiducial)

print("Efficiency: ({0:.1f} +- {1:.1f}) %".format(eff * 100, eff_err * 100))
print("Reco. purity: ({0:.1f} +- {1:.1f}) %"
      .format(p_reco * 100, p_reco_err * 100))
print("Vertex purity: ({0:.1f} +- {1:.1f}) %"
      .format(p_vertex * 100, p_vertex_err * 100))
print("Reco. efficiency x purity: ({0:.1f} +- {1:.1f}) %"
      .format(ep_reco * 100, ep_reco_err * 100))
print("Vertex efficiency x purity: ({0:.1f} +- {1:.1f}) %"
      .format(ep_vertex * 100, ep_vertex_err * 100))
f_energy = ROOT.TFile("plots/f_energy.root", "RECREATE")
e_energy.Write()
f_energy.Close()

pt = ROOT.TPaveText(0.1, 0.91, 0.45, 0.97, "ndc")
pt.AddText("MicroBooNE Preliminary")
pt.SetFillColor(0)
pt.SetBorderSize(0)
pt.SetShadowColor(0)

legend = ROOT.TLegend(0.13, 0.68, 0.84, 0.86)
legend.SetTextSize(16)
legend.AddEntry(e_energy, "#epsilon ({0:.1f} #pm {1:.1f}) %"
                .format(eff * 100, eff_err * 100),
                "lep")

legend.AddEntry(e_energy, "", "")
legend.AddEntry(p_energy, "P_{{reco}} ({0:.1f} #pm {1:.1f}) %"
                .format(p_reco * 100, p_reco_err * 100),
                "lep")

legend.AddEntry(p_dist_energy, "P_{{vertex}} ({0:.1f} #pm {1:.1f}) %"
                .format(p_vertex * 100, p_vertex_err * 100), "lep")
legend.AddEntry(ep_energy, "#epsilon #times P_{{reco}} ({0:.1f} #pm {1:.1f}) %"
                .format(ep_reco * 100, ep_reco_err * 100), "lep")
legend.AddEntry(
    ep_dist_energy, "#epsilon #times P_{{vertex}} ({0:.1f} #pm {1:.1f}) %"
    .format(ep_vertex * 100, ep_vertex_err * 100), "lep")

legend.SetNColumns(2)

c_energy = ROOT.TCanvas("c_energy")
e_energy.Draw("apl")
p_energy.SetMarkerStyle(22)
p_energy.SetLineColor(ROOT.kGreen + 1)
p_energy.SetLineWidth(2)
f_p_energy = ROOT.TFile("plots/f_p_energy.root", "RECREATE")
p_energy.Write()
f_p_energy.Close()
p_energy.Draw("pl same")

p_dist_energy.SetMarkerStyle(22)
p_dist_energy.SetLineColor(ROOT.kGreen + 3)
p_dist_energy.SetLineWidth(2)
f_p_dist_energy = ROOT.TFile("plots/f_p_dist_energy.root", "RECREATE")
p_dist_energy.Write()
f_p_dist_energy.Close()
p_dist_energy.Draw("pl same")

ep_energy.SetMarkerStyle(23)
ep_energy.SetLineColor(ROOT.kBlue + 1)
ep_energy.SetLineWidth(2)
f_ep_energy = ROOT.TFile("plots/f_ep_energy.root", "RECREATE")
ep_energy.Write()
f_ep_energy.Close()
ep_energy.Draw("pl same")

ep_dist_energy.SetMarkerStyle(23)
ep_dist_energy.SetLineColor(ROOT.kBlue + 3)
ep_dist_energy.SetLineWidth(2)
f_ep_dist_energy = ROOT.TFile("plots/f_ep_dist_energy.root", "RECREATE")
ep_dist_energy.Write()
f_ep_dist_energy.Close()
ep_dist_energy.Draw("pl same")

e_energy.SetMarkerStyle(20)
e_energy.SetLineColor(ROOT.kRed + 1)
e_energy.SetLineWidth(2)
c_energy.Update()
e_energy.GetPaintedGraph().GetXaxis().SetRangeUser(0.1, 2)
e_energy.GetPaintedGraph().GetYaxis().SetRangeUser(0, 1.39)
pt.Draw()
legend.Draw()
c_energy.SaveAs("plots/efficiency.pdf")
c_energy.Draw()

c_proton = ROOT.TCanvas("c_proton")
e_proton.Draw("apl")
c_proton.Update()

c_lproton = ROOT.TCanvas("c_lproton")
l_e_proton.Draw("colz")
c_lproton.Update()

c_dist = ROOT.TCanvas("c_dist")
h_dist.Draw()
h_dist_nosce.SetLineColor(ROOT.kRed + 1)
h_dist_nosce.Draw("same")
c_dist.Update()

c_x = ROOT.TCanvas("c_x", "", 500, 500)
h_x.Draw()
c_x.Update()

c_y = ROOT.TCanvas("c_y", "", 500, 500)
h_y.Draw()
c_y.Update()

c_z = ROOT.TCanvas("c_z", "", 500, 500)
h_z.Draw()
c_z.Update()


input()
