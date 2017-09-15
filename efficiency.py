#!/usr/bin/env python3.4

import math
import ROOT
from glob import glob

ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetPalette(87)
ROOT.gStyle.SetNumberContours(99)

nue_cosmic = glob("mc_nue_nofidvol/*/Pandora*.root")[:1]
print(nue_cosmic)
chain_nue = ROOT.TChain("robertoana/pandoratree")
chain_nue_pot = ROOT.TChain("robertoana/pot")

for f in nue_cosmic:
    chain_nue.Add(f)
    chain_nue_pot.Add(f)

# total_pot = 6.6e20
# total_nue_pot = 0
# for i in range(chain_nue_pot.GetEntries()):
#     chain_nue_pot.GetEntry(i)
#     total_nue_pot += chain_nue_pot.pot
# print("Total POT v_e", total_nue_pot)

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

is_fiducial = 0
eNp = 0
yesshower_notrack = 0
noshower_yestrack = 0
noshower_notrack = 0
passed = 0
not_passed = 0
flash_passed = 0
reco_ok = 0
vertex_ok = 0
noflash = 0
entries = chain_nue.GetEntries()
for i in range(entries):
    chain_nue.GetEntry(i)

    protons = 0
    electrons = 0
    photons = 0
    pions = 0
    for i, energy in enumerate(chain_nue.nu_daughters_E):
        if chain_nue.nu_daughters_pdg[i] == 2212:
            if energy - 0.938 > 0.0005:
                protons += 1

        if chain_nue.nu_daughters_pdg[i] == 11:
            if energy > 0.03:
                electrons += 1

        if chain_nue.nu_daughters_pdg[i] == 22:
            #if energy > 0.035:
            photons += 1

        if chain_nue.nu_daughters_pdg[i] == 111:
            #if energy > 0.06:
            pions += 1

    if electrons > 0 and photons == 0 and pions == 0 and protons > 0:
        eNp += 1

        if chain_nue.true_nu_is_fiducial:
            is_fiducial += 1

            p = False
            p_track = False
            p_shower = False

            if protons == chain_nue.nu_matched_tracks:
                p_track = True

            if electrons == chain_nue.nu_matched_showers:
                p_shower = True

            if p_track and p_shower:
                p = True

            proton_energy = sum([chain_nue.nu_daughters_E[i] for i,pdg in enumerate(chain_nue.nu_daughters_pdg) if pdg == 2212])
            electron_energy = max([chain_nue.nu_daughters_E[i] for i,pdg in enumerate(chain_nue.nu_daughters_pdg) if abs(pdg) == 11])

            neutrino_vertex = [chain_nue.vx,chain_nue.vy,chain_nue.vz]

            true_neutrino_vertex = [chain_nue.true_vx_sce,chain_nue.true_vy_sce,chain_nue.true_vz_sce]
            true_neutrino_vertex_nosce = [chain_nue.true_vx,chain_nue.true_vy,chain_nue.true_vz]

            dist = 0

            if chain_nue.passed:
                dist = math.sqrt(sum([(t-r) ** 2 for t, r in zip(neutrino_vertex,true_neutrino_vertex)]))

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

                h_x_diff.Fill(neutrino_vertex[0]-true_neutrino_vertex[0])
                h_y_diff.Fill(neutrino_vertex[1]-true_neutrino_vertex[1])
                h_z_diff.Fill(neutrino_vertex[2]-true_neutrino_vertex[2])

                # print("proton", protons, chain_nue.nu_matched_tracks)
                # print("electrons", electrons, chain_nue.nu_matched_showers)
                p_energy.Fill(p_track and p_shower, chain_nue.nu_E)
                p_dist_energy.Fill(chain_nue.distance < 5, chain_nue.nu_E)

            else:
                not_passed += 1
                if chain_nue.flash_passed > 0:
                    if chain_nue.shower_passed > 0 and chain_nue.track_passed <= 0:
                        yesshower_notrack += 1
                    if chain_nue.track_passed > 0 and chain_nue.shower_passed <= 0:
                        noshower_yestrack += 1
                    if chain_nue.track_passed <= 0 and chain_nue.shower_passed <= 0:
                        noshower_notrack += 1
                if chain_nue.flash_passed == 1 and chain_nue.track_passed == 1 and chain_nue.shower_passed == 1: print(chain_nue.run, chain_nue.subrun, chain_nue.event)
                else:
                    noflash += 1


            if protons == 1:
                e_proton.Fill(chain_nue.passed and p_track and p_shower, proton_energy-0.938)

            if chain_nue.passed and protons == 1 and chain_nue.nu_matched_tracks == 1 and dist < 2:
                l_e_proton.Fill(chain_nue.track_len[0], proton_energy-0.938)

            ep_energy.Fill(chain_nue.passed and p_track and p_shower, chain_nue.nu_E)
            ep_dist_energy.Fill(chain_nue.passed and chain_nue.distance < 5, chain_nue.nu_E)
            e_energy.Fill(chain_nue.passed, chain_nue.nu_E)


print("Entries", entries)
print("1eNp", eNp)
print("1eNp + Is fiducial", is_fiducial)

print("Passed", passed)
print("Not passed", not_passed)

print("No flash {:.1f} %".format(
    noflash / not_passed * 100))
print("Yes shower no track {:.1f} %".format(
    yesshower_notrack / not_passed * 100))
print("No shower yes track {:.1f} %".format(
    noshower_yestrack / not_passed * 100))
print("No shower no track {:.1f} %".format(
    noshower_notrack / not_passed * 100))


eff = passed / is_fiducial
eff_err = math.sqrt((eff * (1 - eff)) / eNp)

p_reco = reco_ok / passed
p_reco_err = math.sqrt((p_reco * (1 - p_reco)) / passed)
p_vertex = vertex_ok / passed
p_vertex_err = math.sqrt((p_vertex * (1 - p_vertex)) / passed)

ep_reco = reco_ok / is_fiducial
ep_reco_err = math.sqrt((ep_reco * (1 - ep_reco)) / is_fiducial)
ep_vertex = vertex_ok / is_fiducial
ep_vertex_err = math.sqrt((ep_vertex * (1 - ep_vertex)) / is_fiducial)

print("Efficiency: ({0:.1f} +- {1:.1f}) %".format(eff * 100, eff_err * 100))
print("Reco. purity: ({0:.1f} +- {1:.1f}) %"
      .format(p_reco * 100, p_reco_err * 100))
print("Vertex purity: ({0:.1f} +- {1:.1f}) %"
      .format(p_vertex * 100, p_vertex_err * 100))
print("Reco. efficiency x purity: ({0:.1f} +- {1:.1f}) %"
      .format(ep_reco * 100, ep_reco_err * 100))
print("Vertex efficiency x purity: ({0:.1f} +- {1:.1f}) %"
      .format(ep_vertex * 100, ep_vertex_err * 100))
f_energy = ROOT.TFile("f_energy.root", "RECREATE")
e_energy.Write()
f_energy.Close()

pt = ROOT.TPaveText(0.1, 0.91, 0.45, 0.97, "ndc")
pt.AddText("MicroBooNE Preliminary")
pt.SetFillColor(0)
pt.SetBorderSize(0)
pt.SetShadowColor(0)

legend = ROOT.TLegend(0.13, 0.68, 0.84, 0.86)
legend.SetTextSize(16)
legend.AddEntry(e_energy, "#epsilon ({0:.1f} \
                #pm {1:.1f}) %".format(eff * 100, eff_err * 100),
                "lep")

legend.AddEntry(e_energy, "", "")
legend.AddEntry(p_energy, "P_{{reco}} \({0:.1f} #pm {1:.1f}) %"
                .format(p_reco * 100, p_reco_err * 100),
                "lep")

legend.AddEntry(p_dist_energy, "P_{{vertex}} ({0:.1f} #pm {1:.1f}) %"
                .format(p_vertex * 100, p_vertex_err * 100), "lep")
legend.AddEntry(ep_energy, "#epsilon #times P_{{reco}} ({0:.1f} #pm {1:.1f}) %"
                .format(ep_reco * 100, ep_reco_err * 100), "lep")
legend.AddEntry(ep_dist_energy, "#epsilon #times P_{{vertex}} \
                ({0:.1f} #pm {1:.1f}) %"
                .format(ep_vertex * 100, ep_vertex_err * 100), "lep")

legend.SetNColumns(2)

c_energy = ROOT.TCanvas("c_energy")
e_energy.Draw("apl")
p_energy.SetMarkerStyle(22)
p_energy.SetLineColor(ROOT.kGreen + 1)
p_energy.SetLineWidth(2)
p_energy.Draw("pl same")

p_dist_energy.SetMarkerStyle(22)
p_dist_energy.SetLineColor(ROOT.kGreen + 3)
p_dist_energy.SetLineWidth(2)
p_dist_energy.Draw("pl same")

ep_energy.SetMarkerStyle(23)
ep_energy.SetLineColor(ROOT.kBlue + 1)
ep_energy.SetLineWidth(2)
ep_energy.Draw("pl same")

ep_dist_energy.SetMarkerStyle(23)
ep_dist_energy.SetLineColor(ROOT.kBlue+3)
ep_dist_energy.SetLineWidth(2)
ep_dist_energy.Draw("pl same")

e_energy.SetMarkerStyle(20)
e_energy.SetLineColor(ROOT.kRed+1)
e_energy.SetLineWidth(2)
c_energy.Update()
e_energy.GetPaintedGraph().GetXaxis().SetRangeUser(0.1,2)
e_energy.GetPaintedGraph().GetYaxis().SetRangeUser(0,1.39)
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
h_dist_nosce.SetLineColor(ROOT.kRed+1)
h_dist_nosce.Draw("same")
c_dist.Update()

c_x = ROOT.TCanvas("c_x","",500,500)
h_x_diff.Draw()
c_x.Update()

c_y = ROOT.TCanvas("c_y","",500,500)
h_y_diff.Draw()
c_y.Update()

c_z = ROOT.TCanvas("c_z", "", 500, 500)
h_z_diff.Draw()
c_z.Update()


input()
