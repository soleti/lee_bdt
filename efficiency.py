#!/usr/local/bin/python3

import math
import ROOT
import random
from glob import glob
from bdt_common import x_start, x_end, y_start, y_end, z_start, z_end
from array import array
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetPalette(87)
ROOT.gStyle.SetNumberContours(99)

# Orange #efac3a
# Red #c46d27
# Green #62b570
# Violet #8779b1
# Azure #609ec6
colors = [ROOT.TColor.GetColor("#e7623d"), ROOT.TColor.GetColor("#c46d27"),
          ROOT.TColor.GetColor("#62b570"), ROOT.TColor.GetColor("#8779b1"),
          ROOT.TColor.GetColor("#609ec6"), ROOT.TColor.GetColor("#c46d27"),
          ROOT.TColor.GetColor("#ffffff"), ROOT.TColor.GetColor("#e7623d")]

pot = 6.6e20

def fix_binning(histo):
    width = histo.GetBinWidth(1)
    for k in range(1, histo.GetNbinsX() + 1):
        bin_width = histo.GetBinWidth(k)
        histo.SetBinError(k, histo.GetBinError(k) / (bin_width / width))
        histo.SetBinContent(k, histo.GetBinContent(k) / (bin_width / width))


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


nue_cosmic = glob("mc_bnb_mcc86_2/*/Pandora*.root") + glob("mc_nue_mcc86_2/*/*.root")
random.shuffle(nue_cosmic)
# nue_cosmic = glob("softmerge.root")
chain_nue = ROOT.TChain("robertoana/pandoratree")
chain_pot = ROOT.TChain("robertoana/pot")

for f in nue_cosmic:
    chain_nue.Add(f)
    chain_pot.Add(f)

print("entries",chain_nue.GetEntries())
e_energy = ROOT.TEfficiency("e_energy",
                            ";#nu_{e} E [GeV];Fraction", 20, 0, 2)
ep_energy = ROOT.TEfficiency("e_energy",
                             ";#nu_{e} energy [GeV];Efficiency #times Purity",
                             20, 0, 2)

e_proton = ROOT.TEfficiency("e_proton",
                            ";p kinetic energy [GeV];#epsilon #times P_{reco}",
                            20, 0, 0.5)

e_nprotons = ROOT.TEfficiency("e_nprotons",
                            ";N. of protons;Fraction",
                            4, 0.5, 4.5)


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

h_x = ROOT.TH1F("h_x", ";#Delta x [cm]; N. Entries / 0.2 cm",
                50, x_start - 40, x_end + 40)
h_y = ROOT.TH1F("h_y", ";#Delta y [cm]; N. Entries / 0.2 cm",
                50, y_start - 40, y_end + 40)
h_z = ROOT.TH1F("h_z", ";#Delta z [cm]; N. Entries / 0.2 cm",
                50, z_start - 40, z_end + 40)
bins = array("f", [0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.6, 0.8, 1])

h_energy = ROOT.TH1F("h_energy",";Reco. energy [GeV]; N. Entries / 0.05 GeV", 10,0.2,1)
fiducial = 0
eNp_events = 0

h_total_electron = ROOT.TH1F(
    "h_total_electron", ";E_{e} [GeV];Efficiency", 20, 0, 0.2)
h_passed_electron = ROOT.TH1F(
    "h_passed_electron", ";E_{e} [GeV];Efficiency", 20, 0, 0.2)
h_total_proton = ROOT.TH1F(
    "h_total_proton", ";E_{p} [GeV];Efficiency", 20, 0, 0.2)
h_passed_proton = ROOT.TH1F(
    "h_passed_proton", ";E_{p} [GeV];Efficiency", 20, 0, 0.2)
h_total_2d = ROOT.TH2F(
    "h_total_2d", ";E_{p} [GeV];E_{e} [GeV]", 40, 0, 1, 40, 0, 1)
h_passed_2d = ROOT.TH2F(
    "h_passed_2d", ";E_{p} [GeV];E_{e} [GeV]", 40, 0, 1, 40, 0, 1)

description = ["Other",
               "Cosmic",
               "Beam Intrinsic #nu_{e}",
               "Beam Intrinsic #nu_{#mu}",
               "Beam Intrinsic NC",
               "Dirt",
               "Data",
               "Cosmic contaminated"]

efficiencies = []
eff_num = [0 for i in range(len(description))]
eff_den = [0 for i in range(len(description))]

for i, d in enumerate(description):
    efficiency = ROOT.TEfficiency("e_%i" % i, ";#nu_{e} energy [GeV];Efficiency", 20, 0, 2)
    efficiencies.append(efficiency)

PROTON_THRESHOLD = 0.040
ELECTRON_THRESHOLD = 0.020

passed = 0
not_passed = 0
flash_passed = 0
reco_ok = 0
vertex_ok = 0
noflash = 0
entries = chain_nue.GetEntries()
cc_events = 0
nue_pot = 0
for i in range(chain_pot.GetEntries()):
    chain_pot.GetEntry(i)
    nue_pot += chain_pot.pot

print("POT", nue_pot, pot, pot/nue_pot)

for i in range(100000):
    chain_nue.GetEntry(i)
    if i % 1000 == 0: print(i)
    protons = 0
    electrons = 0
    photons = 0
    pions = 0
    electron_energy = 0
    proton_energy = 0
    for i, energy in enumerate(chain_nue.nu_daughters_E):
        if chain_nue.nu_daughters_pdg[i] == 2212:
            proton_energy += energy - 0.938
            if energy - 0.938 > PROTON_THRESHOLD:
                protons += 1

        if chain_nue.nu_daughters_pdg[i] == 11:
            electron_energy += energy
            if energy - 0.51e-3 > ELECTRON_THRESHOLD:
                electrons += 1

        if chain_nue.nu_daughters_pdg[i] == 22:
            # if energy > 0.035:
            photons += 1

        if chain_nue.nu_daughters_pdg[i] == 111:
            # if energy > 0.06:
            pions += 1

    eNp = electrons == 1 and photons == 0 and pions == 0 and protons >= 1

    true_neutrino_vertex = [chain_nue.true_vx_sce,
                            chain_nue.true_vy_sce,
                            chain_nue.true_vz_sce]

    if chain_nue.category != 4:
        cc_events += 1

    # if not eNp: print(electrons, photons, protons, pions)

    if 0 < chain_nue.nu_E < 100:
        eNp_events += 1
        h_x.Fill(true_neutrino_vertex[0])
        h_y.Fill(true_neutrino_vertex[1])
        h_z.Fill(true_neutrino_vertex[2])

        if is_fiducial(true_neutrino_vertex):

            p = False
            p_track = False
            p_shower = False

            if protons == chain_nue.nu_matched_tracks:
                p_track = True

            if electrons == chain_nue.nu_matched_showers:
                p_shower = True

            if p_track and p_shower:
                p = True

            proton_energy = sum([chain_nue.nu_daughters_E[i] - 0.938
                                 for i, pdg in
                                 enumerate(chain_nue.nu_daughters_pdg)
                                 if pdg == 2212])

            if electrons > 0:
                electron_energy = max([chain_nue.nu_daughters_E[i] - 0.51e-3
                                        for i, pdg in
                                        enumerate(chain_nue.nu_daughters_pdg)
                                        if abs(pdg) == 11])

            fp_energy = chain_nue.nu_E #electron_energy + proton_energy

            neutrino_vertex = [chain_nue.vx, chain_nue.vy, chain_nue.vz]

            true_neutrino_vertex = [chain_nue.true_vx_sce,
                                    chain_nue.true_vy_sce,
                                    chain_nue.true_vz_sce]
            true_neutrino_vertex_nosce = [chain_nue.true_vx,
                                          chain_nue.true_vy,
                                          chain_nue.true_vz]

            dist = 0

            electron_energies = [chain_nue.nu_daughters_E[i] - 0.51e-3
                                 for i, pdg in
                                 enumerate(chain_nue.nu_daughters_pdg)
                                 if abs(pdg) == 11]

            proton_energies = [chain_nue.nu_daughters_E[i] - 0.938
                               for i, pdg in
                               enumerate(chain_nue.nu_daughters_pdg)
                               if pdg == 2212]

            if eNp and chain_nue.category != 4:
                ep_energy.Fill(chain_nue.passed and chain_nue.category == 2 and p_track and p_shower, fp_energy)
                ep_dist_energy.Fill(chain_nue.passed and chain_nue.distance < 5 and chain_nue.category == 2, fp_energy)
                e_energy.Fill(chain_nue.passed and chain_nue.category == 2, fp_energy)

                if chain_nue.passed and chain_nue.category == 2:
                    p_energy.Fill(p_track and p_shower, fp_energy)
                    p_dist_energy.Fill(chain_nue.distance < 5, fp_energy)

                h_total_2d.Fill(proton_energies[0], electron_energies[0])
                fiducial += 1

            if chain_nue.passed and chain_nue.category == 2:
                if eNp:
                    matched_electron_energies = [e for e, pdg in zip(chain_nue.matched_showers_energy, chain_nue.matched_showers) if abs(pdg) == 11]
                    if matched_electron_energies:# and electron_energies[0] > 0.02:
                        h_passed_electron.Fill(electron_energies[0])

                    matched_protons = len([pdg for pdg in chain_nue.matched_tracks if pdg == 2212])
                    if matched_protons:
                        h_passed_proton.Fill(max(proton_energies))
                    
                    if matched_electron_energies and matched_protons:
                        h_passed_2d.Fill(proton_energies[0], electron_energies[0])

                    h_total_electron.Fill(electron_energies[0])
                    for proton in proton_energies:
                        h_total_proton.Fill(proton)

                    for i_sh in range(chain_nue.n_showers):
                        energy = chain_nue.matched_showers_energy[i_sh]

                    passed += 1
                    if p:
                        reco_ok += 1
                    if chain_nue.distance < 5:
                        vertex_ok += 1
                    dist = math.sqrt(sum([(t - r) ** 2 for t, r in zip(neutrino_vertex, true_neutrino_vertex)]))

                    h_dist.Fill(dist)
                    h_dist_nosce.Fill(math.sqrt(
                        sum([(t - r)**2 for t, r in zip(
                            neutrino_vertex, true_neutrino_vertex_nosce)])
                    ))

                    h_x_diff.Fill(neutrino_vertex[0] - true_neutrino_vertex[0])
                    h_y_diff.Fill(neutrino_vertex[1] - true_neutrino_vertex[1])
                    h_z_diff.Fill(neutrino_vertex[2] - true_neutrino_vertex[2])
            else:
                if eNp:
                    if chain_nue.flash_passed != 1:
                        noflash += 1

            if protons == 1:
                e_proton.Fill(chain_nue.passed and p_track and p_shower,
                              proton_energy)

            if chain_nue.passed and protons == 1 and p_track:
                l_e_proton.Fill(chain_nue.track_len[0], proton_energy)

            if chain_nue.passed:
                e_nprotons.Fill(p_track, protons)

            for i in range(len(description)):
                if chain_nue.category == 4:
                    if chain_nue.passed:
                        eff_num[4] += 1
                    eff_den[4] += 1
                    efficiencies[4].Fill(chain_nue.passed, fp_energy)
                if eNp:
                    if chain_nue.passed and chain_nue.category == 2:
                        eff_num[2] += 1
                    eff_den[2] += 1
                    efficiencies[2].Fill(chain_nue.passed and chain_nue.category == 2, fp_energy)
                if chain_nue.category == 3:
                    if chain_nue.passed:
                        eff_num[3] += 1
                    eff_den[3] += 1
                    efficiencies[3].Fill(chain_nue.passed, fp_energy)



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
ep_vertex = p_vertex * eff
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
#pt = ROOT.TPaveText(0.098, 0.905, 0.576, 0.989, "ndc")

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
# legend.AddEntry(p_energy, "P_{{reco}} ({0:.1f} #pm {1:.1f}) %"
#                 .format(p_reco * 100, p_reco_err * 100),
#                 "lep")

legend.AddEntry(p_dist_energy, "P ({0:.1f} #pm {1:.1f}) %"
                .format(p_vertex * 100, p_vertex_err * 100), "lep")
# legend.AddEntry(ep_energy,"#epsilon #times P_{{reco}} ({0:.1f} #pm {1:.1f})%"
#                 .format(ep_reco * 100, ep_reco_err * 100), "lep")

legend.AddEntry(
    ep_dist_energy, "#epsilon #times P ({0:.1f} #pm {1:.1f}) %"
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
# p_energy.Draw("pl same")

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
# ep_energy.Draw("pl same")

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
e_energy.GetPaintedGraph().GetXaxis().SetRangeUser(0, 2)
e_energy.GetPaintedGraph().GetYaxis().SetRangeUser(0, 1.39)
e_energy.GetPaintedGraph().GetYaxis().SetTitleOffset(0.9)
pt.Draw()
legend.Draw()
c_energy.SaveAs("plots/efficiency.pdf")
c_energy.Draw()

c_proton = ROOT.TCanvas("c_proton")
e_proton.Draw("apl")
c_proton.Update()

c_nproton = ROOT.TCanvas("c_nproton")
e_nprotons.Draw("ab")
c_nproton.Update()

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

c_spectrum = ROOT.TCanvas("c_spectrum")
h_energy.Scale(pot/nue_pot)
print("nu_e events", h_energy.Integral())
fix_binning(h_energy)
h_energy.Draw("hist")
h_energy.SetLineColor(1)
h_energy.SetFillColor(ROOT.kGreen - 2)
pt.Draw()
c_spectrum.Update()


c_electron = ROOT.TCanvas("c_electron")
e_electron = ROOT.TEfficiency(h_passed_electron, h_total_electron)
e_electron.Draw("ap")
e_electron.SetMarkerStyle(20)
e_electron.SetMarkerColor(ROOT.kAzure + 1)
c_electron.Update()
e_electron.GetPaintedGraph().GetYaxis().SetTitleOffset(0.9)
c_electron.Update()

c_proton2 = ROOT.TCanvas("c_proton2")
e_proton = ROOT.TEfficiency(h_passed_proton, h_total_proton)
e_proton.Draw("ap")
e_proton.SetMarkerStyle(20)
e_proton.SetMarkerColor(ROOT.kRed + 1)
c_proton2.Update()
e_proton.GetPaintedGraph().GetYaxis().SetTitleOffset(0.9)
c_proton2.Update()

c_2d = ROOT.TCanvas("c_2d")
e_2d = ROOT.TEfficiency(h_passed_2d, h_total_2d)
e_2d.Draw("colz")
c_2d.Update()
e_2d.GetPaintedHistogram().GetYaxis().SetTitleOffset(0.9)
c_2d.Update()


c_efficiencies = ROOT.TCanvas("c_efficiencies")
efficiencies[1].SetMarkerStyle(20)
efficiencies[1].SetMarkerColor(colors[0])
# efficiencies[1].Draw("ap")
l = ROOT.TLegend(0.1074499, 0.8336842, 0.8925501, 0.9621053)
if eff_den[1] > 0:
    eff = eff_num[1]/eff_den[1] * 100
    l.AddEntry(efficiencies[1], "%s: %.1f%%" % (description[1], eff), "lep")
for i in range(2, len(efficiencies)):
    if description[i] != "Data":
        efficiencies[i].SetMarkerStyle(20)
        if description[i] == "Cosmic contaminated":
            efficiencies[i].SetMarkerStyle(24)
        efficiencies[i].SetMarkerColor(colors[i])
        if eff_den[i] > 0:
            if i == 2:
                efficiencies[i].Draw("ap")
                # efficiencies[i].GetPaintedGraph().GetXaxis().SetRangeUser(0, 2)
                # efficiencies[i].GetPaintedGraph().GetYaxis().SetRangeUser(0, 1.39)
            else:
                efficiencies[i].Draw("p same")
            eff = eff_num[i]/eff_den[i] * 100
            l.AddEntry(efficiencies[i], "%s: %.1f%%" % (description[i], eff), "lep")

l.SetNColumns(2)
l.Draw()
c_efficiencies.SetLeftMargin(0.1346705)
c_efficiencies.SetTopMargin(0.1936842)
c_efficiencies.Range(-0.3872659,-0.02854239,2.48839,0.2548979)
c_efficiencies.Update()
c_efficiencies.SaveAs("plots/efficiencies.pdf")
c_efficiencies.Draw()

input()
