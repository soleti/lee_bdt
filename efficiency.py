#!/usr/local/bin/python3

import math
import ROOT

from glob import glob
from bdt_common import x_start, x_end, y_start, y_end, z_start, z_end, printProgressBar, bins, distance, is_fiducial
from array import array
from proton_energy import length2energy

ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptFit(0)
ROOT.gStyle.SetPalette(87)
ROOT.gStyle.SetNumberContours(99)

scaling = 1
pot_weight = 0.0041654013956019835 * scaling

colors = [ROOT.TColor.GetColor("#e7623d"), ROOT.TColor.GetColor("#c46d27"),
          ROOT.TColor.GetColor("#62b570"), ROOT.TColor.GetColor("#8779b1"),
          ROOT.TColor.GetColor("#609ec6"), ROOT.TColor.GetColor("#c46d27"),
          ROOT.TColor.GetColor("#ffffff"), ROOT.TColor.GetColor("#e7623d")]

pot = 6.6e20


def choose_plane(root_chain):
    total_hits = [0, 0, 0]
    shower_hits = [0, 0, 0]
    track_hits = [0, 0, 0]

    for i_sh in range(root_chain.n_showers):
        for i_plane in range(len(root_chain.shower_nhits[i_sh])):
            total_hits[i_plane] += root_chain.shower_nhits[i_sh][i_plane]
            shower_hits[i_plane] += root_chain.shower_nhits[i_sh][i_plane]

    for i_tr in range(root_chain.n_tracks):
        for i_plane in range(len(root_chain.track_nhits[i_tr])):
            total_hits[i_plane] += root_chain.track_nhits[i_tr][i_plane]
            track_hits[i_plane] += root_chain.track_nhits[i_tr][i_plane]

    product = [t * s for t, s in zip(track_hits, shower_hits)]

    return product.index(max(product))


def fix_binning(histo):
    width = histo.GetBinWidth(1)
    for k in range(1, histo.GetNbinsX() + 1):
        bin_width = histo.GetBinWidth(k)
        histo.SetBinError(k, histo.GetBinError(k) / (bin_width / width))
        histo.SetBinContent(k, histo.GetBinContent(k) / (bin_width / width))

def is_active(point):
    ok_y = point[1] > y_start and point[1] < y_end
    ok_x = point[0] > x_start and point[0] < x_end
    ok_z = point[2] > z_start and point[2] < z_end
    return ok_y and ok_x and ok_z


# glob("mc_nue_ubxsec/*.root")
nue_cosmic = glob("mc_nue_crhit/*.root")  # glob("mc_bnb_ubxsec3/*/*.root")
chain_nue = ROOT.TChain("robertoana/pandoratree")

for f in nue_cosmic:
    chain_nue.Add(f)

print("entries", chain_nue.GetEntries())
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
                            20, 0, 0.2)
p_dist_energy = ROOT.TEfficiency("p_dist_energy",
                                 ";#nu_{e} energy [GeV];Purity",
                                 20, 0, 2)
ep_dist_energy = ROOT.TEfficiency("p_dist_energy",
                                  ";#nu_{e} energy [GeV];Efficiency #times \
                                  Purity",
                                  20, 0, 2)


e_binned = []
for i in range(30):
    e_binned.append(ROOT.TEfficiency("e_energy%i" %
                                     i, ";#nu_{e} energy [GeV];Efficiency", 20, 0, 2))


h_dist = ROOT.TH1F("h_dist", ";Distance [cm];N. Entries / 0.3 cm", 20, 0, 6)
h_dist_nosce = ROOT.TH1F("h_dist_nosce",
                         ";Distance [cm];N. Entries / 0.3 cm",
                         20, 0, 6)

h_x_diff = ROOT.TH1F("h_x_diff",
                     ";#Delta x [cm]; N. Entries / 0.2 cm", 50, -5, 5)
h_y_diff = ROOT.TH1F("h_y_diff",
                     ";#Delta y [cm]; N. Entries / 0.2 cm", 50, -5, 5)
h_z_diff = ROOT.TH1F("h_z_diff",
                     ";#Delta z [cm]; N. Entries / 0.2 cm", 50, -5, 5)

h_energy_binned = ROOT.TH1F(
    "h_energy_binned", ";E_{#nu_{e}}; N. Entries / 0.05 GeV", len(bins) - 1, bins)
h_energy_binned_passed = ROOT.TH1F(
    "h_energy_binned_passed", ";E_{#nu_{e}}; N. Entries / 0.05 GeV", len(bins) - 1, bins)

h_energy = ROOT.TH1F(
    "h_energy", ";Reco. energy [GeV]; N. Entries / 0.05 GeV", 20, 0, 10)
h_energy_passed = ROOT.TH1F(
    "h_energy_passed", ";Reco. energy [GeV]; N. Entries / 0.05 GeV", 20, 0, 10)

eNp_events = 0

h_total_electron = ROOT.TH1F(
    "h_total_electron", ";E_{e} [GeV];Efficiency", 10, 0, 0.2)
h_passed_electron = ROOT.TH1F(
    "h_passed_electron", ";E_{e} [GeV];Efficiency", 10, 0, 0.2)
h_total_proton = ROOT.TH1F(
    "h_total_proton", ";E_{p} [GeV];Efficiency", 20, 0, 2)
h_passed_proton = ROOT.TH1F(
    "h_passed_proton", ";E_{p} [GeV];Efficiency", 20, 0, 2)
h_total_2d = ROOT.TH2F(
    "h_total_2d", ";E_{p} [GeV];E_{e} [GeV]", 40, 0, 1, 40, 0, 1)
h_passed_2d = ROOT.TH2F(
    "h_passed_2d", ";E_{p} [GeV];E_{e} [GeV]", 40, 0, 1, 40, 0, 1)

h_missing_proton_energy = ROOT.TH2F(
    "h_missing_proton_energy", ";Proton E_{k}^{true} [GeV];Proton E_{k}^{reco} [GeV]", 50, 0, 1, 50, 0, 1)
l_missing_proton_energy = [
    ROOT.TH1F("h_%i_p" % i, "", 50, 0, 1) for i in range(50)]
h_electron_energy = ROOT.TH2F(
    "h_electron_energy", ";e^{-} E_{k}^{true} [Gev];e^{-} E_{k}^{reco} [GeV]", 50, 0, 1, 50, 0, 1)
l_electron_energy = [
    ROOT.TH1F("h_%i_e" % i, "", 50, 0, 1) for i in range(50)]
p_proton_energy = ROOT.TProfile(
    "p_proton_energy", ";Proton E_{k} [GeV];E_{reco}/E_{k}", 50, 0, 1, 0, 1)

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
    efficiency = ROOT.TEfficiency(
        "e_%i" % i, ";#nu_{e} energy [GeV];Efficiency", 20, 0, 1)
    efficiencies.append(efficiency)

PROTON_THRESHOLD = 0.020
ELECTRON_THRESHOLD = 0.020

passed = 0
not_passed = 0
flash_passed = 0
reco_ok = 0
vertex_ok = 0
noflash = 0
entries = int(chain_nue.GetEntries() / scaling)
cc_events = 0
fiducial = 0
dirt = 0
intrinsic_numu = 0
intrinsic_NC = 0
nue_cc0pi_np = 0
nue_cc = 0

bnb_weight = 0.02065927860647806
nue_weight = 0.001830225149382802
proton1 = 0
proton2 = 0

protonN = 0
for i in range(entries):
    printProgressBar(i, entries, prefix="Progress:",
                     suffix="Complete", length=20)
    chain_nue.GetEntry(i)
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
    
            if energy - 0.938 > PROTON_THRESHOLD and is_fiducial(p_vertex) and is_fiducial(p_end):
                protons += 1

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
    if chain_nue.category != 4:
        cc_events += chain_nue.bnbweight

    # if not eNp: print(electrons, photons, protons, pions)

    if 0 < chain_nue.nu_E < 100:
        eNp_events += chain_nue.bnbweight

        if not is_fiducial(true_neutrino_vertex_nosce) and chain_nue.category != 0 and chain_nue.category != 6 and chain_nue.category != 1 and chain_nue.category != 7:
            dirt += bnb_weight
        elif chain_nue.nu_pdg == 14 and chain_nue.ccnc == 0 and is_fiducial(true_neutrino_vertex_nosce):
            intrinsic_numu += bnb_weight
        elif chain_nue.ccnc == 1 and is_fiducial(true_neutrino_vertex_nosce):
            intrinsic_NC += bnb_weight
        elif chain_nue.nu_pdg == 12 and eNp and is_fiducial(true_neutrino_vertex_nosce):
            nue_cc0pi_np += nue_weight * chain_nue.bnbweight
        elif chain_nue.nu_pdg == 12 and not eNp and is_fiducial(true_neutrino_vertex_nosce):
            nue_cc += nue_weight * chain_nue.bnbweight

        if is_fiducial(true_neutrino_vertex):

            p_track = False
            p_shower = False

            if protons == chain_nue.nu_matched_tracks:
                p_track = True

            if electrons == chain_nue.nu_matched_showers:
                p_shower = True

            fp_energy = chain_nue.nu_E

            neutrino_vertex = [chain_nue.vx, chain_nue.vy, chain_nue.vz]

            electron_energies = [chain_nue.nu_daughters_E[i] - 0.51e-3
                                 for i, pdg in
                                 enumerate(chain_nue.nu_daughters_pdg)
                                 if abs(pdg) == 11]

            proton_energies = [chain_nue.nu_daughters_E[i] - 0.938
                               for i, pdg in
                               enumerate(chain_nue.nu_daughters_pdg)
                               if pdg == 2212]

            if eNp and chain_nue.ccnc == 0:
                showers_check = True
                if chain_nue.n_tracks == 0 and chain_nue.n_showers == 1:
                    showers_check = False

                contaminated = chain_nue.cosmic_fraction > 0.5 and chain_nue.category == 7
                selected = chain_nue.category == 2 or contaminated

                efficiency_condition = showers_check and chain_nue.passed and selected
                matched_proton_showers = sum(
                    1 for pdg in chain_nue.matched_showers if pdg == 2212)
                matched_proton_tracks = sum(
                    1 for pdg in chain_nue.matched_tracks if pdg == 2212)
                proton_condition = matched_proton_showers + matched_proton_tracks >= protons
                matched_electron_showers = sum(
                    1 for pdg in chain_nue.matched_showers if pdg == 11)
                electron_condition = matched_electron_showers >= electrons
                reco_condition = electron_condition and proton_condition
                ep_energy.Fill(
                    efficiency_condition and reco_condition and chain_nue.distance < 5, fp_energy)
                ep_dist_energy.Fill(
                    efficiency_condition and chain_nue.distance < 5, fp_energy)
                e_energy.Fill(efficiency_condition, fp_energy)
                e_binned[min(int(chain_nue.bnbweight * 10), 29)
                         ].Fill(efficiency_condition, fp_energy)

                h_energy.Fill(chain_nue.nu_E, chain_nue.bnbweight)
                h_energy_binned.Fill(chain_nue.nu_E, chain_nue.bnbweight)
                h_total_2d.Fill(sum(proton_energies), sum(electron_energies))
                fiducial += chain_nue.bnbweight
                h_total_electron.Fill(sum(electron_energies))
                h_total_proton.Fill(sum(proton_energies))
                if efficiency_condition:
                    if protons == 1:
                        proton1 += 1
                    elif protons == 2:
                        proton2 += 1

                    elif protons > 2:
                        protonN += 1
                    h_energy_binned_passed.Fill(
                        chain_nue.nu_E, chain_nue.bnbweight)
                    passed += chain_nue.bnbweight

                    if chain_nue.distance < 5:
                        if reco_condition:
                            reco_ok += chain_nue.bnbweight
                        vertex_ok += chain_nue.bnbweight

                    hit_index = choose_plane(chain_nue)
                    reco_proton_energies = 0
                    reco_electron_energy = 0

                    for i_sh in range(chain_nue.n_showers):
                        if chain_nue.matched_showers[i_sh] == 2212:
                            reco_proton_energies += chain_nue.shower_energy[i_sh][hit_index]
                        if chain_nue.matched_showers[i_sh] == 11:
                            reco_electron_energy += chain_nue.shower_energy[i_sh][hit_index]

                    for i_tr in range(chain_nue.n_tracks):
                        if chain_nue.matched_tracks[i_tr] == 2212:
                            reco_proton_energies += length2energy(
                                chain_nue.track_len[i_tr])
                        if chain_nue.matched_tracks[i_tr] == 11:
                            reco_electron_energy += chain_nue.track_energy_hits[i_tr][hit_index]

                    if 11 in chain_nue.matched_tracks or 11 in chain_nue.matched_showers:
                        h_passed_electron.Fill(sum(electron_energies))

                    p_proton_energy.Fill(
                        sum(proton_energies), reco_proton_energies / sum(proton_energies))

                    if 2212 in chain_nue.matched_tracks or 2212 in chain_nue.matched_showers:
                        h_passed_proton.Fill(sum(proton_energies))

                        h_missing_proton_energy.Fill(
                            sum(proton_energies), reco_proton_energies)
                        if sum(proton_energies) < 1:
                            l_missing_proton_energy[int(
                                sum(proton_energies) / 0.02)].Fill(reco_proton_energies)

                    if 11 in chain_nue.matched_tracks or 11 in chain_nue.matched_showers:
                        h_electron_energy.Fill(
                            sum(electron_energies), reco_electron_energy)
                        if sum(electron_energies) < 1:
                            l_electron_energy[int(
                                sum(electron_energies) / 0.02)].Fill(reco_electron_energy)

                    h_energy_passed.Fill(chain_nue.nu_E, chain_nue.bnbweight)

                    p_energy.Fill(p_track and p_shower, fp_energy)
                    p_dist_energy.Fill(chain_nue.distance < 5, fp_energy)

                    electron_showers = 0
                    for i_sh in range(chain_nue.n_showers):
                        if chain_nue.matched_showers[i_sh] == 11:
                            electron_showers += 1

                    if 2212 in chain_nue.matched_tracks and 11 in chain_nue.matched_showers:
                        h_passed_2d.Fill(sum(proton_energies),
                                         sum(electron_energies))

                    dist = math.sqrt(
                        sum([(t - r) ** 2 for t, r in zip(neutrino_vertex, true_neutrino_vertex)]))
                    dist_nosce = math.sqrt(
                        sum([(t - r)**2 for t, r in zip(neutrino_vertex, true_neutrino_vertex_nosce)]))

                    h_dist.Fill(dist)
                    h_dist_nosce.Fill(dist_nosce)

                    h_x_diff.Fill(neutrino_vertex[0] - true_neutrino_vertex[0])
                    h_y_diff.Fill(neutrino_vertex[1] - true_neutrino_vertex[1])
                    h_z_diff.Fill(neutrino_vertex[2] - true_neutrino_vertex[2])

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
                    efficiencies[2].Fill(
                        chain_nue.passed and chain_nue.category == 2, fp_energy)
                if chain_nue.category == 3:
                    if chain_nue.passed:
                        eff_num[3] += 1
                    eff_den[3] += 1
                    efficiencies[3].Fill(chain_nue.passed, fp_energy)

print(dirt,
      intrinsic_numu * scaling,
      intrinsic_NC * scaling,
      nue_cc0pi_np * scaling,
      nue_cc * scaling)

print("p", proton1, proton2, protonN)
print("Entries", entries)
print("CC", cc_events)
print("1eNp", eNp_events)
print("1eNp + Is fiducial", fiducial)
print("Passed", passed)

eff = passed / fiducial
eff_err = math.sqrt((eff * (1 - eff)) / passed)

p_reco = reco_ok / fiducial
p_reco_err = math.sqrt((p_reco * (1 - p_reco)) / passed)
p_vertex = vertex_ok / passed
p_vertex_err = math.sqrt((p_vertex * (1 - p_vertex)) /
                         passed)

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


pt = ROOT.TPaveText(0.1, 0.91, 0.45, 0.97, "ndc")
pt.AddText("MicroBooNE Preliminary")
pt.SetFillColor(0)
pt.SetBorderSize(0)
pt.SetShadowColor(0)

c_dist = ROOT.TCanvas("c_dist")
h_dist.SetLineWidth(3)
h_dist.SetLineColor(ROOT.kAzure + 1)
h_dist.Draw("hist")
c_dist.Update()

c_electron = ROOT.TCanvas("c_electron")
e_electron = ROOT.TEfficiency(h_passed_electron, h_total_electron)
e_electron.Draw("ap")
e_electron.SetMarkerStyle(20)
e_electron.SetMarkerColor(ROOT.kAzure + 1)
c_electron.Update()
e_electron.GetPaintedGraph().GetYaxis().SetTitleOffset(0.9)
c_electron.Update()

c_proton2 = ROOT.TCanvas("c_proton2")
h_passed_proton.SetBinContent(3, h_passed_proton.GetBinContent(3) / 10)
h_passed_proton.SetBinContent(4, h_passed_proton.GetBinContent(4) / 10)

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
    eff = eff_num[1] / eff_den[1] * 100
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
            eff = eff_num[i] / eff_den[i] * 100
            l.AddEntry(efficiencies[i], "%s: %.1f%%" %
                       (description[i], eff), "lep")

l.SetNColumns(2)
l.Draw()
c_efficiencies.SetLeftMargin(0.1346705)
c_efficiencies.SetTopMargin(0.1936842)
c_efficiencies.Range(-0.3872659, -0.02854239, 2.48839, 0.2548979)
c_efficiencies.Update()
c_efficiencies.SaveAs("plots/efficiencies.pdf")
c_efficiencies.Draw()

c_binned = ROOT.TCanvas("c_binned")
eff_total = e_binned[0]
for e in e_binned[1:]:
    eff_total.Add(e)
eff_total.Draw("AC3")
eff_total.SetFillColor(ROOT.kRed + 1)
eff_total.SetFillStyle(3002)
eff_total.SetLineColor(ROOT.kRed + 1)
eff_total.SetLineWidth(2)
f_energy = ROOT.TFile("plots/f_energy.root", "RECREATE")
eff_total.Write()
f_energy.Close()

ep_dist_energy.SetMarkerSize(0)
ep_dist_energy.SetFillColor(ROOT.kGreen + 2)
ep_dist_energy.SetLineColor(ROOT.kGreen + 2)
ep_dist_energy.SetLineWidth(2)
ep_dist_energy.SetFillStyle(3002)
ep_dist_energy.Draw("C3 same")

ep_energy.SetLineColor(ROOT.kBlue + 1)
ep_energy.SetLineWidth(2)
ep_energy.SetFillStyle(3002)
ep_energy.SetFillColor(ROOT.kBlue + 1)
ep_energy.Draw("C3 same")


legend = ROOT.TLegend(0.13, 0.68, 0.84, 0.86)
legend.SetTextSize(16)
legend.AddEntry(eff_total, "#epsilon ({0:.1f} #pm {1:.1f}) %"
                .format(h_energy_passed.Integral() / h_energy.Integral() * 100, eff_err * 100),
                "fl")

legend.AddEntry(ep_energy, "#epsilon_{{reco}} ({0:.1f} #pm {1:.1f}) %"
                .format(ep_reco * 100, ep_reco_err * 100),
                "fl")

legend.AddEntry(ep_dist_energy, "#epsilon_{{vertex}} ({0:.1f} #pm {1:.1f}) %"
                .format(ep_vertex * 100, ep_vertex_err * 100), "fl")

legend.SetFillStyle(0)
legend.SetNColumns(2)
pt.Draw()
legend.Draw()
c_binned.Update()

c_missing_proton = ROOT.TCanvas("c_missing_proton")
a_missing_proton_energy = array("f", [])
a_errs = array("f", [])
a_bins = array("f", [])
a_bins_errs = array("f", [])

for i, bin in enumerate(l_missing_proton_energy):
    if bin.Integral() > 0:
        a_missing_proton_energy.append(bin.GetMaximumBin() * 0.02 - 0.01)
        a_errs.append(bin.GetStdDev())
    else:
        a_missing_proton_energy.append(0)
        a_errs.append(0)
    a_bins.append(1 / 50 * i + 0.01)
    a_bins_errs.append(0.02 / 2)

g_missing_proton_energy = ROOT.TGraphErrors(
    50, a_bins, a_missing_proton_energy, a_bins_errs, a_errs)
h_missing_proton_energy.Draw("colz")
g_missing_proton_energy.SetMarkerStyle(20)
g_missing_proton_energy.Draw("p same")
f_line = ROOT.TF1("f_line", "[0]*x+[1]", 0, 1)
f_line.SetParNames("m", "q")
f_line.SetParameters(1, 0)
l_p_true_reco = ROOT.TLegend(0.099, 0.913, 0.900, 0.968)
l_p_true_reco.SetNColumns(2)
l_p_true_reco.AddEntry(g_missing_proton_energy, "Most probable values", "lep")
g_missing_proton_energy.Fit(f_line, "", "", 0.06, 0.5)
l_p_true_reco.AddEntry(f_line, "E_{k}^{reco} = %.2f E_{k}^{true} + %.2f GeV" %
                       (f_line.GetParameter(0), f_line.GetParameter(1)), "l")
l_p_true_reco.Draw()
c_missing_proton.Update()

c_electron_energy = ROOT.TCanvas("c_electron_energy")
a_electron_energy = array("f", [])
a_electron_errs = array("f", [])
for i, bin in enumerate(l_electron_energy):
    if bin.Integral() > 0:
        a_electron_energy.append(bin.GetMaximumBin() * 0.02 - 0.01)
        a_electron_errs.append(bin.GetStdDev())
    else:
        a_electron_energy.append(0)
        a_electron_errs.append(0)

g_electron_energy = ROOT.TGraphErrors(
    50, a_bins, a_electron_energy, a_bins_errs, a_electron_errs)
h_electron_energy.Draw("colz")
g_electron_energy.SetMarkerStyle(20)
g_electron_energy.Draw("p same")
f_line2 = f_line.Clone()
f_line2.SetName("f_line2")
l_e_true_reco = ROOT.TLegend(0.099, 0.913, 0.900, 0.968)
l_e_true_reco.SetNColumns(2)
l_e_true_reco.AddEntry(g_electron_energy, "Most probable values", "lep")
g_electron_energy.Fit(f_line2, "", "", 0.02, 1)
l_e_true_reco.AddEntry(f_line2, "E_{k}^{reco} = %.2f E_{k}^{true} + %.2f GeV" %
                       (f_line2.GetParameter(0), f_line2.GetParameter(1)), "l")
l_e_true_reco.Draw()
c_electron_energy.Update()

c_proton_energy = ROOT.TCanvas("c_proton_energy")
p_proton_energy.Draw()
c_proton_energy.Update()

c_energy_binned = ROOT.TCanvas("c_energy_binned")
h_energy_binned.Scale(pot_weight)
# h_energy_binned.Draw("hist")
h_energy_binned_passed.Scale(pot_weight)

h_energy_binned_passed.Draw("hist same")

f_energy_binned = ROOT.TFile("plots/h_energy_binned.root", "RECREATE")
h_energy_binned.Write()
f_energy_binned.Close()
c_energy_binned.Update()

input()
