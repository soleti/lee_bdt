#!/usr/bin/env python3.4

import ROOT
import statistics
import math
from array import array
from glob import glob


def gauss_exp(var, par):
    """
    n:par[0]
    mu:par[1]
    sigma:par[2]
    k:par[3]
    """
    n = par[0]
    mu = par[1]
    sigma = par[2]
    k = par[3]
    x = var[0]

    if (x - mu) / sigma >= -k:
        return n*math.exp(-0.5 * ((x - mu) / sigma)**2)
    else:
        return n*math.exp(k**2 / 2 + k * ((x - mu) / sigma))


ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetNumberContours(99)
ROOT.gStyle.SetOptFit(0)

nue_cosmic = sorted(glob("nue_efficiency/*/Pandora*.root"))
chain = ROOT.TChain("robertoana/pandoratree")

for f in nue_cosmic:
    chain.Add(f)

entries = chain.GetEntries()

total = 0
splitted_event = 0
incomplete_event = 0
perfect_event = 0
wrong_event = 0
track_ok_shower_mis = 0
track_ok_shower_no = 0
shower_ok_track_mis = 0
shower_ok_track_no = 0
track_no_shower_no = 0
flash_not_passed = 0

h_e_true_reco = ROOT.TH2F("h_e_true_reco",
                          ";E_{k}^{true} [GeV];E_{k}^{reco} [GeV]",
                          100, 0, 2, 100, 0, 2)

h_e_res = ROOT.TH1F("h_res",
                    ";(E_{k}^{corr} - E_{k}^{true}) / E_{k}^{true};N. Entries / 0.04",
                    50, -1, 1)

l_true_reco = [[], [], [], [], [], [], [], [], [], []]

for evt in range(entries):
    chain.GetEntry(evt)

    p = 0
    e = 0
    photons = 0
    pions = 0
    electron_energy = 0
    proton_energy = 0
    for i, energy in enumerate(chain.nu_daughters_E):
        if abs(chain.nu_daughters_pdg[i]) == 2212:
            if energy - 0.938 > 0.000005:
                proton_energy += energy - 0.938
                p += 1

        if abs(chain.nu_daughters_pdg[i]) == 11:
            if energy > 0.00003:
                electron_energy += energy
                e += 1

        if chain.nu_daughters_pdg[i] == 22:
            # if energy > 0.035:
            photons += 1

        if chain.nu_daughters_pdg[i] == 111:
            # if energy > 0.06:
            pions += 1

    eNp = e == 1 and photons == 0 and pions == 0 and p > 0

    if eNp and chain.nu_E > 0.1:

        if chain.true_nu_is_fiducial:
            total += 1
            primary_indexes = []
            shower_passed = []
            track_passed = []
            flash_passed = []

            for i in range(chain.n_primaries):
                primary_indexes.append(chain.primary_indexes[i])
                shower_passed.append(chain.shower_passed[i])
                track_passed.append(chain.track_passed[i])
                flash_passed.append(chain.flash_passed[i] + 1)

            if chain.passed:
                candidate_id = primary_indexes.index(chain.chosen_candidate)

                chosen_showers = shower_passed[candidate_id]
                chosen_tracks = track_passed[candidate_id]

                if chosen_showers == e and chosen_tracks == p:
                    perfect_event += 1
                    tot_energy = electron_energy + proton_energy
                    h_e_true_reco.Fill(tot_energy, chain.E)
                    h_e_res.Fill(((chain.E - 2.19184e-02) / 6.35247e-01 -
                                  tot_energy) / tot_energy)
                    if tot_energy < 2:
                        l_true_reco[int(tot_energy / 0.2)].append(chain.E)

                elif chosen_showers > 0 or chosen_tracks > 0:
                    if chosen_showers < e or chosen_tracks < p:
                        incomplete_event += 1
                    elif chosen_showers >= e or chosen_tracks >= p:
                        splitted_event += 1
                else:
                    wrong_event += 1
            else:
                find_track = False
                find_shower = False
                if 1 in flash_passed:
                    for i in range(chain.n_primaries):
                        if track_passed[i] > 0:
                            find_track = True
                            if track_passed[i] > 1:
                                track_ok_shower_mis += 1
                                break
                            else:
                                track_ok_shower_no += 1
                                break
                        if shower_passed[i] > 0:
                            find_shower = True
                            if shower_passed[i] > 1:
                                shower_ok_track_mis += 1
                                break
                            else:
                                shower_ok_track_no += 1
                                break

                    if not find_track and not find_shower:
                        track_no_shower_no += 1

                else:
                    flash_not_passed += 1


print("Passed event, perfect {:.1f} %".format(perfect_event / total * 100))
print("Passed event, incomplete {:.1f} %"
      .format(incomplete_event / total * 100))
print("Passed event, splitted {:.1f} %".format(splitted_event / total * 100))
print("Passed event, wrong {:.1f} %".format(wrong_event / total * 100))
print("Not passed event, flash not passed {:.1f}% "
      .format(flash_not_passed / total * 100))
print("Not passed event, flash passed, track ok, shower misid. {:.1f}% "
      .format(track_ok_shower_mis / total * 100))
print("Not passed event, flash passed, shower ok, track misid. {:.1f}% "
      .format(shower_ok_track_mis / total * 100))
print("Not passed event, flash passed, track ok, no shower {:.1f}% "
      .format(track_ok_shower_no / total * 100))
print("Not passed event, flash passed, shower ok, no track {:.1f}% "
      .format(shower_ok_track_no / total * 100))
print("Not passed event, flash passed, no track, no shower {:.1f}% "
      .format(track_no_shower_no / total * 100))


e_values = array("f", [i * 0.2 + 0.1 for i in range(10)])
e_errs = array("f", [0.1] * 10)

median_values = array("f")
median_errs = array("f")

for i in l_true_reco:
    median_values.append(statistics.median(i))
    median_errs.append(statistics.stdev(i) / math.sqrt(len(i)))

g_e_true_reco = ROOT.TGraphErrors(len(median_values),
                                  e_values, median_values, e_errs, median_errs)

c_e_2d = ROOT.TCanvas("c_e_2d")
h_e_true_reco.Draw("colz")
g_e_true_reco.Draw("ep same")
g_e_true_reco.SetLineWidth(2)
g_e_true_reco.SetMarkerStyle(20)
f_line = ROOT.TF1("f_line", "[0]*x+[1]", 0, 2)
f_line.SetParNames("m", "q")
l_e_true_reco = ROOT.TLegend(0.14, 0.71, 0.49, 0.85)
l_e_true_reco.AddEntry(g_e_true_reco, "Median values", "lep")
g_e_true_reco.Fit(f_line)
l_e_true_reco.AddEntry(f_line, "E_{k}^{reco} = %.2f E_{k}^{true} + %.2f GeV" %
                       (f_line.GetParameter(0), f_line.GetParameter(1)), "l")
h_e_true_reco.GetYaxis().SetTitleOffset(0.8)
l_e_true_reco.Draw()

c_e_2d.Update()

c_e_res = ROOT.TCanvas("c_e_res")

f_gausexp = ROOT.TF1("f_gausexp", gauss_exp, -1, 1, 4)
h_e_res.Draw("ep")
h_e_res.SetMarkerStyle(20)

f_gausexp.SetParLimits(2, 0, 2)
f_gausexp.SetParLimits(3, 0, 1)
f_gausexp.SetParLimits(1, -1, 1)
f_gausexp.SetParameters(250, 0, 0.13, 0.18)

h_e_res.Fit(f_gausexp)
h_e_res.SetLineColor(ROOT.kBlack)
h_e_res.GetYaxis().SetTitleOffset(0.9)

l_e_res = ROOT.TLegend(0.13, 0.76, 0.47, 0.85)
l_e_res.AddEntry(f_gausexp, "GaussExp fit", "l")
l_e_res.AddEntry(f_gausexp, "#mu = %.2f, #sigma = %.2f, k = %.2f" %
        (f_gausexp.GetParameter(1),
        f_gausexp.GetParameter(2),
        f_gausexp.GetParameter(3)), "")
l_e_res.Draw()
l_e_true_reco.AddEntry(g_e_true_reco, "aa", "")

c_e_res.Update()
input()
