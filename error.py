#!/usr/local/bin/python3

import ROOT
from operator import mul
import functools
import math
from bdt_common import manual_cuts, bdt, manual, pre_cuts, variables, spectators, labels, binning, N_UNI, bins, fix_binning, fixed_width_histo, fixed_width_histo_2d

ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetNumberContours(99)
SYS_VARIABLES = ["nu_E"]
MODE = "genie"

if MODE == "flux":
    ROOT.gStyle.SetPalette(ROOT.kMint)
else:
    ROOT.gStyle.SetPalette(ROOT.kLake)


chain = ROOT.TChain("mc_tree")
chain.Add("root_files/mc_file.root")

total_entries = int(chain.GetEntries() / 1)
vars = dict(variables + spectators)
print(total_entries)
h_sys = {}
h_2d = {}
h_cv = {}
h_covariance = {}
h_frac = {}
h_corr = {}
for n, b in vars.items():
    h_uni = []

    for u in range(N_UNI):
        if n == "reco_energy":
            h_uni.append(ROOT.TH1F("h_%s_%i" % (n, u), labels[n], len(bins) - 1, bins))
        else:
            h_uni.append(ROOT.TH1F("h_%s_%i" % (n, u), labels[n], binning[n][0], binning[n][1], binning[n][2]))

    h_sys[n] = h_uni
    if n == "reco_energy":
        h_cv[n] = ROOT.TH1F("h_%s_cv" % n, labels[n], len(bins) - 1, bins)
    else:
        h_cv[n] = ROOT.TH1F("h_%s_cv" % n, labels[n], binning[n][0], binning[n][1], binning[n][2])


for ievt in range(total_entries):
    chain.GetEntry(ievt)

    if ievt % 100 == 0:
        print(ievt)

    if not pre_cuts(chain):
        continue

    if manual:
        apply_manual = manual_cuts(chain)
    else:
        apply_manual = True

    if apply_manual:

        for name, var in variables:
            chain.SetBranchAddress(name, var)

        for name, var in spectators:
            chain.SetBranchAddress(name, var)

        for name, var in vars.items():
            if name not in SYS_VARIABLES:
                continue
            for i, v in enumerate(var):
                if v > -999:
                    h_cv[name].Fill(v, chain.event_weight)

                    for u in range(N_UNI):
                        if MODE == "genie":
                            weight = chain.genie_weights[u]
                        elif MODE == "flux":
                            weight = chain.flux_weights[u]
                        if weight > 100:
                            weight = 1

                        h_sys[name][u].Fill(v, chain.event_weight * weight)


for n, b in vars.items():
    if n == "reco_energy":
        h_2d[n] = ROOT.TH2F("h_%s" % n,
                            labels[n],
                            len(bins) - 1,
                            bins,
                            N_UNI * 2, 0, h_sys[n][0].GetMaximum() * 2)
    else:
        h_2d[n] = ROOT.TH2F("h_%s" % n,
                            labels[n],
                            binning[n][0],
                            binning[n][1],
                            binning[n][2],
                            N_UNI, 0, h_sys[n][0].GetMaximum() * 2)

    if n == "reco_energy":
        h_covariance[n] = ROOT.TH2F("h_cov_%s" % n,
                                    labels[n],
                                    len(bins) - 1,
                                    bins,
                                    len(bins) - 1,
                                    bins)
    else:
        h_covariance[n] = ROOT.TH2F("h_cov_%s" % n,
                                    labels[n],
                                    binning[n][0],
                                    binning[n][1],
                                    binning[n][2],
                                    binning[n][0],
                                    binning[n][1],
                                    binning[n][2])
    h_covariance[n].GetYaxis().SetTitle(h_covariance[n].GetXaxis().GetTitle())

    if n == "reco_energy":
        h_frac[n] = ROOT.TH2F("h_frac_%s" % n,
                              labels[n],
                              len(bins) - 1,
                              bins,
                              len(bins) - 1,
                              bins)
    else:
        h_frac[n] = ROOT.TH2F("h_frac_%s" % n,
                            labels[n],
                            binning[n][0],
                            binning[n][1],
                            binning[n][2],
                            binning[n][0],
                            binning[n][1],
                            binning[n][2])
    h_frac[n].GetYaxis().SetTitle(h_frac[n].GetXaxis().GetTitle())

    if n == "reco_energy":
        h_corr[n] = ROOT.TH2F("h_corr_%s" % n,
                              labels[n],
                              len(bins) - 1,
                              bins,
                              len(bins) - 1,
                              bins)
    else:
        h_corr[n] = ROOT.TH2F("h_corr_%s" % n,
                              labels[n],
                              binning[n][0],
                              binning[n][1],
                              binning[n][2],
                              binning[n][0],
                              binning[n][1],
                              binning[n][2])
    h_corr[n].GetYaxis().SetTitle(h_corr[n].GetXaxis().GetTitle())

sys_err = {}

for n in h_sys:
    sys_err[n] = [0] * (h_cv[n].GetNbinsX() + 1)

    if n == "reco_energy":
        fix_binning(h_cv[n])
        h_cv[n] = fixed_width_histo(h_cv[n])

    for u in range(N_UNI):
        if not h_sys[n][u].Integral() > 0:
            continue

        if n == "reco_energy":
            fix_binning(h_sys[n][u])
            h_sys[n][u] = fixed_width_histo(h_sys[n][u])

        for bin in range(1, h_sys[n][u].GetNbinsX() + 1):
            value = h_sys[n][u].GetBinContent(bin)
            center = h_sys[n][u].GetBinCenter(bin)

            diff = (h_cv[n].GetBinContent(bin) - h_sys[n][u].GetBinContent(bin))

            # if n == "reco_energy":
            #     bin_width = h_sys[n][u].GetBinWidth(bin)
            #     diff = (h_cv[n].GetBinContent(bin) / (bin_width / 0.05) -
            #             h_sys[n][u].GetBinContent(bin) / (bin_width / 0.05))
            #     value /= bin_width / 0.05

            h_2d[n].Fill(center, value)
            sys_err[n][bin] += diff**2



for v in SYS_VARIABLES:
    for i in range(1, h_cv[v].GetNbinsX() + 1):
        for j in range(1, h_cv[v].GetNbinsX() + 1):
            e_ij = 0
            for u in range(N_UNI):
                diff_i = - h_cv[v].GetBinContent(i) + h_sys[v][u].GetBinContent(i)
                diff_j = - h_cv[v].GetBinContent(j) + h_sys[v][u].GetBinContent(j)

                e_ij += diff_i * diff_j

            e_ij /= N_UNI
            f_ij = e_ij/(h_cv[v].GetBinContent(i)*h_cv[v].GetBinContent(j))
            h_frac[v].SetBinContent(i, j, f_ij)
            h_covariance[v].SetBinContent(i, j, e_ij)

for v in SYS_VARIABLES:
    for i in range(1, h_cv[v].GetNbinsX() + 1):
        e_ii = h_covariance[v].GetBinContent(i, i)
        for j in range(1, h_cv[v].GetNbinsX() + 1):
            e_jj = h_covariance[v].GetBinContent(j, j)
            e_ij = h_covariance[v].GetBinContent(i, j)
            h_corr[v].SetBinContent(i, j, e_ij / math.sqrt(e_ii*e_jj))

OBJECTS = []

pt = ROOT.TPaveText(0.057, 0.905, 0.422, 0.951, "ndc")
pt.AddText("MicroBooNE Preliminary")
pt.SetFillColor(0)
pt.SetBorderSize(0)
pt.SetShadowColor(0)

OBJECTS.append(pt)

for v in SYS_VARIABLES:

    for bin in range(1, h_cv[v].GetNbinsX() + 1):
        h_cv[v].SetBinError(bin, math.sqrt(sys_err[v][bin] / N_UNI))

    if v == "reco_energy":
        h_covariance[v] = fixed_width_histo_2d(h_covariance[v])
        h_frac[v] = fixed_width_histo_2d(h_frac[v])
        h_corr[v] = fixed_width_histo_2d(h_corr[v])

    c = ROOT.TCanvas("c_%s" % v)
    h_cv[v].Draw("ep")
    h_2d[v].GetZaxis().SetRangeUser(0, 30)
    h_2d[v].Draw("colz same")
    h_cv[v].Draw("ep same")
    h_cv[v].SetMarkerStyle(0)
    h_cv[v].SetMarkerColor(ROOT.kRed + 1)
    h_cv[v].SetLineColor(ROOT.kRed + 1)
    h_cv[v].GetYaxis().SetTitleOffset(0.9)
    h_cv[v].GetYaxis().SetRangeUser(0.001, h_cv[v].GetMaximum() * 1.5)
    pt.Draw()

    f_cv = ROOT.TFile( "plots/sys/h_%s_%s_sys.root" % (v, MODE), "RECREATE")
    h_cv[v].Write()
    f_cv.Close()

    c.SaveAs("plots/sys/h_%s_%s_err.pdf" % (v, MODE))
    c.Update()

    c_cov = ROOT.TCanvas("c_cov_%s" % v)
    ROOT.gStyle.SetPaintTextFormat(".1f")
    h_covariance[v].Draw("colz text")
    h_covariance[v].GetYaxis().SetTitleOffset(0.9)
    pt.Draw()
    c_cov.SaveAs("plots/sys/h_%s_%s_cov.pdf" % (v, MODE))
    c_cov.Update()

    c_frac = ROOT.TCanvas("c_frac_%s" % v)
    ROOT.gStyle.SetPaintTextFormat(".3f")
    h_frac[v].Draw("colz text")
    h_frac[v].GetYaxis().SetTitleOffset(0.9)
    pt.Draw()
    c_frac.SaveAs("plots/sys/h_%s_%s_frac.pdf" % (v, MODE))
    c_frac.Update()

    c_corr = ROOT.TCanvas("c_corr_%s" % v)
    ROOT.gStyle.SetPaintTextFormat(".3f")
    h_corr[v].Draw("colz text")
    h_corr[v].GetYaxis().SetTitleOffset(0.9)
    pt.Draw()
    c_corr.SaveAs("plots/sys/h_%s_%s_corr.pdf" % (v, MODE))
    c_corr.Update()

    OBJECTS.append(c)
    OBJECTS.append(c_cov)
    OBJECTS.append(c_frac)
    OBJECTS.append(c_corr)

input()
