#!/usr/local/bin/python3

import ROOT
from operator import mul
import functools
import math
from tqdm import tqdm
from bdt_common import pre_cuts, variables, spectators, labels, binning, N_UNI, bins, fix_binning, fixed_width_histo, fixed_width_histo_2d
from bdt_common import bdt_types, load_bdt, load_variables, apply_cuts
from bdt_common import BDT, MANUAL
import sys

if len(sys.argv) > 1:
    mode = sys.argv[1]
else:
    mode = "nue"


if MANUAL:
    folder = "_cuts"
elif BDT:
    folder = "_bdt"
else:
    folder = ""

if mode == "numu":
    folder = "_numu"
    BDT = False
    MANUAL = True
elif mode == "nc":
    folder = "_nc"
    BDT = False
    MANUAL = True

if len(sys.argv) > 3:
    SYS_VARIABLES = sys.argv[3:]
else:
    SYS_VARIABLES = ["reco_energy"]


ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetNumberContours(99)
SYS_MODE = sys.argv[2]

if SYS_MODE == "flux":
    ROOT.gStyle.SetPalette(ROOT.kMint)
else:
    ROOT.gStyle.SetPalette(ROOT.kLake)


chain = ROOT.TChain("mc_tree")
chain.Add("root_files/mc_file.root/mc_tree")
# chain.Add("root_files/nue_file.root/nue_tree")

total_entries = int(chain.GetEntries() / 1)

h_sys = {}
h_2d = {}
h_cv = {}
h_covariance = {}
h_frac = {}
h_corr = {}

vars = load_variables(chain)

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

ROOT.TMVA.Tools.Instance()
reader = ROOT.TMVA.Reader(":".join([
    "!V",
    "!Silent",
    "Color"]))
load_bdt(reader)


for ievt in tqdm(range(total_entries)):
    chain.GetEntry(ievt)

    if not pre_cuts(vars):
        continue

    bdt_values = {}

    for bdt_name in bdt_types:
        bdt_values[bdt_name] = reader.EvaluateMVA("BDT%s" % bdt_name)

    if apply_cuts(bdt_values, vars, BDT, MANUAL, mode):
        for name, var in vars.items():
            if name not in SYS_VARIABLES:
                continue
            for i, v in enumerate(var):
                if v > -999:
                    h_cv[name].Fill(v, chain.event_weight)
                    for u in range(N_UNI):
                        if SYS_MODE == "genie":
                            weight = chain.genie_weights[u]
                        elif SYS_MODE == "flux":
                            weight = chain.flux_weights[u]
                        else:
                            weight = chain.flux_weights[u] * chain.genie_weights[u]
                        if weight > 5:
                            weight = 1
                        h_sys[name][u].Fill(v, chain.event_weight * weight)


for n, b in vars.items():
    if n == "reco_energy":
        h_2d[n] = ROOT.TH2F("h_%s" % n,
                            labels[n],
                            len(bins) - 1,
                            bins,
                            int(h_sys[n][0].GetMaximum()), 0, h_sys[n][0].GetMaximum() * 2)
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
        # fix_binning(h_cv[n])
        h_cv[n] = fixed_width_histo(h_cv[n])

    for u in range(N_UNI):
        if not h_sys[n][u].Integral() > 0:
            continue

        if n == "reco_energy":
            # fix_binning(h_sys[n][u])
            h_sys[n][u] = fixed_width_histo(h_sys[n][u])

        for i_bin in range(1, h_sys[n][u].GetNbinsX() + 1):
            value = h_sys[n][u].GetBinContent(i_bin)
            center = h_sys[n][u].GetBinCenter(i_bin)

            diff = (h_cv[n].GetBinContent(i_bin) - h_sys[n][u].GetBinContent(i_bin))

            h_2d[n].Fill(center, value)
            sys_err[n][i_bin] += diff**2



for v in SYS_VARIABLES:
    for i in range(1, h_cv[v].GetNbinsX() + 1):
        for j in range(1, h_cv[v].GetNbinsX() + 1):
            e_ij = 0
            for u in range(N_UNI):
                diff_i = - h_cv[v].GetBinContent(i) + h_sys[v][u].GetBinContent(i)
                diff_j = - h_cv[v].GetBinContent(j) + h_sys[v][u].GetBinContent(j)

                e_ij += diff_i * diff_j

            e_ij /= N_UNI
            if h_cv[v].GetBinContent(i)*h_cv[v].GetBinContent(j):
                f_ij = e_ij/(h_cv[v].GetBinContent(i)*h_cv[v].GetBinContent(j))
                h_frac[v].SetBinContent(i, j, f_ij)
            h_covariance[v].SetBinContent(i, j, e_ij)

for v in SYS_VARIABLES:
    for i in range(1, h_cv[v].GetNbinsX() + 1):
        e_ii = h_covariance[v].GetBinContent(i, i)
        for j in range(1, h_cv[v].GetNbinsX() + 1):
            e_jj = h_covariance[v].GetBinContent(j, j)
            e_ij = h_covariance[v].GetBinContent(i, j)
            if e_ii*e_jj:
                h_corr[v].SetBinContent(i, j, e_ij / math.sqrt(e_ii*e_jj))

OBJECTS = []

pt = ROOT.TPaveText(0.078, 0.91, 0.7, 0.978, "NDC")
pt.AddText("MicroBooNE Simulation Preliminary")
pt.SetFillColor(0)
pt.SetBorderSize(0)
pt.SetShadowColor(0)

OBJECTS.append(pt)

for v in SYS_VARIABLES:

    for i_bin in range(1, h_cv[v].GetNbinsX() + 1):
        h_cv[v].SetBinError(i_bin, math.sqrt(sys_err[v][i_bin] / N_UNI))

    if v == "reco_energy":
        f_cov = ROOT.TFile("plots%s/sys/h_cov_reco_energy_%s.root" % (folder, SYS_MODE), "RECREATE")
        h_covariance[v].Write()
        f_cov.Close()
        h_covariance[v] = fixed_width_histo_2d(h_covariance[v])
        h_frac[v] = fixed_width_histo_2d(h_frac[v])
        h_corr[v] = fixed_width_histo_2d(h_corr[v])

    c = ROOT.TCanvas("c_%s" % v, "", 900, 44, 700, 645)
    # c.SetTopMargin(0.18)
    h_cv[v].Draw("ep")
    h_2d[v].GetZaxis().SetRangeUser(0, 30)
    h_2d[v].Draw("colz same")
    h_cv[v].Draw("ep same")
    h_cv[v].SetMarkerStyle(0)
    h_cv[v].SetMarkerColor(ROOT.kRed + 1)
    h_cv[v].SetLineColor(ROOT.kRed + 1)
    h_cv[v].GetYaxis().SetTitleOffset(0.9)
    h_cv[v].GetYaxis().SetRangeUser(0.001, h_cv[v].GetMaximum() * 1.2)
    if v == "reco_energy":
        p_15 = ROOT.TPaveText(0.744, 0.0555, 0.794, 0.101, "NDC")
        p_white = ROOT.TPaveText(0.78, 0.0555, 0.84, 0.091, "NDC")
        p_white2 = ROOT.TPaveText(0.875, 0.0545, 0.925, 0.091, "NDC")
        p_3 = ROOT.TPaveText(0.874, 0.0575, 0.924, 0.103, "NDC")

        p_15.AddText("1.5")
        p_15.SetFillStyle(0)
        p_15.SetShadowColor(0)
        p_15.SetBorderSize(0)
        p_15.SetTextFont(42)

        p_white.AddText("  ")
        p_white.SetFillColor(ROOT.kWhite)
        p_white.SetShadowColor(0)
        p_white.SetBorderSize(0)
        p_white.SetTextFont(42)
        p_white.SetTextSize(14)

        p_white2.AddText("  ")
        p_white2.SetFillColor(ROOT.kWhite)
        p_white2.SetShadowColor(0)
        p_white2.SetBorderSize(0)
        p_white2.SetTextFont(42)
        p_white2.SetTextSize(14)

        p_3.AddText("3")
        # p_3.SetFillColor(ROOT.kWhite)
        p_3.SetFillStyle(0)
        p_3.SetShadowColor(0)
        p_3.SetBorderSize(0)
        p_3.SetTextFont(42)
        p_white.Draw()
        p_white2.Draw()
        p_15.Draw()
        p_3.Draw()
    pt.Draw()

    f_cv = ROOT.TFile("plots%s/sys/h_%s_%s_sys.root" % (folder, v, SYS_MODE), "RECREATE")
    h_cv[v].Write()
    f_cv.Close()

    c.SaveAs("plots%s/sys/h_%s_%s_err.pdf" % (folder, v, SYS_MODE))
    c.Update()

    c_cov = ROOT.TCanvas("c_cov_%s" % v, v, 900, 44, 700, 645)
    ROOT.gStyle.SetPaintTextFormat(".3f")
    h_covariance[v].Draw("colz")
    h_covariance[v].GetYaxis().SetTitleOffset(0.9)
    pt.Draw()
    c_cov.SaveAs("plots%s/sys/h_%s_%s_cov.pdf" % (folder, v, SYS_MODE))
    c_cov.Update()

    c_frac = ROOT.TCanvas("c_frac_%s" % v, v, 900, 44, 700, 645)
    ROOT.gStyle.SetPaintTextFormat(".3f")
    h_frac[v].Draw("colz")
    h_frac[v].GetYaxis().SetTitleOffset(0.9)
    pt.Draw()
    if v == "reco_energy":
        p_15_vertical = ROOT.TPaveText(0.045, 0.745, 0.0959, 0.790, "NDC")
        p_15_vertical.AddText("1.5")
        p_15_vertical.SetFillStyle(0)
        p_15_vertical.SetShadowColor(0)
        p_15_vertical.SetBorderSize(0)
        p_15_vertical.SetTextFont(42)
        p_3_vertical = ROOT.TPaveText(0.05, 0.879, 0.0995, 0.924, "NDC")
        p_3_vertical.AddText("3")
        p_3_vertical.SetFillStyle(0)
        p_3_vertical.SetShadowColor(0)
        p_3_vertical.SetBorderSize(0)
        p_3_vertical.SetTextFont(42)
        p_white_vertical = ROOT.TPaveText(0.045, 0.879, 0.095, 0.924, "NDC")
        p_white_vertical.AddText("   ")
        p_white_vertical.SetFillColor(ROOT.kWhite)
        p_white_vertical.SetShadowColor(0)
        p_white_vertical.SetBorderSize(0)
        p_white_vertical.SetTextFont(42)
        p_white_vertical.SetTextSize(14)
        p_white2_vertical = ROOT.TPaveText(0.045, 0.788, 0.095, 0.833, "NDC")
        p_white2_vertical.AddText("   ")
        p_white2_vertical.SetFillColor(ROOT.kWhite)
        p_white2_vertical.SetShadowColor(0)
        p_white2_vertical.SetBorderSize(0)
        p_white2_vertical.SetTextFont(42)
        p_white2_vertical.SetTextSize(14)
        p_white.Draw()
        p_white2.Draw()
        p_white_vertical.Draw()
        p_white2_vertical.Draw()
        p_15.Draw()
        p_15_vertical.Draw()
        p_3_vertical.Draw()
        p_3.Draw()
    c_frac.SetRightMargin(0.2)
    c_frac.SaveAs("plots%s/sys/h_%s_%s_frac.pdf" % (folder, v, SYS_MODE))
    c_frac.Update()

    c_corr = ROOT.TCanvas("c_corr_%s" % v)
    ROOT.gStyle.SetPaintTextFormat(".3f")
    h_corr[v].Draw("colz")
    h_corr[v].GetYaxis().SetTitleOffset(0.9)
    pt.Draw()
    c_corr.SaveAs("plots%s/sys/h_%s_%s_corr.pdf" % (folder, v, SYS_MODE))
    c_corr.Update()

    OBJECTS.append(c)
    OBJECTS.append(c_cov)
    OBJECTS.append(c_frac)
    OBJECTS.append(c_corr)

input()
