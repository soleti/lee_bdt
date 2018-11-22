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
chain.Add("root_files/nue_file.root/nue_tree")

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
                        if weight > 2:
                            weight = 1
                        h_sys[name][u].Fill(v, chain.event_weight * weight)


for n, b in vars.items():
    if n == "reco_energy":
        h_2d[n] = ROOT.TH2F("h_%s" % n,
                            labels[n],
                            len(bins) - 1,
                            bins,
                            int(N_UNI/2), 0, h_sys[n][0].GetMaximum() * 2)
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

pt = ROOT.TPaveText(0.057, 0.905, 0.422, 0.951, "ndc")
pt.AddText("MicroBooNE Preliminary")
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

    f_cv = ROOT.TFile("plots%s/sys/h_%s_%s_sys.root" % (folder, v, SYS_MODE), "RECREATE")
    h_cv[v].Write()
    f_cv.Close()

    c.SaveAs("plots%s/sys/h_%s_%s_err.pdf" % (folder, v, SYS_MODE))
    c.Update()

    c_cov = ROOT.TCanvas("c_cov_%s" % v)
    ROOT.gStyle.SetPaintTextFormat(".3f")
    h_covariance[v].Draw("colz text")
    h_covariance[v].GetYaxis().SetTitleOffset(0.9)
    pt.Draw()
    c_cov.SaveAs("plots%s/sys/h_%s_%s_cov.pdf" % (folder, v, SYS_MODE))
    c_cov.Update()

    c_frac = ROOT.TCanvas("c_frac_%s" % v)
    ROOT.gStyle.SetPaintTextFormat(".3f")
    h_frac[v].Draw("colz text")
    h_frac[v].GetYaxis().SetTitleOffset(0.9)
    pt.Draw()
    c_frac.SaveAs("plots%s/sys/h_%s_%s_frac.pdf" % (folder, v, SYS_MODE))
    c_frac.Update()

    c_corr = ROOT.TCanvas("c_corr_%s" % v)
    ROOT.gStyle.SetPaintTextFormat(".3f")
    h_corr[v].Draw("colz text")
    h_corr[v].GetYaxis().SetTitleOffset(0.9)
    pt.Draw()
    c_corr.SaveAs("plots%s/sys/h_%s_%s_corr.pdf" % (folder, v, SYS_MODE))
    c_corr.Update()

    OBJECTS.append(c)
    OBJECTS.append(c_cov)
    OBJECTS.append(c_frac)
    OBJECTS.append(c_corr)

input()
