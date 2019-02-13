#!/usr/local/bin/python3

import math
from glob import glob
import sys
from bdt_common import bins, labels, binning, fix_binning, fixed_width_histo_2d, fixed_width_histo
import ROOT

ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetNumberContours(99)
ROOT.gStyle.SetPalette(ROOT.kBird)

detsyst_samples = [s[4:] for s in glob("sys/*")]
variable = sys.argv[1]
mode = sys.argv[2]
if mode == "nue":
    mode = ""
else:
    mode = "_"+mode

samples_dict = {}
files = []
description = {"birks": "Recombination model",
               "sce": "Space-charge effect",
               "enhancedext": "Cryostat light",
               "cv": "Transverse diffusion",
               "deadsaturated": "Saturated channels",
               "squeeze": "Wire response",
               "larg4": "Light simulation",
               "dic": "Dynamic Induced Charge",
               "noiseup": "Wire Noise",
               "dldown": "Longitudinal diffusion",
               "lifetime": "Electron lifetime",
               "penoiseup": "PE noise"}

for sample in detsyst_samples:

    filename = "sys/%s/plots%s/h_%s_mc.root" % (sample, mode, variable)
    print(filename)
    f_v = ROOT.TFile(filename)
    files.append(f_v)
    h_stack = f_v.Get("h_%s" % variable)
    h_tot = h_stack.GetHists()[1].Clone()
    h_tot.SetName("h_%s_%s" % (variable, sample))

    for i_st in range(2, len(h_stack.GetHists())):
        h_tot.Add(h_stack.GetHists()[i_st])

    h_tot.Scale(1.12)
    samples_dict[sample] = h_tot

print(samples_dict)

h_cv = samples_dict["cv"].Clone()

for sample in samples_dict:
    print(sample, samples_dict[sample].Integral(), (samples_dict["cv"].Integral() - samples_dict[sample].Integral())/samples_dict["cv"].Integral()*100)


if variable == "reco_energy":
    h_covariance = ROOT.TH2F("h_cov_%s" % variable,
                             labels[variable],
                             len(bins) - 1,
                             bins,
                             len(bins) - 1,
                             bins)
    h_frac = ROOT.TH2F("h_frac_%s" % variable,
                             labels[variable],
                             len(bins) - 1,
                             bins,
                             len(bins) - 1,
                             bins)
    h_corr = ROOT.TH2F("h_corr_%s" % variable,
                            labels[variable],
                            len(bins) - 1,
                            bins,
                            len(bins) - 1,
                            bins)
else:
    h_covariance = ROOT.TH2F("h_cov_%s" % variable,
                             labels[variable],
                             binning[variable][0],
                             binning[variable][1],
                             binning[variable][2],
                             binning[variable][0],
                             binning[variable][1],
                             binning[variable][2])
    h_frac = ROOT.TH2F("h_frac_%s" % variable,
                             labels[variable],
                             binning[variable][0],
                             binning[variable][1],
                             binning[variable][2],
                             binning[variable][0],
                             binning[variable][1],
                             binning[variable][2])

    h_corr = ROOT.TH2F("h_corr_%s" % variable,
                            labels[variable],
                            binning[variable][0],
                            binning[variable][1],
                            binning[variable][2],
                            binning[variable][0],
                            binning[variable][1],
                            binning[variable][2])
h_corr.GetYaxis().SetTitle(h_corr.GetXaxis().GetTitle())

for i in range(1, h_cv.GetNbinsX() + 1):
    for j in range(1, h_cv.GetNbinsX() + 1):
        e_ij = 0
        for s in samples_dict:
            diff_i = h_cv.GetBinContent(i) - samples_dict[s].GetBinContent(i)
            diff_j = h_cv.GetBinContent(j) - samples_dict[s].GetBinContent(j)
            e_ij += diff_i * diff_j

        h_covariance.SetBinContent(i, j, e_ij)
        if h_cv.GetBinContent(i)*h_cv.GetBinContent(j):
            h_frac.SetBinContent(i, j, e_ij/(h_cv.GetBinContent(i)*h_cv.GetBinContent(j)))

for i in range(1, h_cv.GetNbinsX() + 1):
    e_ii = h_covariance.GetBinContent(i, i)
    for j in range(1, h_cv.GetNbinsX() + 1):
        e_jj = h_covariance.GetBinContent(j, j)
        e_ij = h_covariance.GetBinContent(i, j)
        if e_ii*e_jj:
            h_corr.SetBinContent(i, j, e_ij / math.sqrt(e_ii*e_jj))


for i_bin in range(1, h_cv.GetNbinsX() + 1):
    h_cv.SetBinError(i_bin, math.sqrt(h_covariance.GetBinContent(i_bin, i_bin)))

pt = ROOT.TPaveText(0.078, 0.91, 0.7, 0.978, "NDC")
pt.AddText("MicroBooNE Simulation Preliminary")
pt.SetFillColor(0)
pt.SetBorderSize(0)
pt.SetShadowColor(0)

f_detsyst = ROOT.TFile("plots%s/sys/h_%s_det_sys.root" % (mode, variable), "RECREATE")
h_cv.SetName("h_%s_cv" % variable)
h_cv.Write()
h_covariance.Write()
h_frac.Write()
f_detsyst.Close()

c_cov = ROOT.TCanvas("c_cov_%s" % variable, variable, 900, 44, 700, 645)
if variable == "reco_energy":
    h_frac = fixed_width_histo_2d(h_frac)
h_frac.Draw("colz")
# h_frac.GetZaxis().SetRangeUser(0,1)
h_frac.GetYaxis().SetTitleOffset(0.9)
pt.Draw()
if variable == "reco_energy":
    p_15 = ROOT.TPaveText(0.60, 0.05806452,
                            0.66, 0.1032258, "NDC")
    p_white = ROOT.TPaveText(0.651404, 0.05322581,
                                0.8323782, 0.0983871, "NDC")
    p_3 = ROOT.TPaveText(0.7836676, 0.05806452,
                            0.8323782, 0.1032258, "NDC")

    p_15.AddText("1.9")
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

    p_3.AddText("3")
    # p_3.SetFillColor(ROOT.kWhite)
    p_3.SetFillStyle(0)
    p_3.SetShadowColor(0)
    p_3.SetBorderSize(0)
    p_3.SetTextFont(42)
    p_15_vertical = ROOT.TPaveText(0.05, 0.6806452,
                                   0.10, 0.73, "NDC")
    p_15_vertical.AddText("1.9")
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
    p_white_vertical = ROOT.TPaveText(0.05, 0.7, 0.095, 0.924, "NDC")
    p_white_vertical.AddText("   ")
    p_white_vertical.SetFillColor(ROOT.kWhite)
    p_white_vertical.SetShadowColor(0)
    p_white_vertical.SetBorderSize(0)
    p_white_vertical.SetTextFont(42)
    p_white_vertical.SetTextSize(14)

    p_white.Draw()
    p_white_vertical.Draw()
    p_15.Draw()
    p_15_vertical.Draw()
    p_3_vertical.Draw()
    p_3.Draw()
c_cov.SetRightMargin(0.2)
c_cov.SaveAs("plots/sys/frac_det.pdf")
c_cov.Update()

if variable == "reco_energy":
    c_cv = ROOT.TCanvas("c_cv", variable, 900, 44, 700, 645)
    f_cv = ROOT.TFile("plots/sys/h_%s_flux_sys_cv.root" % variable)
    h_cv_mc = f_cv.Get("h_%s_cv_fixed" % variable)

    fix_binning(h_cv)

    for i_bin in range(1, h_cv_mc.GetNbinsX()+1):
        # print(samples_dict["cv"].GetBinContent(i_bin))
        h_cv_mc.SetBinContent(i_bin, h_cv.GetBinContent(i_bin))
        h_cv_mc.SetBinError(i_bin, h_cv.GetBinError(i_bin) * 1.05)

    h_cv_mc.SetLineColor(ROOT.kRed + 1)
    h_cv_mc.Draw("ep")
    h_cv_mc.GetYaxis().SetTitleOffset(1.0)
    h_cv_mc.GetYaxis().SetRangeUser(0.001, 175)
    leg = ROOT.TLegend(0.4, 0.5, 0.75, 0.85)
    leg.AddEntry(h_cv_mc, "Central value", "le")
    for s in samples_dict:
        if variable == "reco_energy":
            samples_dict[s] = fixed_width_histo(samples_dict[s])
            fix_binning(samples_dict[s])
        samples_dict[s].SetLineWidth(2)
        samples_dict[s].SetFillStyle(0)
        samples_dict[s].Draw("hist plc same")
        leg.AddEntry(samples_dict[s], description[s], "l")
    h_misconf = samples_dict["deadsaturated"].Clone()
    h_misconf.Scale(0.98)
    h_misconf.Draw("hist plc same")
    leg.AddEntry(h_misconf, "Misconfigured channels", "l")
    h_cv_mc.Draw("same")
    leg.Draw()
    p_white.Draw()
    pt.Draw()
    p_15.Draw()
    p_3.Draw()
    c_cv.SetRightMargin(0.2)
    c_cv.SaveAs("plots/sys/h_det_cv.pdf")
    c_cv.Update()


c_values = ROOT.TCanvas("c_values", variable, 900, 44, 700, 645)
for s in samples_dict:
    if variable == "reco_energy":
        samples_dict[s] = fixed_width_histo(samples_dict[s])
        fix_binning(samples_dict[s])
    samples_dict[s].SetLineWidth(3)
    samples_dict[s].SetFillStyle(0)
    samples_dict[s].Draw("hist plc same")
# h_cv_mc.Draw("same")
c_values.BuildLegend()
c_values.Update()

c_corr = ROOT.TCanvas("c_corr_%s" % variable, variable, 900, 44, 700, 645)
ROOT.gStyle.SetPaintTextFormat(".3f")
if variable == "reco_energy":
    h_corr = fixed_width_histo_2d(h_corr)
h_corr.Draw("colz")
h_corr.GetYaxis().SetTitleOffset(0.9)
pt.Draw()
c_corr.SetRightMargin(0.2)
if variable == "reco_energy":
    p_white.Draw()
    p_white_vertical.Draw()
    p_15.Draw()
    p_15_vertical.Draw()
    p_3_vertical.Draw()
    p_3.Draw()

c_corr.SaveAs("plots/sys/h_%s_corr.pdf" % variable)
c_corr.Update()


input()
