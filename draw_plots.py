#!/usr/local/bin/python3

import pickle
import math
import numpy as np
from root_numpy import hist2array
import ROOT

from bdt_common import variables, spectators, bins, total_data_bnb_pot
from bdt_common import description, total_pot, sigmaCalc

ROOT.gStyle.SetOptStat(0)

SYS_ERR = 0.1

DRAW_SUBTRACTION = False
DRAW_EXT = False
DRAW_LEE = False
DRAW_DATA = True
DRAW_COSMIC = False

OBJECTS = []
VARIABLES = variables + spectators
RECO_ENERGY = list(dict(VARIABLES)).index("reco_energy")

histograms_bnb = []
histograms_bnbext = []
histograms_mc = []
histograms_lee = []

histograms = [histograms_bnb, histograms_bnbext, histograms_mc, histograms_lee]


if DRAW_DATA:
    POST_SCALING = 1
else:
    POST_SCALING = 6.6e20 / total_data_bnb_pot

total_pot *= POST_SCALING


def set_axis(histogram, y_max = 0):
    histogram.GetYaxis().SetTitleSize(0.06)
    histogram.GetYaxis().SetTitleOffset(0.8)
    histogram.SetMinimum(0.01)
    if y_max != 0:
        histogram.SetMaximum(y_max)
    else:
        histogram.SetMaximum(histogram.GetMaximum() * 1.3)


def fix_binning(histogram, width=0.05):
    for bin_i in range(1, histogram.GetNbinsX() + 1):
        bin_width = histogram.GetBinWidth(bin_i)
        histogram.SetBinError(bin_i, histogram.GetBinError(
            bin_i) / (bin_width / width))
        histogram.SetBinContent(bin_i, histogram.GetBinContent(
            bin_i) / (bin_width / width))


def draw_top():
    pad_top = ROOT.TPad("pad_top", "", 0, 0.3, 1, 1)
    pad_top.SetBottomMargin(0)
    pad_top.Range(-22.21825, -2.003018, 202.5403, 2073.676)
    pad_top.SetTopMargin(0.2411347)
    pad_top.Draw()
    pad_top.cd()
    OBJECTS.append(pad_top)


def draw_ratio(num, den):
    pad_bottom = ROOT.TPad("pad_bottom", "", 0, 0, 1, 0.3)
    pad_bottom.Range(-22.5, -0.6346511, 202.5, 1.99)

    pad_bottom.SetFrameBorderMode(0)
    pad_bottom.SetFrameBorderMode(0)
    pad_bottom.SetBorderMode(0)
    pad_bottom.SetTopMargin(0)
    pad_bottom.SetBottomMargin(0.275614)
    pad_bottom.Draw()
    pad_bottom.cd()
    OBJECTS.append(pad_bottom)
    h_ratio = num.Clone()
    h_ratio_sys = num.Clone()
    den_sys = den.Clone()
    den_sys.SetName("den_sys%i" % i)
    for k in range(1, den_sys.GetNbinsX() + 1):
        stat_err = den_sys.GetBinError(k)
        sys_err = SYS_ERR * den_sys.GetBinContent(k)
        den_sys.SetBinError(k, math.sqrt(stat_err**2 + sys_err**2))

    OBJECTS.append(h_ratio)
    OBJECTS.append(h_ratio_sys)

    h_ratio.SetName("h_ratio%i" % i)
    h_ratio_sys.SetName("h_ratio_sys%i" % i)

    h_ratio.GetYaxis().SetRangeUser(0.01, 1.99)
    h_ratio.Divide(den)
    h_ratio_sys.Divide(den_sys)

    h_ratio.GetXaxis().SetLabelFont(42)
    h_ratio.GetXaxis().SetLabelSize(0.13)
    h_ratio.GetXaxis().SetTitleSize(0.13)
    h_ratio.GetXaxis().SetTitleOffset(0.91)
    h_ratio.GetYaxis().SetTitle("Data / MC")
    h_ratio.GetYaxis().SetNdivisions(509)
    h_ratio.GetYaxis().SetLabelFont(42)
    h_ratio.GetYaxis().SetLabelSize(0.13)
    h_ratio.GetYaxis().SetTitleSize(0.13)
    h_ratio.GetYaxis().SetTitleOffset(0.36)

    h_ratio.Draw("hist")
    h_ratio.SetMarkerSize(0)
    h_ratio_2 = h_ratio.Clone()
    h_ratio_2.SetFillStyle(3002)
    h_ratio_2.SetFillColor(1)
    h_ratio_2.SetName("h_ratio_2%i" % i)
    OBJECTS.append(h_ratio_2)
    h_ratio_2.Draw("e2 same")
    h_ratio_sys.Draw("e1 same")
    h_ratio_sys.SetMarkerSize(0)

    line = ROOT.TLine(h_ratio.GetXaxis().GetXmin(), 1,
                      h_ratio.GetXaxis().GetXmax(), 1)
    line.SetLineWidth(2)
    line.SetLineStyle(2)
    line.Draw()
    OBJECTS.append(line)


legends = []
l_errs = []
samples = ["bnb", "bnbext", "mc", "lee"]

for name, var in VARIABLES:

    for histos, s in zip(histograms, samples):
        f = ROOT.TFile("plots/h_%s_%s.root" % (name, s))
        h = f.Get("h_%s" % name)
        OBJECTS.append(f)
        histos.append(h)

    legend = ROOT.TLegend(0.094555, 0.78502, 0.89234, 0.97919, "", "brNDC")
    legend.SetTextSize(16)
    legend.SetTextFont(63)
    legend.SetHeader("MicroBooNE Preliminary %.1e POT" % total_pot)
    legend.SetTextFont(43)
    legend.SetNColumns(2)
    l_err = ROOT.TLegend(0.4713467, 0.5990783, 0.7421203,
                         0.6935484, "", "brNDC")
    l_err.SetTextSize(18)
    l_err.SetTextFont(43)

    l_errs.append(l_err)
    legends.append(legend)

for i, h in enumerate(histograms_mc):
    for d, histo in zip(description, h.GetHists()):
        if histo.Integral():
            legends[i].AddEntry(
                histo,
                "{}: {:.0f} events".format(d, histo.Integral() * POST_SCALING), "f")

    if DRAW_DATA:
        legends[i].AddEntry(histograms_bnb[i],
                            "Data BNB: {:.0f} events"
                            .format(histograms_bnb[i].Integral()
                                    * POST_SCALING), "lep")
    if DRAW_LEE:
        histograms_lee[i].SetLineColor(ROOT.kBlack)
        histograms_lee[i].SetFillColor(ROOT.kGreen - 2)
        histograms_lee[i].SetFillStyle(3002)
        if i != RECO_ENERGY:
            legends[i].AddEntry(histograms_lee[i],
                                "Low-energy excess: {:.0f} events"
                                .format(histograms_lee[i].Integral()
                                        * POST_SCALING), "f")


for h in histograms_mc:
    h.GetHists()[4].SetFillStyle(3002)


if DRAW_SUBTRACTION:
    for h_bnb, h_ext, leg in zip(histograms_bnb, histograms_bnbext, legends):
        for i in range(1, h_bnb.GetNbinsX() + 1):
            bnb = h_bnb.GetBinContent(i)
            ext = h_ext.GetBinContent(i)

            h_bnb.SetBinContent(i, bnb - ext)
            if bnb > 0:
                err1 = h_bnb.GetBinError(i)
                err2 = h_ext.GetBinError(i)
                h_bnb.SetBinError(i, math.sqrt(err1**2 + err2**2))

        leg.AddEntry(
            h_bnb, "Data BNB - BNB EXT: {:.0f} events".format(h_bnb.Integral() * POST_SCALING), "lep")


h_lee = ROOT.TH1F("h_lee", "", len(bins) - 1, bins)

scaling = [0, 0, 0, 6.0920944819073988, 3.6447414342239273, 3.2123920194399913,
           2.6504659907742409, 3.2558450032216988, 2.5826310533377432,
           2, 1, 1, 1, 1]

h_lee.SetFillColor(ROOT.kGreen - 2)
h_lee.SetFillStyle(3002)
h_lee.SetLineColor(1)
for i, scale in enumerate(scaling):
    if scaling[i] - 1 > 0:
        h_lee.SetBinContent(i + 1, histograms_mc[RECO_ENERGY].GetHists()[0].
                            GetBinContent(i + 1) * (scale - 1))

histograms_lee[RECO_ENERGY] = h_lee

if DRAW_LEE:
    legends[RECO_ENERGY].AddEntry(histograms_lee[RECO_ENERGY],
                                  "Low-energy excess: {:.0f} events"
                                  .format(histograms_lee[RECO_ENERGY].Integral() * POST_SCALING),
                                  "f")

c_energy = ROOT.TCanvas("c_energy")
OBJECTS.append(c_energy)
signal_spectrum = histograms_mc[RECO_ENERGY].GetHists()[0].Clone()
signal_spectrum.Scale(POST_SCALING)
print("Spectrum integral", signal_spectrum.Integral())
signal_spectrum.Draw("hist")
fix_binning(signal_spectrum)

pt = ROOT.TPaveText(0.098, 0.905, 0.576, 0.989, "ndc")
pt.AddText("MicroBooNE Preliminary - 6.6e20 POT")
pt.SetFillColor(0)
pt.SetBorderSize(0)
pt.SetShadowColor(0)
pt.Draw()


OBJECTS.append(pt)
c_energy.Update()

for i in range(len(VARIABLES)):
    if histograms_mc[i].GetHists()[0].Integral() > 0:

        c = ROOT.TCanvas("c%i" % i, "", 900, 44, 700, 645)

        histograms_bnbext[i].SetLineColor(ROOT.kBlack)
        histograms_bnbext[i].SetMarkerStyle(20)

        histograms_bnb[i].SetLineColor(1)
        histograms_bnb[i].SetMarkerStyle(20)

        # if not DRAW_SUBTRACTION and DRAW_EXT:
        #     if i == RECO_ENERGY:
        #         fix_binning(histograms_cosmic[i])
        #     histograms_mc[i].Add(histograms_bnbext[i])
        #
        # if not DRAW_EXT:
        #     if i == RECO_ENERGY:
        #         fix_binning(histograms_bnbext[i])
        #     histograms_mc[i].Add(histograms_cosmic[i])

        h_mc_err = histograms_mc[i].GetHists()[0].Clone()
        h_mc_err.SetName("h_mc_err%i" % i)
        h_mc_err.Reset()

        h_mc_err_nobinning = h_mc_err.Clone()
        h_mc_err_nobinning.SetName("h_mc_err_nobin%i" % i)

        for j in range(histograms_mc[i].GetNhists()):
            histograms_mc[i].GetHists()[j].Scale(POST_SCALING)

            h_mc_err_nobinning.Add(histograms_mc[i].GetHists()[j])

            if i == RECO_ENERGY:
                fix_binning(histograms_mc[i].GetHists()[j])

            h_mc_err.Add(histograms_mc[i].GetHists()[j])

        for k in range(1, h_mc_err_nobinning.GetNbinsX() + 1):
            new_err = h_mc_err_nobinning.GetBinContent(k)
            h_mc_err_nobinning.SetBinError(k, math.sqrt(new_err))

        h_mc_err_sys = h_mc_err.Clone()
        h_mc_err_sys.SetName("h_mc_err_sys%i" % i)
        for k in range(1, h_mc_err.GetNbinsX() + 1):
            stat_err = h_mc_err.GetBinError(k)
            sys_err = SYS_ERR * h_mc_err.GetBinContent(k)
            h_mc_err_sys.SetBinError(k, math.sqrt(stat_err**2 + sys_err**2))

        if DRAW_LEE:
            histograms_lee[i].Scale(POST_SCALING)
            if i == RECO_ENERGY:
                fix_binning(histograms_lee[i])
                print("Sigma stat", sigmaCalc(
                    histograms_lee[i], h_mc_err_nobinning))

            histograms_mc[i].Add(histograms_lee[i])

        if DRAW_DATA:
            draw_top()
        else:
            c.SetTopMargin(0.2274194)

        histograms_mc[i].Draw("hist")

        if DRAW_DATA:
            set_axis(histograms_mc[i], histograms_bnb[i].GetMaximum() * 1.35)
        else:
            set_axis(histograms_mc[i])

        if i == RECO_ENERGY:
            new_max = h_mc_err.GetMaximum() + histograms_lee[i].GetMaximum()
            histograms_mc[i].SetMaximum(new_max * 1.3)

        h_mc_err.SetFillStyle(3002)
        h_mc_err.SetFillColor(1)
        h_mc_err_sys.SetFillStyle(3002)
        h_mc_err_sys.SetFillColor(1)
        l_errs[i].AddEntry(h_mc_err, "Stat. uncertainties", "lf")
        l_errs[i].AddEntry(h_mc_err_sys, "Stat. #oplus 10% sys. uncertainties", "le")
        l_errs[i].Draw("same")
        h_mc_err.Draw("e2 same")

        h_mc_err_sys.Draw("e1 same")

        if DRAW_DATA:
            if i == RECO_ENERGY:
                fix_binning(histograms_bnb[i])

            p_chi2 = ROOT.TPaveText(0.65, 0.62, 0.86, 0.73, "NDC")
            p_chi2.AddText("#chi^{2} / ndf = %.2f" %
                           histograms_bnb[i].Chi2Test(h_mc_err_sys, "UW"))
            histograms_bnb[i].Draw("ep same")
            # p_chi2.Draw("same")
            p_chi2.SetFillStyle(0)
            p_chi2.SetBorderSize(0)
            OBJECTS.append(p_chi2)

        OBJECTS.append(h_mc_err)
        OBJECTS.append(h_mc_err_sys)

        legends[i].Draw("same")
        c.cd()

        if DRAW_DATA:
            if i == RECO_ENERGY:
                print("Data/MC ratio: ",
                      histograms_bnb[i].Integral() / h_mc_err_sys.Integral())
            draw_ratio(histograms_bnb[i], h_mc_err)

        c.Update()
        c.SaveAs("plots/%s.pdf" % histograms_bnb[i].GetName())

        OBJECTS.append(c)

h_true_e = ROOT.THStack(
    "h_true_e", ";Reco. energy [GeV]; N. Entries / 0.05 GeV")
with open("a_e_true_reco.bin", "rb") as f:
    a_e_true_reco = pickle.load(f)

for j in range(histograms_mc[RECO_ENERGY].GetNhists()):
    h_clone = histograms_mc[RECO_ENERGY].GetHists()[j].Clone()
    h_fixed = ROOT.TH1F("h_fixed%i" % j, "", len(bins) - 1, 0, len(bins) - 1)
    a_e_reco = hist2array(h_clone)
    a_true = np.dot(a_e_true_reco, a_e_reco)

    for i in range(1, h_clone.GetNbinsX() + 1):
        h_clone.SetBinContent(i, a_true[i - 1])
        h_clone.SetBinError(i, math.sqrt(a_true[i - 1]))
    fix_binning(h_clone)
    h_clone.Scale(POST_SCALING)
    ax = h_fixed.GetXaxis()
    for i in range(1, h_clone.GetNbinsX() + 1):
        h_fixed.SetBinContent(i, h_clone.GetBinContent(i))
        h_fixed.GetXaxis().SetBinLabel(i, "")

    h_fixed.SetLineColor(ROOT.kBlack)
    h_fixed.SetFillColor(h_clone.GetFillColor())
    h_fixed.SetFillStyle(h_clone.GetFillStyle())

    OBJECTS.append(h_clone)
    h_true_e.Add(h_fixed)

c_true = ROOT.TCanvas("c_true")
h_true_e.Draw("hist")
legends[RECO_ENERGY].Draw()
c_true.SetTopMargin(0.2274194)

ax = ROOT.TGaxis(0, 0, len(bins) - 1, 0, 0, len(bins) - 1, 515, "")
for i, i_bin in enumerate(bins):
    ax.ChangeLabel(i + 1, -1, -1, -1, -1, -1,
                   "{0}".format(str(round(i_bin, 2) if i_bin % 1 else int(i_bin))))
ax.SetLabelFont(42)
ax.SetLabelSize(0.05)
ax.Draw()
c_true.Update()

if DRAW_COSMIC:
    legend_cosmic = ROOT.TLegend(0.099, 0.909, 0.900, 0.987, "", "brNDC")
    legend_cosmic.AddEntry(histograms_bnbext[RECO_ENERGY], "Data EXT: {:.0f} events".format(
        histograms_bnbext[RECO_ENERGY].Integral() * POST_SCALING), "lep")
    legend_cosmic.AddEntry(histograms_mc[RECO_ENERGY].GetHists()[
                           1], "CORSIKA in-time Monte Carlo: integral normalized", "f")
    legend_cosmic.SetNColumns(2)

    for i in range(len(VARIABLES)):
        h_ext = histograms_bnbext[i].Clone()
        h_intime = histograms_mc[i].GetHists()[1].Clone()
        if i == RECO_ENERGY:
            fix_binning(h_ext)
        h_ext.Scale(POST_SCALING)

        if h_intime.Integral() > 0:
            c_cosmic = ROOT.TCanvas("c{}_c".format(i), "", 900, 44, 700, 645)
            OBJECTS.append(h_ext)
            OBJECTS.append(h_intime)
            draw_top()

            h_intime.Draw("hist")
            h_ext.Draw("ep same")
            set_axis(h_intime)
            legend_cosmic.Draw()
            if i == RECO_ENERGY:
                print(h_ext.Integral() / h_intime.Integral(), i)

            c_cosmic.cd()
            draw_ratio(h_ext, h_intime)
            c_cosmic.cd()

            c_cosmic.Update()
            c_cosmic.SaveAs("plots/%s_cosmic.pdf" % h_intime.GetName())
            OBJECTS.append(c_cosmic)

input()
