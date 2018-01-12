#!/usr/local/bin/python3

import ROOT
import math
import numpy as np
from bdt_common import variables, spectators, bins, colors
from bdt_common import description, total_pot, sigmaCalc, total_data_bnb_pot
from root_numpy import hist2array
import pickle
from array import array

ROOT.gStyle.SetOptStat(0)

pads = []
h_ratios = []
lines = []

variables = variables + spectators

reco_energy = list(dict(variables)).index("reco_energy")

histograms_bnb = []
histograms_bnbext = []
histograms_mc = []
histograms_lee = []

histograms = [histograms_bnb, histograms_bnbext, histograms_mc, histograms_lee]

post_scaling = 1#6.6e20 / total_data_bnb_pot
total_pot *= post_scaling

draw_subtraction = False
draw_ext = False
draw_lee = False
draw_data = True


def set_axis(histo):
    histo.GetYaxis().SetTitleSize(0.06)
    histo.GetYaxis().SetTitleOffset(0.8)
    histo.SetMinimum(0.1)
    histo.SetMaximum(histo.GetMaximum() * 1.3)


def fix_binning(histo, width=0.05):
    for k in range(1, histo.GetNbinsX() + 1):
        bin_width = histo.GetBinWidth(k)
        histo.SetBinError(k, histo.GetBinError(k) / (bin_width / width))
        histo.SetBinContent(k, histo.GetBinContent(k) / (bin_width / width))


def draw_top():
    pad_top = ROOT.TPad("pad_top", "", 0, 0.3, 1, 1)
    pad_top.SetBottomMargin(0)
    pad_top.Range(-22.21825, -2.003018, 202.5403, 2073.676)
    pad_top.SetTopMargin(0.2411347)
    pad_top.Draw()
    pad_top.cd()
    pads.append(pad_top)


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
    pads.append(pad_bottom)
    h_ratio = num.Clone()
    h_ratios.append(h_ratio)

    h_ratio.SetName("h_ratio%i" % i)

    h_ratio.GetYaxis().SetRangeUser(0.01, 1.99)
    h_ratio.Divide(den)

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
    h_ratio.Draw("ep")
    line = ROOT.TLine(h_ratio.GetXaxis().GetXmin(), 1,
                      h_ratio.GetXaxis().GetXmax(), 1)
    line.SetLineWidth(2)
    line.SetLineStyle(2)
    line.Draw()
    lines.append(line)


files = []
legends = []

samples = ["bnb", "bnbext", "mc", "lee"]

for name, var in variables:

    for i, s in enumerate(samples):
        f = ROOT.TFile("plots/h_%s_%s.root" % (name, s))
        h = f.Get("h_%s" % name)
        files.append(f)
        histograms[i].append(h)

    legend = ROOT.TLegend(0.094555, 0.78502, 0.89234, 0.97919, "", "brNDC")
    legend.SetTextSize(16)
    legend.SetTextFont(63)
    legend.SetHeader("MicroBooNE Preliminary %.1e POT" % total_pot)
    legend.SetTextFont(43)
    legend.SetNColumns(2)

    legends.append(legend)

for i, h in enumerate(histograms_mc):
    for j in range(h.GetNhists()):
        if h.GetHists()[j].Integral():
            legends[i].AddEntry(
                h.GetHists()[j],
                "{}: {:.0f} events".format(description[j],
                                           h.GetHists()[j].
                                           Integral() * post_scaling), "f")

    # if draw_ext:
    #     legends[i].AddEntry(histograms_bnbext[i],
    #                         "Data EXT: {:.0f} events"
    #                         .format(histograms_bnbext[i].Integral()
    #                                 * post_scaling), "f")
    # else:
    #     legends[i].AddEntry(histograms_cosmic[i],
    #                         "Cosmic in-time: {:.0f} events"
    #                         .format(histograms_cosmic[i].Integral()
    #                                 * post_scaling), "f")
    if draw_data:
        legends[i].AddEntry(histograms_bnb[i],
                            "Data BNB: {:.0f} events"
                            .format(histograms_bnb[i].Integral()
                                    * post_scaling), "lep")
    if draw_lee:
        histograms_lee[i].SetLineColor(ROOT.kBlack)
        histograms_lee[i].SetFillColor(ROOT.kGreen - 2)
        histograms_lee[i].SetFillStyle(3002)
        if i != reco_energy:
            legends[i].AddEntry(histograms_lee[i],
                                "Low-energy excess: {:.0f} events"
                                .format(histograms_lee[i].Integral()
                                        * post_scaling), "f")


for h in histograms_mc:
    h.GetHists()[3].SetFillStyle(3002)


if draw_subtraction:
    for i in range(len(histograms_bnb)):
        for j in range(1, histograms_bnb[i].GetNbinsX() + 1):
            bnb = histograms_bnb[i].GetBinContent(j)
            ext = histograms_bnbext[i].GetBinContent(j)

            histograms_bnb[i].SetBinContent(j, bnb - ext)
            if bnb > 0:
                err1 = histograms_bnb[i].GetBinError(j)
                err2 = histograms_bnbext[i].GetBinError(j)
                histograms_bnb[i].SetBinError(j, math.sqrt(err1**2 + err2**2))

        legends[i].AddEntry(histograms_bnb[i], "Data BNB - BNB EXT: {:.0f} events"
                           .format(histograms_bnb[i].Integral() * post_scaling),
                           "lep")


legend_cosmic = ROOT.TLegend(0.099, 0.909, 0.900, 0.987, "", "brNDC")
legend_cosmic.AddEntry(histograms_bnbext[reco_energy],
                       "Data EXT: {:.0f} events"
                       .format(histograms_bnbext[reco_energy].Integral()
                               * post_scaling), "lep")
legend_cosmic.AddEntry(histograms_mc[reco_energy].GetHists()[1],
                       "CORSIKA in-time Monte Carlo: integral normalized", "f")
legend_cosmic.SetNColumns(2)

canvases = []
h_errs = []
h_errs_sys = []

canvases_cosmic = []
paves = []

h_lee = ROOT.TH1F("h_lee", "", len(bins) - 1, bins)

scaling = [0,0,0,6.0920944819073988, 3.6447414342239273, 3.2123920194399913,
           2.6504659907742409, 3.2558450032216988, 2.5826310533377432,
           2, 1, 1, 1, 1]

h_lee.SetFillColor(ROOT.kGreen - 2)
h_lee.SetFillStyle(3002)
h_lee.SetLineColor(1)
for i in range(len(scaling)):
    if scaling[i] - 1 > 0:
        h_lee.SetBinContent(i + 1, histograms_mc[reco_energy].GetHists()[0].
                            GetBinContent(i + 1) * (scaling[i] - 1))

histograms_lee[reco_energy] = h_lee

legends[reco_energy].AddEntry(histograms_lee[reco_energy],
                              "Low-energy excess: {:.0f} events"
                              .format(histograms_lee[reco_energy].Integral()
                                      * post_scaling), "f")

c_energy = ROOT.TCanvas("c_energy")
canvases.append(c_energy)
signal_spectrum = histograms_mc[reco_energy].GetHists()[0].Clone()
signal_spectrum.Scale(post_scaling)
print("Spectrum integral", signal_spectrum.Integral())
signal_spectrum.Draw("hist")
fix_binning(signal_spectrum)

pt = ROOT.TPaveText(0.098, 0.905, 0.576, 0.989, "ndc")
pt.AddText("MicroBooNE Preliminary - 6.6e20 POT")
pt.SetFillColor(0)
pt.SetBorderSize(0)
pt.SetShadowColor(0)
pt.Draw()
canvases.append(pt)
c_energy.Update()

h_true_e = ROOT.THStack("h_true_e", ";#nu_{e} energy [GeV]; N. Entries / 0.05 GeV")
with open("a_e_true_reco.bin", "rb") as f:
    a_e_true_reco = pickle.load(f)

clones = []
for j in range(histograms_mc[reco_energy].GetNhists()):
    h_clone = histograms_mc[reco_energy].GetHists()[j].Clone()
    h_fixed = ROOT.TH1F("h_fixed%i" % j, "", len(bins) - 1, 0, len(bins) - 1)
    a_e_reco = hist2array(h_clone)
    a_true = np.dot(a_e_true_reco, a_e_reco)

    for i in range(1, h_clone.GetNbinsX()+1):
        h_clone.SetBinContent(i, a_e_reco[i - 1])
        h_clone.SetBinError(i, math.sqrt(a_e_reco[i - 1]))
    fix_binning(h_clone)
    h_clone.Scale(post_scaling)
    ax = h_fixed.GetXaxis()
    for i in range(1, h_clone.GetNbinsX() + 1):
        h_fixed.SetBinContent(i, h_clone.GetBinContent(i))
        h_fixed.GetXaxis().SetBinLabel(i, "")

    h_fixed.SetLineColor(ROOT.kBlack)
    h_fixed.SetFillColor(h_clone.GetFillColor())
    h_fixed.SetFillStyle(h_clone.GetFillStyle())

    clones.append(h_clone)
    h_true_e.Add(h_fixed)

for i in range(len(variables)):
    if histograms_mc[i].GetHists()[0].Integral() > 0:

        c = ROOT.TCanvas("c%i" % i, "", 900, 44, 700, 645)

        histograms_bnbext[i].SetLineColor(ROOT.kBlack)
        histograms_bnbext[i].SetMarkerStyle(20)

        histograms_bnb[i].SetLineColor(1)
        histograms_bnb[i].SetMarkerStyle(20)

        # if not draw_subtraction and draw_ext:
        #     if i == reco_energy:
        #         fix_binning(histograms_cosmic[i])
        #     histograms_mc[i].Add(histograms_bnbext[i])
        #
        # if not draw_ext:
        #     if i == reco_energy:
        #         fix_binning(histograms_bnbext[i])
        #     histograms_mc[i].Add(histograms_cosmic[i])

        h_mc_err = histograms_mc[i].GetHists()[0].Clone()
        h_mc_err.SetName("h_mc_err%i" % i)
        h_mc_err.Reset()

        h_mc_err_nobinning = h_mc_err.Clone()
        h_mc_err_nobinning.SetName("h_mc_err_nobin%i" % i)

        for j in range(histograms_mc[i].GetNhists()):
            histograms_mc[i].GetHists()[j].Scale(post_scaling)

            h_mc_err_nobinning.Add(histograms_mc[i].GetHists()[j])

            if i == reco_energy:
                fix_binning(histograms_mc[i].GetHists()[j])

            h_mc_err.Add(histograms_mc[i].GetHists()[j])

        for k in range(1, h_mc_err_nobinning.GetNbinsX() + 1):
            new_err = h_mc_err_nobinning.GetBinContent(k)
            h_mc_err_nobinning.SetBinError(k, math.sqrt(new_err))

        h_mc_err_sys = h_mc_err.Clone()
        h_mc_err_sys.SetName("h_mc_err_sys%i" % i)
        for k in range(1, h_mc_err.GetNbinsX() + 1):
            new_err = h_mc_err.GetBinContent(k)
            h_mc_err_sys.SetBinError(k, math.sqrt(new_err))

        h_mc_err_stat = h_mc_err_sys.Clone()
        h_mc_err_stat.SetName("h_mc_err_stat%i" % i)
        for k in range(1, h_mc_err_stat.GetNbinsX() + 1):
            sys_err = (0.2 * h_mc_err_sys.GetBinContent(k))**2
            new_err = h_mc_err_sys.GetBinError(k)**2 + sys_err
            h_mc_err_stat.SetBinError(k, math.sqrt(new_err))

        if draw_lee:
            histograms_lee[i].Scale(post_scaling)
            if i == reco_energy:
                fix_binning(histograms_lee[i])
                print("Sigma stat", sigmaCalc(histograms_lee[i], h_mc_err_nobinning))
                print("Sigma 20 sys", sigmaCalc(histograms_lee[i], h_mc_err_nobinning, 0.2))

            histograms_mc[i].Add(histograms_lee[i])

        if draw_data:
            draw_top()
        else:
            c.SetTopMargin(0.2274194)

        histograms_mc[i].Draw("hist")
        set_axis(histograms_mc[i])

        if i == reco_energy:
            new_max = h_mc_err.GetMaximum() + histograms_lee[i].GetMaximum()
            histograms_mc[i].SetMaximum(new_max * 1.3)

        h_mc_err.SetFillStyle(3002)
        h_mc_err.SetFillColor(1)
        h_mc_err_sys.SetFillStyle(3002)
        h_mc_err_sys.SetFillColor(1)
        if i == reco_energy:
            h_mc_err_sys.Draw("e1 same")
        else:
            h_mc_err.Draw("e1 same")

        if draw_data:
            if i == reco_energy:
                fix_binning(histograms_bnb[i])

            # histograms_bnb[i].Scale(post_scaling)
            p_chi2 = ROOT.TPaveText(0.65, 0.62, 0.86, 0.73, "NDC")
            p_chi2.AddText("#chi^{2} / ndf = %.2f" % histograms_bnb[i].Chi2Test(h_mc_err_sys, "UW"))
            histograms_bnb[i].Draw("ep same")
            p_chi2.Draw("same")
            p_chi2.SetFillStyle(0)
            p_chi2.SetBorderSize(0)
            paves.append(p_chi2)

        h_errs.append(h_mc_err)
        h_errs_sys.append(h_mc_err_sys)
        h_errs_sys.append(h_mc_err_stat)

        legends[i].Draw("same")
        c.cd()

        if draw_data:
            if i == reco_energy:
                print("Data/MC ratio: ", histograms_bnb[i].Integral()/h_mc_err_sys.Integral())
            draw_ratio(histograms_bnb[i], h_mc_err_sys)
        #histograms_bnb[i].Draw("same")
        c.Update()
        c.SaveAs("plots/%s.pdf" % histograms_bnb[i].GetName())

        canvases.append(c)

c_true = ROOT.TCanvas("c_true")
h_true_e.Draw("hist")
legends[reco_energy].Draw()
c_true.SetTopMargin(0.2274194)

ax = ROOT.TGaxis(0, 0, len(bins) - 1, 0, 0, len(bins) - 1 , 515, "")
for i, bin in enumerate(bins):
    ax.ChangeLabel(i+1, -1, -1, -1, -1, -1, "{0}".format(str(round(bin, 2) if bin % 1 else int(bin))))
ax.SetLabelFont(42)
ax.SetLabelSize(0.05)
ax.Draw()
c_true.Update()

h_cosmics = []
draw_cosmic = False
if draw_cosmic:
    for i in range(len(variables)):
        # if i == reco_energy:
        h_ext = histograms_bnbext[i].Clone()
        h_intime = histograms_mc[i].GetHists()[1].Clone()
        if i == reco_energy:
            fix_binning(h_ext)
        h_ext.Scale(post_scaling)

        if h_intime.Integral() > 0:
            c_cosmic = ROOT.TCanvas("c{}_canvas".format(i), "", 900, 44, 700, 645)
            h_cosmics.append(h_ext)
            h_cosmics.append(h_intime)
            draw_top()

            h_intime.Draw("hist")
            h_ext.Draw("ep same")
            set_axis(h_intime)
            legend_cosmic.Draw()
            if i == reco_energy:
                print(h_ext.Integral() / h_intime.Integral(), i)

            c_cosmic.cd()
            draw_ratio(h_ext, h_intime)
            c_cosmic.cd()

            c_cosmic.Update()
            c_cosmic.SaveAs("plots/%s_cosmic.pdf" % h_intime.GetName())
            canvases_cosmic.append(c_cosmic)

input()
