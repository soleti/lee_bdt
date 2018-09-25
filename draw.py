#!/usr/local/bin/python3

import pickle
import math
import numpy as np
from array import array
# from root_numpy import hist2array
import ROOT
import os.path

from bdt_common import variables, spectators, bins, bins2, total_data_bnb_pot, labels
from bdt_common import description, total_pot, fix_binning, sigma_calc_matrix

ROOT.gStyle.SetOptStat(0)


SYS_ERR = 0.1

DRAW_POT = True
DRAW_SUBTRACTION = False
DRAW_LEE = False
DRAW_DATA = True
DRAW_COSMIC = False
DRAW_SYS = True
OBJECTS = []
VARIABLES = variables + spectators
RECO_ENERGY = list(dict(VARIABLES)).index("reco_energy")
TRUE_ENERGY = list(dict(VARIABLES)).index("nu_E")
# RECO_ENERGY = TRUE_ENERGY
histograms_bnb = []
histograms_bnbext = []
histograms_mc = []
histograms_lee = []

histograms = [histograms_bnb, histograms_mc, histograms_lee]


if DRAW_DATA:
    POST_SCALING = 1
else:
    POST_SCALING = 1#6.6e20 / total_data_bnb_pot

total_pot *= POST_SCALING


def set_axis(histogram, y_max=0):
    histogram.GetYaxis().SetTitleSize(0.06)
    histogram.GetYaxis().SetTitleOffset(0.75)
    histogram.SetMinimum(0.01)
    if y_max != 0:
        histogram.SetMaximum(y_max)
    else:
        histogram.SetMaximum(histogram.GetMaximum() * 1.3)

def draw_top():
    pad_top = ROOT.TPad("pad_top", "", 0, 0.3, 1, 1)
    pad_top.SetBottomMargin(0)
    pad_top.Range(-22.21825, -2.003018, 202.5403, 2073.676)
    pad_top.SetTopMargin(0.2411347)
    pad_top.Draw()
    pad_top.cd()
    OBJECTS.append(pad_top)


def draw_ratio(num, den, den_sys=None):
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

    OBJECTS.append(h_ratio)
    OBJECTS.append(h_ratio_sys)

    h_ratio.SetName("h_ratio%i" % i)
    h_ratio_sys.SetName("h_ratio_sys%i" % i)

    h_ratio.GetYaxis().SetRangeUser(0.01, 1.99)
    h_ratio.Divide(den)
    if DRAW_SYS:
        h_ratio_sys.Divide(den_sys)

    h_ratio.GetXaxis().SetLabelFont(42)
    h_ratio.GetXaxis().SetLabelSize(0.1)
    h_ratio.GetXaxis().SetTitleSize(0.13)
    h_ratio.GetXaxis().SetTitleOffset(0.91)
    h_ratio.GetYaxis().SetTitle("Ratio")
    h_ratio.GetYaxis().SetNdivisions(509)
    h_ratio.GetYaxis().SetLabelFont(42)
    h_ratio.GetYaxis().SetLabelSize(0.1)
    h_ratio.GetYaxis().SetTitleSize(0.14)
    h_ratio.GetYaxis().SetTitleOffset(0.3)

    h_ratio.Draw("hist")
    # h_ratio.Fit("pol1")
    h_ratio.SetMarkerSize(0)
    h_ratio_2 = h_ratio.Clone()
    h_ratio_2.SetFillStyle(3002)
    h_ratio_2.SetFillColor(1)
    h_ratio_2.SetName("h_ratio_2%i" % i)
    OBJECTS.append(h_ratio_2)
    h_ratio_2.Draw("e2 same")
    if DRAW_SYS:
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
samples = ["bnb", "mc", "lee"]

for name, var in VARIABLES:

    for histos, s in zip(histograms, samples):
        f = ROOT.TFile("plots/h_%s_%s.root" % (name, s))
        h = f.Get("h_%s" % name)
        OBJECTS.append(f)
        histos.append(h)

    legend = ROOT.TLegend(0.085, 0.7747, 0.891117, 0.97926, "", "brNDC")
    legend.SetTextSize(16)
    legend.SetTextFont(63)
    legend.SetHeader("MicroBooNE Preliminary %.1e POT" % total_pot)
    legend.SetTextFont(43)
    legend.SetNColumns(2)

    l_err = ROOT.TLegend(0.4713467, 0.5990783, 0.7421203,
                         0.6935484, "", "brNDC")
    if not DRAW_DATA:
        l_err = ROOT.TLegend(0.4756, 0.6403, 0.7464, 0.7338)

    l_err.SetTextSize(18)
    l_err.SetTextFont(43)

    l_errs.append(l_err)
    legends.append(legend)

for i, h in enumerate(histograms_mc):

    if DRAW_DATA:
        legends[i].AddEntry(histograms_bnb[i],
                            "Data beam-on: {:.0f} entries"
                            .format(histograms_bnb[i].Integral()
                                    * POST_SCALING), "lep")

    for d, histo in zip(description, h.GetHists()):
        if histo.Integral():

            n_events = histo.Integral() * POST_SCALING
            legends[i].AddEntry(
                histo,
                "{}: {:.1f} entries".format(d, n_events), "f")

    if DRAW_LEE:
        histograms_lee[i].SetLineWidth(0)
        histograms_lee[i].SetFillColor(ROOT.kGreen - 10)
        if i != RECO_ENERGY:
            legends[i].AddEntry(histograms_lee[i],
                                "Low-energy excess: {:.1f} entries"
                                .format(histograms_lee[i].Integral()
                                        * POST_SCALING), "f")


for h in histograms_mc:
    h.GetHists()[4].SetFillColor(ROOT.TColor.GetColor("#fca4d0"))
    h.GetHists()[0].SetFillStyle(3004)


# if DRAW_SUBTRACTION:
#     for h_bnb, h_ext, leg in zip(histograms_bnb, histograms_bnbext, legends):
#         for i in range(1, h_bnb.GetNbinsX() + 1):
#             bnb = h_bnb.GetBinContent(i)
#             ext = h_ext.GetBinContent(i)

#             h_bnb.SetBinContent(i, bnb - ext)
#             if bnb > 0:
#                 err1 = h_bnb.GetBinError(i)
#                 err2 = h_ext.GetBinError(i)
#                 h_bnb.SetBinError(i, math.sqrt(err1**2 + err2**2))

#         leg.AddEntry(
#             h_bnb,
#             "Data BNB - BNB EXT: {:.0f} events".format(h_bnb.Integral() * POST_SCALING),
#             "lep")


h_lee = ROOT.TH1F("h_lee", "", len(bins) - 1, bins)

scaling = [6.0920944819073988, 3.6447414342239273, 3.2123920194399913,
           2.6504659907742409, 3.2558450032216988, 2.5826310533377432,
           2, 1, 1, 1]

h_lee.SetFillColor(ROOT.kGreen - 10)
# h_lee.SetFillStyle(3001)
h_lee.SetLineWidth(0)
for i, scale in enumerate(scaling):
    if scaling[i] - 1 > 0:
        h_lee.SetBinContent(i + 1, histograms_mc[RECO_ENERGY].GetHists()[-1].
                            GetBinContent(i + 1) * (scale - 1))

# histograms_lee[RECO_ENERGY] = h_lee

if DRAW_LEE:
    legends[RECO_ENERGY].AddEntry(histograms_lee[RECO_ENERGY],
                                  "Low-energy excess: {:.1f} entries"
                                  .format(histograms_lee[RECO_ENERGY].Integral() * POST_SCALING),
                                  "f")

# c_energy = ROOT.TCanvas("c_energy")
# OBJECTS.append(c_energy)
# signal_spectrum = histograms_mc[RECO_ENERGY].GetHists()[0].Clone()
# signal_spectrum_nuecc = histograms_mc[RECO_ENERGY].GetHists()[1].Clone()
# spectrum_stack = ROOT.THStack("spectrum_stack", ";E_{deposited} [GeV]; N. Entries / 0.05 GeV")

# signal_spectrum.Scale(POST_SCALING)
# signal_spectrum_nuecc.Scale(POST_SCALING)
# legend_spectrum = ROOT.TLegend(0.6733, 0.92, 0.913, 0.9747, "", "brNDC")

# # signal_spectrum.Add(signal_spectrum_nuecc)
# legend_spectrum.AddEntry(
#     signal_spectrum, "Beam intrinsic #nu_{e}: %.1f" % signal_spectrum.Integral(), "f")
# signal_spectrum_integral = signal_spectrum.Integral()
# fix_binning(signal_spectrum)
# fix_binning(signal_spectrum_nuecc)

# spectrum_stack.Add(signal_spectrum_nuecc)
# spectrum_stack.Add(signal_spectrum)

# spectrum_stack.Draw("hist")
# #signal_spectrum.Draw("hist")
# legend_spectrum.Draw()

# pt = ROOT.TPaveText(0.098, 0.905, 0.576, 0.989, "ndc")
# pt.AddText("MicroBooNE Preliminary %.1e POT" % total_pot)
# f_spectrum = ROOT.TFile("plots/f_spectrum.root", "RECREATE")
# signal_spectrum.Write()
# f_spectrum.Close()
# pt.SetFillColor(0)
# pt.SetBorderSize(0)
# pt.SetShadowColor(0)
# pt.Draw()
# OBJECTS.append(signal_spectrum)
# OBJECTS.append(spectrum_stack)
# OBJECTS.append(pt)
# c_energy.Update()

sigma = array("f", [])
sigma5 = array("f", [])
sigma20 = array("f", [])
sigma_err = array("f", [])
sigma5_err = array("f", [])
sigma20_err = array("f", [])
pots_err = array("f", [])
pots = array("f", [])
step = 100

# for i in range(len(VARIABLES)):
for i in range(RECO_ENERGY, RECO_ENERGY+1):
    if histograms_mc[i].GetHists()[0].Integral() > -10:

        c = ROOT.TCanvas("c%i" % i, "", 900, 44, 700, 645)

        histograms_bnb[i].SetLineColor(1)
        histograms_bnb[i].SetMarkerStyle(20)

        h_mc_err = histograms_mc[i].GetHists()[0].Clone()
        h_mc_err.SetName("h_mc_err%i" % i)
        h_mc_err.Reset()

        h_mc_err_nobinning = h_mc_err.Clone()
        h_mc_err_nobinning.SetName("h_mc_err_nobin%i" % i)

        h_sig = histograms_mc[i].GetHists()[-1].Clone()

        for j in range(histograms_mc[i].GetNhists()):
            histograms_mc[i].GetHists()[j].Scale(POST_SCALING)
            h_mc_err_nobinning.Add(histograms_mc[i].GetHists()[j])

            if i == RECO_ENERGY:
                fix_binning(histograms_mc[i].GetHists()[j])

            h_mc_err.Add(histograms_mc[i].GetHists()[j])

        h_bkg = h_mc_err_nobinning.Clone()
        h_bkg.Add(h_sig, -1)


        if i == RECO_ENERGY:
            reco_err = h_mc_err
            print(h_sig.Integral(), h_bkg.Integral()+h_sig.Integral())
            print("Purity", h_sig.Integral()/(h_bkg.Integral() + h_sig.Integral()))

        if DRAW_LEE:

            histograms_lee[i].Scale(POST_SCALING)
            if i == RECO_ENERGY:
                if DRAW_POT:
                    sig_err = []
                    bkg_err = []
                    sig = []
                    bkg = []

                    for i_bin in range(1, histograms_lee[i].GetNbinsX()+1):
                        sig.append(histograms_lee[i].GetBinContent(i_bin))
                        bkg.append(h_mc_err_nobinning.GetBinContent(i_bin))
                        sig_err.append(histograms_lee[i].GetBinError(i_bin))
                        bkg_err.append(h_mc_err_nobinning.GetBinError(i_bin))

                    bkg[0] = 0.01
                    print(sig, bkg)
                    sig = np.array(sig)
                    bkg = np.array(bkg)
                    sig_err = np.array(sig_err)
                    bkg_err = np.array(bkg_err)
                    print("Significance: ", sigma_calc_matrix(sig, bkg, 30.41))

                    # for i_pot in range(step):
                    #     pot = total_pot + (1.32e21 - total_pot)/step * i_pot
                    #     scale = pot/total_pot
                    #     pots_err.append(0)
                    #     significance = sigma_calc_matrix(sig, bkg, scale)
                    #     significance20 = sigma_calc_matrix(sig, bkg, scale, 0.2)
                    #     significance5 = sigma_calc_matrix(sig, bkg, scale, 0.05)
                    #     # sigma.append(significance[0])
                    #     # sigma_err.append(significance[1])
                    #     # sigma20.append(significance20[0])
                    #     # sigma20_err.append(significance[1])
                    #     # sigma5.append(significance5[0])
                    #     # sigma5_err.append(significance5[1])
                    #     pots.append(pot)

                fix_binning(histograms_lee[i])
            # histograms_mc[i].RecursiveRemove(histograms_mc[i].GetHists()[0])
            # histograms_mc[i].RecursiveRemove(histograms_mc[i].GetHists()[0])
            # histograms_mc[i].RecursiveRemove(histograms_mc[i].GetHists()[0])
            # histograms_mc[i].RecursiveRemove(histograms_mc[i].GetHists()[0])
            # histograms_mc[i].RecursiveRemove(histograms_mc[i].GetHists()[0])
            # histograms_mc[i].RecursiveRemove(histograms_mc[i].GetHists()[0])
            # histograms_mc[i].RecursiveRemove(histograms_mc[i].GetHists()[0])
            # histograms_mc[i].RecursiveRemove(histograms_mc[i].GetHists()[0])
            # histograms_mc[i].RecursiveRemove(histograms_mc[i].GetHists()[0])
            h_mc_err.Add(histograms_lee[i])
            histograms_mc[i].Add(histograms_lee[i])

        if DRAW_DATA:
            draw_top()
        else:
            c.SetTopMargin(0.2274194)
        if DRAW_SUBTRACTION:
            histograms_bnb[i].Add(histograms_mc[i].GetHists()[0], -1)
            h_mc_err.Add(histograms_mc[i].GetHists()[0], -1)
            histograms_mc[i].RecursiveRemove(histograms_mc[i].GetHists()[0])
        histograms_mc[i].Draw("hist")

        if DRAW_DATA:
            set_axis(histograms_mc[i], histograms_bnb[i].GetMaximum() * 1.35)
        else:
            set_axis(histograms_mc[i])

        if i == RECO_ENERGY:
            new_max = h_mc_err.GetMaximum() + histograms_lee[i].GetMaximum()
            # histograms_mc[i].SetMaximum(new_max * 1.3)

        h_mc_err.SetFillStyle(3002)
        h_mc_err.SetFillColor(1)

        max_hist = histograms_bnb[i].GetMaximum() * 1.35

        if VARIABLES[i][0] == "dedx":
            l1 = ROOT.TLine(1, 0, 1, max_hist)
            l1.SetLineStyle(2)
            l1.SetLineWidth(3)
            l2 = ROOT.TLine(3.2, 0, 3.2, max_hist)
            l2.SetLineStyle(2)
            l2.SetLineWidth(3)
            l1.Draw()
            l2.Draw()
            OBJECTS.append(l1)
            OBJECTS.append(l2)

        if VARIABLES[i][0] == "track_distance":
            l1 = ROOT.TLine(5, 0, 5, max_hist)
            l1.SetLineStyle(2)
            l1.SetLineWidth(3)
            l1.Draw()
            OBJECTS.append(l1)

        if VARIABLES[i][0] == "shower_distance":
            l1 = ROOT.TLine(5, 0, 5, max_hist)
            l1.SetLineStyle(2)
            l1.SetLineWidth(3)
            l1.Draw()
            OBJECTS.append(l1)

        if VARIABLES[i][0] == "shower_open_angle":
            l1 = ROOT.TLine(1, 0, 1, max_hist)
            l1.SetLineStyle(2)
            l1.SetLineWidth(3)
            l2 = ROOT.TLine(20, 0, 20, max_hist)
            l2.SetLineStyle(2)
            l2.SetLineWidth(3)
            l1.Draw()
            l2.Draw()
            OBJECTS.append(l1)
            OBJECTS.append(l2)

        if VARIABLES[i][0] == "track_shower_angle":
            l1 = ROOT.TLine(-0.9, 0, -0.9, max_hist)
            l1.SetLineStyle(2)
            l1.SetLineWidth(3)
            l1.Draw()
            OBJECTS.append(l1)

        if VARIABLES[i][0] == "total_hits_y":
            l1 = ROOT.TLine(50, 0, 50, max_hist)
            l1.SetLineStyle(2)
            l1.SetLineWidth(3)
            l1.Draw()
            OBJECTS.append(l1)

        if VARIABLES[i][0] == "shower_energy":
            l1 = ROOT.TLine(0.050, 0, 0.050, max_hist)
            l1.SetLineStyle(2)
            l1.SetLineWidth(3)
            l1.Draw()
            OBJECTS.append(l1)

        if VARIABLES[i][0] == "track_length":
            l1 = ROOT.TLine(80, 0, 80, max_hist)
            l1.SetLineStyle(2)
            l1.SetLineWidth(3)
            l1.Draw()
            OBJECTS.append(l1)

        if VARIABLES[i][0] == "dqdx_bdt":
            l1 = ROOT.TLine(0.1, 0, 0.1, max_hist)
            l1.SetLineStyle(2)
            l1.SetLineWidth(3)
            l1.Draw()
            OBJECTS.append(l1)

        h_mc_err.Draw("e2 same")
        h_mc_err_clone = h_mc_err.Clone()
        h_mc_err_clone.SetLineWidth(2)
        h_mc_err_clone.SetFillStyle(0)
        h_mc_err_clone.Draw("hist same")
        OBJECTS.append(h_mc_err_clone)
        if DRAW_SYS:
            legends[i].AddEntry(h_mc_err, "Stat. #oplus sys. uncertainty", "lef")
        else:
            legends[i].AddEntry(h_mc_err, "Stat. uncertainty", "lf")

        if DRAW_SYS:
            h_mc_err_sys = h_mc_err.Clone()
            fname_flux = "plots/sys/h_%s_flux_sys.root" % VARIABLES[i][0]
            fname_genie = "plots/sys/h_%s_genie_sys.root" % VARIABLES[i][0]
            OBJECTS.append(h_mc_err_sys)
            if DRAW_SYS and os.path.isfile(fname_flux) and os.path.isfile(fname_genie):
                f_flux = ROOT.TFile(fname_flux)
                h_flux = f_flux.Get("h_%s_cv" % VARIABLES[i][0])
                f_genie = ROOT.TFile(fname_genie)
                h_genie = f_genie.Get("h_%s_cv" % VARIABLES[i][0])
                print(h_flux, h_genie)
                OBJECTS.append(h_flux)
                OBJECTS.append(h_genie)

                for k in range(1, h_mc_err_sys.GetNbinsX() + 1):
                    stat_err = h_mc_err_sys.GetBinError(k)
                    flux_err = h_flux.GetBinError(k)
                    genie_err = h_genie.GetBinError(k)
                    h_mc_err_sys.SetBinError(
                        k, math.sqrt(flux_err**2 + genie_err**2 + stat_err**2))
                f_genie.Close()
                f_flux.Close()
            p_chi2 = ROOT.TPaveText(0.679, 0.604, 0.872, 0.715, "NDC")
            p_chi2.AddText("#chi^{2} prob. = %.2f" %
                           histograms_bnb[i].Chi2Test(h_mc_err_sys, "UW"))
            ks = histograms_bnb[i].KolmogorovTest(h_mc_err_sys)
            p_chi2.AddText("K-S prob. = %.2f" % ks)
            p_chi2.SetTextAlign(11)
            if i == RECO_ENERGY:
                reco_chi2 = p_chi2
                h_reco_sys = h_mc_err_sys
            h_mc_err_sys.Draw("e1p same")
            OBJECTS.append(h_mc_err_sys)


        if DRAW_DATA:

            if i == RECO_ENERGY:
                ratio = histograms_bnb[i].Integral() / (h_mc_err_nobinning.Integral())
                fix_binning(histograms_bnb[i])

            histograms_bnb[i].Draw("e1p same")
            p_chi2.Draw("same")
            p_chi2.SetFillStyle(0)
            p_chi2.SetBorderSize(0)
            OBJECTS.append(p_chi2)

        OBJECTS.append(h_mc_err)

        legends[i].Draw("same")
        c.cd()
        # histograms_mc[i].GetHists()[-1].Draw("hist")

        if DRAW_DATA and h_mc_err_sys.Integral() > 0:
            if i == RECO_ENERGY:
                print("Data/(MC+EXT) ratio: ", ratio)

            draw_ratio(histograms_bnb[i], h_mc_err, h_mc_err_sys)
        c.Update()
        c.SaveAs("plots/%s.pdf" % histograms_bnb[i].GetName())

        OBJECTS.append(c)


# *******************************
# START FIXED BIN WIDTH PLOT
# *******************************

h_true_e = ROOT.THStack(
    "h_true_e", ";E_{deposited} [GeV]; N. Entries / 0.05 GeV")

h_mc_fixed = ROOT.TH1F("h_mc_fixed", "", len(bins) - 1, bins2)

for j in range(histograms_mc[RECO_ENERGY].GetNhists()):

    h_clone = histograms_mc[RECO_ENERGY].GetHists()[j].Clone()
    h_fixed = ROOT.TH1F("h_fixed%i" % j, "", len(bins) - 1, bins2)

    for i in range(1, h_clone.GetNbinsX() + 1):
        h_fixed.SetBinContent(i, h_clone.GetBinContent(i))
        h_fixed.SetBinError(i, h_clone.GetBinError(i))
        h_mc_fixed.SetBinContent(i, h_mc_fixed.GetBinContent(i) + h_clone.GetBinContent(i))


    if j != 0 and not DRAW_SUBTRACTION:
        h_fixed.SetLineWidth(0)

    h_fixed.SetLineColor(1)
    h_fixed.SetFillColor(h_clone.GetFillColor())
    h_fixed.SetFillStyle(h_clone.GetFillStyle())

    h_true_e.Add(h_fixed)


h_mc_fixed_sys = h_mc_fixed.Clone()
for i in range(1, h_mc_fixed.GetNbinsX() + 1):
    h_mc_fixed.SetBinError(i, reco_err.GetBinError(i))
    h_mc_fixed_sys.SetBinError(i, h_reco_sys.GetBinError(i))

h_clone_data = histograms_bnb[RECO_ENERGY].Clone()
h_fixed_data = ROOT.TH1F("h_fixed_data", labels["reco_energy"], len(bins) - 1, bins2)

for i in range(1, h_clone_data.GetNbinsX() + 1):
    h_fixed_data.SetBinContent(i, h_clone_data.GetBinContent(i))
    h_fixed_data.SetBinError(i, h_clone_data.GetBinError(i))
    # h_fixed_data.GetXaxis().SetBinLabel(i, "")

h_fixed_data.SetLineColor(ROOT.kBlack)
h_fixed_data.SetMarkerStyle(20)
c_fixed = ROOT.TCanvas("c_true", "", 900, 44, 700, 645)

# c_fixed.SetTopMargin(0.245)
# c_fixed.SetBottomMargin(0.105)
draw_top()
h_true_e.Draw("hist")
h_true_e.GetYaxis().SetTitleOffset(0.95)
h_mc_fixed.SetLineWidth(2)
h_mc_fixed.SetLineColor(1)
h_mc_fixed.Draw("e2 same")
h_mc_fixed_sys.SetLineColor(1)
h_mc_fixed_sys.Draw("e1 same")
h_mc_fixed_clone = h_mc_fixed.Clone()
h_mc_fixed_clone.Draw("hist same")

h_mc_fixed.SetFillColor(ROOT.kBlack)
h_mc_fixed.SetFillStyle(3002)
h_true_e.SetMaximum(max(h_true_e.GetMaximum(), h_fixed_data.GetMaximum()) * 1.3)
h_true_e.SetMinimum(0.001)

if DRAW_DATA:
    h_fixed_data.Draw("e1p same")
    # reco_chi2.Draw()
    p_datamc = ROOT.TPaveText(0.563, 0.58947, 0.882, 0.703, "NDC")
    p_datamc.SetFillStyle(0)
    p_datamc.SetShadowColor(0)
    p_datamc.SetBorderSize(0)
    p_datamc.AddText("Data / (MC + EXT) = %.2f" % ratio)
    # p_datamc.Draw()

legends[RECO_ENERGY].Draw()

c_fixed.cd()
draw_ratio(h_fixed_data, h_mc_fixed, h_mc_fixed_sys)


p_15 = ROOT.TPaveText(0.77,0.17,0.84,0.282, "NDC")
p_15.AddText("1.5")
p_15.SetFillStyle(0)
p_15.SetShadowColor(0)
p_15.SetBorderSize(0)
p_15.SetTextFont(42)
p_15.Draw()

p_white = ROOT.TPaveText(0.83,0.17,0.87,0.27, "NDC")
p_white.AddText("  ")
p_white.SetFillColor(ROOT.kWhite)
p_white.SetShadowColor(0)
p_white.SetBorderSize(0)
p_white.SetTextFont(42)
p_white.SetTextSize(14)
p_white.Draw()


p_3 = ROOT.TPaveText(0.865,0.17,0.935,0.282, "NDC")
p_3.AddText("3")
p_3.SetFillStyle(0)
p_3.SetShadowColor(0)
p_3.SetBorderSize(0)
p_3.SetTextFont(42)
p_3.Draw()
# h_true_e.GetXaxis().ChangeLabel(2,  -1, -1, -1, -1, -1, "aaaaa")
# ax = ROOT.TGaxis(0, 0, 1.7, 0, 0, 1.7, 515, "")
# for i, i_bin in enumerate(bins):
#     ax.ChangeLabel(i + 1, -1, -1, -1, -1, -1,
#                    "{0}".format(str(round(i_bin, 3) if i_bin % 1 else int(i_bin))))
# ax.SetLabelFont(42)
# ax.SetLabelSize(0.04)
# ax.Draw()
# pt3.Draw()
c_fixed.Update()
c_fixed.SaveAs("plots/h_fixed_energy.pdf")
# *******************************
# END FIXED BIN WIDTH PLOT
# *******************************

if DRAW_LEE and DRAW_POT and False:
    c_pot = ROOT.TCanvas("c_pot")
    g_pot = ROOT.TGraphErrors(step, pots, sigma, pots_err, sigma_err)
    g_pot5 = ROOT.TGraphErrors(step, pots, sigma5, pots_err, sigma_err)
    g_pot20 = ROOT.TGraphErrors(step, pots, sigma20, pots_err, sigma_err)

    g_pot.Draw("AC4")
    ROOT.TGaxis().SetMaxDigits(2)
    run1_label = ROOT.TPaveText(0.4484, 0.00442, 0.5257, 0.067, "NDC")
    run1_label.AddText("6.6e20")
    run1_label.SetLineWidth(0)
    run1_label.SetFillStyle(0)
    run1_label.Draw()
    run1_label.SetBorderSize(0)
    run1_label.SetTextFont(62)
    run1_line = ROOT.TLine(6.6e20, 0, 6.6e20, 8)
    run1_line.SetLineWidth(2)
    run1_line.Draw()

    five_line = ROOT.TLine(total_pot, 5, 1.32e21, 5)
    three_line = ROOT.TLine(total_pot, 3, 1.32e21, 3)
    five_line.SetLineWidth(2)
    three_line.SetLineWidth(2)
    five_line.Draw()
    three_line.Draw()
    three_line.SetLineStyle(7)
    five_line.SetLineStyle(7)

    g_pot.GetXaxis().SetTitle("POT")
    g_pot.GetYaxis().SetTitle("Significance [#sigma]")
    g_pot.GetYaxis().SetRangeUser(0, 8)
    g_pot.GetXaxis().SetRangeUser(total_pot, 1.32e21)

    g_pot.SetLineColor(ROOT.kOrange + 1)
    g_pot.SetFillColor(ROOT.kOrange + 1)
    g_pot.SetFillStyle(3002)

    g_pot.SetLineWidth(2)
    g_pot5.SetLineColor(ROOT.kRed + 1)
    g_pot5.SetFillColor(ROOT.kRed + 1)
    g_pot5.SetFillStyle(3002)

    g_pot5.SetLineWidth(2)
    g_pot5.Draw("C4 SAME")
    g_pot20.SetLineColor(ROOT.kGreen + 1)
    g_pot20.SetFillColor(ROOT.kGreen + 1)
    g_pot20.SetFillStyle(3002)

    g_pot20.SetLineWidth(2)
    g_pot20.Draw("C4 SAME")
    l_pot = ROOT.TLegend(0.55, 0.69, 0.86, 0.86)
    l_pot.AddEntry(g_pot, "Stat. only", "lf")
    l_pot.AddEntry(g_pot5, "Stat. only + 5% sys.", "lf")
    l_pot.AddEntry(g_pot20, "Stat. only + 20% sys", "lf")

    l_pot.Draw()
    c_pot.Update()
    c_pot.SaveAs("plots/pots.pdf")

DRAW_NORMALIZED = False
if DRAW_NORMALIZED:

    for i in range(len(VARIABLES)):
        c_norm = ROOT.TCanvas("c%i_norm" % i, "", 900, 44, 700, 645)


        h_signal = histograms_lee[i]
        print(h_signal.GetName())
        if h_signal.Integral() <= 0:
            continue
        h_signal.SetLineColor(ROOT.TColor.GetColor("#1e7a2e"))
        h_signal.SetLineWidth(3)
        h_signal.SetFillStyle(0)
        h_signal.Scale(1 / h_signal.Integral())


        h_neutrino_bkg = histograms_mc[i].GetHists()[1].Clone()
        h_neutrino_bkg.Add(histograms_mc[i].GetHists()[2])
        h_neutrino_bkg.Add(histograms_mc[i].GetHists()[3])
        h_neutrino_bkg.Add(histograms_mc[i].GetHists()[4])
        if h_neutrino_bkg.Integral() > 0:
            h_neutrino_bkg.Scale(1 / h_neutrino_bkg.Integral())
        h_neutrino_bkg.SetFillStyle(0)
        h_neutrino_bkg.SetLineWidth(3)
        h_neutrino_bkg.SetLineColor(ROOT.kBlue + 1)

        h_cosmic_bkg = histograms_mc[i].GetHists()[5].Clone()
        # h_cosmic_bkg.Add(histograms_mc[i].GetHists()[6])
        # h_cosmic_bkg.Add(histograms_mc[i].GetHists()[7])
        if h_cosmic_bkg.Integral() > 0:
            h_cosmic_bkg.Scale(1 / h_cosmic_bkg.Integral())
        h_cosmic_bkg.SetFillStyle(0)
        h_cosmic_bkg.SetLineWidth(3)
        h_cosmic_bkg.SetLineColor(ROOT.kRed + 1)
        max_hist = max(h_signal.GetMaximum(),
                       h_neutrino_bkg.GetMaximum(),
                       h_cosmic_bkg.GetMaximum()) * 1.2

        h_signal.SetMaximum(max_hist)
        h_signal.SetMinimum(0.00001)

        l_norm = ROOT.TLegend(0.59, 0.73, 0.86, 0.82)
        l_norm.AddEntry(h_signal, "#nu_{e} CC0#pi-Np", "f")
        l_norm.AddEntry(h_neutrino_bkg, "Neutrino background", "f")
        l_norm.AddEntry(h_cosmic_bkg, "Cosmic background", "f")
        h_signal.Draw("hist")
        h_neutrino_bkg.Draw("hist same")
        h_cosmic_bkg.Draw("hist same")
        l_norm.Draw()

    # track_distance = chain.track_distance < 5
    # track_res = chain.track_res_std < 2
    # shower_angle = not (5 < chain.shower_angle < 45)
    # dqdx = chain.dqdx_bdt_max > 0.1  # and chain.dqdx_bdt > 0.1
    # corrected_energy = 0.55 < chain.total_shower_energy / \

        if VARIABLES[i][0] == "dedx":
            l1 = ROOT.TLine(1, 0, 1, max_hist)
            l1.SetLineStyle(2)
            l1.SetLineWidth(3)
            l2 = ROOT.TLine(3.2, 0, 3.2, max_hist)
            l2.SetLineStyle(2)
            l2.SetLineWidth(3)
            l1.Draw()
            l2.Draw()
            OBJECTS.append(l1)
            OBJECTS.append(l2)

        if VARIABLES[i][0] == "track_distance":
            l1 = ROOT.TLine(5, 0, 5, max_hist)
            l1.SetLineStyle(2)
            l1.SetLineWidth(3)
            l1.Draw()
            OBJECTS.append(l1)

        if VARIABLES[i][0] == "shower_distance":
            l1 = ROOT.TLine(5, 0, 5, max_hist)
            l1.SetLineStyle(2)
            l1.SetLineWidth(3)
            l1.Draw()
            OBJECTS.append(l1)

        if VARIABLES[i][0] == "shower_open_angle":
            l1 = ROOT.TLine(1, 0, 1, max_hist)
            l1.SetLineStyle(2)
            l1.SetLineWidth(3)
            l2 = ROOT.TLine(20, 0, 20, max_hist)
            l2.SetLineStyle(2)
            l2.SetLineWidth(3)
            l1.Draw()
            l2.Draw()
            OBJECTS.append(l1)
            OBJECTS.append(l2)

        if VARIABLES[i][0] == "track_shower_angle":
            l1 = ROOT.TLine(-0.9, 0, -0.9, max_hist)
            l1.SetLineStyle(2)
            l1.SetLineWidth(3)
            l1.Draw()
            OBJECTS.append(l1)

        if VARIABLES[i][0] == "total_hits_y":
            l1 = ROOT.TLine(50, 0, 50, max_hist)
            l1.SetLineStyle(2)
            l1.SetLineWidth(3)
            l1.Draw()
            OBJECTS.append(l1)

        if VARIABLES[i][0] == "shower_energy":
            l1 = ROOT.TLine(0.050, 0, 0.050, max_hist)
            l1.SetLineStyle(2)
            l1.SetLineWidth(3)
            l1.Draw()
            OBJECTS.append(l1)

        if VARIABLES[i][0] == "track_length":
            l1 = ROOT.TLine(80, 0, 80, max_hist)
            l1.SetLineStyle(2)
            l1.SetLineWidth(3)
            l1.Draw()
            OBJECTS.append(l1)


        if VARIABLES[i][0] == "dqdx_bdt":
            l1 = ROOT.TLine(0.1, 0, 0.1, max_hist)
            l1.SetLineStyle(2)
            l1.SetLineWidth(3)
            l1.Draw()
            OBJECTS.append(l1)

        pt = ROOT.TPaveText(0.136, 0.906, 0.501, 0.953, "ndc")
        pt.AddText("MicroBooNE Preliminary")
        pt.SetFillColor(0)
        pt.SetBorderSize(0)
        pt.SetShadowColor(0)
        pt.Draw()

        pt2 = ROOT.TPaveText(0.59, 0.83, 0.80, 0.87, "ndc")
        pt2.AddText("Area normalized")
        pt2.SetFillColor(0)
        pt2.SetBorderSize(0)
        pt2.SetShadowColor(0)
        pt2.Draw()

        c_norm.SetLeftMargin(0.14)
        c_norm.SetRightMargin(0.06)
        c_norm.Update()
        c_norm.SaveAs("plots/integral/%s_norm.pdf" % histograms_bnb[i].GetName())

        OBJECTS.append(pt)
        OBJECTS.append(pt2)

        OBJECTS.append(l_norm)
        OBJECTS.append(c_norm)
        OBJECTS.append(h_signal)
        OBJECTS.append(h_neutrino_bkg)
        OBJECTS.append(h_cosmic_bkg)

if DRAW_COSMIC:
    legend_cosmic = ROOT.TLegend(0.099, 0.909, 0.900, 0.987, "", "brNDC")
    legend_cosmic.AddEntry(histograms_bnbext[RECO_ENERGY], "Data EXT: {:.0f} entries".format(
        histograms_bnbext[RECO_ENERGY].Integral() * POST_SCALING), "lep")
    legend_cosmic.AddEntry(histograms_mc[RECO_ENERGY].GetHists()[1],
                           "CORSIKA in-time Monte Carlo: integral normalized", "f")
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
