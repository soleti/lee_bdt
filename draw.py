#!/usr/local/bin/python3

import math
import pickle
import numpy as np
import ROOT
import os.path
import sys

from bdt_common import variables, spectators, bins, bins2, total_data_bnb_pot, labels
from bdt_common import description, total_pot, fix_binning, sigma_calc_matrix, BDT, MANUAL
from bdt_common import save_histo_sbnfit, binning, inv_pdgs, pdg_colors

ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetNumberContours(99)

# if len(sys.argv) > 1:
#     plot = sys.argv[1]
# else:
#     plot = ""


DRAW_NORMALIZED = False
DRAW_POT = True
DRAW_SUBTRACTION = False
DRAW_LEE = True
DRAW_DATA = True
DRAW_SYS = True
OBJECTS = []
VARIABLES = variables + spectators
RECO_ENERGY = list(dict(VARIABLES)).index("reco_energy")
POST_SCALING = 1

def plot_pdg(variable, mode):
    if variable in ("nu_E", "E_dep"):
        DRAW_DATA = False
    else:
        DRAW_DATA = True
    if mode == "selection":
        folder = ""
    elif mode == "bdt":
        folder = "_bdt"
    elif mode == "cuts":
        folder = "_cuts"
    elif mode == "numu":
        folder = "_numu"
    elif mode == "nc":
        folder = "_nc"

    v = variable
    h_pdgs = pickle.load(open("plots%s/pdg_plots.p" % folder, "rb"))
    h_stack = ROOT.THStack("h_stack_%s" % v, labels[v])
    OBJECTS.append(h_stack)
    h_tot_mc = ROOT.TH1F("h_tot_mc_%s" % v,
                            labels[v],
                            binning[v][0],
                            binning[v][1],
                            binning[v][2])
    h_tot_mc.SetTitle("Stat. uncertainty")

    l_pdg = ROOT.TLegend(0.09, 0.7747, 0.904, 0.984)
    l_pdg.SetNColumns(3)
    l_pdg.SetTextSize(16)
    l_pdg.SetTextFont(63)
    l_pdg.SetHeader("MicroBooNE Preliminary %.1e POT" %
                    total_data_bnb_pot)
    l_pdg.SetTextFont(43)
    for pdg in h_pdgs[v]:
        h_pdgs[v][pdg].SetLineWidth(0)
        h_pdgs[v][pdg].SetMarkerStyle(0)
        h_pdgs[v][pdg].SetMarkerSize(0)
        h_pdgs[v][pdg].SetFillColor(
            ROOT.TColor.GetColor(pdg_colors[pdg]))
        if pdg != 2147483648:
            h_stack.Add(h_pdgs[v][pdg])
            h_tot_mc.Add(h_pdgs[v][pdg])

    for pdg in h_pdgs[v]:
        integral = h_pdgs[v][pdg].Integral()
        if integral > 0 and pdg != 2147483648 and h_tot_mc.Integral() > 0:
            l_pdg.AddEntry(h_pdgs[v][pdg],
                            "%s: %.1f%%" % (
                inv_pdgs[pdg], (integral / h_tot_mc.Integral()) * 100),
                "f")

    h_tot_mc_clone = h_tot_mc.Clone()
    h_tot_mc.SetFillStyle(3002)
    h_tot_mc.SetFillColor(1)
    c = ROOT.TCanvas("c_%s" % v, v, 900, 44, 700, 645)

    if DRAW_DATA and not mode in ("bdt", "cuts", "numu", "nc"):
        draw_top(OBJECTS, margin_top=0.25)
    else:
        c.SetTopMargin(0.23)

    h_stack.Draw("hist")
    if DRAW_SYS:
        h_mc_err_sys = h_tot_mc.Clone()
        fname_flux = "plots%s/sys/h_%s_flux_sys.root" % (folder, v)
        fname_genie = "plots%s/sys/h_%s_genie_sys.root" % (folder, v)
        if mode == "selection":
            fname_detsys = "plots%s/sys/h_%s_det_sys.root" % (folder, v)

        if v == "track_energy_length":
            fname_flux = fname_flux.replace("track", "total_track")
            fname_genie = fname_genie.replace("track", "total_track")
        OBJECTS.append(h_mc_err_sys)
        if DRAW_SYS:
            if os.path.isfile(fname_flux) and os.path.isfile(fname_genie):
                var_name = v
                if v == "track_energy_length":
                    var_name = v.replace("track", "total_track")
                f_flux = ROOT.TFile(fname_flux)
                h_flux = f_flux.Get("h_%s_cv" % var_name)
                f_genie = ROOT.TFile(fname_genie)
                h_genie = f_genie.Get("h_%s_cv" % var_name)
                if mode == "selection" and os.path.isfile(fname_detsys):
                    f_detsys = ROOT.TFile(fname_detsys)
                    h_detsys = f_detsys.Get("h_%s_cv" % var_name)
                    OBJECTS.append(h_detsys)
                else:
                    print("Det sys files not available")

                OBJECTS.append(h_flux)
                OBJECTS.append(h_genie)

                for k in range(1, h_mc_err_sys.GetNbinsX() + 1):
                    stat_err = h_mc_err_sys.GetBinError(k)
                    flux_err = h_flux.GetBinError(k)
                    genie_err = h_genie.GetBinError(k)
                    if mode == "selection" and os.path.isfile(fname_detsys):
                        detsys_err = h_detsys.GetBinError(k) * 1.2
                    else:
                        detsys_err = 0

                    h_mc_err_sys.SetBinError(k, math.sqrt(flux_err**2 + genie_err**2 + detsys_err**2 + stat_err**2))
                    # h_mc_err_sys.SetBinError(k, math.sqrt(stat_err**2))
                f_genie.Close()
                f_flux.Close()
            else:
                print("GENIE or FLUX sys files not available")

        h_mc_err_sys.SetLineWidth(0)
        h_mc_err_sys.SetLineColor(1)
        h_mc_err_sys.Draw("e2 same")
        OBJECTS.append(h_mc_err_sys)

    if DRAW_SYS:
        l_pdg.AddEntry(h_mc_err_sys, "Sys. uncertainty", "f")
    if DRAW_DATA:
        f = ROOT.TFile("plots%s/h_%s_bnb.root" % (folder, v))
        h = f.Get("h_%s" % v)
        h.SetLineColor(1)
        h.SetTitle("Data beam-on - beam-off")
        h.SetMarkerStyle(20)
        h_tot_mc.SetLineColor(1)
        h_tot_mc.SetLineWidth(2)
        h.Add(h_pdgs[v][2147483648], -1)
        set_axis(h_stack, max(h_stack.GetHistogram().GetMaximum(), h.GetMaximum()) * 1.3)
        h.Draw("e1p same")
        l_pdg.AddEntry(h, "Data (beam-on - beam-off)", "lep")
    else:
        set_axis(h_stack)

    OBJECTS.append(f)
    OBJECTS.append(h_tot_mc)
    if not DRAW_SYS:
        h_tot_mc.Draw("e2 same")
    h_tot_mc_clone.SetLineWidth(2)
    h_tot_mc_clone.SetLineColor(1)
    h_tot_mc_clone.Draw("hist same")
    l_pdg.Draw()

    if DRAW_SYS and DRAW_DATA:
        chi2 = h.Chi2Test(h_mc_err_sys, "UW")
        print("chi2", chi2)
        ks = h.KolmogorovTest(h_mc_err_sys)
    else:
        chi2 = h.Chi2Test(h_tot_mc, "UW")
        ks = h.KolmogorovTest(h_tot_mc)

    p_test = ROOT.TPaveText(0.656, 0.587, 0.849, 0.695, "NDC")
    p_test.AddText("#chi^{2} prob. = %.2f" % chi2)
    p_test.AddText("K-S prob. = %.2f" % ks)
    p_test.SetFillStyle(0)
    p_test.SetBorderSize(0)
    p_test.SetTextAlign(11)
    p_test.Draw()

    OBJECTS.append(p_test)
    OBJECTS.append(h_tot_mc_clone)
    h_stack.GetYaxis().SetTitleOffset(0.9)
    c.cd()
    if DRAW_DATA and mode not in ("bdt", "cuts", "numu", "nc"):
        if DRAW_SYS:
            draw_ratio(h, h_tot_mc, h_mc_err_sys, OBJECTS)
        else:
            draw_ratio(h, h_tot_mc, h_tot_mc, OBJECTS)

    c.Update()
    c.SaveAs("plots%s/pdg/%s_pdg.pdf" % (folder, h_stack.GetName()))
    OBJECTS.append(c)
    OBJECTS.append(l_pdg)

    integral_err = 0
    integral = 0
    for i_bin in range(1, h_mc_err_sys.GetNbinsX()+1):
        integral += h_mc_err_sys.GetBinContent(i_bin)
        integral_err += h_mc_err_sys.GetBinContent(i_bin) + h_mc_err_sys.GetBinError(i_bin)
    print("Sys uncert. %.1f %%" % ((1 - integral / integral_err) * 100))

    return c, l_pdg


def set_axis(histogram, y_max=0):
    h = histogram.GetHistogram()
    h.SetTitleSize(0.06)
    h.SetTitleOffset(0.75)
    histogram.SetHistogram(h)
    histogram.SetMinimum(0.01)
    if y_max != 0:
        histogram.SetMaximum(y_max)
    else:
        histogram.SetMaximum(histogram.GetMaximum() * 1.3)
    return h


def draw_top(OBJECTS=[], margin_top=0.33):
    pad_top = ROOT.TPad("pad_top", "", 0, 0.3, 1, 1)
    pad_top.SetBottomMargin(0)
    pad_top.Range(-22.21825, -2.003018, 202.5403, 2073.676)
    pad_top.SetTopMargin(margin_top)
    pad_top.Draw()
    pad_top.cd()
    OBJECTS.append(pad_top)


def draw_ratio(num, den, den_sys=None, OBJECTS=[]):
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
    h_ratio_sys = den.Clone()

    OBJECTS.append(h_ratio)
    OBJECTS.append(h_ratio_sys)

    h_ratio.GetYaxis().SetRangeUser(0.5, 1.5)
    h_ratio.Divide(den)
    if den_sys:
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

    h_ratio.Draw("e1")
    if den_sys:
        h_ratio_sys.Draw("e2 same")
    h_ratio_sys.SetMarkerSize(0)

    line = ROOT.TLine(h_ratio.GetXaxis().GetXmin(), 1,
                      h_ratio.GetXaxis().GetXmax(), 1)
    line.SetLineWidth(2)
    line.SetLineStyle(2)
    line.Draw()
    OBJECTS.append(line)


def plot_energy(histograms, reco_err, reco_sys, mode="selection"):
    c_fixed = ROOT.TCanvas("c_energy", "", 900, 44, 700, 645)
    c_fixed.cd()
    if DRAW_DATA and not mode in ("bdt", "cuts", "numu", "nc"):
        draw_top()
    else:
        c_fixed.SetTopMargin(0.23)

    h_true_e = ROOT.THStack(
        "h_true_e", ";E_{deposited} [GeV]; N. Entries / 0.05 GeV")

    h_mc_fixed = ROOT.TH1F("h_mc_fixed", "", len(bins) - 1, bins2)

    for j in range(histograms["mc"].GetNhists()):

        h_clone = histograms["mc"].GetHists()[j].Clone()
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
        h_mc_fixed_sys.SetBinError(i, reco_sys.GetBinError(i))

    h_clone_data = histograms["bnb"].Clone()
    h_fixed_data = ROOT.TH1F("h_fixed_data", labels["reco_energy"], len(bins) - 1, bins2)

    for i in range(1, h_clone_data.GetNbinsX() + 1):
        h_fixed_data.SetBinContent(i, h_clone_data.GetBinContent(i))
        h_fixed_data.SetBinError(i, h_clone_data.GetBinError(i))

    h_fixed_data.SetLineColor(ROOT.kBlack)
    h_fixed_data.SetMarkerStyle(20)

    h_mc_fixed.SetLineWidth(2)
    h_mc_fixed.SetLineColor(1)
    h_mc_fixed_sys.SetLineColor(1)
    h_mc_fixed_sys.SetFillStyle(3002)
    h_mc_fixed_sys.SetFillColor(ROOT.kBlack)

    h_mc_fixed.SetFillColor(ROOT.kBlack)
    h_mc_fixed.SetFillStyle(3002)
    h_true_e.SetMaximum(max(h_true_e.GetMaximum(), h_fixed_data.GetMaximum()) * 1.3)
    h_true_e.SetMinimum(0.001)


    h_true_e.Draw("hist")
    if DRAW_SYS:
        h_mc_fixed_sys.Draw("e2 same")
    else:
        h_mc_fixed.Draw("e2 same")

    h_mc_fixed_clone = h_mc_fixed.Clone()
    h_mc_fixed_clone.SetFillStyle(0)
    h_mc_fixed_clone.Draw("hist same")
    if DRAW_DATA:
        h_fixed_data.Draw("e1p same")

    h_true_e.GetYaxis().SetTitleOffset(0.95)

    OBJECTS.append(h_true_e)
    OBJECTS.append(h_mc_fixed)
    OBJECTS.append(h_mc_fixed_clone)
    OBJECTS.append(h_mc_fixed_sys)
    OBJECTS.append(h_fixed_data)
    c_fixed.cd()
    if DRAW_DATA and mode not in ("bdt", "cuts", "numu", "nc"):
        draw_ratio(h_fixed_data, h_mc_fixed, h_mc_fixed_sys, OBJECTS)
        p_15 = ROOT.TPaveText(0.77, 0.17, 0.84, 0.282, "NDC")
        p_white = ROOT.TPaveText(0.83, 0.17, 0.87, 0.27, "NDC")
        p_3 = ROOT.TPaveText(0.865, 0.17, 0.935, 0.282, "NDC")
    else:
        p_15 = ROOT.TPaveText(0.78, 0.026, 0.83, 0.129, "NDC")
        p_white = ROOT.TPaveText(0.83, 0.06, 0.88, 0.091, "NDC")
        p_3 = ROOT.TPaveText(0.8624642, 0.05645161,
                             0.9326648, 0.1032258, "NDC")
    p_15.AddText("1.5")
    p_15.SetFillStyle(0)
    p_15.SetShadowColor(0)
    p_15.SetBorderSize(0)
    p_15.SetTextFont(42)
    p_15.Draw()

    p_white.AddText("  ")
    p_white.SetFillColor(ROOT.kWhite)
    p_white.SetShadowColor(0)
    p_white.SetBorderSize(0)
    p_white.SetTextFont(42)
    p_white.SetTextSize(14)
    p_white.Draw()

    p_3.AddText("3")
    p_3.SetFillStyle(0)
    p_3.SetShadowColor(0)
    p_3.SetBorderSize(0)
    p_3.SetTextFont(42)
    p_3.Draw()
    c_fixed.cd()
    return c_fixed

def plot_normalized(name):
    histograms = {}
    samples = ["bnb", "mc", "lee"]

    for s in samples:
        f = ROOT.TFile("plots/h_%s_%s.root" % (name, s))
        h = f.Get("h_%s" % name)
        OBJECTS.append(f)
        histograms[s] = h


    c_norm = ROOT.TCanvas("c%s_norm" % name, "", 900, 44, 700, 645)

    h_signal = histograms["mc"].GetHists()[8].Clone()
    h_signal.SetLineColor(ROOT.TColor.GetColor("#62b570"))
    h_signal.SetLineWidth(3)
    h_signal.SetFillStyle(0)
    h_signal.Scale(1 / h_signal.Integral())
    h_neutrino_bkg = histograms["mc"].GetHists()[1].Clone()
    h_neutrino_bkg.Add(histograms["mc"].GetHists()[2])
    h_neutrino_bkg.Add(histograms["mc"].GetHists()[3])
    h_neutrino_bkg.Add(histograms["mc"].GetHists()[4])
    if h_neutrino_bkg.Integral() > 0:
        h_neutrino_bkg.Scale(1 / h_neutrino_bkg.Integral())
    h_neutrino_bkg.SetFillStyle(0)
    h_neutrino_bkg.SetLineWidth(3)
    h_neutrino_bkg.SetLineColor(ROOT.kBlue + 1)

    h_cosmic_bkg = histograms["mc"].GetHists()[0].Clone()
    h_cosmic_bkg.Add(histograms["mc"].GetHists()[5])
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

    l_norm = ROOT.TLegend(0.135, 0.867, 0.891117, 1, "", "brNDC")
    l_norm.SetTextFont(63)
    l_norm.SetHeader("MicroBooNE Preliminary Simulation")
    l_norm.SetTextSize(16)
    l_norm.SetTextFont(43)
    l_norm.SetNColumns(3)
    header = l_norm.GetListOfPrimitives().First()
    header.SetTextSize(20)
    l_norm.AddEntry(h_signal, "#nu_{e} CC0#pi-Np", "f")
    l_norm.AddEntry(h_neutrino_bkg, "Neutrino background", "f")
    l_norm.AddEntry(h_cosmic_bkg, "Cosmic background", "f")
    c_norm.SetTopMargin(0.145)
    h_signal.Draw("hist")
    h_neutrino_bkg.Draw("hist same")
    h_cosmic_bkg.Draw("hist same")
    OBJECTS.append(h_signal)
    OBJECTS.append(h_neutrino_bkg)
    OBJECTS.append(h_cosmic_bkg)

    l_norm.Draw()

    pt2 = ROOT.TPaveText(0.7, 0.73, 0.9, 0.83, "ndc")
    pt2.AddText("Area normalized")
    pt2.SetFillColor(0)
    pt2.SetBorderSize(0)
    pt2.SetShadowColor(0)
    pt2.Draw()
    OBJECTS.append(pt2)

    c_norm.SetLeftMargin(0.14)
    c_norm.SetRightMargin(0.06)
    c_norm.Update()

    return c_norm, l_norm


def plot_variable(name, mode="selection"):
    if name in ("nu_E", "E_dep", "is_signal"):
        DRAW_DATA = False
    else:
        DRAW_DATA = True
    histograms = {}

    if mode == "selection":
        folder = ""
    elif mode == "bdt":
        folder = "_bdt"
    elif mode == "cuts":
        folder = "_cuts"
    elif mode == "numu":
        folder = "_numu"
    elif mode == "nc":
        folder = "_nc"

    samples = ["bnb", "mc", "lee"]

    for s in samples:
        f = ROOT.TFile("plots%s/h_%s_%s.root" % (folder, name, s))
        h = f.Get("h_%s" % name)
        OBJECTS.append(f)
        histograms[s] = h

    legend = ROOT.TLegend(0.085, 0.7847, 0.891117, 0.97926, "", "brNDC")
    legend.SetTextSize(16)
    legend.SetTextFont(63)
    legend.SetHeader("MicroBooNE Preliminary %.1e POT" % total_pot)
    legend.SetTextFont(43)
    legend.SetNColumns(2)

    if DRAW_DATA:
        legend.AddEntry(histograms["bnb"],
                        "Data beam-on: {:.0f} entries"
                        .format(histograms["bnb"].Integral() * POST_SCALING), "lep")

    for d, histo in zip(description, histograms["mc"].GetHists()):
        if histo.Integral():
            n_events = histo.Integral() * POST_SCALING
            legend.AddEntry(histo, "{}: {:.1f} entries".format(d, n_events), "f")

    if DRAW_LEE:
        histograms["lee"].Scale(1.5)
        histograms["lee"].SetLineWidth(0)
        histograms["lee"].SetFillColor(ROOT.kGreen - 10)
        legend.AddEntry(histograms["lee"],
                        "Low-energy excess: {:.1f} entries"
                        .format(histograms["lee"].Integral() * POST_SCALING),
                        "f")

    histograms["mc"].GetHists()[4].SetFillColor(ROOT.TColor.GetColor("#fca4d0"))
    histograms["mc"].GetHists()[0].SetFillStyle(3004)

    c = ROOT.TCanvas("c%s" % name, name, 900, 44, 700, 645)

    histograms["bnb"].SetLineColor(1)
    histograms["bnb"].SetMarkerStyle(20)

    h_mc_err = histograms["mc"].GetHists()[0].Clone()
    h_mc_err.SetName("h_mc_err%s" % name)
    h_mc_err.Reset()

    h_mc_err_nobinning = h_mc_err.Clone()
    h_mc_err_nobinning.SetName("h_mc_err_nobin%s" % name)

    h_sig = histograms["mc"].GetHists()[-1].Clone()

    for j in range(histograms["mc"].GetNhists()):
        histograms["mc"].GetHists()[j].Scale(POST_SCALING)
        h_mc_err_nobinning.Add(histograms["mc"].GetHists()[j])

        if name == "reco_energy":
            fix_binning(histograms["mc"].GetHists()[j])

        h_mc_err.Add(histograms["mc"].GetHists()[j])

    if name == "reco_energy" and mode not in ("numu", "selection", "nc"):
        print("Purity", h_sig.Integral()/h_mc_err_nobinning.Integral())

    if DRAW_LEE:
        if name == "reco_energy" and mode not in ("numu", "selection", "nc"):
            sig = []
            bkg = []
            nu_e = []
            for i_bin in range(1, histograms["lee"].GetNbinsX()+1):
                sig.append(histograms["lee"].GetBinContent(i_bin))
                bkg.append(max(0.001, h_mc_err_nobinning.GetBinContent(i_bin)))
                nu_e.append(h_sig.GetBinContent(i_bin))

            sig = np.array(sig)
            bkg = np.array(bkg)
            nu_e = np.array(nu_e)
            print("Significance @ 4.4e19 POT stat: ",
                sigma_calc_matrix(sig, bkg, 1, False))
            print("Significance @ 1.3e21 POT stat: ",
                sigma_calc_matrix(sig, bkg, 30.4, False))
            if DRAW_SYS:
                print("Significance @ 1.3e21 POT sys: ", sigma_calc_matrix(sig, bkg, 30.4, True, mode))
            fix_binning(histograms["lee"])

        h_mc_err.Add(histograms["lee"])
        histograms["mc"].Add(histograms["lee"])

    if DRAW_DATA and mode not in ("bdt", "cuts", "numu", "nc"):
        draw_top()
    else:
        c.SetTopMargin(0.2274194)
    if DRAW_SUBTRACTION:
        histograms["bnb"].Add(histograms["mc"].GetHists()[0], -1)
        h_mc_err.Add(histograms["mc"].GetHists()[0], -1)
        histograms["mc"].RecursiveRemove(histograms["mc"].GetHists()[0])
    histograms["mc"].Draw("hist")
    histograms["mc"].GetHistogram().GetXaxis().SetTitleOffset(0.8)

    if DRAW_DATA:
        set_axis(histograms["mc"], max(histograms["mc"].GetHistogram().GetMaximum(), histograms["bnb"].GetMaximum()) * 1.35)
    else:
        set_axis(histograms["mc"])


    h_mc_err.SetFillStyle(3002)
    h_mc_err.SetFillColor(1)

    h_mc_err_clone = h_mc_err.Clone()
    h_mc_err_clone.SetLineWidth(2)
    h_mc_err_clone.SetFillStyle(0)
    h_mc_err_clone.Draw("hist same")
    OBJECTS.append(h_mc_err_clone)
    h_mc_err_sys = h_mc_err.Clone()
    h_mc_err_sys_nobinning = h_mc_err_nobinning.Clone()
    if DRAW_SYS:
        h_mc_err_sys.SetLineWidth(0)
        legend.AddEntry(h_mc_err_sys, "Sys. uncertainty", "f")
    else:
        legend.AddEntry(h_mc_err, "Stat. uncertainty", "lf")
    if DRAW_SYS:
        fname_flux = "plots%s/sys/h_%s_flux_sys.root" % (folder, name)
        fname_genie = "plots%s/sys/h_%s_genie_sys.root" % (folder, name)
        fname_detsys = "plots%s/sys/h_%s_det_sys.root" % (folder, name)

        if "_before" in name:
            fname_flux = fname_flux.replace("_before", "")
            fname_genie = fname_genie.replace("_before", "")
        if name == "track_energy_length":
            fname_flux = fname_flux.replace("track", "total_track")
            fname_genie = fname_genie.replace("track", "total_track")

        OBJECTS.append(h_mc_err_sys)
        print(os.path.isfile(fname_flux), os.path.isfile(fname_genie))
        if os.path.isfile(fname_flux) and os.path.isfile(fname_genie):
            f_flux = ROOT.TFile(fname_flux)
            if name == "reco_energy":
                fixed = "_fixed"
            else:
                fixed = ""
            var_name = name.replace("_before", "")
            if name == "track_energy_length":
                var_name = var_name.replace("track", "total_track")

            h_flux = f_flux.Get("h_%s_cv%s" % (var_name, fixed))
            f_genie = ROOT.TFile(fname_genie)
            h_genie = f_genie.Get("h_%s_cv%s" % (var_name, fixed))
            h_flux_nobinning = h_flux.Clone()
            h_genie_nobinning = h_genie.Clone()

            if mode == "selection":
                f_detsys = ROOT.TFile(fname_detsys)
                h_detsys = f_detsys.Get("h_%s_cv" % var_name)
                h_detsys_nobinning = h_detsys.Clone()

            if name == "reco_energy":
                fix_binning(h_flux)
                fix_binning(h_genie)
                if mode == "selection":
                    fix_binning(h_detsys)
            OBJECTS.append(h_flux)
            OBJECTS.append(h_genie)
            if mode == "selection":
                OBJECTS.append(h_detsys)

            for k in range(1, h_mc_err_sys.GetNbinsX() + 1):
                stat_err = h_mc_err_sys.GetBinError(k)
                flux_err = h_flux.GetBinError(k)
                genie_err = h_genie.GetBinError(k)
                if mode == "selection":
                    detsys_err = h_detsys.GetBinError(k) * 1.07
                    detsys_err_nobinning = h_detsys_nobinning.GetBinError(k) * 1.07
                else:
                    detsys_err, detsys_err_nobinning = h_mc_err_sys.GetBinContent(k) * 0.2, h_mc_err_sys_nobinning.GetBinContent(k) * 0.2

                stat_err_nobinning = h_mc_err_sys_nobinning.GetBinError(k)
                flux_err_nobinning = h_flux_nobinning.GetBinError(k)
                genie_err_nobinning = h_genie_nobinning.GetBinError(k)

                h_mc_err_sys.SetBinError(k, math.sqrt(flux_err**2 + genie_err**2 + stat_err**2 + detsys_err**2))
                h_mc_err_sys_nobinning.SetBinError(k, math.sqrt(stat_err_nobinning**2+flux_err_nobinning**2+genie_err_nobinning**2+detsys_err_nobinning**2))
                # h_mc_err_sys.SetBinError(k, math.sqrt(detsys_err**2))
                # h_mc_err_sys_nobinning.SetBinError(k, math.sqrt(detsys_err_nobinning**2))

            f_genie.Close()
            f_flux.Close()
            if mode == "selection":
                f_detsys.Close()
        else:
            print("GENIE or FLUX or DET sys files not available")
        h_mc_err_sys.Draw("e2 same")
    else:
        h_mc_err.Draw("e2 same")

    OBJECTS.append(h_mc_err_sys)
    p_chi2 = ROOT.TPaveText(0.679, 0.654, 0.872, 0.735, "NDC")

    if DRAW_DATA:
        if name == "reco_energy":
            ratio = histograms["bnb"].Integral() / h_mc_err_nobinning.Integral()
            if DRAW_SYS:
                chi2 = histograms["bnb"].Chi2Test(h_mc_err_sys_nobinning, "UW")
                ks = histograms["bnb"].KolmogorovTest(h_mc_err_sys_nobinning)
            else:
                chi2  = histograms["bnb"].Chi2Test(h_mc_err_nobinning, "UW")
                ks = histograms["bnb"].KolmogorovTest(h_mc_err_nobinning)
            fix_binning(histograms["bnb"])
        else:
            if DRAW_SYS:
                chi2 = histograms["bnb"].Chi2Test(h_mc_err_sys, "UW")
                ks = histograms["bnb"].KolmogorovTest(h_mc_err_sys)
            else:
                chi2 = histograms["bnb"].Chi2Test(h_mc_err, "UW")
                ks = histograms["bnb"].KolmogorovTest(h_mc_err)
        p_chi2.AddText("#chi^{2} prob. = %.2f" % chi2)
        p_chi2.AddText("K-S prob. = %.2f" % ks)
        p_chi2.SetTextAlign(11)



    if name == "reco_energy":
        reco_chi2 = p_chi2
        h_reco_sys = h_mc_err_sys
    if DRAW_DATA:
        histograms["bnb"].Draw("e1p same")

    OBJECTS.append(h_mc_err)

    c.cd()

    if DRAW_DATA and h_mc_err_sys.Integral() > 0:
        if name == "RECO_ENERGY":
            print("Data/(MC+EXT) ratio: ", ratio)
        if mode not in ("bdt", "cuts", "numu", "nc"):
            draw_ratio(histograms["bnb"], h_mc_err, h_mc_err_sys, OBJECTS)

    c.Update()
    histograms["mc"].GetYaxis().SetTitleOffset(0.9)
    c.cd()
    if name == "reco_energy":
        c = plot_energy(histograms, h_mc_err, h_mc_err_sys, mode)
    if DRAW_DATA:
        p_chi2.Draw("same")
        p_chi2.SetFillStyle(0)
        p_chi2.SetBorderSize(0)
        OBJECTS.append(p_chi2)
    else:
        c.SetBottomMargin(0.15)
    OBJECTS.append(c)
    legend.Draw("same")
    c.SaveAs("plots%s/%s.pdf" % (folder, histograms["bnb"].GetName()))
    integral_err = 0
    integral = 0
    for i_bin in range(1, h_mc_err_nobinning.GetNbinsX()+1):
        integral += h_mc_err_nobinning.GetBinContent(i_bin)
        integral_err += h_mc_err_nobinning.GetBinContent(i_bin) + h_mc_err_sys_nobinning.GetBinError(i_bin)
    print(integral, integral_err, "Sys uncert. %.1f %%" % (((integral_err / integral) - 1) * 100))
    return c, legend


# histograms_bnb = []
# histograms_bnbext = []
# histograms_mc = []
# histograms_lee = []

# histograms = [histograms_bnb, histograms_mc, histograms_lee]


# if DRAW_DATA:
#     POST_SCALING = 1
# else:
#     POST_SCALING = 1#6.6e20 / total_data_bnb_pot

# total_pot *= POST_SCALING

# legends = []
# l_errs = []
# samples = ["bnb", "mc", "lee"]

# for name, var in VARIABLES:

#     for histos, s in zip(histograms, samples):
#         f = ROOT.TFile("plots%s/h_%s_%s.root" % (folder, name, s))
#         h = f.Get("h_%s" % name)
#         OBJECTS.append(f)
#         histos.append(h)

#     legend = ROOT.TLegend(0.085, 0.7747, 0.891117, 0.97926, "", "brNDC")
#     legend.SetTextSize(16)
#     legend.SetTextFont(63)
#     legend.SetHeader("MicroBooNE Preliminary %.1e POT" % total_pot)
#     legend.SetTextFont(43)
#     legend.SetNColumns(2)

#     l_err = ROOT.TLegend(0.4713467, 0.5990783, 0.7421203,
#                          0.6935484, "", "brNDC")
#     if not DRAW_DATA:
#         l_err = ROOT.TLegend(0.4756, 0.6403, 0.7464, 0.7338)

#     l_err.SetTextSize(18)
#     l_err.SetTextFont(43)

#     l_errs.append(l_err)
#     legends.append(legend)

# for i, h in enumerate(histograms_mc):

#     if DRAW_DATA:
#         legends[i].AddEntry(histograms_bnb[i],
#                             "Data beam-on: {:.0f} entries"
#                             .format(histograms_bnb[i].Integral()
#                                     * POST_SCALING), "lep")

#     for d, histo in zip(description, h.GetHists()):
#         if histo.Integral():

#             n_events = histo.Integral() * POST_SCALING
#             legends[i].AddEntry(histo,
#                                 "{}: {:.1f} entries".format(d, n_events), "f")

#     if DRAW_LEE:
#         histograms_lee[i].Scale(1.5)
#         histograms_lee[i].SetLineWidth(0)
#         histograms_lee[i].SetFillColor(ROOT.kGreen - 10)
#         if i != RECO_ENERGY:
#             legends[i].AddEntry(histograms_lee[i],
#                                 "Low-energy excess: {:.1f} entries"
#                                 .format(histograms_lee[i].Integral()
#                                         * POST_SCALING), "f")


# for h in histograms_mc:
#     h.GetHists()[4].SetFillColor(ROOT.TColor.GetColor("#fca4d0"))
#     h.GetHists()[0].SetFillStyle(3004)


# if DRAW_LEE:
#     legends[RECO_ENERGY].AddEntry(histograms_lee[RECO_ENERGY],
#                                   "Low-energy excess: {:.1f} entries"
#                                   .format(histograms_lee[RECO_ENERGY].Integral() * POST_SCALING),
#                                   "f")

# pots_err = array("f", [])
# pots = array("f", [])
# step = 100

# if plot:
#     to_plot = range(VARIABLE, VARIABLE+1)
# else:
#     to_plot = range(len(VARIABLES))

# for i in to_plot:
#     # if ("track" in VARIABLES[i][0] or "shower" in VARIABLES[i][0]) and not "total" in VARIABLES[i][0] and not "n_" in VARIABLES[i][0]:
#     #     continue
#     if histograms_mc[i].GetHists()[0].Integral() > -10:

#         c = ROOT.TCanvas("c%i" % i, "", 900, 44, 700, 645)

#         histograms_bnb[i].SetLineColor(1)
#         histograms_bnb[i].SetMarkerStyle(20)

#         h_mc_err = histograms_mc[i].GetHists()[0].Clone()
#         h_mc_err.SetName("h_mc_err%i" % i)
#         h_mc_err.Reset()

#         h_mc_err_nobinning = h_mc_err.Clone()
#         h_mc_err_nobinning.SetName("h_mc_err_nobin%i" % i)

#         h_sig = histograms_mc[i].GetHists()[-1].Clone()

#         for j in range(histograms_mc[i].GetNhists()):
#             histograms_mc[i].GetHists()[j].Scale(POST_SCALING)
#             # if j == 2 or j == 1:
#             #     histograms_mc[i].GetHists()[j].Scale(0.5)
#             #     histograms_mc[i].GetHists()[j].Scale(0.5)

#             h_mc_err_nobinning.Add(histograms_mc[i].GetHists()[j])

#             if i == RECO_ENERGY:
#                 fix_binning(histograms_mc[i].GetHists()[j])

#             h_mc_err.Add(histograms_mc[i].GetHists()[j])

#         if i == RECO_ENERGY:
#             reco_err = h_mc_err
#             print("Purity", h_sig.Integral()/h_mc_err_nobinning.Integral())

#         if DRAW_LEE:
#             if i == RECO_ENERGY:
#                 sig = []
#                 bkg = []
#                 nu_e = []
#                 for i_bin in range(1, histograms_lee[i].GetNbinsX()+1):
#                     sig.append(histograms_lee[i].GetBinContent(i_bin))
#                     bkg.append(max(0.001, h_mc_err_nobinning.GetBinContent(i_bin)))
#                     nu_e.append(h_sig.GetBinContent(i_bin))

#                 sig = array("f", sig)
#                 bkg = array("f", bkg)
#                 nu_e = array("f",nu_e)
#                 print("Significance @ 4.4e19 POT stat: ", sigma_calc_matrix(sig, bkg, 1, False))
#                 print("Significance @ 1.3e21 POT stat: ", sigma_calc_matrix(sig, bkg, 30.4, False))
#                 # print("Significance @ 1.3e21 POT sys: ", sigma_calc_matrix(sig, bkg, 30.4, True))

#                 # h_eff_red = ROOT.TH2F("h_eff_red", ";Efficiency increase;Background rejection increase", 40, 1, 5, 40, 1, 5)
#                 # for eff in range(10, 50):
#                 #     for red in range(10, 50):
#                 #         bkg_new = (bkg-nu_e)/(red/10) + nu_e*eff/10
#                 #         sig_new = sig*eff/10
#                 #         if eff == 23 and red == 18:
#                 #             print("Sigma", sigma_calc_matrix(sig_new, bkg_new, 30.4, True))
#                 #         h_eff_red.Fill(eff/10+0.001, red/10+0.001, sigma_calc_matrix(sig_new, bkg_new, 30.4, True))
#                 # OBJECTS.append(h_eff_red)

#                 # c_eff_red = ROOT.TCanvas("c_eff_red")
#                 # h_eff_red.Draw("colz")
#                 # h_eff_red.GetZaxis().SetRangeUser(0.9,7)
#                 # c_eff_red.Update()
#                 # c_eff_red.SaveAs("plots/h_2d.pdf")
#                 # OBJECTS.append(c_eff_red)

#                 f = ROOT.TFile("root_files/sbnfit.root", "RECREATE")
#                 save_histo_sbnfit(histograms_lee[i], "nu_uBooNE_nue_leesignal", 30.4).Write()
#                 save_histo_sbnfit(h_mc_err_nobinning, "nu_uBooNE_nue_intrinsic", 30.4).Write()
#                 f.Close()
#                 fix_binning(histograms_lee[i])
#             # histograms_mc[i].RecursiveRemove(histograms_mc[i].GetHists()[0])
#             # histograms_mc[i].RecursiveRemove(histograms_mc[i].GetHists()[0])
#             # histograms_mc[i].RecursiveRemove(histograms_mc[i].GetHists()[0])
#             # histograms_mc[i].RecursiveRemove(histograms_mc[i].GetHists()[0])
#             # histograms_mc[i].RecursiveRemove(histograms_mc[i].GetHists()[0])
#             # histograms_mc[i].RecursiveRemove(histograms_mc[i].GetHists()[0])
#             # histograms_mc[i].RecursiveRemove(histograms_mc[i].GetHists()[0])
#             # histograms_mc[i].RecursiveRemove(histograms_mc[i].GetHists()[0])
#             # histograms_mc[i].RecursiveRemove(histograms_mc[i].GetHists()[0])
#             h_mc_err.Add(histograms_lee[i])
#             histograms_mc[i].Add(histograms_lee[i])

#         if DRAW_DATA and not (BDT or MANUAL):
#             draw_top()
#         else:
#             c.SetTopMargin(0.2274194)
#         if DRAW_SUBTRACTION:
#             histograms_bnb[i].Add(histograms_mc[i].GetHists()[0], -1)
#             h_mc_err.Add(histograms_mc[i].GetHists()[0], -1)
#             histograms_mc[i].RecursiveRemove(histograms_mc[i].GetHists()[0])
#         histograms_mc[i].Draw("hist")

#         if DRAW_DATA:
#             set_axis(histograms_mc[i], histograms_bnb[i].GetMaximum() * 1.35)
#         else:
#             set_axis(histograms_mc[i])

#         if i == RECO_ENERGY:
#             new_max = h_mc_err.GetMaximum() + histograms_lee[i].GetMaximum()
#             # histograms_mc[i].SetMaximum(new_max * 1.3)

#         h_mc_err.SetFillStyle(3002)
#         h_mc_err.SetFillColor(1)

#         max_hist = histograms_bnb[i].GetMaximum() * 1.35

#         h_mc_err_clone = h_mc_err.Clone()
#         h_mc_err_clone.SetLineWidth(2)
#         h_mc_err_clone.SetFillStyle(0)
#         h_mc_err_clone.Draw("hist same")
#         OBJECTS.append(h_mc_err_clone)
#         h_mc_err_sys = h_mc_err.Clone()
#         if DRAW_SYS:
#             h_mc_err_sys.SetLineWidth(0)
#             legends[i].AddEntry(h_mc_err_sys, "Sys. uncertainty", "f")
#         else:
#             legends[i].AddEntry(h_mc_err, "Stat. uncertainty", "lf")
#         if DRAW_SYS:
#             fname_flux = "plots/sys/h_%s_flux_sys.root" % VARIABLES[i][0]
#             fname_genie = "plots/sys/h_%s_genie_sys.root" % VARIABLES[i][0]
#             if "_before" in VARIABLES[i][0]:
#                 fname_flux = fname_flux.replace("_before", "")
#                 fname_genie = fname_genie.replace("_before", "")
#             if VARIABLES[i][0] == "track_energy_length":
#                 fname_flux = fname_flux.replace("track", "total_track")
#                 fname_genie = fname_genie.replace("track", "total_track")

#             OBJECTS.append(h_mc_err_sys)
#             if os.path.isfile(fname_flux) and os.path.isfile(fname_genie):
#                 f_flux = ROOT.TFile(fname_flux)
#                 if i == RECO_ENERGY:
#                     fixed = "_fixed"
#                 else:
#                     fixed = ""
#                 var_name = VARIABLES[i][0].replace("_before", "")
#                 if VARIABLES[i][0] == "track_energy_length":
#                     var_name = var_name.replace("track", "total_track")

#                 h_flux = f_flux.Get("h_%s_cv%s" % (var_name, fixed))
#                 f_genie = ROOT.TFile(fname_genie)
#                 h_genie = f_genie.Get("h_%s_cv%s" % (var_name, fixed))
#                 if i == RECO_ENERGY:
#                     fix_binning(h_flux)
#                     fix_binning(h_genie)
#                 OBJECTS.append(h_flux)
#                 OBJECTS.append(h_genie)

#                 for k in range(1, h_mc_err_sys.GetNbinsX() + 1):
#                     stat_err = h_mc_err_sys.GetBinError(k)
#                     flux_err = h_flux.GetBinError(k)
#                     genie_err = h_genie.GetBinError(k)
#                     h_mc_err_sys.SetBinError(
#                         k, math.sqrt(flux_err**2 + genie_err**2 + stat_err**2))
#                 f_genie.Close()
#                 f_flux.Close()
#             else:
#                 print("GENIE or FLUX sys files not available")
#             h_mc_err_sys.Draw("e2 same")
#         else:
#             h_mc_err.Draw("e2 same")

#         OBJECTS.append(h_mc_err_sys)

#         if DRAW_DATA:
#             if i == RECO_ENERGY:
#                 ratio = histograms_bnb[i].Integral() / (h_mc_err_nobinning.Integral())

#                 fix_binning(histograms_bnb[i])

#         p_chi2 = ROOT.TPaveText(0.679, 0.604, 0.872, 0.715, "NDC")
#         p_chi2.AddText("#chi^{2} prob. = %.2f" %
#                        histograms_bnb[i].Chi2Test(h_mc_err_sys, "UW"))
#         ks = histograms_bnb[i].KolmogorovTest(h_mc_err_sys)
#         p_chi2.AddText("K-S prob. = %.2f" % ks)
#         p_chi2.SetTextAlign(11)
#         if i == RECO_ENERGY:
#             reco_chi2 = p_chi2
#             h_reco_sys = h_mc_err_sys
#         if DRAW_DATA:
#             histograms_bnb[i].Draw("e1p same")
#             p_chi2.Draw("same")
#             p_chi2.SetFillStyle(0)
#             p_chi2.SetBorderSize(0)
#             OBJECTS.append(p_chi2)

#         OBJECTS.append(h_mc_err)

#         legends[i].Draw("same")
#         c.cd()
#         # histograms_mc[i].GetHists()[-1].Draw("hist")

#         if DRAW_DATA and h_mc_err_sys.Integral() > 0:
#             # if i == RECO_ENERGY:
#             #     print("Data/(MC+EXT) ratio: ", ratio)
#             if not (BDT or MANUAL):
#                 draw_ratio(histograms_bnb[i], h_mc_err, h_mc_err_sys, OBJECTS)
#         c.Update()
#         c.SaveAs("plots%s/%s.pdf" % (folder, histograms_bnb[i].GetName()))

#         OBJECTS.append(c)


# # *******************************
# # START FIXED BIN WIDTH PLOT
# # *******************************
# if not plot or plot == "reco_energy":
#     h_true_e = ROOT.THStack(
#         "h_true_e", ";E_{deposited} [GeV]; N. Entries / 0.05 GeV")

#     h_mc_fixed = ROOT.TH1F("h_mc_fixed", "", len(bins) - 1, bins2)

#     for j in range(histograms_mc[RECO_ENERGY].GetNhists()):

#         h_clone = histograms_mc[RECO_ENERGY].GetHists()[j].Clone()
#         h_fixed = ROOT.TH1F("h_fixed%i" % j, "", len(bins) - 1, bins2)

#         for i in range(1, h_clone.GetNbinsX() + 1):
#             h_fixed.SetBinContent(i, h_clone.GetBinContent(i))
#             h_fixed.SetBinError(i, h_clone.GetBinError(i))
#             h_mc_fixed.SetBinContent(i, h_mc_fixed.GetBinContent(i) + h_clone.GetBinContent(i))

#         if j != 0 and not DRAW_SUBTRACTION:
#             h_fixed.SetLineWidth(0)

#         h_fixed.SetLineColor(1)
#         h_fixed.SetFillColor(h_clone.GetFillColor())
#         h_fixed.SetFillStyle(h_clone.GetFillStyle())

#         h_true_e.Add(h_fixed)


#     h_mc_fixed_sys = h_mc_fixed.Clone()
#     for i in range(1, h_mc_fixed.GetNbinsX() + 1):
#         h_mc_fixed.SetBinError(i, reco_err.GetBinError(i))
#         h_mc_fixed_sys.SetBinError(i, h_reco_sys.GetBinError(i))

#     h_clone_data = histograms_bnb[RECO_ENERGY].Clone()
#     h_fixed_data = ROOT.TH1F("h_fixed_data", labels["reco_energy"], len(bins) - 1, bins2)

#     for i in range(1, h_clone_data.GetNbinsX() + 1):
#         h_fixed_data.SetBinContent(i, h_clone_data.GetBinContent(i))
#         h_fixed_data.SetBinError(i, h_clone_data.GetBinError(i))
#         # h_fixed_data.GetXaxis().SetBinLabel(i, "")

#     h_fixed_data.SetLineColor(ROOT.kBlack)
#     h_fixed_data.SetMarkerStyle(20)
#     c_fixed = ROOT.TCanvas("c_true", "", 900, 44, 700, 645)

#     # c_fixed.SetTopMargin(0.245)
#     # c_fixed.SetBottomMargin(0.105)
#     if DRAW_DATA and not (BDT or MANUAL):
#         draw_top()
#     else:
#         c_fixed.SetTopMargin(0.23)

#     h_true_e.Draw("hist")
#     h_true_e.GetYaxis().SetTitleOffset(0.95)
#     h_mc_fixed.SetLineWidth(2)
#     h_mc_fixed.SetLineColor(1)
#     h_mc_fixed_sys.SetLineColor(1)
#     h_mc_fixed_sys.SetFillStyle(3002)
#     h_mc_fixed_sys.SetFillColor(ROOT.kBlack)
#     if DRAW_SYS:
#         h_mc_fixed_sys.Draw("e2 same")
#     else:
#         h_mc_fixed.Draw("e2 same")

#     h_mc_fixed_clone = h_mc_fixed.Clone()
#     h_mc_fixed_clone.Draw("hist same")

#     h_mc_fixed.SetFillColor(ROOT.kBlack)
#     h_mc_fixed.SetFillStyle(3002)
#     h_true_e.SetMaximum(max(h_true_e.GetMaximum(), h_fixed_data.GetMaximum()) * 1.3)
#     h_true_e.SetMinimum(0.001)

#     if DRAW_DATA:
#         h_fixed_data.Draw("e1p same")
#         reco_chi2.Draw()
#         p_datamc = ROOT.TPaveText(0.563, 0.58947, 0.882, 0.703, "NDC")
#         p_datamc.SetFillStyle(0)
#         p_datamc.SetShadowColor(0)
#         p_datamc.SetBorderSize(0)
#         # p_datamc.AddText("Data / (MC + EXT) = %.2f" % ratio)
#         # p_datamc.Draw()

#     legends[RECO_ENERGY].Draw()

#     c_fixed.cd()

#     if DRAW_DATA and not (BDT or MANUAL):
#         draw_ratio(h_fixed_data, h_mc_fixed, h_mc_fixed_sys, OBJECTS)
#         p_15 = ROOT.TPaveText(0.77,0.17,0.84,0.282, "NDC")
#         p_white = ROOT.TPaveText(0.83,0.17,0.87,0.27, "NDC")
#         p_3 = ROOT.TPaveText(0.865,0.17,0.935,0.282, "NDC")
#     else:
#         p_15 = ROOT.TPaveText(0.78, 0.026, 0.83, 0.129, "NDC")
#         p_white = ROOT.TPaveText(0.83,0.06,0.88,0.091, "NDC")
#         p_3 = ROOT.TPaveText(0.8624642, 0.05645161, 0.9326648, 0.1032258, "NDC")

#     p_15.AddText("1.5")
#     p_15.SetFillStyle(0)
#     p_15.SetShadowColor(0)
#     p_15.SetBorderSize(0)
#     p_15.SetTextFont(42)
#     p_15.Draw()

#     p_white.AddText("  ")
#     p_white.SetFillColor(ROOT.kWhite)
#     p_white.SetShadowColor(0)
#     p_white.SetBorderSize(0)
#     p_white.SetTextFont(42)
#     p_white.SetTextSize(14)
#     p_white.Draw()

#     p_3.AddText("3")
#     p_3.SetFillStyle(0)
#     p_3.SetShadowColor(0)
#     p_3.SetBorderSize(0)
#     p_3.SetTextFont(42)
#     p_3.Draw()
#     # h_true_e.GetXaxis().ChangeLabel(2,  -1, -1, -1, -1, -1, "aaaaa")
#     # ax = ROOT.TGaxis(0, 0, 1.7, 0, 0, 1.7, 515, "")
#     # for i, i_bin in enumerate(bins):
#     #     ax.ChangeLabel(i + 1, -1, -1, -1, -1, -1,
#     #                    "{0}".format(str(round(i_bin, 3) if i_bin % 1 else int(i_bin))))
#     # ax.SetLabelFont(42)
#     # ax.SetLabelSize(0.04)
#     # ax.Draw()
#     # pt3.Draw()
#     c_fixed.Update()
#     c_fixed.SaveAs("plots%s/h_fixed_energy.pdf" % folder)
#     # *******************************
#     # END FIXED BIN WIDTH PLOT
#     # *******************************

# if DRAW_NORMALIZED:

#     if plot:
#         to_plot = range(VARIABLE, VARIABLE+1)
#     else:
#         to_plot = range(len(VARIABLES))
#     for i in to_plot:
#         c_norm = ROOT.TCanvas("c%i_norm" % i, "", 900, 44, 700, 645)

#         h_signal = histograms_lee[i]
#         if h_signal.Integral() <= 0:
#             continue
#         h_signal.SetLineColor(ROOT.TColor.GetColor("#1e7a2e"))
#         h_signal.SetLineWidth(3)
#         h_signal.SetFillStyle(0)
#         h_signal.Scale(1 / h_signal.Integral())


#         h_neutrino_bkg = histograms_mc[i].GetHists()[1].Clone()
#         h_neutrino_bkg.Add(histograms_mc[i].GetHists()[2])
#         h_neutrino_bkg.Add(histograms_mc[i].GetHists()[3])
#         h_neutrino_bkg.Add(histograms_mc[i].GetHists()[4])
#         if h_neutrino_bkg.Integral() > 0:
#             h_neutrino_bkg.Scale(1 / h_neutrino_bkg.Integral())
#         h_neutrino_bkg.SetFillStyle(0)
#         h_neutrino_bkg.SetLineWidth(3)
#         h_neutrino_bkg.SetLineColor(ROOT.kBlue + 1)

#         h_cosmic_bkg = histograms_mc[i].GetHists()[5].Clone()
#         # h_cosmic_bkg.Add(histograms_mc[i].GetHists()[6])
#         # h_cosmic_bkg.Add(histograms_mc[i].GetHists()[7])
#         if h_cosmic_bkg.Integral() > 0:
#             h_cosmic_bkg.Scale(1 / h_cosmic_bkg.Integral())
#         h_cosmic_bkg.SetFillStyle(0)
#         h_cosmic_bkg.SetLineWidth(3)
#         h_cosmic_bkg.SetLineColor(ROOT.kRed + 1)
#         max_hist = max(h_signal.GetMaximum(),
#                        h_neutrino_bkg.GetMaximum(),
#                        h_cosmic_bkg.GetMaximum()) * 1.2

#         h_signal.SetMaximum(max_hist)
#         h_signal.SetMinimum(0.00001)

#         l_norm = ROOT.TLegend(0.59, 0.73, 0.86, 0.82)
#         l_norm.AddEntry(h_signal, "Low-energy excess", "f")
#         l_norm.AddEntry(h_neutrino_bkg, "Neutrino background", "f")
#         l_norm.AddEntry(h_cosmic_bkg, "Cosmic background", "f")
#         h_signal.Draw("hist")
#         h_neutrino_bkg.Draw("hist same")
#         h_cosmic_bkg.Draw("hist same")
#         l_norm.Draw()

#     # track_distance = chain.track_distance < 5
#     # track_res = chain.track_res_std < 2
#     # shower_angle = not (5 < chain.shower_angle < 45)
#     # dqdx = chain.dqdx_bdt_max > 0.1  # and chain.dqdx_bdt > 0.1
#     # corrected_energy = 0.55 < chain.total_shower_energy / \

#         if VARIABLES[i][0] == "shower_dqdx":
#             l1 = ROOT.TLine(1/3.85e-5, 0, 1/3.85e-5, max_hist)
#             l1.SetLineStyle(2)
#             l1.SetLineWidth(3)
#             l2 = ROOT.TLine(3.2/3.85e-5, 0, 3.2/3.85e-5, max_hist)
#             l2.SetLineStyle(2)
#             l2.SetLineWidth(3)
#             l1.Draw()
#             l2.Draw()
#             print("Integral ", h_signal.Integral(
#                 h_signal.FindBin(1.001/3.85e-5), h_signal.FindBin(3.19/3.85e-5)))
#             OBJECTS.append(l1)
#             OBJECTS.append(l2)

#         if VARIABLES[i][0] == "track_distance":
#             l1 = ROOT.TLine(5, 0, 5, max_hist)
#             l1.SetLineStyle(2)
#             l1.SetLineWidth(3)
#             l1.Draw()
#             print("Integral ", h_signal.Integral(0, h_signal.FindBin(4.99)))
#             OBJECTS.append(l1)

#         if VARIABLES[i][0] == "shower_distance":
#             l1 = ROOT.TLine(5, 0, 5, max_hist)
#             l1.SetLineStyle(2)
#             l1.SetLineWidth(3)
#             l1.Draw()
#             h_signal.FindBin(5)
#             print("Integral ", h_signal.Integral(0, h_signal.FindBin(4.99)))
#             OBJECTS.append(l1)

#         if VARIABLES[i][0] == "shower_open_angle":
#             l1 = ROOT.TLine(1, 0, 1, max_hist)
#             l1.SetLineStyle(2)
#             l1.SetLineWidth(3)
#             l2 = ROOT.TLine(20, 0, 20, max_hist)
#             l2.SetLineStyle(2)
#             l2.SetLineWidth(3)
#             l1.Draw()
#             l2.Draw()
#             print("Integral ", h_signal.Integral(h_signal.FindBin(1.01), h_signal.FindBin(18.9)))
#             OBJECTS.append(l1)
#             OBJECTS.append(l2)

#         if VARIABLES[i][0] == "track_shower_angle":
#             l1 = ROOT.TLine(-0.9, 0, -0.9, max_hist)
#             l1.SetLineStyle(2)
#             l1.SetLineWidth(3)
#             l1.Draw()
#             print("Integral ", h_signal.Integral(h_signal.FindBin(-0.89), h_signal.FindBin(0.999)))
#             OBJECTS.append(l1)

#         if VARIABLES[i][0] == "total_hits_y":
#             l1 = ROOT.TLine(100, 0, 100, max_hist)
#             l1.SetLineStyle(2)
#             l1.SetLineWidth(3)
#             l1.Draw()
#             print("Integral ", h_signal.Integral(h_signal.FindBin(101), h_signal.FindBin(10000)))
#             OBJECTS.append(l1)

#         if VARIABLES[i][0] == "shower_energy":
#             l1 = ROOT.TLine(0.050, 0, 0.050, max_hist)
#             l1.SetLineStyle(2)
#             l1.SetLineWidth(3)
#             l1.Draw()
#             print("Integral ", h_signal.Integral(h_signal.FindBin(0.051), h_signal.FindBin(1)))
#             OBJECTS.append(l1)

#         if VARIABLES[i][0] == "track_length":
#             l1 = ROOT.TLine(80, 0, 80, max_hist)
#             l1.SetLineStyle(2)
#             l1.SetLineWidth(3)
#             l1.Draw()
#             print("Integral ", h_signal.Integral(h_signal.FindBin(0), h_signal.FindBin(79)))
#             OBJECTS.append(l1)

#         if VARIABLES[i][0] == "track_pidchipr":
#             l1 = ROOT.TLine(80, 0, 80, max_hist)
#             l1.SetLineStyle(2)
#             l1.SetLineWidth(3)
#             l1.Draw()
#             print("Integral ", h_signal.Integral(h_signal.FindBin(0), h_signal.FindBin(79)))

#             OBJECTS.append(l1)
#         if VARIABLES[i][0] == "shower_pidchipi":
#             l1 = ROOT.TLine(12, 0, 12, max_hist)
#             l1.SetLineStyle(2)
#             l1.SetLineWidth(3)
#             l1.Draw()
#             print("Integral ", h_signal.Integral(
#                 h_signal.FindBin(12.1), h_signal.FindBin(200)))

#             OBJECTS.append(l1)

#         pt = ROOT.TPaveText(0.136, 0.906, 0.501, 0.953, "ndc")
#         pt.AddText("MicroBooNE Preliminary")
#         pt.SetFillColor(0)
#         pt.SetBorderSize(0)
#         pt.SetShadowColor(0)
#         pt.Draw()

#         pt2 = ROOT.TPaveText(0.59, 0.83, 0.80, 0.87, "ndc")
#         pt2.AddText("Area normalized")
#         pt2.SetFillColor(0)
#         pt2.SetBorderSize(0)
#         pt2.SetShadowColor(0)
#         pt2.Draw()

#         c_norm.SetLeftMargin(0.14)
#         c_norm.SetRightMargin(0.06)
#         c_norm.Update()
#         c_norm.SaveAs("plots/integral/%s_norm.pdf" % histograms_bnb[i].GetName())

#         OBJECTS.append(pt)
#         OBJECTS.append(pt2)

#         OBJECTS.append(l_norm)
#         OBJECTS.append(c_norm)
#         OBJECTS.append(h_signal)
#         OBJECTS.append(h_neutrino_bkg)
#         OBJECTS.append(h_cosmic_bkg)


# input()
