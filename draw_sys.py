#!/usr/local/bin/python3

import ROOT

from bdt_common import variables, spectators, description, total_pot, fix_binning

ROOT.gStyle.SetOptStat(0)

SYS_COLORS = [ROOT.kBlack, ROOT.kBlue + 1, ROOT.kGreen + 1, ROOT.kRed + 1]
DESCRIPTIONS = ["Standard MC + DIC",
                "Standard MC + LArG4 fix",
                "Standard MC + data SCE"]
OBJECTS = []
SYS_ERR = 0.1

histograms_bnb = []
histograms_mc = []
histograms_mc_dic = []
histograms_mc_larg4 = []
histograms_mc_sce = []

VARIABLES = variables + spectators
RECO_ENERGY = list(dict(VARIABLES)).index("reco_energy")

histograms = [histograms_bnb,
              histograms_mc,
              histograms_mc_dic,
              histograms_mc_sce,
              histograms_mc_larg4]


def draw_histo(histo, color):
    histo.SetFillColor(color)
    histo.SetLineColor(color)
    histo.SetFillStyle(3002)
    histo.SetLineWidth(2)

    histo.Draw("e2 same")
    OBJECTS.append(histo)

    histo_clone = histo.Clone()
    histo_clone.SetLineWidth(2)
    histo_clone.SetFillStyle(0)
    histo_clone.Draw("hist same")
    OBJECTS.append(histo_clone)


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

    OBJECTS.append(h_ratio)

    h_ratio.GetYaxis().SetRangeUser(0.01, 1.99)
    h_ratio.Divide(den[0])

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
    h_ratio.SetMarkerSize(0)
    h_ratio_2 = h_ratio.Clone()
    h_ratio_2.SetFillStyle(3002)
    h_ratio_2.SetFillColor(1)
    OBJECTS.append(h_ratio_2)
    h_ratio_2.Draw("e2 same")

    for i, h in enumerate(den[1:]):
        h_ratio_sys = num.Clone()
        h_ratio_sys.Divide(h)
        h_ratio_sys.SetMarkerSize(0)
        h_ratio_sys.Draw("hist same")
        h_ratio_sys_2 = h_ratio_sys.Clone()
        h_ratio_sys_2.SetFillStyle(3002)
        h_ratio_sys_2.SetFillColor(SYS_COLORS[i + 1])
        h_ratio_sys_2.Draw("e2 same")
        h_ratio_sys.SetLineColor(SYS_COLORS[i+1])
        OBJECTS.append(h_ratio_sys)
        OBJECTS.append(h_ratio_sys_2)


    line = ROOT.TLine(h_ratio.GetXaxis().GetXmin(), 1,
                      h_ratio.GetXaxis().GetXmax(), 1)
    line.SetLineWidth(2)
    line.SetLineStyle(2)
    line.Draw()
    OBJECTS.append(line)


legends = []
l_errs = []

cat = ["bnb", "mc", "mc", "mc", "mc"]
dirs = ["plots_cv", "plots_cv", "plots_with_dic", "plots_with_larg4", "plots_with_sce"]

for name, var in VARIABLES:

    for i, histos in enumerate(histograms):

        f = ROOT.TFile("%s/h_%s_%s.root" % (dirs[i], name, cat[i]))
        h = f.Get("h_%s" % name)

        h.SetName("h_%s_%s" % (name, cat[i]))
        if cat[i] == "bnb":
            h.SetDirectory(0)
        histos.append(h)

        f.Close()

    legend = ROOT.TLegend(0.085, 0.7747, 0.891117, 0.97926, "", "brNDC")
    legend.SetTextSize(16)
    legend.SetTextFont(63)
    legend.SetHeader("MicroBooNE Preliminary %.1e POT" % total_pot)
    legend.SetTextFont(43)
    legend.SetNColumns(2)
    legends.append(legend)

for i, h in enumerate(histograms_mc):
    if h:
        on_off = histograms_bnb[i].Integral() - h.GetHists()[0].Integral()
        legends[i].AddEntry(histograms_bnb[i], 
                            "Data on - off: {:.0f} events".format(on_off), "lep")

for i in range(len(VARIABLES)):
# for i in range(RECO_ENERGY, RECO_ENERGY+1):
    if not histograms_mc[i] or histograms_mc[i].GetHists()[0].Integral() <= 0:
        continue 
        
    c = ROOT.TCanvas("c%i" % i, "", 900, 44, 700, 645)

    histograms_bnb[i].SetLineColor(1)
    histograms_bnb[i].SetMarkerStyle(20)

    draw_top()

    histograms_bnb[i].Add(histograms_mc[i].GetHists()[0], -1)
    h_sum = histograms_mc[i].GetStack().Last().Clone()
    h_sum.Add(histograms_mc[i].GetHists()[0], -1)

    legends[i].AddEntry(h_sum, "Standard MC: %.1f events" %
                        h_sum.Integral(), "lf")

    if i == RECO_ENERGY:
        ratio = histograms_bnb[i].Integral() / (h_sum.Integral())
        fix_binning(h_sum)
        fix_binning(histograms_bnb[i])

    h_sum.Draw("e2")
    h_sum.SetFillColor(ROOT.kBlack)
    h_sum.SetLineColor(ROOT.kBlack)
    h_sum.SetLineWidth(2)
    h_sum.SetFillStyle(3002)
    OBJECTS.append(h_sum)

    h_sum_clone = h_sum.Clone()
    h_sum_clone.SetLineWidth(2)
    h_sum_clone.SetFillStyle(0)
    h_sum_clone.Draw("hist same")
    OBJECTS.append(h_sum_clone)

    set_axis(h_sum, h_sum.GetMaximum() * 1.35)

    h_sums = [h_sum]

    for j, h_sys in enumerate(histograms[2:]):
        h_sys_sum = h_sys[i].GetStack().Last().Clone()
        h_sys_sum.Add(h_sys[i].GetHists()[0], -1)
        legends[i].AddEntry(h_sys_sum, "%s: %.1f events" %
                            (DESCRIPTIONS[j], h_sys_sum.Integral()), "lf")

        if i == RECO_ENERGY:
            fix_binning(h_sys_sum)

        draw_histo(h_sys_sum, SYS_COLORS[j + 1])
        h_sums.append(h_sys_sum)

    histograms_bnb[i].Draw("ep same")
    legends[i].Draw("same")
    c.cd()

    draw_ratio(histograms_bnb[i], h_sums)

    c.Update()
    c.SaveAs("plots_sys/%s.pdf" % histograms_bnb[i].GetName()[:-4])

    OBJECTS.append(c)

input()
