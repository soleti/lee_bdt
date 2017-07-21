#!/usr/bin/env python3.4

import ROOT
import math
from bdt_common import binning, labels, variables, spectators, bdt_cut
ROOT.gStyle.SetOptStat(0)
f_data = ROOT.TFile("bnbext_file.root")
t_data = f_data.Get("bnbext_tree")

ROOT.TMVA.Tools.Instance()
reader = ROOT.TMVA.Reader(":".join([
    "!V",
    "!Silent",
    "Color",]))


for name, var in variables:
    t_data.SetBranchAddress(name, var)

for name, var in variables:
    reader.AddVariable(name, var)

for name, var in spectators:
    reader.AddSpectator(name, var)

reader.BookMVA("BDT method","dataset/weights/TMVAClassification_BDT.weights.xml")
description = ["Other", "Cosmic", "Beam Intrinsic #nu_{e}", "Beam Intrinsic #nu_{#mu}", "Beam Intrinsic NC", "Dirt", "Cosmic contaminated"]

variables_dict = dict(variables)

histograms = []

for i,n in enumerate(variables_dict.keys()):
    #h = ROOT.TH1F(n,labels[i],binning[i][0],binning[i][1],binning[i][2])
    h = ROOT.TH1F("h_"+n,labels[n],binning[n][0],binning[n][1],binning[n][2])
    histograms.append(h)

histo_dict = dict(zip(variables_dict.keys(),histograms))


for i in range(t_data.GetEntries()):
    t_data.GetEntry(i)
    BDT_response = reader.EvaluateMVA("BDT method")

    if BDT_response > bdt_cut:

        for name, var in variables:
            histo_dict[name].Fill(var[0], t_data.event_weight)


histograms_data = []
data_files = []
for h in histograms:
    data_files.append(ROOT.TFile("%s_data.root" % h.GetName()))

for i,h in enumerate(histograms):
    h_data = data_files[i].Get(h.GetName())
    histograms_data.append(h_data)


histograms_mc = []
for h in histograms:
    mc_file = ROOT.TFile("%s.root" % h.GetName())
    histograms_mc.append(ROOT.gDirectory.Get(h.GetName()))
    mc_file.Close()

legend = ROOT.TLegend(0.09455587,0.7850208,0.8223496,0.9791956,"","brNDC")
legend.SetTextSize(16)

for j in range(histograms_mc[0].GetNhists()):
    legend.AddEntry(histograms_mc[0].GetHists()[j], "%s: %.0f events" % (description[j], histograms_mc[0].GetHists()[j].Integral()), "f")

legend.AddEntry(histograms_data[0], "Data BNB - BNB EXT shape normalized", "lep")
legend.SetNColumns(2)
pt = ROOT.TPaveText(0.09,0.91,0.60,0.97)
pt.AddText("MicroBooNE Preliminary 6.6e20 POT")
pt.SetFillColor(0)
pt.SetBorderSize(0)
pt.SetShadowColor(0)

canvases = []
h_errs = []
h_ratios = []
pads = []
lines = []
for i in range(len(histograms)):
    c = ROOT.TCanvas("c%i" % i,"",900,44,700,645)

    h_mc_err = histograms_mc[i].GetHists()[0].Clone()
    h_mc_err.SetName("h_mc_err%i" % i)
    for j in range(1,histograms_mc[i].GetNhists()):
        h_mc_err.Add(histograms_mc[i].GetHists()[j])

    pad_top = ROOT.TPad("pad_top","",0, 0.3,1,1)
    pad_top.SetBottomMargin(0)
    pad_top.Range(-22.21825,-2.003018,202.5403,2073.676)
    pad_top.SetTopMargin(0.2411347)
    pad_top.Draw()
    pad_top.cd()
    pads.append(pad_top)
    pt.Draw()

    histograms_mc[i].Draw("hist")
    histograms_mc[i].GetYaxis().SetTitleSize(0.06)
    histograms_mc[i].GetYaxis().SetTitleOffset(0.8)
    histograms_mc[i].SetMinimum(0.1)
    histograms_mc[i].SetMaximum(h_mc_err.GetMaximum()*1.3)

    integral = sum([histograms_mc[i].GetHists()[j].Integral() for j in range(histograms_mc[i].GetNhists())])

    for j in range(histograms_data[i].GetNbinsX()):
        histograms_data[i].SetBinContent(j, histograms_data[i].GetBinContent(j)-histograms[i].GetBinContent(j))
        if histograms_data[i].GetBinContent(j) > 0:
            histograms_data[i].SetBinError(j, math.sqrt(histograms_data[i].GetBinError(j)**2+histograms[i].GetBinError(j)**2))


    histograms_data[i].Scale(integral/histograms_data[i].Integral())
    histograms_data[i].SetLineColor(1)
    histograms_data[i].SetMarkerStyle(20)

    h_mc_err.SetFillStyle(3002)
    h_mc_err.SetFillColor(1)
    h_mc_err.Draw("e2 same")
    histograms_data[i].Draw("ep same")
    h_errs.append(h_mc_err)
    legend.Draw("same")
    c.cd()

    pad_bottom = ROOT.TPad("pad_bottom","",0,0,1,0.3)
    pad_bottom.Range(-22.5,-0.6346511,202.5,1.99)

    pad_bottom.SetFrameBorderMode(0)
    pad_bottom.SetFrameBorderMode(0)
    pad_bottom.SetBorderMode(0)
    pad_bottom.SetTopMargin(0)
    pad_bottom.SetBottomMargin(0.245614)
    pad_bottom.Draw()
    pad_bottom.cd()
    pads.append(pad_bottom)
    h_ratio = histograms_data[i].Clone()
    h_ratios.append(h_ratio)

    h_ratio.SetName("h_ratio%i" % i)

    h_ratio.GetYaxis().SetRangeUser(0.01,1.99)
    h_ratio.Divide(h_mc_err)

    h_ratio.GetXaxis().SetLabelFont(42);
    h_ratio.GetXaxis().SetLabelSize(0.13);
    h_ratio.GetXaxis().SetTitleSize(0.13);
    h_ratio.GetXaxis().SetTitleOffset(0.91);
    h_ratio.GetYaxis().SetTitle("Data/MC");
    h_ratio.GetYaxis().SetNdivisions(509);
    h_ratio.GetYaxis().SetLabelFont(42);
    h_ratio.GetYaxis().SetLabelSize(0.13);
    h_ratio.GetYaxis().SetTitleSize(0.13);
    h_ratio.GetYaxis().SetTitleOffset(0.4);
    h_ratio.Draw("ep")
    line = ROOT.TLine(h_ratio.GetXaxis().GetXmin(), 1, h_ratio.GetXaxis().GetXmax(), 1)
    line.SetLineWidth(2)
    line.SetLineStyle(2)
    line.Draw()
    lines.append(line)
    c.cd()
    c.Update()
    c.SaveAs("%s.pdf" % histograms[i].GetName())
    canvases.append(c)

input()
