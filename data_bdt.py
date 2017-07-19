#!/usr/bin/env python3.4

import ROOT
import math
from bdt_common import binning, labels, variables, spectators, bdt_cut

f_data = ROOT.TFile("bnb_file.root")
t_data = f_data.Get("bnb_tree")

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


histograms_mc = []
for h in histograms:
    mc_file = ROOT.TFile("%s.root" % h.GetName())
    histograms_mc.append(ROOT.gDirectory.Get(h.GetName()))
    mc_file.Close()

legend = ROOT.TLegend(0.455,0.53,0.70,0.85)
legend.SetTextSize(16)

for j in range(histograms_mc[0].GetNhists()):
    legend.AddEntry(histograms_mc[0].GetHists()[j], "%s: %.0f events" % (description[j], histograms_mc[0].GetHists()[j].Integral()), "f")

legend.AddEntry(histograms[0], "Data BNB - BNB EXT shape normalized", "lep")

pt = ROOT.TPaveText(0.09,0.91,0.60,0.97, "ndc")
pt.AddText("MicroBooNE Preliminary 6.6e20 POT")
pt.SetFillColor(0)
pt.SetBorderSize(0)
pt.SetShadowColor(0)

canvases = []*len(histograms)

for i in range(len(histograms)):
    c = ROOT.TCanvas("c%i" % i)
    histograms_mc[i].Draw("hist")
    histograms_mc[i].GetYaxis().SetRangeUser(0,histograms[i].GetMaximum()*1.4)

    integral = sum([histograms_mc[i].GetHists()[j].Integral() for j in range(histograms_mc[i].GetNhists())])
    histograms[i].Scale(integral/histograms[i].Integral())
    histograms[i].SetLineColor(1)
    histograms[i].SetMarkerStyle(20)
    histograms[i].Draw("ep same")
    legend.Draw("same")

    pt.Draw()
    c.Update()
    c.SaveAs("%s.pdf" % histograms[i].GetName())
    canvases.append(c)

input()
