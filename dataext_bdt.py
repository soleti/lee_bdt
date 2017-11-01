#!/usr/local/bin/python3

import ROOT
import math

from bdt_common import binning, labels, variables, spectators, bdt_cut, manual_cuts
from bdt_common import description, total_pot, sigmaCalc

from glob import glob
from array import array

reco_energy_index = 22
post_scaling = 66 / 5
total_pot *= post_scaling

draw_subtraction = False
draw_ext = False
draw_lee = True
draw_data = False
ROOT.gStyle.SetOptStat(0)

f_data = ROOT.TFile("bnbext_file.root")
t_data = f_data.Get("bnbext_tree")
data_bnb = glob("data_files_bnb_6_42_energy/*/*.root")

chain_data_bnb_pot = ROOT.TChain("robertoana/pot")
for f in data_bnb:
    chain_data_bnb_pot.Add(f)
total_data_bnb_pot = 0

for i in range(chain_data_bnb_pot.GetEntries()):
    chain_data_bnb_pot.GetEntry(i)
    total_data_bnb_pot += chain_data_bnb_pot.pot
total_data_bnb_pot *= 1e12

ROOT.TMVA.Tools.Instance()
reader = ROOT.TMVA.Reader(":".join([
    "!V",
    "!Silent",
    "Color"]))

for name, var in variables:
    t_data.SetBranchAddress(name, var)
for name, var in spectators:
    t_data.SetBranchAddress(name, var)

for name, var in variables:
    reader.AddVariable(name, var)

for name, var in spectators:
    reader.AddSpectator(name, var)

reader.BookMVA("BDT method",
               "dataset/weights/TMVAClassification_BDT.weights.xml")

variables_dict = dict(variables + spectators)

histograms = []

for i, n in enumerate(variables_dict.keys()):
    print(i, labels[n])

    if n != "reco_energy":
        h = ROOT.TH1F("h_%s" % n, labels[n],
                      binning[n][0], binning[n][1], binning[n][2])
    else:
        bins = array("f", [0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.6, 0.8, 1])
        h = ROOT.TH1F("h_%s" % n, labels[n], len(bins) - 1, bins)
    h.SetLineColor(1)
    h.SetMarkerStyle(20)
    histograms.append(h)

histo_dict = dict(zip(variables_dict.keys(), histograms))

h_bdt = ROOT.TH1F("h_bdt", ";BDT response; N. Entries / 0.05", 40, -1, 1)

for i in range(t_data.GetEntries()):
    t_data.GetEntry(i)
    BDT_response = reader.EvaluateMVA("BDT method")
    h_bdt.Fill(BDT_response, t_data.event_weight)

    if BDT_response > bdt_cut and manual_cuts(t_data):

        for name, var in variables:
            histo_dict[name].Fill(var[0], t_data.event_weight)
        for name, var in spectators:
            histo_dict[name].Fill(var[0], t_data.event_weight)

f_bdt = ROOT.TFile("plots/h_bdt_dataext.root", "RECREATE")
h_bdt.Write()
f_bdt.Close()

histograms_data = []
data_files = []
for h in histograms:
    data_files.append(ROOT.TFile("plots/%s_bnb.root" % h.GetName()))

histograms_cosmic = []
cosmic_files = []
for h in histograms:
    cosmic_files.append(ROOT.TFile("plots/%s_cosmic_mc.root" % h.GetName()))

for i, h in enumerate(histograms):
    h_data = data_files[i].Get(h.GetName())
    h_data.SetLineColor(1)
    h_data.SetMarkerStyle(20)
    h_cosmic = cosmic_files[i].Get(h.GetName())
    histograms_data.append(h_data)
    histograms_cosmic.append(h_cosmic)

histograms_mc = []
legends = []

for h in histograms:
    mc_file = ROOT.TFile("plots/%s.root" % h.GetName())
    histograms_mc.append(ROOT.gDirectory.Get(h.GetName()))
    mc_file.Close()

    legend = ROOT.TLegend(0.09455587, 0.7850208, 0.8923496, 0.9791956, "", "brNDC")
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
                                           Integral()*post_scaling), "f")

    if draw_ext:
        legends[i].AddEntry(histograms[i], "Data EXT: {:.0f} events"
                            .format(histograms[i].Integral()*post_scaling), "f")
    else:
        legends[i].AddEntry(histograms_cosmic[i], "Cosmic in-time: {:.0f} events"
                            .format(histograms_cosmic[i].Integral()*post_scaling), "f")
    if draw_data:
        legends[i].AddEntry(histograms[i], "Data BNB: {:.0f} events"
                            .format(histograms_data[i].Integral()*post_scaling), "lep")



for i in range(len(histograms)):
    histograms_mc[i].GetHists()[2].SetFillStyle(3001)


if draw_subtraction:
    for i in range(len(histograms)):
        for j in range(histograms_data[i].GetNbinsX()):
            bnb = histograms_data[i].GetBinContent(j)
            ext = histograms[i].GetBinContent(j)

            histograms_data[i].SetBinContent(j, bnb - ext)

            if bnb > 0:
                err1 = histograms_data[i].GetBinError(j)
                err2 = histograms[i].GetBinError(j)
                histograms_data[i].SetBinError(j, math.sqrt(err1**2 + err2**2))

    legend[i].AddEntry(histograms_data[i], "Data BNB - BNB EXT: {:.0f} events"
                    .format(histograms_data[i].Integral()*post_scaling), "lep")


legend_cosmic = ROOT.TLegend(0.099, 0.909, 0.900, 0.987, "", "brNDC")
legend_cosmic.AddEntry(histograms[reco_energy_index], "Data EXT: {:.0f} events"
                       .format(histograms[reco_energy_index].Integral()*post_scaling), "lep")
legend_cosmic.AddEntry(histograms_cosmic[0],
                       "CORSIKA in-time Monte Carlo: integral normalized", "f")
legend_cosmic.SetNColumns(2)
canvases = []
h_errs = []
h_errs_sys = []
h_ratios = []
pads = []
lines = []
canvases_cosmic = []
bins = array("f", [0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.6, 0.8, 1])

h_lee = ROOT.TH1F("h_lee", "", len(bins) - 1, bins)

scaling = [6.0920944819073988, 3.6447414342239273, 3.2123920194399913,
           2.6504659907742409, 3.2558450032216988, 2.5826310533377432,
           2, 1, 1, 1, 1]


# scaling = [6.0920944819073988, 3.6447414342239273, 3.2123920194399913,
#            2.6504659907742409, 3.2558450032216988, 2.5826310533377432,
#            2.6614353575699727, 1.4145769564088304, 1.0206172427887652,
#            0.9972444259255292, 0.79323702430381904, 0.63892043872491167,
#            0.61676413081900316, 0.3541651442224471, 0.28310400773433003,
#            0.94342108559739024]

# scaling = [5.015008606, 4.755966764, 4.240843625, 3.494299576, 2.596148682, 1.715699034, 1.051751175, 0.8966301403, 0.9134166508, 1.010169572, 1.075023103]

h_lee.SetFillColor(ROOT.kGreen - 2)
h_lee.SetFillStyle(3005)
h_lee.SetLineColor(1)

for i in range(len(scaling)):
    if scaling[i] - 1 > 0:
        h_lee.SetBinContent(i+1, histograms_mc[reco_energy_index].GetHists()[3].GetBinContent(i+1)*(scaling[i] - 1))

if draw_lee:
    legends[reco_energy_index].AddEntry(h_lee, "Low-energy excess: {:.0f}"
                    .format(h_lee.Integral()*post_scaling), "f")

rebin = 1


for i in range(len(histograms)):

    histograms[i].SetLineColor(ROOT.kBlack)
    histograms_cosmic[i].SetLineColor(ROOT.kBlack)
    histograms_cosmic[i].SetFillColor(ROOT.kOrange + 1)
    if i == reco_energy_index:
        h_ext = histograms[i].Clone()
        h_intime = histograms_cosmic[i].Clone()
        # histograms[i].Scale(post_scaling)
        # histograms_cosmic[i].Scale(post_scaling)

        c_cosmic = ROOT.TCanvas("c{}_canvas".format(i), "", 900, 44, 700, 645)
        # print(histograms[i].Integral() / histograms_cosmic[i].Integral())
        # if histograms_cosmic[i].Integral() > 0:
        #     histograms_cosmic[i].Scale(
        #         histograms[i].Integral() / histograms_cosmic[i].Integral())
        for k in range(1, h_intime.GetNbinsX()+1):
            bin_width = h_intime.GetBinWidth(k)
            h_intime.SetBinError(k, h_intime.GetBinError(k)/(bin_width/0.05))
            h_intime.SetBinContent(k, h_intime.GetBinContent(k)/(bin_width/0.05))
            h_ext.SetBinError(k, h_ext.GetBinError(k)/(bin_width/0.05))
            h_ext.SetBinContent(k, h_ext.GetBinContent(k)/(bin_width/0.05))

        pad_top = ROOT.TPad("pad_top_cosmic%i" % i, "", 0, 0.3, 1, 1)
        pad_top.SetBottomMargin(0)
        pad_top.Range(-22.21825, -2.003018, 202.5403, 2073.676)
        pad_top.SetTopMargin(0.2411347)
        pad_top.Draw()
        pad_top.cd()
        pads.append(pad_top)
        h_intime.Draw("hist")

        h_ext.SetMarkerStyle(20)
        h_ext.Draw("ep same")
        legend_cosmic.Draw()

        upper_limit = h_intime.GetMaximum() * 1.3
        h_intime.GetYaxis().SetRangeUser(0.01, upper_limit)
        c_cosmic.cd()
        pad_bottom = ROOT.TPad("pad_bottom_cosmic%i" % i, "", 0, 0, 1, 0.3)
        pad_bottom.Range(-22.5, -0.6346511, 202.5, 1.99)

        pad_bottom.SetFrameBorderMode(0)
        pad_bottom.SetFrameBorderMode(0)
        pad_bottom.SetBorderMode(0)
        pad_bottom.SetTopMargin(0)
        pad_bottom.SetBottomMargin(0.245614)
        pad_bottom.Draw()
        pad_bottom.cd()
        pads.append(pad_bottom)
        h_ratio = h_ext.Clone()
        h_ratios.append(h_ratio)

        h_ratio.SetName("h_ratio_cosmic%i" % i)

        h_ratio.GetYaxis().SetRangeUser(0.01, 1.99)
        h_ratio.Divide(h_intime)

        h_ratio.GetXaxis().SetLabelFont(42)
        h_ratio.GetXaxis().SetLabelSize(0.13)
        h_ratio.GetXaxis().SetTitleSize(0.13)
        h_ratio.GetXaxis().SetTitleOffset(0.91)
        h_ratio.GetYaxis().SetTitle("EXT / MC")
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
        c_cosmic.cd()

        c_cosmic.Update()
        c_cosmic.SaveAs("plots/%s_cosmic.pdf" % histograms[i].GetName())
        canvases_cosmic.append(c_cosmic)

    histograms_mc[i].GetHists()[2].SetFillStyle(3001)

    c = ROOT.TCanvas("c%i" % i, "", 900, 44, 700, 645)

    if not draw_subtraction and draw_ext:
        histograms_mc[i].Add(histograms[i])
    if not draw_ext:
        histograms_mc[i].Add(histograms_cosmic[i])

    #histograms_mc[i].GetHists()[0].Rebin(rebin)
    histograms_mc[i].GetHists()[0].Scale(post_scaling)
    if i == reco_energy_index:

        for k in range(1,histograms_mc[i].GetHists()[0].GetNbinsX()+1):
            bin_width = histograms_mc[i].GetHists()[0].GetBinWidth(k)
            histograms_mc[i].GetHists()[0].SetBinError(k, histograms_mc[i].GetHists()[0].GetBinError(k)/(bin_width/0.05))
            histograms_mc[i].GetHists()[0].SetBinContent(k, histograms_mc[i].GetHists()[0].GetBinContent(k)/(bin_width/0.05))

    h_mc_err = histograms_mc[i].GetHists()[0].Clone()
    h_mc_err.SetName("h_mc_err%i" % i)

    for j in range(1, histograms_mc[i].GetNhists()):
        #histograms_mc[i].GetHists()[j].Rebin(rebin)
        histograms_mc[i].GetHists()[j].Scale(post_scaling)

        if i == reco_energy_index:
            for k in range(1,histograms_mc[i].GetHists()[j].GetNbinsX()+1):
                bin_width = histograms_mc[i].GetHists()[j].GetBinWidth(k)
                histograms_mc[i].GetHists()[j].SetBinError(k, histograms_mc[i].GetHists()[j].GetBinError(k)/(bin_width/0.05))
                histograms_mc[i].GetHists()[j].SetBinContent(k, histograms_mc[i].GetHists()[j].GetBinContent(k)/(bin_width/0.05))


        h_mc_err.Add(histograms_mc[i].GetHists()[j])

    h_mc_err_sys = h_mc_err.Clone()
    h_mc_err_sys.SetName("h_mc_err_sys%i" % i)
    for k in range(1, h_mc_err.GetNbinsX()+1):
        #new_err = h_mc_err_sys.GetBinError(k)**2 + (0.15*h_mc_err_sys.GetBinContent(k))**2
        new_err = h_mc_err.GetBinContent(k)
        h_mc_err_sys.SetBinError(k, math.sqrt(new_err))

    h_mc_err_stat = h_mc_err_sys.Clone()
    h_mc_err_stat.SetName("h_mc_err_stat%i" % i)
    for k in range(1, h_mc_err_stat.GetNbinsX()+1):
        new_err = h_mc_err_sys.GetBinError(k)**2 + (0.2*h_mc_err_sys.GetBinContent(k))**2
        h_mc_err_stat.SetBinError(k, math.sqrt(new_err))
    if i == reco_energy_index and draw_lee:
        #h_lee.Rebin(rebin)
        h_lee.Scale(post_scaling)

        # for j in range(1,histograms_mc[i].GetHists()[0].GetNbinsX()+1):
        #     unscaled_nue = histograms_mc[i].GetHists()[3].GetBinContent(j) / 0.05
        #     nue = histograms_mc[i].GetHists()[3].GetBinContent(j)
        #     unscaled = 0
        #     for k in range(1,len(histograms_mc[i].GetHists())):
        #         if k != 3:
        #             unscaled += histograms_mc[i].GetHists()[k].GetBinContent(j) / 1.3
        #         print(j, k, histograms_mc[i].GetHists()[k].GetBinContent(j))
        #     print(j, h_mc_err.GetBinError(j))
        #     print(j, nue)
        #     print(j, math.sqrt(nue))
        #     print(j, math.sqrt(unscaled_nue)/unscaled_nue * nue)
        for k in range(1, h_lee.GetNbinsX()+1):
            bin_width = h_lee.GetBinWidth(k)
            h_lee.SetBinError(k, h_lee.GetBinError(k)/(bin_width/0.05))
            h_lee.SetBinContent(k, h_lee.GetBinContent(k)/(bin_width/0.05))
        # for isys in range(30):
        print("Sigma stat", sigmaCalc(h_lee, h_mc_err_sys))
        print("Sigma 20 sys", sigmaCalc(h_lee, h_mc_err_sys, 0.2))

        histograms_mc[i].Add(h_lee)

    if draw_data:
        pad_top = ROOT.TPad("pad_top", "", 0, 0.3, 1, 1)
        pad_top.SetBottomMargin(0)
        pad_top.Range(-22.21825, -2.003018, 202.5403, 2073.676)
        pad_top.SetTopMargin(0.2411347)
        pad_top.Draw()
        pad_top.cd()
        pads.append(pad_top)
    else:
        c.SetTopMargin(0.2274194)

    histograms_mc[i].Draw("hist")
    histograms_mc[i].GetYaxis().SetTitleSize(0.06)
    histograms_mc[i].GetYaxis().SetTitleOffset(0.8)
    histograms_mc[i].SetMinimum(0.1)
    histograms_mc[i].SetMaximum(h_mc_err.GetMaximum() * 1.3)
    histograms_data[i].SetLineColor(1)
    histograms_data[i].SetMarkerStyle(20)

    h_mc_err.SetFillStyle(3002)
    h_mc_err.SetFillColor(1)
    h_mc_err_sys.SetFillStyle(3002)
    h_mc_err_sys.SetFillColor(1)
    if i == reco_energy_index:
        h_mc_err_sys.Draw("e1 same")
        #h_mc_err_stat.Draw("e1 same")
    else:
        h_mc_err.Draw("e1 same")

    #histograms_data[i].Rebin(rebin)
    if draw_data:
        if i == reco_energy_index:
            for k in range(1,histograms_data[i].GetNbinsX()+1):
                bin_width = histograms_data[i].GetBinWidth(k)
                histograms_data[i].SetBinError(k, histograms_data[i].GetBinError(k)/(bin_width/0.05))
                histograms_data[i].SetBinContent(k, histograms_data[i].GetBinContent(k)/(bin_width/0.05))

        histograms_data[i].Scale(post_scaling)
        histograms_data[i].Draw("ep same")

    h_errs.append(h_mc_err)
    h_errs_sys.append(h_mc_err_sys)
    h_errs_sys.append(h_mc_err_stat)

    legends[i].Draw("same")
    c.cd()

    if draw_data:
        pad_bottom = ROOT.TPad("pad_bottom", "", 0, 0, 1, 0.3)
        pad_bottom.Range(-22.5, -0.6346511, 202.5, 1.99)

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

        h_ratio.GetYaxis().SetRangeUser(0.01, 1.99)
        h_ratio.Divide(h_mc_err_sys)

        h_ratio.GetXaxis().SetLabelFont(42)
        h_ratio.GetXaxis().SetLabelSize(0.13)
        h_ratio.GetXaxis().SetTitleSize(0.13)
        h_ratio.GetXaxis().SetTitleOffset(0.91)
        h_ratio.GetYaxis().SetTitle("BNB / (MC + EXT)")
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
        c.cd()

    c.Update()
    c.SaveAs("plots/%s.pdf" % histograms[i].GetName(), "Q")
    canvases.append(c)
input()
