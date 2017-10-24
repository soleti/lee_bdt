#!/usr/local/bin/python3

import ROOT
import math

from bdt_common import binning, labels, variables, spectators, bdt_cut
from bdt_common import description, total_pot, sigmaCalc

from glob import glob

draw_subtraction = False
draw_lee = True
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


variables_dict = dict(variables+spectators)

histograms = []

for i, n in enumerate(variables_dict.keys()):
    print(i, labels[n])

    h = ROOT.TH1F("h_" + n, labels[n],
                  binning[n][0], binning[n][1], binning[n][2])
    h.SetLineColor(1)
    h.SetMarkerStyle(20)
    histograms.append(h)

histo_dict = dict(zip(variables_dict.keys(), histograms))

h_bdt = ROOT.TH1F("h_bdt", ";BDT response; N. Entries / 0.05", 40, -1, 1)

passed_events = open("dataext_passed.txt", "w")

for i in range(t_data.GetEntries()):
    t_data.GetEntry(i)
    BDT_response = reader.EvaluateMVA("BDT method")
    h_bdt.Fill(BDT_response, t_data.event_weight)

    if BDT_response > bdt_cut:
        print("{} {} {} {}".format(int(t_data.run),
                                   int(t_data.subrun),
                                   int(t_data.event),
                                   t_data.event_weight * 2),
              file=passed_events)
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
    data_files.append(ROOT.TFile("plots/%s_data.root" % h.GetName()))

histograms_cosmic = []
cosmic_files = []
for h in histograms:
    cosmic_files.append(ROOT.TFile("plots/%s_cosmic.root" % h.GetName()))

for i, h in enumerate(histograms):
    h_data = data_files[i].Get(h.GetName())
    h_data.SetLineColor(1)
    h_data.SetMarkerStyle(20)
    h_cosmic = cosmic_files[i].Get(h.GetName())
    histograms_data.append(h_data)
    histograms_cosmic.append(h_cosmic)

histograms_mc = []
for h in histograms:
    mc_file = ROOT.TFile("plots/%s.root" % h.GetName())
    histograms_mc.append(ROOT.gDirectory.Get(h.GetName()))
    mc_file.Close()

legend = ROOT.TLegend(0.09455587, 0.7850208, 0.8923496, 0.9791956, "", "brNDC")
legend.SetTextSize(16)
legend.SetTextFont(63)
legend.SetHeader("MicroBooNE Preliminary %.1e POT" % total_pot)
legend.SetTextFont(43)

for j in range(histograms_mc[24].GetNhists()):
    if histograms_mc[24].GetHists()[j].Integral():
        legend.AddEntry(
            histograms_mc[24].GetHists()[j],
            "{}: {:.0f} events".format(description[j],
                                       histograms_mc[24].GetHists()[j].
                                       Integral(1, 8)),
            "f")

# legend.AddEntry(histograms_data[0], "Data BNB: {:.0f} events"
#                 .format(histograms_data[24].Integral()), "lep")
#
# legend.AddEntry(histograms[0], "Data EXT: {:.0f} events"
#                 .format(histograms[24].Integral()), "f")

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

    legend.AddEntry(histograms_data[0], "Data BNB - BNB EXT: {:.0f} events"
                    .format(histograms_data[0].Integral()), "lep")

legend.SetNColumns(2)

legend_cosmic = ROOT.TLegend(0.099, 0.909, 0.900, 0.987, "", "brNDC")
legend_cosmic.AddEntry(histograms[24], "Data EXT: {:.0f} events"
                       .format(histograms[24].Integral()), "lep")
legend_cosmic.AddEntry(histograms_cosmic[0],
                       "CORSIKA in-time Monte Carlo: integral normalized", "f")
legend_cosmic.SetNColumns(2)
canvases = []
h_errs = []
h_ratios = []
pads = []
lines = []
canvases_cosmic = []

h_lee = ROOT.TH1F("h_lee", "", binning["reco_energy"][0], binning["reco_energy"][1], binning["reco_energy"][2])

scaling = [6.0920944819073988, 3.6447414342239273, 3.2123920194399913, 2.6504659907742409, 3.2558450032216988, 2.5826310533377432, 2.6614353575699727, 1.4145769564088304, 1.0206172427887652, 0.9972444259255292, 0.79323702430381904, 0.63892043872491167, 0.61676413081900316, 0.3541651442224471, 0.28310400773433003, 0.94342108559739024]
#scaling = [5.015008606, 4.755966764, 4.240843625, 3.494299576, 2.596148682, 1.715699034, 1.051751175, 0.8966301403, 0.9134166508, 1.010169572, 1.075023103]
h_lee.SetFillColor(ROOT.kGreen - 2)
h_lee.SetFillStyle(3005)
h_lee.SetLineColor(1)

for i in range(len(scaling)):
    if scaling[i] - 1 > 0:
        h_lee.SetBinContent(i+1, histograms_mc[24].GetHists()[3].GetBinContent(i+1)*(scaling[i] - 1))

if draw_lee:
    legend.AddEntry(h_lee, "Low-energy excess: {:.0f}".format(h_lee.Integral()), "f")

for i in range(len(histograms)):

    histograms[i].SetLineColor(ROOT.kBlack)

    # c_cosmic = ROOT.TCanvas("c{}_canvas".format(i), "", 900, 44, 700, 645)
    # if histograms_cosmic[i].Integral() > 0:
    #     histograms_cosmic[i].Scale(
    #         histograms[i].Integral() / histograms_cosmic[i].Integral())
    # histograms_cosmic[i].SetLineColor(ROOT.kBlack)
    # histograms_cosmic[i].SetFillColor(ROOT.kRed - 3)
    # histograms_cosmic[i].Draw("hist")
    #
    # histograms[i].SetMarkerStyle(20)
    # histograms[i].Draw("ep same")
    # legend_cosmic.Draw()
    #
    # upper_limit = histograms_cosmic[i].GetMaximum() * 1.3
    # histograms_cosmic[i].GetYaxis().SetRangeUser(0.01, upper_limit)
    #
    # c_cosmic.Update()
    # c_cosmic.SaveAs("plots/%s_cosmic.pdf" % histograms[i].GetName())
    # canvases_cosmic.append(c_cosmic)

    histograms_mc[i].GetHists()[2].SetFillStyle(3001)

    c = ROOT.TCanvas("c%i" % i, "", 900, 44, 700, 645)

    # if not draw_subtraction:
    #     histograms_mc[i].Add(histograms[i])

    h_mc_err = histograms_mc[i].GetHists()[0].Clone()
    h_mc_err.SetName("h_mc_err%i" % i)
    for j in range(1, histograms_mc[i].GetNhists()):
        h_mc_err.Add(histograms_mc[i].GetHists()[j])

    if i == 24 and draw_lee:
        #for s in range(100):
        print("Sigma", sigmaCalc(h_lee, h_mc_err))
        histograms_mc[i].Add(h_lee)

    # pad_top = ROOT.TPad("pad_top", "", 0, 0.3, 1, 1)
    # pad_top.SetBottomMargin(0)
    # pad_top.Range(-22.21825, -2.003018, 202.5403, 2073.676)
    # pad_top.SetTopMargin(0.2411347)
    # pad_top.Draw()
    # pad_top.cd()
    # pads.append(pad_top)
    histograms_mc[i].Draw("hist")
    histograms_mc[i].GetYaxis().SetTitleSize(0.06)
    histograms_mc[i].GetYaxis().SetTitleOffset(0.8)
    histograms_mc[i].SetMinimum(0.1)
    histograms_mc[i].SetMaximum(h_mc_err.GetMaximum() * 1.3)

    histograms_data[i].SetLineColor(1)
    histograms_data[i].SetMarkerStyle(20)

    h_mc_err.SetFillStyle(3002)
    h_mc_err.SetFillColor(1)
    h_mc_err.Draw("e2 same")
    #histograms_data[i].Draw("ep same")
    h_errs.append(h_mc_err)
    legend.Draw("same")
    c.cd()
    #
    # pad_bottom = ROOT.TPad("pad_bottom", "", 0, 0, 1, 0.3)
    # pad_bottom.Range(-22.5, -0.6346511, 202.5, 1.99)
    #
    # pad_bottom.SetFrameBorderMode(0)
    # pad_bottom.SetFrameBorderMode(0)
    # pad_bottom.SetBorderMode(0)
    # pad_bottom.SetTopMargin(0)
    # pad_bottom.SetBottomMargin(0.245614)
    # pad_bottom.Draw()
    # pad_bottom.cd()
    # pads.append(pad_bottom)
    # h_ratio = histograms_data[i].Clone()
    # h_ratios.append(h_ratio)
    #
    # h_ratio.SetName("h_ratio%i" % i)
    #
    # h_ratio.GetYaxis().SetRangeUser(0.01, 1.99)
    # h_ratio.Divide(h_mc_err)
    #
    # h_ratio.GetXaxis().SetLabelFont(42)
    # h_ratio.GetXaxis().SetLabelSize(0.13)
    # h_ratio.GetXaxis().SetTitleSize(0.13)
    # h_ratio.GetXaxis().SetTitleOffset(0.91)
    # h_ratio.GetYaxis().SetTitle("BNB / (MC + EXT)")
    # h_ratio.GetYaxis().SetNdivisions(509)
    # h_ratio.GetYaxis().SetLabelFont(42)
    # h_ratio.GetYaxis().SetLabelSize(0.13)
    # h_ratio.GetYaxis().SetTitleSize(0.13)
    # h_ratio.GetYaxis().SetTitleOffset(0.36)
    # h_ratio.Draw("ep")
    # line = ROOT.TLine(h_ratio.GetXaxis().GetXmin(), 1,
    #                   h_ratio.GetXaxis().GetXmax(), 1)
    # line.SetLineWidth(2)
    # line.SetLineStyle(2)
    # line.Draw()
    # lines.append(line)
    # c.cd()
    c.Update()
    c.SaveAs("plots/%s.pdf" % histograms[i].GetName())
    canvases.append(c)
input()
