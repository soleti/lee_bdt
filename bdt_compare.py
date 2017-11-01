#!/usr/local/bin/python3

import ROOT
from bdt_common import description, bdt_cut, total_pot

ROOT.gStyle.SetOptStat(0)

c_bdt = ROOT.TCanvas("c_bdt", "BDT response", 900, 44, 700, 645)
f_bdt_data = ROOT.TFile("plots/h_bdt_bnb.root")
h_bdt_data = f_bdt_data.Get("h_bdt_bnb")
f_bdt_mc = ROOT.TFile("plots/h_bdt.root")
h_bdt_mc = f_bdt_mc.Get("h_bdt")

f_bdt_dataext = ROOT.TFile("plots/h_bdt_bnbext.root")
h_bdt_dataext = f_bdt_dataext.Get("h_bdt_bnbext")


f_bdt_intime = ROOT.TFile("plots/h_bdt_cosmic_mc.root")
h_bdt_intime = f_bdt_intime.Get("h_bdt_cosmic_mc")
# for i in range(h_bdt_data.GetNbinsX()):
#     h_bdt_data.SetBinContent(i,h_bdt_data.GetBinContent(i)-h_bdt_dataext.GetBinContent(i))

# h_bdt_data.Scale(h_bdt_mc.Integral()/h_bdt_data.Integral())

h_bdt_mc.GetHists()[3].SetFillStyle(3001)
h_bdt_dataext.SetLineColor(1)

legend = ROOT.TLegend(0.09455587, 0.7850208, 0.8923496, 0.9791956, "", "brNDC")
legend.SetTextSize(16)
legend.SetTextFont(63)
legend.SetHeader("MicroBooNE Preliminary %.1e POT" % total_pot)
legend.SetTextFont(43)
for j in range(h_bdt_mc.GetNhists()):
    if h_bdt_mc.GetHists()[j].Integral():
        legend.AddEntry(h_bdt_mc.GetHists()[j],
                        "%s: %.0f events" %
                        (description[j], h_bdt_mc.GetHists()[j].Integral()),
                        "f")

print(h_bdt_dataext.Integral() / h_bdt_intime.Integral())
h_bdt_intime.SetFillColor(ROOT.kOrange + 1)
h_bdt_intime.SetLineColor(ROOT.kBlack)

h_bdt_mc.Add(h_bdt_intime)
h_mc_err = h_bdt_mc.GetHists()[0].Clone()
h_mc_err.SetName("h_bdt_mc_err")

for i in range(1, h_bdt_mc.GetNhists()):
    h_mc_err.Add(h_bdt_mc.GetHists()[i])

h_mc_err.SetFillStyle(3002)
h_mc_err.SetFillColor(1)

legend.AddEntry(h_bdt_intime, "Cosmic in-time: %.0f events" %
                h_bdt_intime.Integral(), "f")

legend.AddEntry(h_bdt_data, "Data BNB: %.0f events" %
                h_bdt_data.Integral(), "lep")
#
# legend.AddEntry(h_bdt_dataext, "Data EXT: %.0f events" %
#                (h_bdt_dataext.Integral()), "f")
legend.SetNColumns(2)

pad_top = ROOT.TPad("pad_top", "", 0, 0.3, 1, 1)
pad_top.SetBottomMargin(0)
pad_top.Range(-22.21825, -2.003018, 202.5403, 2073.676)
pad_top.SetTopMargin(0.2411347)
pad_top.Draw()
pad_top.cd()

h_bdt_mc.Draw("hist")
h_mc_err.Draw("e2 same")

h_bdt_data.SetMarkerStyle(20)
h_bdt_data.SetLineColor(1)
h_bdt_mc.GetYaxis().SetTitleOffset(0.8)
h_bdt_mc.GetYaxis().SetTitleSize(0.06)

h_bdt_data.Draw("ep same")
legend.Draw()

line_bdt = ROOT.TLine(bdt_cut, 0, bdt_cut, h_mc_err.GetMaximum() * 1.6)
line_bdt.SetLineWidth(2)
line_bdt.SetLineStyle(2)
line_bdt.Draw()

pad_top.SetLogy()
c_bdt.cd()
pad_bottom = ROOT.TPad("pad_bottom", "", 0, 0, 1, 0.3)
pad_bottom.Range(-22.5, -0.6346511, 202.5, 1.99)

pad_bottom.SetFrameBorderMode(0)
pad_bottom.SetFrameBorderMode(0)
pad_bottom.SetBorderMode(0)
pad_bottom.SetTopMargin(0)
pad_bottom.SetBottomMargin(0.245614)
pad_bottom.Draw()
pad_bottom.cd()
h_ratio = h_bdt_data.Clone()

h_ratio.SetName("h_ratio")
h_ratio.GetYaxis().SetRangeUser(0.01, 1.99)
h_ratio.Divide(h_mc_err)
h_ratio.GetXaxis().SetLabelFont(42)
h_ratio.GetXaxis().SetLabelSize(0.13)
h_ratio.GetXaxis().SetTitleSize(0.13)
h_ratio.GetXaxis().SetTitleOffset(0.91)
h_ratio.GetYaxis().SetTitle("BNB/(MC+EXT)")
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


c_bdt.Update()
c_bdt.SaveAs("plots/bdt.pdf")

input()
