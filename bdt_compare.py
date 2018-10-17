#!/usr/local/bin/python3

import sys
import ROOT
from bdt_common import description, bdt_cut_cosmic, total_pot

if len(sys.argv) > 1:
    bdt_type = sys.argv[1]
else:
    bdt_type = ""

ROOT.gStyle.SetOptStat(0)

f_bdt_data = ROOT.TFile("plots/h_bdt_bnb.root")
h_bdt_data = f_bdt_data.Get("h_bdt_%s_bnb" % bdt_type)

f_bdt_mc = ROOT.TFile("plots/h_bdt.root")
h_bdt_mc = f_bdt_mc.Get("h_bdt_%s" % bdt_type)

f_bdt_lee = ROOT.TFile("plots/h_bdt_lee.root")
h_bdt_lee = f_bdt_lee.Get("h_bdt_%s_lee" % bdt_type)

h_bdt_lee.SetFillColor(ROOT.kGreen - 10)
h_bdt_mc.Add(h_bdt_lee)

h_bdt_mc.GetHists()[4].SetFillStyle(3001)
h_bdt_mc.GetHists()[0].SetFillStyle(3004)

c_bdt = ROOT.TCanvas("c_bdt", "BDT response", 900, 44, 700, 645)


legend = ROOT.TLegend(0.09455587, 0.7850208, 0.8923496, 0.9791956, "", "brNDC")
legend.SetTextSize(16)
legend.SetTextFont(63)
legend.SetHeader("MicroBooNE Preliminary %.1e POT" % total_pot)
legend.SetTextFont(43)
for j in range(h_bdt_mc.GetNhists()-1):
    if h_bdt_mc.GetHists()[j].Integral():
        legend.AddEntry(h_bdt_mc.GetHists()[j],
                        "%s: %.1f events" %
                        (description[j], h_bdt_mc.GetHists()[j].Integral()),
                        "f")

h_mc_err = h_bdt_mc.GetHists()[0].Clone()
h_mc_err.Clear()
h_mc_err.SetName("h_bdt_mc_err")

for i in range(1, h_bdt_mc.GetNhists()):
    h_mc_err.Add(h_bdt_mc.GetHists()[i])

h_mc_err.SetFillStyle(3002)
h_mc_err.SetFillColor(1)

legend.AddEntry(h_bdt_data, "Data BNB: %.0f events" %
                h_bdt_data.Integral(), "lep")

legend.SetNColumns(2)

h_bdt_mc.Draw("hist")
h_bdt_mc.GetYaxis().SetRangeUser(0.5, h_bdt_mc.GetMaximum() * 2)
h_bdt_mc.GetXaxis().SetRangeUser(-0.5, 0.5)
h_bdt_mc.GetYaxis().SetTitleOffset(0.9)
h_mc_err.Draw("e2 same")
h_mc_err_clone = h_mc_err.Clone()
h_mc_err_clone.SetLineColor(1)
h_mc_err_clone.SetFillStyle(0)

h_mc_err_clone.Draw("hist same")

h_bdt_data.SetMarkerStyle(20)
h_bdt_data.SetLineColor(1)

h_bdt_data.Draw("ep same")
legend.AddEntry(h_mc_err, "Stat. error", "lf")
legend.Draw()

c_bdt.SetTopMargin(0.22)

c_bdt.SetLogy()
c_bdt.Update()
c_bdt.SaveAs("plots/bdt.pdf")

input()

f_bdt_data.Close()
f_bdt_lee.Close()
f_bdt_mc.Close()
