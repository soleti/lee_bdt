#!/usr/local/bin/python3

import ROOT
from bdt_common import colors, total_pot

ROOT.gStyle.SetOptStat(0)

OBJECTS = []

def apply_style(histo, color, style):
    histo.SetFillColor(color)
    histo.SetLineColor(color)
    histo.SetFillStyle(style)

def draw_top(*histos):
    pad_top = ROOT.TPad("pad_top", "", 0, 0.3, 1, 1)
    pad_top.SetBottomMargin(0)
    pad_top.SetTopMargin(0.2411347)
    pad_top.Draw()
    pad_top.cd()
    for h in histos:
        h.GetYaxis().SetRangeUser(0.001, 0.5)
        h.Draw("hist same")
        h.Draw("e1 same")
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

    h_ratio.SetName("h_ratio")

    h_ratio.GetYaxis().SetRangeUser(0.01, 1.99)
    h_ratio.Divide(den)

    h_ratio.SetFillStyle(0)
    h_ratio.SetLineColor(1)
    h_ratio.GetXaxis().SetLabelFont(42)
    h_ratio.GetXaxis().SetLabelSize(0.13)
    h_ratio.GetXaxis().SetTitleSize(0.13)
    h_ratio.GetXaxis().SetTitleOffset(0.91)
    h_ratio.GetYaxis().SetTitle("Tune 1 / Tune 2")
    h_ratio.GetYaxis().SetNdivisions(509)
    h_ratio.GetYaxis().SetLabelFont(42)
    h_ratio.GetYaxis().SetLabelSize(0.13)
    h_ratio.GetYaxis().SetTitleSize(0.13)
    h_ratio.GetYaxis().SetTitleOffset(0.36)

    h_ratio.Draw("hist")
    h_ratio.SetMarkerSize(0)
    h_ratio_2 = h_ratio.Clone()
    h_ratio_2.SetFillStyle(3002)
    h_ratio_2.SetFillColor(1)
    h_ratio_2.SetName("h_ratio_2")
    OBJECTS.append(h_ratio_2)
    h_ratio_2.Draw("e1 same")


    line = ROOT.TLine(h_ratio.GetXaxis().GetXmin(), 1,
                      h_ratio.GetXaxis().GetXmax(), 1)
    line.SetLineWidth(2)
    line.SetLineStyle(2)
    line.Draw()
    OBJECTS.append(line)


f_spectrum_normal = ROOT.TFile("genie/spectrum_normal.root")
spectrum_normal = f_spectrum_normal.Get("h_reco_energy_nu_e")
f_spectrum_tune2 = ROOT.TFile("genie/spectrum_tune2.root")
spectrum_tune2 = f_spectrum_tune2.Get("h_reco_energy_nu_e")

apply_style(spectrum_normal, ROOT.kGreen + 1, 3004)
apply_style(spectrum_tune2, ROOT.kAzure + 1, 3005)

f_spectrum_normal_nocuts = ROOT.TFile("genie/spectrum_normal_nocuts.root")
spectrum_normal_nocuts = f_spectrum_normal_nocuts.Get("h_reco_energy_nu_e")
f_spectrum_tune2_nocuts = ROOT.TFile("genie/spectrum_tune2_nocuts.root")
spectrum_tune2_nocuts = f_spectrum_tune2_nocuts.Get("h_reco_energy_nu_e")

apply_style(spectrum_normal_nocuts, ROOT.kGreen + 1, 3004)
apply_style(spectrum_tune2_nocuts, ROOT.kAzure + 1, 3005)

c = ROOT.TCanvas("c" , "", 900, 44, 700, 645)
c.cd()
legend = ROOT.TLegend(0.091, 0.7767, 0.34, 0.905, "", "brNDC")
legend.SetTextSize(16)
legend.SetTextFont(63)
legend.SetHeader("MicroBooNE Preliminary 4.9e19 POT")
legend.SetTextFont(43)
legend.AddEntry(spectrum_normal, "#nu_{e} CC0#pi GENIE tune 1: 3.2 events", "f")
legend.AddEntry(spectrum_tune2, "#nu_{e} CC0#pi GENIE tune 2: 2.9 events", "f")
draw_top(spectrum_normal, spectrum_tune2)
legend.Draw()
c.cd()

draw_ratio(spectrum_normal, spectrum_tune2)
c.Update()

c_nocuts = ROOT.TCanvas("c_nocuts", "", 900, 44, 700, 645)
c_nocuts.cd()
legend2 = ROOT.TLegend(0.091, 0.7767, 0.34, 0.905, "", "brNDC")
legend2.SetTextSize(16)
legend2.SetTextFont(63)
legend2.SetHeader("MicroBooNE Preliminary 6.6e20 POT")
legend2.SetTextFont(43)
legend2.AddEntry(spectrum_normal, "#nu_{e} CC0#pi-Np GENIE tune 1: 9.0 events", "f")
legend2.AddEntry(spectrum_tune2, "#nu_{e} CC0#pi-Np GENIE tune 2: 8.7 events", "f")
draw_top(spectrum_normal_nocuts, spectrum_tune2_nocuts)
legend2.Draw()
c_nocuts.cd()

draw_ratio(spectrum_normal_nocuts, spectrum_tune2_nocuts)
c_nocuts.Update()
input()
