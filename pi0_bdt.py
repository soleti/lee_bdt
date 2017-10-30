#!/usr/local/bin/python3

import ROOT

from bdt_common import binning, labels, variables, spectators, bdt_cut, manual_cuts
from bdt_common import x_start, x_end, y_start, y_end, z_start, z_end
from array import array
f_data = ROOT.TFile("pi0_file.root")
t_data = f_data.Get("pi0_tree")
ROOT.gStyle.SetOptStat(0)

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
    if n != "reco_energy":
        h = ROOT.TH1F("h_%s" % n, labels[n],
                      binning[n][0], binning[n][1], binning[n][2])
    else:
        bins = array("f", [0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.6, 0.8, 1])
        h = ROOT.TH1F("h_%s" % n, labels[n],len(bins)-1, bins)
    histograms.append(h)

histo_dict = dict(zip(variables_dict.keys(), histograms))

h_bdt = ROOT.TH1F("h_bdt_pi0", "BDT response; N. Entries / 0.05", 40, -1, 1)
bins = array("f", [0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.6, 0.8, 1])

h_reco_before = ROOT.TH1F("h_before",";Reco. energy [GeV]; N. Entries / 0.05 GeV", len(bins)-1, bins)
h_reco_after = ROOT.TH1F("h_after",";Reco. energy [GeV]; N. Entries / 0.05 GeV", len(bins)-1, bins)

print(t_data.GetEntries())

for i in range(t_data.GetEntries()):
    t_data.GetEntry(i)
    BDT_response = reader.EvaluateMVA("BDT method")
    h_bdt.Fill(BDT_response, t_data.event_weight)

    if BDT_response > bdt_cut and manual_cuts(t_data):
        for name, var in variables:
            histo_dict[name].Fill(var[0], t_data.event_weight)
        for name, var in spectators:
            histo_dict[name].Fill(var[0], t_data.event_weight)
        h_reco_after.Fill(t_data.reco_energy, t_data.event_weight)

    h_reco_before.Fill(t_data.reco_energy, t_data.event_weight)


f_bdt = ROOT.TFile("plots/h_bdt_pi0.root", "RECREATE")
h_bdt.Write()
f_bdt.Close()
for k in range(1,h_reco_before.GetNbinsX()+1):
    bin_width = h_reco_before.GetBinWidth(k)
    h_reco_before.SetBinError(k, h_reco_before.GetBinError(k)/(bin_width/0.05))
    h_reco_before.SetBinContent(k, h_reco_before.GetBinContent(k)/(bin_width/0.05))
    h_reco_after.SetBinError(k, h_reco_after.GetBinError(k)/(bin_width/0.05))
    h_reco_after.SetBinContent(k, h_reco_after.GetBinContent(k)/(bin_width/0.05))

h_reco_after.Scale(45/5)
h_reco_before.Scale(66/5)

legend_before = ROOT.TLegend(0.62, 0.76, 0.87, 0.86, "", "brNDC")
legend_before.AddEntry(h_reco_before, "NC #pi^{0} events: %.0f" % h_reco_before.Integral(), "f")
legend_after = ROOT.TLegend(0.62, 0.76, 0.87, 0.86, "", "brNDC")
legend_after.AddEntry(h_reco_before, "NC #pi^{0} events: %.0f" % h_reco_after.Integral(), "f")
h_reco_before.SetLineColor(ROOT.kBlack)
h_reco_after.SetLineColor(ROOT.kBlack)
h_reco_before.SetFillColor(ROOT.kBlue-7)
h_reco_after.SetFillColor(ROOT.kBlue-7)
c = ROOT.TCanvas("c")
h_reco_before.Draw("hist")
legend_before.Draw()
c.Update()

c2 = ROOT.TCanvas("c2")
h_reco_after.Draw("hist")
legend_after.Draw("hist")
c2.Update()

for h in histograms:
    f = ROOT.TFile("plots/%s_pi0.root" % h.GetName(), "RECREATE")
    h.Write()
    f.Close()

input()
