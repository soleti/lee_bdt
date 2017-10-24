#!/usr/local/bin/python3

import ROOT

categories = ["other", "electron", "muon", "hadron", "photon", "pion"]

colors = [ROOT.kGray, ROOT.kRed - 4, ROOT.kAzure + 1,
          ROOT.kOrange + 1, ROOT.kGreen + 1, ROOT.kCyan - 9]

files = []
total = 0
h_stack = ROOT.THStack("h_stack", ";dE/dx [cm];a.u.")
l_datamc = ROOT.TLegend(0.54, 0.51, 0.84, 0.85)

for i, c in enumerate(categories):
    f = ROOT.TFile("f_dedx_hits_%s.root" % c)
    files.append(f)
    h = files[i].Get("h_dedx_hits_%s" % c)
    h.SetLineColor(ROOT.kBlack)
    h.SetFillColor(colors[i])
    print(h.Integral())
    l_datamc.AddEntry(h, "%s: %.1f %%" % (h.GetTitle(),
                                          h.Integral() * 100 / 197306.0), "f")
    total += h.Integral()
    h_stack.Add(h)
print(total)
f_data = ROOT.TFile("f_dedx_hits_databnb.root")
h_data = f_data.Get("h_dedx_hits_data")

f_dataext = ROOT.TFile("f_dedx_hits_dataext.root")
h_dataext = f_dataext.Get("h_dedx_hits_data")

h_data.Add(h_dataext, -1)
h_data.Scale(total / h_data.Integral())


c_compare = ROOT.TCanvas("c_compare")

l_datamc.AddEntry(h_data, "Data BNB - EXT", "lep")
h_stack.Draw("hist")
h_data.Draw("ep same")
l_datamc.Draw()

h_data.SetLineColor(ROOT.kBlack)
h_data.SetMarkerStyle(20)

c_compare.Update()

input()
