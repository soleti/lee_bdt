#!/usr/local/bin/python3

import ROOT

f_mc_electron = ROOT.TFile("f_dedx_hits_electron.root")
h_mc_electron = f_mc_electron.Get("h_dedx_hits_electron")

f_mc_photon = ROOT.TFile("f_dedx_hits_photon.root")
h_mc_photon = f_mc_photon.Get("h_dedx_hits_photon")

f_mc_other = ROOT.TFile("f_dedx_hits_other.root")
h_mc_other = f_mc_other.Get("h_dedx_hits_other")

f_data = ROOT.TFile("f_dedx_hits_databnb.root")
h_data = f_data.Get("h_dedx_hits_data")

f_dataext = ROOT.TFile("f_dedx_hits_dataext.root")
h_dataext = f_dataext.Get("h_dedx_hits_data")


print(h_data.Integral())
print(h_dataext.Integral())

h_data.Add(h_dataext, -1)
h_data.Scale(1 / h_data.Integral())

mc = ROOT.TObjArray(3)
mc.Add(h_mc_electron)
mc.Add(h_mc_photon)
mc.Add(h_mc_other)

fit = ROOT.TFractionFitter(h_data, mc)
fit.Constrain(1, 0.0, 1.0)

fit.Fit()

c_compare = ROOT.TCanvas("c_compare")


h_fit = fit.GetPlot()
electron = ROOT.Double(0)
electron_err = ROOT.Double(0)
photon = ROOT.Double(0)
photon_err = ROOT.Double(0)
other = ROOT.Double(0)
other_err = ROOT.Double(0)

fit.GetResult(0, electron, electron_err)
fit.GetResult(1, photon, photon_err)
fit.GetResult(2, other, other_err)

h_mc_electron.Scale(electron)
h_mc_photon.Scale(photon)
h_mc_other.Scale(other)

h_mc_electron.SetLineColor(ROOT.kBlack)
h_mc_photon.SetLineColor(ROOT.kBlack)
h_mc_other.SetLineColor(ROOT.kBlack)

h_mc_electron.SetFillColor(ROOT.kAzure + 1)
h_mc_photon.SetFillColor(ROOT.kRed - 4)
h_mc_other.SetFillColor(ROOT.kGray)


h_stack = ROOT.THStack("h_stack", ";dE/dx [cm];a.u.")

h_stack.Add(h_mc_other)
h_stack.Add(h_mc_photon)
h_stack.Add(h_mc_electron)
l_datamc = ROOT.TLegend(0.37, 0.65, 0.84, 0.85)
l_datamc.AddEntry(h_mc_electron, "Electron component", "f")
l_datamc.AddEntry(h_mc_photon, "Photon component", "f")
l_datamc.AddEntry(h_mc_other, "Other component", "f")
l_datamc.AddEntry(h_data, "Data", "lep")
h_stack.Draw("hist")
h_data.Draw("ep same")
l_datamc.Draw()

h_data.SetLineColor(ROOT.kBlack)
h_data.SetMarkerStyle(20)

c_compare.Update()

input()
