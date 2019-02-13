#!/usr/local/bin/python3

import ROOT

f_mc_cuts = ROOT.TFile("plots_cuts/h_reco_energy_mc.root")
h_mc_cuts = f_mc_cuts.Get("h_reco_energy")
h_signal_cuts = h_mc_cuts.GetHists()[-1].Clone()
h_tot_cuts = h_mc_cuts.GetHists()[0].Clone()
for h in h_mc_cuts.GetHists()[1:]:
    h_tot_cuts.Add(h)

f_mc_precuts = ROOT.TFile("plots/h_reco_energy_mc.root")
h_mc_precuts = f_mc_precuts.Get("h_reco_energy")
h_signal_precuts = h_mc_precuts.GetHists()[-1].Clone()
h_tot_precuts = h_mc_precuts.GetHists()[0].Clone()
for h in h_mc_precuts.GetHists()[1:]:
    h_tot_precuts.Add(h)

f_mc_bdt = ROOT.TFile("plots_bdt/h_reco_energy_mc.root")
h_mc_bdt = f_mc_bdt.Get("h_reco_energy")
h_signal_bdt = h_mc_bdt.GetHists()[-1].Clone()
h_tot_bdt = h_mc_bdt.GetHists()[0].Clone()
for h in h_mc_bdt.GetHists()[1:]:
    h_tot_bdt.Add(h)

c_purity = ROOT.TCanvas("c_purity")
h_tot_cuts.GetYaxis().SetTitle("Purity")
pur_e_cuts = ROOT.TEfficiency(h_signal_cuts, h_tot_cuts)
pur_e_cuts.SetLineColor(ROOT.kOrange+1)
pur_e_cuts.SetLineWidth(2)

h_tot_bdt.GetYaxis().SetTitle("Purity")
pur_e_bdt = ROOT.TEfficiency(h_signal_bdt, h_tot_bdt)
pur_e_bdt.SetLineColor(ROOT.kMagenta+1)
pur_e_bdt.SetLineWidth(2)

h_tot_precuts.GetYaxis().SetTitle("Purity")
pur_e_precuts = ROOT.TEfficiency(h_signal_precuts, h_tot_precuts)
pur_e_precuts.SetLineColor(ROOT.kGreen+1)
pur_e_precuts.SetLineWidth(2)
pur_e_precuts.Draw("AP")
pur_e_cuts.Draw("P SAME")
pur_e_bdt.Draw("P SAME")
c_purity.Update()
pur_e_precuts.GetPaintedGraph().SetMinimum(0)
pur_e_precuts.GetPaintedGraph().SetMaximum(1)

c_purity.SetLeftMargin(0.13)
leg = ROOT.TLegend(0.16, 0.74, 0.46, 0.84)
# leg.AddEntry(eff_e, "Topology and flash requirements: %.1f %%" %
#                 (h_selected[0].Integral()/h_tot[0].Integral()*100), "le")
# leg.AddEntry(eff_e_numu, "#nu_{#mu} rejection: %.1f %%" % (
#     h_selected_numu[0].Integral()/h_tot[0].Integral()*100), "le")
leg.AddEntry(pur_e_precuts, "Quality precuts: %.1f%%" %
             (h_signal_precuts.Integral()/h_tot_precuts.Integral()*100), "le")
leg.AddEntry(pur_e_cuts, "Rectangular cuts: %.1f%%" %
             (h_signal_cuts.Integral()/h_tot_cuts.Integral()*100), "le")
leg.AddEntry(pur_e_bdt, "BDT: %.1f%%" %
             (h_signal_bdt.Integral()/h_tot_bdt.Integral()*100), "le")
leg.Draw()
c_purity.Update()
c_purity.SaveAs("plots/purity.pdf")

input()
