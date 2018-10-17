#!/usr/local/bin/python3

import ROOT
import os
import sys
import math
import pickle
from bdt_common import labels, binning, total_data_bnb_pot, pdg_colors, pdgs, inv_pdgs, draw_top, draw_ratio

if len(sys.argv) > 1:
    plot = sys.argv[1]
else:
    plot = ""

ROOT.gStyle.SetOptStat(0)

DRAW_SYS = True

h_pdgs = pickle.load(open("plots/pdg_plots.p", "rb"))

if __name__ == "__main__":
    OBJECTS = []

    for v in h_pdgs:
        if plot and plot != v:
            continue
        if h_pdgs[v] and ("track" in v or "shower" in v):
            h_stack = ROOT.THStack("h_stack_%s" % v, labels[v])
            OBJECTS.append(h_stack)
            h_tot_mc = ROOT.TH1F("h_tot_mc_%s" % v,
                                 labels[v],
                                 binning[v][0],
                                 binning[v][1],
                                 binning[v][2])
            h_tot_mc.SetTitle("Stat. uncertainty")

            l_pdg = ROOT.TLegend(0.09, 0.7747, 0.904, 0.984)
            l_pdg.SetNColumns(3)
            l_pdg.SetTextSize(16)
            l_pdg.SetTextFont(63)
            l_pdg.SetHeader("MicroBooNE Preliminary %.1e POT" %
                            total_data_bnb_pot)
            l_pdg.SetTextFont(43)
            for pdg in h_pdgs[v]:
                h_pdgs[v][pdg].SetLineWidth(0)
                h_pdgs[v][pdg].SetMarkerStyle(0)
                h_pdgs[v][pdg].SetMarkerSize(0)
                h_pdgs[v][pdg].SetFillColor(
                    ROOT.TColor.GetColor(pdg_colors[pdg]))
                if pdg != 2147483648:
                    h_stack.Add(h_pdgs[v][pdg])
                    h_tot_mc.Add(h_pdgs[v][pdg])

            for pdg in h_pdgs[v]:
                integral = h_pdgs[v][pdg].Integral()
                if integral > 0 and pdg != 2147483648:
                    l_pdg.AddEntry(h_pdgs[v][pdg],
                                    "%s: %.1f%%" % (
                        inv_pdgs[pdg], (integral / h_tot_mc.Integral()) * 100),
                        "f")

            h_tot_mc_clone = h_tot_mc.Clone()
            h_tot_mc.SetFillStyle(3002)
            h_tot_mc.SetFillColor(1)
            c = ROOT.TCanvas("c_%s" % v, v, 900, 44, 700, 645)
            # c.SetTopMargin(0.1978947)
            # c.Range(-0.4035779,-33.09252,6.713775,294.3855)
            draw_top(OBJECTS)

            h_stack.Draw("hist")

            if DRAW_SYS:
                h_mc_err_sys = h_tot_mc.Clone()
                fname_flux = "plots/sys/h_%s_flux_sys.root" % v
                fname_genie = "plots/sys/h_%s_genie_sys.root" % v
                if v == "track_energy_length":
                    fname_flux = fname_flux.replace("track", "total_track")
                    fname_genie = fname_genie.replace("track", "total_track")
                OBJECTS.append(h_mc_err_sys)
                if DRAW_SYS and os.path.isfile(fname_flux) and os.path.isfile(fname_genie):
                    var_name = v
                    if v == "track_energy_length":
                        var_name = v.replace("track", "total_track")
                    f_flux = ROOT.TFile(fname_flux)
                    h_flux = f_flux.Get("h_%s_cv" % var_name)
                    f_genie = ROOT.TFile(fname_genie)
                    h_genie = f_genie.Get("h_%s_cv" % var_name)

                    OBJECTS.append(h_flux)
                    OBJECTS.append(h_genie)

                    for k in range(1, h_mc_err_sys.GetNbinsX() + 1):
                        stat_err = h_mc_err_sys.GetBinError(k)
                        flux_err = h_flux.GetBinError(k)
                        genie_err = h_genie.GetBinError(k)
                        h_mc_err_sys.SetBinError(k, math.sqrt(
                            flux_err**2 + genie_err**2 + stat_err**2))
                    f_genie.Close()
                    f_flux.Close()

                h_mc_err_sys.SetLineColor(1)
                h_mc_err_sys.Draw("e2 same")
                OBJECTS.append(h_mc_err_sys)

            f = ROOT.TFile("plots/h_%s_bnb.root" % v)
            h = f.Get("h_%s" % v)
            h.SetLineColor(1)
            h.SetTitle("Data beam-on - beam-off")
            h.SetMarkerStyle(20)
            h_tot_mc.SetLineColor(1)
            h_tot_mc.SetLineWidth(2)
            h.Add(h_pdgs[v][2147483648], -1)
            h.Draw("e1p same")
            l_pdg.AddEntry(h, "Data (beam-on - beam-off)", "lep")
            h_mc_err_sys.SetLineWidth(0)
            if DRAW_SYS:
                l_pdg.AddEntry(h_mc_err_sys, "Sys. uncertainty", "f")

            OBJECTS.append(f)
            OBJECTS.append(h_tot_mc)
            if not DRAW_SYS:
                h_tot_mc.Draw("e2 same")
            h_tot_mc_clone.SetLineWidth(2)
            h_tot_mc_clone.SetLineColor(1)
            h_tot_mc_clone.Draw("hist same")
            l_pdg.Draw()
            if DRAW_SYS:
                chi2 = h.Chi2Test(h_mc_err_sys, "WW")
                ks = h.KolmogorovTest(h_mc_err_sys)
            else:
                chi2 = h.Chi2Test(h_tot_mc, "WW")
                ks = h.KolmogorovTest(h_tot_mc)

            p_test = ROOT.TPaveText(0.656, 0.587, 0.849, 0.695, "NDC")
            p_test.AddText("#chi^{2} prob. = %.2f" % chi2)
            p_test.AddText("K-S prob. = %.2f" % ks)
            p_test.SetFillStyle(0)
            p_test.SetBorderSize(0)
            p_test.SetTextAlign(11)
            p_test.Draw()

            OBJECTS.append(p_test)
            OBJECTS.append(h_tot_mc_clone)
            h_stack.SetMaximum(
                max(h_tot_mc.GetMaximum(), h.GetMaximum()) * 1.2)
            h_stack.GetYaxis().SetTitleOffset(0.9)
            c.cd()
            draw_ratio(h, h_tot_mc, h_mc_err_sys, OBJECTS)
            c.Update()
            c.SaveAs("plots/pdg/%s_pdg.pdf" % h_stack.GetName())
            OBJECTS.append(c)
            OBJECTS.append(l_pdg)

    input()
