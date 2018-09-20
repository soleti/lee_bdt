#!/usr/local/bin/python3

import ROOT
import collections
from bdt_common import bdt_cut, binning, labels, variables, spectators, fill_histos_data
from bdt_common import manual_cuts, bdt, manual, bins, colors, inv_interactions, interactions
from bdt_common import pre_cuts, rectangular_cut, fix_binning, total_data_bnb_pot, description2
from array import array
import math
import collections
import sys
import os.path


ROOT.gStyle.SetPalette(ROOT.kCMYK)
print("BDT: ", bdt, "Manual: ", manual)

DRAW_SYS = True

pdg_colors = {
    2147483648: "#ffffff",
    22: '#e6194b',
    11: '#f58231',
    2112: '#ffe119',
    13: '#bfef45',
    211: '#3cb44b',
    2212: '#42d4f4',
    999999: '#4363db',
    99999: '#a9a9a9'
}

pdgs = [
    ("Data off-beam", 2147483648),
    ("#gamma", 22),
    ("e^{#pm}", 11),
    ("n", 2112),
    ("#mu^{#pm}", 13),
    ("#pi^{#pm}", 211),
    ("p", 2212),
    ("Overlay cosmic", 999999),
    ("Other", 99999),
]

pdgs = collections.OrderedDict(pdgs)

inv_pdgs = {v: k for k, v in pdgs.items()}


h_pdgs = {}
vars = dict(variables + spectators)

for n, b in vars.items():
    pdg_dict = {}
    for pdg in pdgs:
        pdg_dict[pdgs[pdg]] = ROOT.TH1F("h_%s_%s" % (n, pdg), pdg+labels[n], binning[n][0], binning[n][1], binning[n][2])
    h_pdgs[n] = pdg_dict

h_int = {}
for name in interactions:
    h_int[interactions[name]] = ROOT.TH1F(name, ";E_{corr};N. Entries / 0.05 GeV", len(bins) - 1, bins)

h_angle_energy_bkg = ROOT.TH2F("h_angle_energy_bkg",
                               ";Shower opening angle [#circ];Shower energy [GeV]",
                               45, 0, 45, 20, 0.05, 0.5)
h_angle_energy_sig = ROOT.TH2F("h_angle_energy_sig",
                               ";Shower opening angle [#circ];Shower energy [GeV]",
                               45, 0, 45, 20, 0.05, 0.5)

def fill_histos(chain, histo_dict, h_bdts, option=""):
    ROOT.TMVA.Tools.Instance()
    reader = ROOT.TMVA.Reader(":".join([
        "!V",
        "!Silent",
        "Color"]))

    for name, var in variables:
        chain.SetBranchAddress(name, var)

    for name, var in spectators:
        chain.SetBranchAddress(name, var)

    for name, var in variables:
        reader.AddVariable(name, var)

    for name, var in spectators:
        reader.AddSpectator(name, var)

    reader.BookMVA("BDT method",
                   "dataset/weights/TMVAClassification_BDT.weights.xml")
    # reader.BookMVA("Likelihood method",
    #                "dataset/weights/TMVAClassification_Likelihood.weights.xml")
    # reader.BookMVA("Cuts method",
    #                 "dataset/weights/TMVAClassification_Cuts.weights.xml")

    correction_factor = 1
    h_energy_sig = ROOT.TH1F("h_energy_sig_%s" % option, ";E_{#nu_{e}};N. Entries / 0.05 GeV", len(bins) - 1, bins)
    h_energy_bkg = ROOT.TH1F("h_energy_bkg_%s" % option, ";E_{#nu_{e}};N. Entries / 0.05 GeV", len(bins) - 1, bins)
    h_reco_sig = ROOT.TH1F("h_reco_sig_%s" % option,
                           ";Reco. energy [GeV];N. Entries / 0.05 GeV", len(bins) - 1, bins)
    h_reco_bkg = ROOT.TH1F("h_reco_bkg_%s" % option,
                           ";Reco. energy [GeV];N. Entries / 0.05 GeV", len(bins) - 1, bins)


    passed_events = 0
    entries = int(chain.GetEntries() / 1)
    for i in range(entries):
        chain.GetEntry(i)
        category = int(chain.category)

        BDT_response = reader.EvaluateMVA("BDT method")
        # likelihood_response = reader.EvaluateMVA("Likelihood method")
        # cuts_response = reader.EvaluateMVA("Cuts method", rectangular_cut)

        if pre_cuts(chain):
            h_bdts[category].Fill(BDT_response, chain.event_weight)
        else:
            continue

        if bdt:
            apply_bdt = BDT_response > bdt_cut
        else:
            apply_bdt = True

        if manual:
            apply_manual = manual_cuts(chain)
        else:
            apply_manual = True
        if apply_bdt and apply_manual:
            passed_events += chain.event_weight
            if not chain.true_nu_is_fidvol and category != 0 and category != 6 and category != 1 and category != 7:
                category = 5

            corr = 1
            # if category == 2 and not manual:
            #     corr = 1.2

            corrected_energy_cali = ((chain.total_shower_energy_cali + 1.36881e-02) /
                                     7.69908e-01) + chain.total_track_energy_length

            corrected_energy_hits = ((chain.total_shower_energy_cali + 1.36881e-02) /
                                     7.69908e-01) + (chain.total_track_energy + 3.57033e-02) / 7.70870e-01

            if chain.true_nu_is_fidvol:
                h_energy_sig.Fill(chain.nu_E, chain.event_weight * corr)

            h_reco_sig.Fill(corrected_energy_cali, chain.event_weight * corr)

            if (option == "mc" or option == "bnbext" or (option == "nue" and category != 2)):
                if chain.true_nu_is_fidvol:
                    h_energy_bkg.Fill(chain.nu_E, chain.event_weight * corr)
                h_reco_bkg.Fill(corrected_energy_hits, chain.event_weight * corr)


            var_dict = dict(variables + spectators)
            for name, var in variables:
                for i, v in enumerate(var):
                    if v > -999:
                        if "track" in name:
                            pdg_code = int(var_dict["track_pdg"][i])
                            if option != "bnbext" and abs(pdg_code) == 2147483648:
                                pdg_code = int(var_dict["shower_pdg"][0])
                        elif "shower" in name:
                            pdg_code = int(var_dict["shower_pdg"][i])
                            if option != "bnbext" and abs(pdg_code) == 2147483648:
                                pdg_code = int(var_dict["track_pdg"][0])

                        if option != "bnbext" and abs(pdg_code) == 2147483648:
                            pdg_code = 99999

                        if abs(pdg_code) in h_pdgs[name]:
                            h_pdgs[name][abs(pdg_code)].Fill(v, chain.event_weight)
                        else:
                            h_pdgs[name][99999].Fill(v, chain.event_weight)

                        histo_dict[name][category].Fill(v, chain.event_weight)

            for name, var in spectators:
                for i, v in enumerate(var):
                    if v > -999:
                        if "track" in name:
                            pdg_code = int(var_dict["track_pdg"][i])
                            if option != "bnbext" and abs(pdg_code) == 2147483648:
                                pdg_code = int(var_dict["shower_pdg"][0])
                        elif "shower" in name:
                            pdg_code = int(var_dict["shower_pdg"][i])
                            if option != "bnbext" and abs(pdg_code) == 2147483648:
                                pdg_code = int(var_dict["track_pdg"][0])


                        if option != "bnbext" and abs(pdg_code) == 2147483648:
                            pdg_code = 99999
                        if abs(pdg_code) in h_pdgs[name]:
                            h_pdgs[name][abs(pdg_code)].Fill(v, chain.event_weight)
                        else:
                            h_pdgs[name][99999].Fill(v, chain.event_weight)

                        histo_dict[name][category].Fill(v, chain.event_weight)

            # if chain.category != 2 and 0 < chain.reco_energy < 2:
            #     print("nu_mu",
            #           int(chain.run), int(chain.subrun), int(chain.event),
            #           int(chain.category),
            #           inv_interactions[int(chain.interaction_type)])

    if manual or bdt:
        f_energy_binned_selected = ROOT.TFile("plots/h_energy_%s_after.root" % option, "RECREATE")
        h_energy_sig.Write()
        h_energy_bkg.Write()
        h_reco_sig.Write()
        h_reco_bkg.Write()
        f_energy_binned_selected.Close()
    else:
        f_energy_binned_selected = ROOT.TFile("plots/h_energy_%s.root" % option, "RECREATE")
        h_energy_sig.Write()
        h_energy_bkg.Write()
        h_reco_sig.Write()
        h_reco_bkg.Write()
        f_energy_binned_selected.Close()

    return passed_events



print("LEE events", fill_histos_data("lee", bdt, manual))
print("Data events", fill_histos_data("bnb", bdt, manual))


mc_chain = ROOT.TChain("mc_tree")
nue_chain = ROOT.TChain("nue_tree")
bnbext_chain = ROOT.TChain("bnbext_tree")

mc_chain.Add("root_files_overlay/mc_file.root")
nue_chain.Add("root_files_overlay/nue_file.root")
bnbext_chain.Add("root_files/bnbext_file.root")

kinds = []
categories = ["intime", "cosmic", "nu_e", "nu_mu", "nc", "dirt", "data",
              "mixed", "cc0pi"]
stacked_histos = []

variables_dict = dict(variables + spectators)

for i, n in enumerate(variables_dict.keys()):
    histos = []

    h_stack = ROOT.THStack("h_" + n, labels[n])

    for c in categories:
        if n != "reco_energy" and n != "nu_E":
            h = ROOT.TH1F("h_%s_%s" % (n, c), labels[n],
                          binning[n][0], binning[n][1], binning[n][2])
        else:
            h = ROOT.TH1F("h_%s_%s" % (n, c), labels[n], len(bins) - 1, bins)
        histos.append(h)

    h_stack.Add(histos[0])
    h_stack.Add(histos[3])
    h_stack.Add(histos[4])
    h_stack.Add(histos[5])
    h_stack.Add(histos[7])
    h_stack.Add(histos[1])
    h_stack.Add(histos[6])
    h_stack.Add(histos[8])
    h_stack.Add(histos[2])
    stacked_histos.append(h_stack)
    kinds.append(histos)


histo_dict = dict(zip(variables_dict.keys(), kinds))

h_bdt_stack = ROOT.THStack("h_bdt", ";BDT response; N. Entries / 0.05")
h_bdts = []

for i, c in enumerate(categories):
    h = ROOT.TH1F("h_bdt_%s" % c, "BDT response; N. Entries / 0.05", 20, -1, 1)
    h.SetLineColor(1)
    h.SetFillColor(colors[i])
    if i != 0:
        h.SetLineWidth(0)
    h_bdts.append(h)

h_bdt_stack.Add(h_bdts[0])
h_bdt_stack.Add(h_bdts[3])
h_bdt_stack.Add(h_bdts[4])
h_bdt_stack.Add(h_bdts[5])
h_bdt_stack.Add(h_bdts[7])
h_bdt_stack.Add(h_bdts[1])
h_bdt_stack.Add(h_bdts[6])
h_bdt_stack.Add(h_bdts[8])
h_bdt_stack.Add(h_bdts[2])
print("nu_e events", fill_histos(nue_chain, histo_dict, h_bdts, "nue"))
print("BNB + cosmic events", fill_histos(mc_chain, histo_dict, h_bdts, "mc"))
print("EXT events", fill_histos(bnbext_chain, histo_dict, h_bdts, "bnbext"))

h_interactions = ROOT.THStack(
    "h_interactions", ";E_{corr};N. Entries / 0.05 GeV")
l_interactions = ROOT.TLegend(0.13467, 0.7747, 0.9255, 0.8968, "", "brNDC")

l_interactions.SetNColumns(3)
for interaction in h_int:
    histo = h_int[interaction]
    if histo.Integral() > 0:
        histo.SetLineColor(1)
        h_interactions.Add(histo)

h_true_e = ROOT.THStack(
    "h_true_e", ";E_{corr} [GeV]; N. Entries / 0.05 GeV")

for j in range(h_interactions.GetNhists()):
    h_clone = h_interactions.GetHists()[j].Clone()
    h_clone.Scale(1)
    fix_binning(h_clone)

    h_fixed = ROOT.TH1F("h_fixed%i" % j, "", len(bins) - 1, 0, len(bins) - 1)

    for i in range(1, h_clone.GetNbinsX() + 1):
        h_fixed.SetBinContent(i, h_clone.GetBinContent(i))
        h_fixed.SetBinError(i, h_clone.GetBinError(i))
        h_fixed.GetXaxis().SetBinLabel(i, "")

    h_fixed.SetLineColor(ROOT.kBlack)
    l_interactions.AddEntry(h_fixed, "%s: %.1f events" % (
        h_clone.GetName(), h_interactions.GetHists()[j].Integral()), "f")

    h_true_e.Add(h_fixed)

f_bdt = ROOT.TFile("plots/h_bdt.root", "RECREATE")
h_bdt_stack.Write()
f_bdt.Close()

for i, histos in enumerate(kinds):
    for j, h in enumerate(histos):
        h.SetLineColor(1)
        if j != 0:
            h.SetLineWidth(0)
        h.SetFillColor(colors[j])

for h in stacked_histos:
    f = ROOT.TFile("plots/%s_mc.root" % h.GetName(), "RECREATE")
    h.Write()
    f.Close()

OBJECTS = []
for v in h_pdgs:
    if h_pdgs[v] and ("track" in v or "shower" in v):
        h_stack = ROOT.THStack("h_stack_%s" % v, labels[v])
        OBJECTS.append(h_stack)
        h_tot_mc = ROOT.TH1F("h_tot_mc_%s" % v,
                             labels[v],
                             binning[v][0],
                             binning[v][1],
                             binning[v][2])
        h_tot_mc.SetTitle("Stat. uncertainty")

        l_pdg = ROOT.TLegend(0.09025788,0.8105263,0.9040115,0.9852632)
        l_pdg.SetNColumns(3)
        l_pdg.SetTextSize(16)
        l_pdg.SetTextFont(63)
        l_pdg.SetHeader("MicroBooNE Preliminary %.1e POT" % total_data_bnb_pot)
        l_pdg.SetTextFont(43)
        for pdg in h_pdgs[v]:
            h_pdgs[v][pdg].SetLineWidth(0)
            h_pdgs[v][pdg].SetMarkerStyle(0)
            h_pdgs[v][pdg].SetMarkerSize(0)
            h_pdgs[v][pdg].SetFillColor(ROOT.TColor.GetColor(pdg_colors[pdg]))
            if pdg != 2147483648:
                h_stack.Add(h_pdgs[v][pdg])
                h_tot_mc.Add(h_pdgs[v][pdg])

        for pdg in h_pdgs[v]:
            integral = h_pdgs[v][pdg].Integral()
            if integral > 0 and pdg != 2147483648:
                l_pdg.AddEntry(h_pdgs[v][pdg],
                               "%s: %.1f%%" % (inv_pdgs[pdg], (integral / h_tot_mc.Integral()) * 100),
                               "f")

        h_tot_mc_clone = h_tot_mc.Clone()
        h_tot_mc.SetFillStyle(3002)
        h_tot_mc.SetFillColor(1)
        c = ROOT.TCanvas("c_%s" % v)
        c.SetTopMargin(0.1978947)
        c.Range(-0.4035779,-33.09252,6.713775,294.3855)

        h_stack.Draw("hist")

        if DRAW_SYS:
            h_mc_err_sys = h_tot_mc.Clone()
            fname_flux = "plots/sys/h_%s_flux_sys.root" % v
            fname_genie = "plots/sys/h_%s_genie_sys.root" % v
            OBJECTS.append(h_mc_err_sys)
            if DRAW_SYS and os.path.isfile(fname_flux) and os.path.isfile(fname_genie):
                f_flux = ROOT.TFile(fname_flux)
                h_flux = f_flux.Get("h_%s_cv" % v)
                f_genie = ROOT.TFile(fname_genie)
                h_genie = f_genie.Get("h_%s_cv" % v)

                OBJECTS.append(h_flux)
                OBJECTS.append(h_genie)

                for k in range(1, h_mc_err_sys.GetNbinsX() + 1):
                    stat_err = h_mc_err_sys.GetBinError(k)
                    flux_err = h_flux.GetBinError(k)
                    genie_err = h_genie.GetBinError(k)
                    h_mc_err_sys.SetBinError(
                        k, math.sqrt(flux_err**2 + genie_err**2 + stat_err**2))
                f_genie.Close()
                f_flux.Close()

            h_mc_err_sys.SetLineColor(1)
            h_mc_err_sys.Draw("e1p same")
            OBJECTS.append(h_mc_err_sys)

        f = ROOT.TFile("plots/h_%s_bnb.root" % v)
        h = f.Get("h_%s" % v)
        h.SetLineColor(1)
        h.SetTitle("Data beam-on - beam-off")
        h.SetMarkerStyle(20)
        h_tot_mc.SetLineColor(1)
        h_tot_mc.SetLineWidth(2)
        h.Add(h_pdgs[v][2147483648], -1)
        h.Draw("ep same")
        l_pdg.AddEntry(h, "Data (beam-on - beam-off)", "lep")

        OBJECTS.append(f)
        OBJECTS.append(h_tot_mc)
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

        p_test = ROOT.TPaveText(0.677, 0.642, 0.871, 0.751, "NDC")
        p_test.AddText("#chi^{2} prob. = %.2f" % chi2)
        p_test.AddText("K-S prob. = %.2f" % ks)
        p_test.SetFillStyle(0)
        p_test.SetBorderSize(0)
        p_test.SetTextAlign(11)
        p_test.Draw()

        OBJECTS.append(p_test)
        OBJECTS.append(h_tot_mc_clone)
        h_stack.SetMaximum(max(h_tot_mc.GetMaximum(), h.GetMaximum()) * 1.2)
        h_stack.GetYaxis().SetTitleOffset(0.9)
        c.Update()
        c.SaveAs("plots/pdg/%s_pdg.pdf" % h_stack.GetName())
        OBJECTS.append(c)
        OBJECTS.append(l_pdg)
input()
