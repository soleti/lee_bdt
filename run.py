#!/usr/local/bin/python3

import ROOT

from bdt_common import bdt_cut, binning, labels, variables, spectators, fill_histos_data
from bdt_common import manual_cuts, bdt, manual, bins, colors, inv_interactions, interactions
from bdt_common import pre_cuts, rectangular_cut, fix_binning, total_data_bnb_pot
from array import array

ROOT.gStyle.SetPalette(ROOT.kBird)
print("BDT: ", bdt, "Manual: ", manual)

h_int = {}
for name in interactions:
    h_int[interactions[name]] = ROOT.TH1F(name, ";E_{corr};N. Entries / 0.05 GeV", len(bins) - 1, bins)


h_angle_energy_bkg = ROOT.TH2F("h_angle_energy_bkg", ";Shower #theta [#circ];Shower energy [GeV]", 20, 0, 180, 20, 0.05, 1)
h_angle_energy_sig = ROOT.TH2F("h_angle_energy_sig", ";Shower #theta [#circ];Shower energy [GeV]", 20, 0, 180, 20, 0.05, 1)

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

    nue_file = open("nue_passed.txt", "w")
    numu_file = open("numu_passed.txt", "w")

    passed_events = 0

    for i in range(chain.GetEntries()):
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
            if category == 2 and not manual:
                corr = 1/1.2
            corrected_energy = chain.total_shower_energy / 0.78 + 0.02 + chain.total_track_energy_length


            if category == 2:
                h_angle_energy_sig.Fill(chain.shower_theta, chain.shower_energy, chain.event_weight * corr)
            elif category != 0 and category != 1:
                h_angle_energy_bkg.Fill(chain.shower_theta, chain.shower_energy, chain.event_weight * corr)
            if category != 0 and category != 6 and category != 1 and category != 7 and category == 8:
                h_int[chain.interaction_type].Fill(corrected_energy, chain.event_weight * corr)


            if category == 2 and option == "nue":
                print(int(chain.run), int(chain.subrun), int(chain.event), file=nue_file)

                if chain.true_nu_is_fidvol:
                    h_energy_sig.Fill(chain.nu_E, chain.event_weight * corr)

                h_reco_sig.Fill(corrected_energy, chain.event_weight * corr)

            if (option == "mc" or option == "bnbext" or (option == "nue" and category != 2)):
                if chain.true_nu_is_fidvol:
                    h_energy_bkg.Fill(chain.nu_E, chain.event_weight * corr)
                h_reco_bkg.Fill(corrected_energy, chain.event_weight * corr)

            for name, var in variables:
                histo_dict[name][category].Fill(var[0], chain.event_weight)

            for name, var in spectators:

                if name == "reco_energy":
                    histo_dict[name][category].Fill(corrected_energy, chain.event_weight * corr)
                else:
                    histo_dict[name][category].Fill(var[0], chain.event_weight * corr)

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
        f_energy_binned_selected = ROOT.TFile( "plots/h_energy_%s.root" % option, "RECREATE")
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

mc_chain.Add("mc_file.root")
nue_chain.Add("nue_file.root")
bnbext_chain.Add("bnbext_file.root")

kinds = []
categories = ["intime", "cosmic", "nu_e", "nu_mu", "nc", "dirt", "data",
              "mixed", "cc0pi"]
stacked_histos = []

variables_dict = dict(variables + spectators)

for i, n in enumerate(variables_dict.keys()):
    histos = []

    h_stack = ROOT.THStack("h_" + n, labels[n])

    for c in categories:
        if n != "reco_energy":
            h = ROOT.TH1F("h_%s_%s" % (n, c), labels[n],
                          binning[n][0], binning[n][1], binning[n][2])
        else:
            h = ROOT.TH1F("h_%s_%s" % (n, c), labels[n], len(bins) - 1, bins)
        histos.append(h)

    h_stack.Add(histos[2])
    h_stack.Add(histos[8])
    h_stack.Add(histos[3])
    h_stack.Add(histos[4])
    h_stack.Add(histos[5])
    h_stack.Add(histos[7])
    h_stack.Add(histos[1])
    h_stack.Add(histos[0])
    h_stack.Add(histos[6])

    stacked_histos.append(h_stack)
    kinds.append(histos)



histo_dict = dict(zip(variables_dict.keys(), kinds))

h_bdt_stack = ROOT.THStack("h_bdt", ";BDT response; N. Entries / 0.05")
h_bdts = []

for i, c in enumerate(categories):
    h = ROOT.TH1F("h_bdt_%s" % c, "BDT response; N. Entries / 0.05", 40, -1, 1)
    h.SetLineColor(1)
    h.SetFillColor(colors[i])
    h.SetLineWidth(2)
    h_bdts.append(h)

h_bdt_stack.Add(h_bdts[2])
h_bdt_stack.Add(h_bdts[8])
h_bdt_stack.Add(h_bdts[3])
h_bdt_stack.Add(h_bdts[4])
h_bdt_stack.Add(h_bdts[5])
h_bdt_stack.Add(h_bdts[7])
h_bdt_stack.Add(h_bdts[1])
h_bdt_stack.Add(h_bdts[0])
h_bdt_stack.Add(h_bdts[6])
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


c_interactions = ROOT.TCanvas("c_interactions", "")
c_interactions.SetLeftMargin(0.14)
h_true_e.Draw("pfc hist")
l_interactions.Draw("same")
ax = ROOT.TGaxis(0, 0, len(bins) - 1, 0, 0, len(bins) - 1, 515, "")
for i, i_bin in enumerate(bins):
    ax.ChangeLabel(i + 1, -1, -1, -1, -1, -1,
                   "{0}".format(str(round(i_bin, 3) if i_bin % 1 else int(i_bin))))
ax.SetLabelFont(42)
ax.SetLabelSize(0.04)
ax.Draw()
pt = ROOT.TPaveText(0.117, 0.89, 0.801, 0.975, "ndc")
pt.AddText("MicroBooNE Preliminary %.1e POT - #nu_{e} CC" % total_data_bnb_pot)
pt.SetFillColor(0)
pt.SetBorderSize(0)
pt.SetShadowColor(0)
pt.Draw()
c_interactions.SetTopMargin(0.2411347)
c_interactions.Update()


f_bdt = ROOT.TFile("plots/h_bdt.root", "RECREATE")
h_bdt_stack.Write()
f_bdt.Close()

for i, histos in enumerate(kinds):
    for j, h in enumerate(histos):
        h.SetLineColor(1)
        h.SetLineWidth(2)
        h.SetFillColor(colors[j])

for h in stacked_histos:
    f = ROOT.TFile("plots/%s_mc.root" % h.GetName(), "RECREATE")
    h.Write()
    f.Close()

f_data = ROOT.TFile("bnb_file.root")
t_data = f_data.Get("bnb_tree")

for name, var in variables:
    t_data.SetBranchAddress(name, var)

for name, var in spectators:
    t_data.SetBranchAddress(name, var)

ROOT.gStyle.SetOptStat(0)

# c_angle = ROOT.TCanvas("c_angle")
# c_angle.Divide(2)
# c_angle.cd(1)
# h_angle_energy_bkg.SetMinimum(-0.001)
# h_angle_energy_bkg.Draw("col")
# c_angle.cd(2)
# h_angle_energy_sig.SetMinimum(-0.001)
# h_angle_energy_sig.Draw("col")
# c_angle.Update()
input()
