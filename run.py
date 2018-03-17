#!/usr/local/bin/python3

import ROOT

from bdt_common import bdt_cut, binning, labels, variables, spectators, fill_histos_data
from bdt_common import sigmaCalc, manual_cuts, bdt, manual, bins, colors
from bdt_common import pre_cuts, rectangular_cut

print("BDT: ", bdt, "Manual: ", manual)


def fill_histos(chain, histo_dict, h_bdts):
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
            for name, var in variables:
                histo_dict[name][category].Fill(var[0], chain.event_weight)

            for name, var in spectators:
                if name == "reco_energy":
                    # histo_dict[name][category].Fill(correction_factor * (0.917*chain.total_track_energy_length + chain.total_shower_energy_y), chain.event_weight)
                    histo_dict[name][category].Fill(chain.reco_energy, chain.event_weight)

                elif name == "total_shower_energy_y":
                    histo_dict[name][category].Fill(var[0], chain.event_weight)
                elif name == "total_track_energy_length":
                    histo_dict[name][category].Fill(var[0], chain.event_weight)
                else:
                    histo_dict[name][category].Fill(var[0], chain.event_weight)

                
            # if chain.category != 2 and 0 < chain.reco_energy < 2:
            #     print("nu_mu",
            #           int(chain.run), int(chain.subrun), int(chain.event),
            #           int(chain.category),
            #           inv_interactions[int(chain.interaction_type)])

print("LEE events", fill_histos_data("lee", bdt, manual))


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

fill_histos(mc_chain, histo_dict, h_bdts)
fill_histos(nue_chain, histo_dict, h_bdts)
fill_histos(bnbext_chain, histo_dict, h_bdts)

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

# variables_dict = dict(variables + spectators)
# histo2D = []
# for i, n in enumerate(variables_dict.keys()):
#     h = ROOT.TH2F("h%s" % n, ";%s;Reco. energy" % n, 
#                   100, binning[n][1], binning[n][2], 1, 0.375, 0.475)
#     histo2D.append(h)

# histo2d_dict = dict(zip(variables_dict.keys(), histo2D))

# for i in range(t_data.GetEntries()):
#     t_data.GetEntry(i)
#     if pre_cuts(t_data) and manual_cuts(t_data):
#         for name, var in variables:
#             histo2d_dict[name].Fill(
#                 var[0], t_data.reco_energy, t_data.event_weight)
#         for name, var in spectators:
#             if spectators != "reco_energy":
#                 histo2d_dict[name].Fill(
#                     var[0], t_data.reco_energy, t_data.event_weight)

# OBJECTS = []
# for h in histo2d_dict:
#     c = ROOT.TCanvas(h)
#     histo2d_dict[h].SetMarkerStyle(20)

#     histo2d_dict[h].Draw()
#     c.Update()
#     OBJECTS.append(c)

print("Data events", fill_histos_data("bnb", bdt, manual))

# fill_histos_data("bnbext", bdt, manual)
input()
