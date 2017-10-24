#!/usr/local/bin/python3

import ROOT
from array import array
import math
from bdt_common import bdt_cut, binning, labels, find_interaction, variables, spectators, sigmaCalc
from bdt_common import interactions, x_start, x_end, y_start, y_end, z_start, z_end

fout = ROOT.TFile("test.root")

tout = ROOT.TChain("dataset/TestTree")
ttrain = ROOT.TChain("dataset/TrainTree")
ttrain.Add("test.root")
tout.Add("test.root")
colors = [ROOT.kGray + 2, ROOT.kRed - 3, ROOT.kGreen - 2, ROOT.kBlue - 5,
          ROOT.kBlue - 9, ROOT.kOrange + 3, ROOT.kWhite, ROOT.kRed - 3]

kinds = []
categories = ["other", "cosmic", "nu_e", "nu_mu", "nc", "dirt", "data",
              "mixed"]
stacked_histos = []

for name, var in variables:
    tout.SetBranchAddress(name, var)
    ttrain.SetBranchAddress(name, var)

for name, var in spectators:
    tout.SetBranchAddress(name, var)
    ttrain.SetBranchAddress(name, var)

variables_dict = dict(variables+spectators)

for i, n in enumerate(variables_dict.keys()):
    # h = ROOT.TH1F(n,labels[i],binning[i][0],binning[i][1],binning[i][2])
    histos = []

    h_stack = ROOT.THStack("h_" + n, labels[n])

    for c in categories:
        h = ROOT.TH1F("h_%s_%s" % (n, c), labels[n],
                      binning[n][0], binning[n][1], binning[n][2])
        histos.append(h)

    h_stack.Add(histos[0])
    h_stack.Add(histos[1])
    h_stack.Add(histos[7])
    h_stack.Add(histos[2])
    h_stack.Add(histos[3])
    h_stack.Add(histos[4])
    h_stack.Add(histos[5])
    h_stack.Add(histos[6])

    stacked_histos.append(h_stack)
    kinds.append(histos)

histo_dict = dict(zip(variables_dict.keys(), kinds))
bkg_types = [0] * 2000

h_bdt_stack = ROOT.THStack("h_bdt", ";BDT response; N. Entries / 0.05")
h_bdts = []
for i, c in enumerate(categories):
    h = ROOT.TH1F("h_bdt_%s" % c, "BDT response; N. Entries / 0.05", 40, -1, 1)
    h.SetLineColor(1)
    h.SetFillColor(colors[i])
    h.SetLineWidth(2)
    h_bdts.append(h)

h_bdt_stack.Add(h_bdts[0])
h_bdt_stack.Add(h_bdts[1])
h_bdt_stack.Add(h_bdts[7])
h_bdt_stack.Add(h_bdts[2])
h_bdt_stack.Add(h_bdts[3])
h_bdt_stack.Add(h_bdts[4])
h_bdt_stack.Add(h_bdts[5])
h_bdt_stack.Add(h_bdts[6])

passed_events = open("mc_passed.txt", "w")

nue_selected = 0
nue_bdt = 0

h_bkg = ROOT.TH1F("h_bkg","",16, 0.2, 1)
h_sig = ROOT.TH1F("h_sig","",16, 0.2, 1)
h_lee = ROOT.TH1F("h_lee","",16, 0.2, 1)

scaling = [6.0920944819073988, 3.6447414342239273, 3.2123920194399913, 2.6504659907742409, 3.2558450032216988, 2.5826310533377432, 2.6614353575699727, 1.4145769564088304, 1.0206172427887652, 0.9972444259255292, 0.79323702430381904, 0.63892043872491167, 0.61676413081900316, 0.3541651442224471, 0.28310400773433003, 0.94342108559739024]
best_sigma = 0
chosen_cut = 0
h_chosen_lee = ROOT.TH1F("h_chosen_lee", "", 16, 0.2, 1)
h_chosen_bkg = ROOT.TH1F("h_chosen_bkg", "", 16, 0.2, 1)

for i in range(200):
    bdt_cut = 0.2+i*0.004
    h_bkg.Reset()
    h_sig.Reset()
    h_lee.Reset()
    for i in range(tout.GetEntries()):
        tout.GetEntry(i)
        if tout.BDT > bdt_cut:

            if tout.category == 2:
                h_sig.Fill(tout.reco_energy, tout.event_weight)
            else:
                h_bkg.Fill(tout.reco_energy, tout.event_weight)

    for i in range(ttrain.GetEntries()):
        ttrain.GetEntry(i)
        if ttrain.BDT > bdt_cut:
            if ttrain.category == 2:
                h_sig.Fill(ttrain.reco_energy, ttrain.event_weight)
            else:
                h_bkg.Fill(ttrain.reco_energy, ttrain.event_weight)

    for i in range(len(scaling)):
        if scaling[i] - 1 > 0:
            h_lee.SetBinContent(i+1, h_sig.GetBinContent(i+1)*(scaling[i] - 1))

    h_bkg.Add(h_sig)
    f_cosmic_in_time = ROOT.TFile("bdt_cosmic/cut%.3f.root" % bdt_cut)
    h_cosmic_in_time = f_cosmic_in_time.Get("h_reco")
    h_bkg.Add(h_cosmic_in_time)

    if sigmaCalc(h_lee, h_bkg) > best_sigma:
        best_sigma = sigmaCalc(h_lee, h_bkg)
        h_chosen_lee = h_lee.Clone()
        h_chosen_bkg = h_bkg.Clone()
        chosen_cut = bdt_cut
    f_cosmic_in_time.Close()
print(best_sigma, chosen_cut)

for i in range(ttrain.GetEntries()):
    ttrain.GetEntry(i)
    category = int(ttrain.category)
    h_bdts[category].Fill(ttrain.BDT, ttrain.event_weight)

    if ttrain.BDT > chosen_cut:
        for name, var in variables:
            histo_dict[name][category].Fill(var[0], ttrain.event_weight)

        if ttrain.reco_energy > 0.2 and (ttrain.category == 3 or ttrain.category == 4):
            bkg_types[int(ttrain.interaction_type)] += ttrain.event_weight


for i in range(tout.GetEntries()):
    tout.GetEntry(i)
    category = int(tout.category)

    h_bdts[category].Fill(tout.BDT, tout.event_weight)

    if tout.category == 2:
        nue_selected += 2

    if tout.BDT > chosen_cut:

        if tout.category == 2:
            nue_bdt += 2

        # Store interaction types of background events
        if tout.reco_energy > 0.2 and (tout.category == 3 or tout.category == 4):
            bkg_types[int(tout.interaction_type)] += tout.event_weight

        print("{} {} {} {}".format(int(tout.run),
                                   int(tout.subrun),
                                   int(tout.event),
                                   tout.event_weight), file=passed_events)

        for name, var in variables:
            histo_dict[name][category].Fill(var[0], tout.event_weight)

        for name, var in spectators:
            histo_dict[name][category].Fill(var[0], tout.event_weight)


print(nue_selected, nue_bdt)
passed_events.close()

f_bdt = ROOT.TFile("plots/h_bdt.root", "RECREATE")
h_bdt_stack.Write()
f_bdt.Close()

for i, histos in enumerate(kinds):
    for j, h in enumerate(histos):
        h.SetLineColor(1)
        h.SetLineWidth(2)
        h.SetFillColor(colors[j])

for h in stacked_histos:
    f = ROOT.TFile("plots/%s.root" % h.GetName(), "RECREATE")
    h.Write()
    f.Close()

# Save interaction types of background events
event_file = open("events.txt", "w")
for i, interaction in enumerate(bkg_types):
    if int(interaction) > 0:
        print(find_interaction(interactions, i), int(interaction),
              file=event_file)
event_file.close()

# c = ROOT.TCanvas("c")
# hstack = ROOT.THStack("h","")
# hstack.Add(h_chosen_bkg)
# hstack.Add(h_chosen_lee)
# h_chosen_lee.SetLineColor(ROOT.kRed)
# hstack.Draw("hist")
# c.Update()
input()
