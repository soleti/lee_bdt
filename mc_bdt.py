#!/usr/bin/env python3.4

import ROOT
from array import array
import math
from bdt_common import bdt_cut, binning, labels, find_interaction, variables, interactions

fout = ROOT.TFile("test.root")

tout = ROOT.TChain("dataset/TestTree")
tout.Add("test.root")
colors = [ROOT.kGray+2, ROOT.kRed - 3, ROOT.kGreen - 2, ROOT.kBlue - 5, ROOT.kBlue - 9, ROOT.kOrange+3, ROOT.kRed - 3 ]

kinds = []
categories = ["other", "cosmic", "nu_e", "nu_mu", "nc", "dirt"]
stacked_histos = []

for name, var in variables:
    tout.SetBranchAddress(name, var)

variables_dict = dict(variables)

for i,n in enumerate(variables_dict.keys()):
    #h = ROOT.TH1F(n,labels[i],binning[i][0],binning[i][1],binning[i][2])
    histos = []

    h_stack = ROOT.THStack("h_"+n,labels[n])

    for c in categories:
        h = ROOT.TH1F("h_%s_%s" % (n, c), labels[n],binning[n][0],binning[n][1],binning[n][2])
        histos.append(h)
        h_stack.Add(h)

    stacked_histos.append(h_stack)
    kinds.append(histos)

histo_dict = dict(zip(variables_dict.keys(), kinds))
bkg_types = [0]*2000

for i in range(tout.GetEntries()):
    tout.GetEntry(i)
    if tout.BDT > bdt_cut:

        category = int(tout.category)
        if int(tout.category) == 7: category = 1

        # Store interaction types of background events
        # if tout.reco_energy > 0.1 and tout.category != 1 and tout.category != 2:
        #     bkg_types[int(tout.interaction_type)] += tout.event_weight*2
        #h_energies[int(tout.category)].Fill(tout.reco_energy, tout.event_weight*2)

        for name, var in variables:
            histo_dict[name][category].Fill(var[0],tout.event_weight*2)



for i,histos in enumerate(kinds):
    for j,h in enumerate(histos):
        h.SetLineColor(1)
        h.SetLineWidth(2)
        h.SetFillColor(colors[j])

for h in stacked_histos:
    f = ROOT.TFile("%s.root" % h.GetName(),"RECREATE")
    h.Write()
    f.Close()

# Save interaction types of background events
event_file = open("events.txt","w")
for i,interaction in enumerate(bkg_types):
    if int(interaction)>0:
        print(find_interaction(interactions,i),int(interaction),file=event_file)
event_file.close()
