#!/usr/bin/env python3.4

import ROOT
from glob import glob

bnb_cosmic = glob("nu_files/*/*.root")
nue_cosmic = glob("nue_files/*/*.root")

f = ROOT.TFile("kin_file.root")
t = ROOT.TTree()
t = f.Get("kin_tree")

ROOT.TMVA.Tools.Instance()
fout = ROOT.TFile("test.root","RECREATE")
factory = ROOT.TMVA.Factory("TMVAClassification", fout,
                            ":".join([
                                "!V",
                                "!Silent",
                                "Color",
                                "DrawProgressBar",
                                "Transformations=I;D;P;G,D",
                                "AnalysisType=Classification"]))
dataloader = ROOT.TMVA.DataLoader("dataset");

dataloader.AddVariable("reco_energy","F")
dataloader.AddVariable("track_length","F")
dataloader.AddVariable("track_theta","F")
dataloader.AddVariable("track_phi","F")
dataloader.AddVariable("shower_energy","F")
dataloader.AddVariable("shower_theta","F")
dataloader.AddVariable("shower_phi","F")
dataloader.AddVariable("shower_z","F")
dataloader.AddVariable("track_z","F")
dataloader.AddVariable("pt","F")
dataloader.AddVariable("n_tracks","F")
dataloader.AddVariable("n_showers","F")
dataloader.AddVariable("track_shower_angle","F")
dataloader.AddVariable("shower_distance","F")
dataloader.AddSpectator("category","F")
dataloader.AddSpectator("event_weight","F")
dataloader.AddSpectator("event","F")
dataloader.AddSpectator("run","F")
dataloader.AddSpectator("subrun","F")

sigCut = ROOT.TCut("is_signal > 0.5")
bgCut = ROOT.TCut("is_signal <= 0.5")

dataloader.AddSignalTree(t)
dataloader.AddBackgroundTree(t)
dataloader.PrepareTrainingAndTestTree(sigCut, bgCut, "SplitMode=Random:NormMode=NumEvents:!V")

method = factory.BookMethod(dataloader, ROOT.TMVA.Types.kBDT, "BDT",
                   "!H:!V:NTrees=850:MinNodeSize=2.5%:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=20")

factory.TrainAllMethods()
factory.TestAllMethods()
factory.EvaluateAllMethods()

tout = fout.Get("dataset/TestTree")
cout = ROOT.TCanvas("cout")
h_sig = ROOT.TH1F("h_sig",";Score;N. Entries / 0.025",80,-1,1)
h_bkg = ROOT.TH1F("h_bkg",";Score;N. Entries / 0.025",80,-1,1)
h_bkg.SetLineColor(ROOT.kRed+1)
tout.Draw("BDT>>h_bkg","classID==1")
tout.Draw("BDT>>h_sig","classID==0","same")
h_bkg.Scale(1/h_bkg.Integral())
h_sig.Scale(1/h_sig.Integral())
h_bkg.Draw("hist")
h_sig.Draw("hist same")
cout.Update()

h_energy = ROOT.THStack("h_energy",";Reco. energy [GeV];N.Entries / 0.1 GeV")

h_cosmic = ROOT.TH1F("h_cosmic",";Reco. energy [GeV];N.Entries / 0.1 GeV",20,0,2)
h_dirt = ROOT.TH1F("h_dirt",";Reco. energy [GeV];N.Entries / 0.1 GeV",20,0,2)
h_nu_e = ROOT.TH1F("h_nu_e",";Reco. energy [GeV];N.Entries / 0.1 GeV",20,0,2)
h_nu_mu = ROOT.TH1F("h_nu_mu",";Reco. energy [GeV];N.Entries / 0.1 GeV",20,0,2)
h_other = ROOT.TH1F("h_other",";Reco. energy [GeV];N.Entries / 0.1 GeV",20,0,2)
h_nc = ROOT.TH1F("h_nc",";Reco. energy [GeV];N.Entries / 0.1 GeV",20,0,2)

h_energies = [h_other,h_cosmic,h_nu_e,h_nu_mu,h_nc,h_dirt]
colors = [ROOT.kGray+2, ROOT.kRed - 3, ROOT.kGreen - 2, ROOT.kBlue - 5, ROOT.kBlue - 9, ROOT.kOrange+3]
l_energy = ROOT.TLegend(0.48,0.55,0.84,0.84)
description = ["Other", "Cosmic", "Beam Intrinsic #nu_{e}", "Beam Intrinsic #nu_{#mu}", "Beam Intrinsic NC", "Dirt"]
i = 0
event_file = open("events.txt","w")

for entry in tout:
    if entry.BDT > -0.0068:
        if entry.event_weight > 1:
            dataset = "bnb"
        else:
            dataset = "nu_e"

        print(dataset,int(entry.run),int(entry.subrun),int(entry.event),int(entry.category),file=event_file)
        if entry.category == 2:
            h_nu_e.Fill(entry.reco_energy, entry.event_weight*2)
        else:
            h_energies[int(entry.category)].Fill(entry.reco_energy, entry.event_weight*2)

event_file.close()

for i,h in enumerate(h_energies):
    h.SetLineColor(1)
    h.SetLineWidth(2)
    h.SetFillColor(colors[i])
    h_energy.Add(h)
    l_energy.AddEntry(h, "%s: %.0f events" % (description[i], h.Integral(3,20)), "f")

cenergy = ROOT.TCanvas("cenergy","")

h_energy.Draw("hist")
h_energy.GetXaxis().SetRangeUser(0.2,2)
l_energy.Draw()
# Loop over tout and fill the energy histogram with the right weight
cenergy.Update()
cenergy.SaveAs("plots/energy.pdf")

input()
