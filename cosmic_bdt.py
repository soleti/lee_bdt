#!/usr/local/bin/python3

import ROOT
from bdt_common import variables, spectators


f_data = ROOT.TFile("cosmic_mc_file.root")
t_data = f_data.Get("cosmic_mc_tree")

ROOT.TMVA.Tools.Instance()
reader = ROOT.TMVA.Reader(":".join([
    "!V",
    "!Silent",
    "Color"]))

for name, var in variables:
    t_data.SetBranchAddress(name, var)

for name, var in spectators:
    t_data.SetBranchAddress(name, var)

for name, var in variables:
    reader.AddVariable(name, var)

for name, var in spectators:
    reader.AddSpectator(name, var)

reader.BookMVA("BDT method",
               "dataset/weights/TMVAClassification_BDT.weights.xml")

h_reco = ROOT.TH1F("h_reco", "", 16, 0.2, 1)

for i in range(300):
    cut = -0.2 + i * 0.004
    h_reco.Reset()
    for i in range(t_data.GetEntries()):
        t_data.GetEntry(i)
        BDT_response = reader.EvaluateMVA("BDT method")
        if BDT_response > cut:
            h_reco.Fill(t_data.reco_energy, t_data.event_weight)

    f_bdt = ROOT.TFile("bdt_cosmic/cut%.3f.root" % cut, "RECREATE")
    h_reco.Write()
    f_bdt.Close()
