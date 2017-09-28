#!/usr/local/bin/python3

import ROOT
from bdt_common import binning, labels, variables, spectators, bdt_cut

f_data = ROOT.TFile("bnb_file.root")
t_data = f_data.Get("bnb_tree")

ROOT.TMVA.Tools.Instance()
reader = ROOT.TMVA.Reader(":".join([
    "!V",
    "!Silent",
    "Color"]))

for name, var in variables:
    t_data.SetBranchAddress(name, var)

for name, var in variables:
    reader.AddVariable(name, var)

for name, var in spectators:
    reader.AddSpectator(name, var)

reader.BookMVA("BDT method",
               "dataset/weights/TMVAClassification_BDT.weights.xml")


variables_dict = dict(variables)

histograms = []

for i, n in enumerate(variables_dict.keys()):
    h = ROOT.TH1F("h_" + n, labels[n], binning[n][0], binning[n][1],
                  binning[n][2])
    histograms.append(h)

histo_dict = dict(zip(variables_dict.keys(), histograms))

h_bdt = ROOT.TH1F("h_bdt", ";BDT response; N. Entries / 0.05", 40, -1, 1)

passed_events = open("data_passed.txt", "w")

for i in range(t_data.GetEntries()):
    t_data.GetEntry(i)
    BDT_response = reader.EvaluateMVA("BDT method")
    h_bdt.Fill(BDT_response, t_data.event_weight)

    if BDT_response > bdt_cut:
        print("{} {} {} {}".format(int(tout.run),
                                   int(tout.subrun),
                                   int(tout.event),
                                   tout.event_weight * 2), file = passed_events)

        for name, var in variables:
            histo_dict[name].Fill(var[0], t_data.event_weight)

f_bdt = ROOT.TFile("plots/h_bdt_data.root", "RECREATE")
h_bdt.Write()
f_bdt.Close()

for h in histograms:
    f = ROOT.TFile("plots/%s_data.root" % h.GetName(), "RECREATE")
    h.Write()
    f.Close()
