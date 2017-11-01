#!/usr/local/bin/python3

import ROOT
from bdt_common import fill_histos, variables, spectators

fill_histos("cosmic_mc")

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
#
# h_bdt = ROOT.TH1F("h_bdt_intime", "BDT response; N. Entries / 0.05", 40, -1, 1)
#
#
# for i in range(t_data.GetEntries()):
#     t_data.GetEntry(i)
#     BDT_response = reader.EvaluateMVA("BDT method")
#     h_bdt.Fill(BDT_response, t_data.event_weight)
#
#
#     if BDT_response > bdt_cut and manual_cuts(t_data):
#         for name, var in variables:
#             histo_dict[name].Fill(var[0], t_data.event_weight)
#         for name, var in spectators:
#             histo_dict[name].Fill(var[0], t_data.event_weight)
#
#
#
# f_bdt = ROOT.TFile("plots/h_bdt_intime.root", "RECREATE")
# h_bdt.Write()
# f_bdt.Close()
#
# f_2d = ROOT.TFile("2d_mc.root", "RECREATE")
# h_xy_track_start.Write()
# h_xz_track_start.Write()
# h_yz_track_start.Write()
# f_2d.Close()
#
# for h in histograms:
#     f = ROOT.TFile("plots/%s_cosmic.root" % h.GetName(), "RECREATE")
#     h.Write()
#     f.Close()
