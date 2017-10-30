#!/usr/local/bin/python3

import ROOT

from bdt_common import binning, labels, variables, spectators, bdt_cut, manual_cuts
from bdt_common import x_start, x_end, y_start, y_end, z_start, z_end
from array import array
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


variables_dict = dict(variables + spectators)

histograms = []

for i, n in enumerate(variables_dict.keys()):
    if n != "reco_energy":
        h = ROOT.TH1F("h_%s" % n, labels[n],
                      binning[n][0], binning[n][1], binning[n][2])
    else:
        bins = array("f", [0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.6, 0.8, 1])
        h = ROOT.TH1F("h_%s" % n, labels[n],len(bins)-1, bins)
    histograms.append(h)

histo_dict = dict(zip(variables_dict.keys(), histograms))

h_xy_track_start = ROOT.TH2F("h_xy_track_start", "x [cm]; y [cm]",
                             10, x_start, x_end, 10, y_start, y_end)
h_yz_track_start = ROOT.TH2F("h_yz_track_start", "x [cm]; y [cm]",
                             10, y_start, y_end, 10, z_start, z_end)
h_xz_track_start = ROOT.TH2F("h_xz_track_start", "x [cm]; y [cm]",
                             10, x_start, x_end, 10, z_start, z_end)
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

h_bdt = ROOT.TH1F("h_bdt_intime", "BDT response; N. Entries / 0.05", 40, -1, 1)


for i in range(t_data.GetEntries()):
    t_data.GetEntry(i)
    BDT_response = reader.EvaluateMVA("BDT method")
    h_bdt.Fill(BDT_response, t_data.event_weight)


    if BDT_response > bdt_cut and manual_cuts(t_data):
        for name, var in variables:
            histo_dict[name].Fill(var[0], t_data.event_weight)
        for name, var in spectators:
            histo_dict[name].Fill(var[0], t_data.event_weight)


    h_xy_track_start.Fill(variables_dict["track_start_x"][0],
                          variables_dict["track_start_y"][0],
                          t_data.event_weight)
    h_yz_track_start.Fill(variables_dict["track_start_y"][0],
                          variables_dict["track_start_z"][0],
                          t_data.event_weight)
    h_xz_track_start.Fill(variables_dict["track_start_x"][0],
                          variables_dict["track_start_z"][0],
                          t_data.event_weight)

f_bdt = ROOT.TFile("plots/h_bdt_intime.root", "RECREATE")
h_bdt.Write()
f_bdt.Close()

f_2d = ROOT.TFile("2d_mc.root", "RECREATE")
h_xy_track_start.Write()
h_xz_track_start.Write()
h_yz_track_start.Write()
f_2d.Close()

for h in histograms:
    f = ROOT.TFile("plots/%s_cosmic.root" % h.GetName(), "RECREATE")
    h.Write()
    f.Close()
