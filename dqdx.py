#!/usr/local/bin/python3

import ROOT
from glob import glob
import math
from array import array
from bdt_common import printProgressBar, is_1eNp, PROTON_MASS, is_fiducial

def train_BDT(name, variables, signal_tree, bkg_tree):
    ROOT.TMVA.Tools.Instance()
    fout = ROOT.TFile("%s.root" % name, "RECREATE")
    factory = ROOT.TMVA.Factory("TMVAClassification", fout,
                                ":".join([
                                    "!V",
                                    "!Silent",
                                    "Color",
                                    "DrawProgressBar",
                                    "Transformations=I;D;P;G,D",
                                    "AnalysisType=Classification"]))
    dataloader = ROOT.TMVA.DataLoader("dataset")
    for v in variables:
        dataloader.AddVariable(v, "F")
    
    dataloader.AddSignalTree(signal_tree)
    dataloader.AddBackgroundTree(bkg_tree)
    sigCut = ROOT.TCut("signal > 0.5")
    bgCut = ROOT.TCut("signal <= 0.5")
    dataloader.PrepareTrainingAndTestTree(sigCut, bgCut,
                                        ":".join([
                                            "SplitMode=Alternate",
                                            "NormMode=NumEvents:!V"]))

    method_bdt = factory.BookMethod(dataloader, ROOT.TMVA.Types.kBDT, "%s BDT" % name,
                                    ":".join([
                                        "!H:!V:NTrees=600",
                                        "MinNodeSize=5%",
                                        "MaxDepth=5",
                                        "BoostType=AdaBoost",
                                        "AdaBoostBeta=0.5",
                                        "UseBaggedBoost",
                                        "BaggedSampleFraction=0.5",
                                        "SeparationType=GiniIndex",
                                        "nCuts=20"]))

    factory.TrainAllMethods()
    factory.TestAllMethods()
    factory.EvaluateAllMethods()
    fout.Close()

ROOT.gStyle.SetOptTitle(1)
ROOT.gStyle.SetTitleFillColor(ROOT.kWhite)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetPalette(87)
ROOT.gStyle.SetNumberContours(99)
ROOT.gStyle.SetStatW(0.16)
ROOT.gStyle.SetOptFit(0)

nue_cosmic = glob("mc_bnb_cali/*/*.root")
c = ROOT.TChain("robertoana/pandoratree")

for f in nue_cosmic:
    c.Add(f)

h_dqdx_length_proton = ROOT.TH2F(
    "h_dqdx_length_proton", ";dQ/dx [a.u.];Track length [cm]", 40, 0, 1000, 40, 0, 80)

h_dqdx_length_muon = ROOT.TH2F(
    "h_dqdx_length_muon", ";dQ/dx [a.u.];Track length [cm]", 40, 0, 1000, 40, 0, 80)

h_dqdx_length_electron = ROOT.TH2F(
    "h_dqdx_length_electron", ";dE/dx [a.u.];Track length [cm]", 40, 0, 6, 40, 0, 100)

h_dqdx_length_other = ROOT.TH2F(
    "h_dqdx_length_other", ";dE/dx [a.u.];Track length [cm]", 40, 0, 6, 40, 0, 100)

dqdx_ntuple = ROOT.TNtuple("dqdx_ntuple", "dqdx_ntuple", "dqdx:len:signal")
dedx_ntuple = ROOT.TNtuple("dedx_ntuple", "dedx_ntuple", "dedx:nhits:signal")

entries = int(c.GetEntries() / 5)

for i in range(entries):
    c.GetEntry(i)
    printProgressBar(i, entries, prefix='Progress:',
                     suffix='Complete', length=20)

    if not c.passed:
        continue

    for i_sh in range(c.n_showers):
        if c.shower_dQdx[i_sh][2] < 0:
            continue

        if abs(c.matched_showers[i_sh]) == 13:
            h_dqdx_length_electron.Fill(
                c.shower_dEdx[i_sh][2],
                c.shower_nhits[i_sh][2])
            dedx_ntuple.Fill(c.shower_dEdx[i_sh][2],
                             c.shower_nhits[i_sh][2],
                             1)

        elif abs(c.matched_showers[i_sh]) == 2212:
            h_dqdx_length_other.Fill(
                c.shower_dEdx[i_sh][2],
                c.shower_nhits[i_sh][2])
            dedx_ntuple.Fill(c.shower_dEdx[i_sh][2],
                             c.shower_nhits[i_sh][2],
                             0)

    for i_tr in range(c.n_tracks):
        track_start = [c.track_start_x[i_tr],
                       c.track_start_y[i_tr],
                       c.track_start_z[i_tr]]
        track_end = [c.track_end_x[i_tr],
                       c.track_end_y[i_tr],
                       c.track_end_z[i_tr]]
        if not is_fiducial(track_start) or not is_fiducial(track_end):
            continue
        if c.track_dQdx[i_tr][2] < 0:
            continue
        if abs(c.matched_tracks[i_tr]) == 2212:
            h_dqdx_length_proton.Fill(c.track_dQdx[i_tr][2], c.track_len[i_tr])
            dqdx_ntuple.Fill(c.track_dQdx[i_tr][2], c.track_len[i_tr], 1)
        elif abs(c.matched_tracks[i_tr]) == 13:
            h_dqdx_length_muon.Fill(c.track_dQdx[i_tr][2], c.track_len[i_tr])
            dqdx_ntuple.Fill(c.track_dQdx[i_tr][2], c.track_len[i_tr], 0)


train_BDT("dqdx", ["dqdx", "len"], dqdx_ntuple, dqdx_ntuple)
train_BDT("dedx", ["dedx", "nhits"], dedx_ntuple, dedx_ntuple)

c_dqdx_length = ROOT.TCanvas("c_dqdx_length", "")
h_dqdx_length_muon.Scale(1 / h_dqdx_length_muon.Integral())
h_dqdx_length_muon.SetLineColor(ROOT.kRed - 4)
h_dqdx_length_muon.SetMarkerStyle(6)
h_dqdx_length_proton.SetMarkerStyle(6)

h_dqdx_length_proton.SetLineColor(ROOT.kAzure + 1)
h_dqdx_length_proton.Scale(1 / h_dqdx_length_proton.Integral())
h_dqdx_length_muon.Draw("BOX")
h_dqdx_length_proton.Draw("BOX SAME")
c_dqdx_length.Update()


c_dqdx_length_shower = ROOT.TCanvas("c_dqdx_length_shower", "")
h_dqdx_length_electron.Scale(1 / h_dqdx_length_electron.Integral())
h_dqdx_length_electron.SetLineColor(ROOT.kRed - 4)
h_dqdx_length_electron.SetMarkerStyle(6)
h_dqdx_length_other.SetMarkerStyle(6)

h_dqdx_length_other.SetLineColor(ROOT.kAzure + 1)
h_dqdx_length_other.Scale(1 / h_dqdx_length_other.Integral())
h_dqdx_length_electron.Draw("BOX")
h_dqdx_length_other.Draw("BOX SAME")
c_dqdx_length_shower.Update()


input()

    
