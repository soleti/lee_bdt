#!/usr/bin/env python3.4

import ROOT
from bdt_common import binning, labels, variables, spectators, bdt_cut


f_input = ROOT.TFile("kin_file.root")
t = ROOT.TTree()
t = f_input.Get("kin_tree")

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


dataloader.AddVariable("track_length","F")
dataloader.AddVariable("track_theta","F")
dataloader.AddVariable("track_phi","F")
dataloader.AddVariable("shower_energy","F")
dataloader.AddVariable("shower_theta","F")
dataloader.AddVariable("shower_phi","F")
dataloader.AddVariable("pt","F")
dataloader.AddVariable("n_tracks","F")
dataloader.AddVariable("n_showers","F")
dataloader.AddVariable("track_shower_angle","F")
dataloader.AddVariable("proton_score","F")
dataloader.AddVariable("shower_distance","F")
dataloader.AddVariable("track_distance","F")
dataloader.AddVariable("track_start_y","F")
dataloader.AddVariable("track_end_y","F")
dataloader.AddVariable("track_start_x","F")
dataloader.AddVariable("track_end_x","F")
dataloader.AddVariable("track_start_z","F")
dataloader.AddVariable("track_end_z","F")
dataloader.AddVariable("shower_start_y","F")
dataloader.AddVariable("shower_start_x","F")
dataloader.AddVariable("shower_start_z","F")
dataloader.AddVariable("reco_energy","F")

dataloader.AddSpectator("category","F")
dataloader.AddSpectator("event_weight","F")
dataloader.AddSpectator("event","F")
dataloader.AddSpectator("run","F")
dataloader.AddSpectator("subrun","F")
dataloader.AddSpectator("interaction_type","F")

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
