#!/usr/local/bin/python3

import ROOT
from bdt_common import variables, spectators, manual, bdt, bins, binning


f_input = ROOT.TFile("root_files/mc_file.root")
t = ROOT.TTree()
t = f_input.Get("mc_tree")

f_cosmic = ROOT.TFile("root_files/bnbext_file.root")
t_cosmic = ROOT.TTree()
t_cosmic = f_cosmic.Get("bnbext_tree")

f_nue = ROOT.TFile("root_files/nue_file.root")
t_nue = ROOT.TTree()
t_nue = f_nue.Get("nue_tree")

f_lee = ROOT.TFile("root_files/lee_file.root")
t_lee = ROOT.TTree()
t_lee = f_lee.Get("lee_tree")

ROOT.TMVA.Tools.Instance()
fout = ROOT.TFile("root_files/test.root", "RECREATE")
factory = ROOT.TMVA.Factory("TMVAClassification", fout,
                            ":".join([
                                "!V",
                                "!Silent",
                                "Color",
                                "DrawProgressBar",
                                "Transformations=I;D;P;G,D",
                                "AnalysisType=Classification"]))
dataloader = ROOT.TMVA.DataLoader("dataset")

for name, var in variables:
    dataloader.AddVariable(name, "F")


for name, var in spectators:
    dataloader.AddSpectator(name, "F")

# sigCut = ROOT.TCut("shower_track_d < 5 && track_distance < 5 && shower_distance < 5 && total_shower_energy > 0.01 && total_track_energy_length > 0 && numu_score < 17 && track_hits > 5 && shower_hits_y > 5 && total_hits_y > 0 && total_hits_u > 0 && total_hits_v > 0 && category == 2")
# bgCut = ROOT.TCut("shower_track_d < 5 && track_distance < 5 && shower_distance < 5 && total_shower_energy > 0.01  && total_track_energy_length > 0 && category != 2 && numu_score < 17 && track_hits > 5 && shower_hits_y > 5 && total_hits_y > 0 && total_hits_u > 0 && total_hits_v > 0")

sigCut = ROOT.TCut("category == 10 && shower_dedx_cali > -999 && track_distance > -999 && track_distance < 6 && shower_distance > -999 && numu_score < 17 && track_hits > 5 && shower_hits > 5")
bgCut = ROOT.TCut("category != 2 && shower_dedx_cali > -999 && track_distance > -999 && track_distance < 6 && shower_distance > -999 && numu_score < 17 && track_hits > 5 && shower_hits > 5")

dataloader.AddBackgroundTree(t_nue)
dataloader.AddSignalTree(t_lee)
dataloader.AddBackgroundTree(t)
dataloader.AddBackgroundTree(t_cosmic)

dataloader.PrepareTrainingAndTestTree(sigCut, bgCut,
                                      ":".join([
                                        #   "SplitMode=Alternate",
                                          "NormMode=NumEvents:!V"]))


method_bdt = factory.BookMethod(dataloader, ROOT.TMVA.Types.kBDT, "BDT",
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

var_list = [var[0] for var in variables]
range_min = ["CutRangeMin[%i]=%.2f" % (var_list.index(var), binning[var][1]) for var in var_list]
range_max = ["CutRangeMax[%i]=%.2f" % (var_list.index(var), binning[var][2]) for var in var_list]

cuts_min = ":".join(range_min)
cuts_max = ":".join(range_max)

# method_cuts = factory.BookMethod(dataloader, ROOT.TMVA.Types.kCuts, "Cuts", cuts_min + ":" + cuts_max)

factory.TrainAllMethods()
factory.TestAllMethods()
factory.EvaluateAllMethods()
