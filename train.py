#!/usr/local/bin/python3

import ROOT
from bdt_common import variables, spectators, manual, bdt, bins, binning


f_input = ROOT.TFile("mc_file.root")
t = ROOT.TTree()
t = f_input.Get("mc_tree")

f_cosmic = ROOT.TFile("bnbext_file.root")
t_cosmic = ROOT.TTree()
t_cosmic = f_cosmic.Get("bnbext_tree")

f_nue = ROOT.TFile("nue_file.root")
t_nue = ROOT.TTree()
t_nue = f_nue.Get("nue_tree")

f_lee = ROOT.TFile("lee_file.root")
t_lee = ROOT.TTree()
t_lee = f_lee.Get("lee_tree")

ROOT.TMVA.Tools.Instance()
fout = ROOT.TFile("test.root", "RECREATE")
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


sigCut = ROOT.TCut("shower_track_d < 5 && track_distance < 5 && shower_distance < 5 && dedx > 1 && dedx < 3.4 && shower_energy > 0.03 && total_shower_energy > 0.01 && shower_energy_y > 0.01 && total_track_energy_length > 0 && numu_score < 17 && track_hits > 5 && shower_hits_y > 5 && total_hits_y > 0 && total_hits_u > 0 && total_hits_v > 0 && category == 2")
bgCut = ROOT.TCut("shower_track_d < 5 && track_distance < 5 && shower_distance < 5 && dedx > 1 && dedx < 3.4 && shower_energy > 0.03 && total_shower_energy > 0.01 && shower_energy_y > 0.01 && total_track_energy_length > 0 && category != 2 && numu_score < 17 && track_hits > 5 && shower_hits_y > 5 && total_hits_y > 0 && total_hits_u > 0 && total_hits_v > 0")

# dataloader.AddSignalTree(t_lee)
dataloader.AddBackgroundTree(t_nue)
dataloader.AddSignalTree(t_nue)
dataloader.AddBackgroundTree(t)
dataloader.AddBackgroundTree(t_cosmic)

dataloader.PrepareTrainingAndTestTree(sigCut, bgCut,
                                      ":".join([
                                          "SplitMode=Alternate",
                                          "NormMode=NumEvents:!V"]))

# if not manual:
#     method_cuts = factory.BookMethod(dataloader, ROOT.TMVA.Types.kCuts, "Cuts",
#                                     ":".join([
#                                         "!H:!V",
#                                         "FitMethod=GA",
#                                         "EffSel",
#                                         "CutRangeMin[0]=0:CutRangeMax[0]=6",
#                                         "CutRangeMin[1]=0:CutRangeMax[1]=1",
#                                         "VarProp[1]=FSmart",
#                                         "CutRangeMin[2]=0:CutRangeMax[2]=10",
#                                         "CutRangeMin[3]=0:CutRangeMax[3]=10",
#                                         "CutRangeMin[5]=0:CutRangeMax[5]=180",
#                                         "CutRangeMin[6]=-180:CutRangeMax[6]=180",
#                                         "CutRangeMin[7]=-180:CutRangeMax[7]=180",
#                                         "CutRangeMin[8]=0:CutRangeMax[8]=180"]))

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

# method_likelihood2 = factory.BookMethod(dataloader, ROOT.TMVA.Types.kPDERS, "PDERS",
#                                        ":".join(["V", "H"]))

# method_likelihood = factory.BookMethod(dataloader, ROOT.TMVA.Types.kLikelihood, "Likelihood",
#                                        ":".join(["V", "H"]))

var_list = [var[0] for var in variables]
range_min = ["CutRangeMin[%i]=%.2f" % (var_list.index(var), binning[var][1]) for var in var_list]
range_max = ["CutRangeMax[%i]=%.2f" % (var_list.index(var), binning[var][2]) for var in var_list]

cuts_min = ":".join(range_min)
cuts_max = ":".join(range_max)
# print(cuts_min, cuts_max)

# method_cuts = factory.BookMethod(dataloader, ROOT.TMVA.Types.kCuts, "Cuts", cuts_min + ":" + cuts_max)

factory.TrainAllMethods()
factory.TestAllMethods()
factory.EvaluateAllMethods()
