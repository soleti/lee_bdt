#!/usr/local/bin/python3

import sys
import ROOT
from bdt_common import variables, spectators, MANUAL, BDT, bins, binning

if len(sys.argv) > 1:
    bdt_type = sys.argv[1]
else:
    bdt_type = ""

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
fout = ROOT.TFile("root_files/tmva_%s.root" % bdt_type, "RECREATE")
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


pre_cuts = "total_shower_energy_cali > 0.01 && \
            total_track_energy_length > 0 && \
            shower_dqdx[shower_id] > 10000 && \
            shower_dqdx[shower_id] < 80000 && \
            shower_theta > -999 && \
            track_length > -999 && \
            shower_energy > -999 && \
            track_start_y > -999 && \
            shower_start_y > -999 && \
            track_res_mean > -999 && \
            track_theta > -999 && \
            track_pidchipr > -999 && \
            track_distance > -999 && \
            shower_distance > -999 && \
            track_likelihood > -999 && \
            numu_score == 0 && \
            track_hits > 5 && \
            shower_hits > 5"

sigCut = ROOT.TCut("(category == 2 || category == 10) && " + pre_cuts)

bgCutNeutrino = ROOT.TCut("(category == 3 || category == 4 || category == 7) && " + pre_cuts)
bgCutNuMu = ROOT.TCut("category == 3 && " + pre_cuts)
bgCutNC = ROOT.TCut("category == 4 && " + pre_cuts)
bgCutCosmic = ROOT.TCut("(category == 0 || category == 1 || category == 7) && " + pre_cuts)

bgCut = ROOT.TCut("(category == 3 || category == 4 || category == 7 || category == 0 || category == 1 || category == 5) && " + pre_cuts)


dataloader.AddBackgroundTree(t_nue)
# dataloader.AddSignalTree(t_lee)
dataloader.AddSignalTree(t_nue)
dataloader.AddBackgroundTree(t)
dataloader.AddBackgroundTree(t_cosmic)

if bdt_type == "cosmic":
    dataloader.PrepareTrainingAndTestTree(sigCut, bgCutCosmic,
                                        ":".join([
                                            "SplitMode=Alternate",
                                            "NormMode=NumEvents:!V"]))
elif bdt_type == "neutrino":
    dataloader.PrepareTrainingAndTestTree(sigCut, bgCutNeutrino,
                                        ":".join([
                                            "SplitMode=Alternate",
                                            "NormMode=NumEvents:!V"]))
elif bdt_type == "nc":
    dataloader.PrepareTrainingAndTestTree(sigCut, bgCutNC,
                                          ":".join([
                                              "SplitMode=Alternate",
                                              "NormMode=NumEvents:!V"]))
elif bdt_type == "numu":
    dataloader.PrepareTrainingAndTestTree(sigCut, bgCutNuMu,
                                          ":".join([
                                              "SplitMode=Alternate",
                                              "NormMode=NumEvents:!V"]))
elif bdt_type == "":
    dataloader.PrepareTrainingAndTestTree(sigCut, bgCut,
                                          ":".join([
                                              "SplitMode=Alternate",
                                              "NormMode=NumEvents:!V"]))

method_bdt = factory.BookMethod(dataloader, ROOT.TMVA.Types.kBDT, "BDT%s" % bdt_type,
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

# method_svm = factory.BookMethod(dataloader, ROOT.TMVA.Types.kBDT, "SVM%s" % bdt_type)

var_list = [var[0] for var in variables]
range_min = ["CutRangeMin[%i]=%.2f" % (var_list.index(var), binning[var][1]) for var in var_list]
range_max = ["CutRangeMax[%i]=%.2f" % (var_list.index(var), binning[var][2]) for var in var_list]

cuts_min = ":".join(range_min)
cuts_max = ":".join(range_max)

# method_cuts = factory.BookMethod(dataloader, ROOT.TMVA.Types.kCuts, "Cuts", cuts_min + ":" + cuts_max)

factory.TrainAllMethods()
factory.TestAllMethods()
factory.EvaluateAllMethods()
