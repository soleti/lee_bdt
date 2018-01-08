#!/usr/local/bin/python3

import ROOT
from bdt_common import variables, spectators, fill_histos, manual, bdt


f_input = ROOT.TFile("mc_file.root")
t = ROOT.TTree()
t = f_input.Get("mc_tree")

f_cosmic = ROOT.TFile("cosmic_mc_file.root")
t_cosmic = ROOT.TTree()
t_cosmic = f_cosmic.Get("cosmic_mc_tree")

f_nue = ROOT.TFile("nue_file.root")
t_nue = ROOT.TTree()
t_nue = f_nue.Get("nue_tree")

# f_pi0 = ROOT.TFile("pi0_file.root")
# t_pi0 = ROOT.TTree()
# t_pi0 = f_nue.Get("pi0_tree")

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

sigCut = ROOT.TCut("is_signal > 0.5 && reco_energy > 0.2 && reco_energy < 0.6 && shower_energy > 0.1")
bgCut = ROOT.TCut("(is_signal <= 0.5 && reco_energy > 0.2 && reco_energy < 0.6 && shower_energy > 0.1)")#" || (is_signal > 0.5 && (reco_energy < 0.2 || reco_energy > 0.6))")


dataloader.AddSignalTree(t_nue)
dataloader.AddBackgroundTree(t_nue)
dataloader.AddBackgroundTree(t)
dataloader.AddBackgroundTree(t_cosmic)

dataloader.PrepareTrainingAndTestTree(sigCut, bgCut,
                                      ":".join([
                                          "SplitMode=Alternate",
                                          "NormMode=NumEvents:!V"]))

method_cuts = factory.BookMethod(dataloader, ROOT.TMVA.Types.kCuts, "Cuts",
                                 ":".join([
                                     "!H:!V",
                                     "FitMethod=GA",
                                     "EffSel",
                                     "CutRangeMin[0]=0:CutRangeMax[0]=6",
                                     "CutRangeMin[1]=0:CutRangeMax[1]=1",
                                     "VarProp[1]=FSmart",
                                     "CutRangeMin[2]=0:CutRangeMax[2]=10",
                                     "CutRangeMin[3]=0:CutRangeMax[3]=10"]))

method_bdt = factory.BookMethod(dataloader, ROOT.TMVA.Types.kBDT, "BDT",
                                ":".join([
                                    "!H:!V:NTrees=50",
                                    "MinNodeSize=2.5%",
                                    "MaxDepth=3",
                                    "BoostType=AdaBoost",
                                    "AdaBoostBeta=0.5",
                                    "UseBaggedBoost",
                                    "BaggedSampleFraction=0.5",
                                    "SeparationType=GiniIndex",
                                    "nCuts=20"]))

method_likelihood = factory.BookMethod(dataloader, ROOT.TMVA.Types.kLikelihood,
                                       "Likelihood",
                                       ":".join(["!H:!V",
                                                 "TransformOutput",
                                                 "PDFInterpol=Spline2",
                                                 "NSmoothSig[0]=20",
                                                 "NSmoothBkg[0]=20",
                                                 "NSmooth=1",
                                                 "NAvEvtPerBin=50"]))

factory.TrainAllMethods()
factory.TestAllMethods()
factory.EvaluateAllMethods()

fill_histos("bnb", bdt, manual)
fill_histos("bnbext", bdt, manual)
fill_histos("lee", bdt, manual)
