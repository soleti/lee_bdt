import ROOT
from bdt_common import bdt_types
from array import array

def train_SVM(name, variables):
    ROOT.TMVA.Tools.Instance()
    fout = ROOT.TFile("root_files/%s.root" % name, "RECREATE")
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

    f_mc = ROOT.TFile("plots/bdt_ntuple_mc.root")
    mc_tree = f_mc.Get("bdt_ntuple")
    f_nue = ROOT.TFile("plots/bdt_ntuple_nue.root")
    nue_tree = f_nue.Get("bdt_ntuple")
    f_ext = ROOT.TFile("plots/bdt_ntuple_bnbext.root")
    ext_tree = f_ext.Get("bdt_ntuple")

    dataloader.AddSignalTree(nue_tree)
    dataloader.AddBackgroundTree(nue_tree)
    dataloader.AddBackgroundTree(mc_tree)
    dataloader.AddBackgroundTree(ext_tree)

    sigCut = ROOT.TCut("category == 10 || category == 2")
    bgCut = ROOT.TCut("category != 10 && category != 2")

    dataloader.PrepareTrainingAndTestTree(sigCut, bgCut,
                                          ":".join([
                                              "SplitMode=Alternate",
                                              "NormMode=NumEvents:!V"]))

    method_svm = factory.BookMethod(dataloader, ROOT.TMVA.Types.kSVM, "SVM", "")
                                    # ":".join(["!H:!V:NTrees=600",
                                    #           "MinNodeSize=5%",
                                    #           "MaxDepth=5",
                                    #           "BoostType=AdaBoost",
                                    #           "AdaBoostBeta=0.5",
                                    #           "UseBaggedBoost",
                                    #           "BaggedSampleFraction=0.5",
                                    #           "SeparationType=GiniIndex",
                                    #           "nCuts=20"]))

    factory.TrainAllMethods()
    factory.TestAllMethods()
    factory.EvaluateAllMethods()
    fout.Close()


train_SVM("svm",  ["cosmic", "numu", "nc"])
