#!/usr/local/bin/python3

import sys
import pickle
from array import array

from tqdm import tqdm
import ROOT

from bdt_common import binning, labels, variables, spectators
from bdt_common import bins, colors, interactions, load_variables
from bdt_common import pre_cuts, fix_binning
from bdt_common import bdt_types, apply_cuts, load_bdt
from bdt_common import pdgs, BDT, MANUAL

if len(sys.argv) > 1:
    mode = sys.argv[1]

if mode in ("numu", "nc"):
    BDT, MANUAL = False, True

if mode == "nue":
    if len(sys.argv) > 2:
        if sys.argv[2] == "cuts":
            BDT, MANUAL = False, True
        elif sys.argv[2] == "bdt":
            BDT, MANUAL = True, False
        else:
            BDT, MANUAL = False, False


ROOT.gStyle.SetOptStat(0)

DRAW_SYS = True
DRAW_PLOTS = True


def fill_histos_data(tree_name, mode="nue"):
    f_data = ROOT.TFile("root_files/%s_file.root" % tree_name)
    t_data = f_data.Get("%s_tree" % tree_name)

    ROOT.TMVA.Tools.Instance()
    reader = ROOT.TMVA.Reader(":".join([
        "!V",
        "!Silent",
        "Color"]))

    variables_dict = load_variables(t_data)
    load_bdt(reader)

    histograms = []

    for i, n in enumerate(variables_dict.keys()):
        if n != "reco_energy" and n != "nu_E":
            h = ROOT.TH1F("h_%s" % n, labels[n],
                          binning[n][0], binning[n][1], binning[n][2])
        else:
            h = ROOT.TH1F("h_%s" % n, labels[n], len(bins) - 1, bins)
        histograms.append(h)

    histo_dict = dict(zip(variables_dict.keys(), histograms))

    h_bdts = {}
    for bdt_name in bdt_types:
        h_bdts[bdt_name] = ROOT.TH1F("h_bdt_%s_%s" % (bdt_name, tree_name),
                                     "BDT %s response; N. Entries / 0.04" % bdt_name,
                                     50, -1, 1)

    passed_events = 0
    entries = int(t_data.GetEntries() / 1)

    # bdt_ntuple = ROOT.TNtuple(
    #     "bdt_ntuple", "bdt_ntuple", "single"+":".join(bdt_types) + ":category")
    f_clone = ROOT.TFile("root_files/filtered_%s.root" % tree_name, "RECREATE")
    clone = t_data.CloneTree(0)

    for i_evt in tqdm(range(entries)):
        t_data.GetEntry(i_evt)
        bdt_values = {}

        for bdt_name in bdt_types:
            bdt_values[bdt_name] = reader.EvaluateMVA("BDT%s" % bdt_name)

        numu = mode == "numu"

        if pre_cuts(variables_dict, numu):
            for bdt_name in bdt_types:
                h_bdts[bdt_name].Fill(bdt_values[bdt_name], t_data.event_weight)

            if apply_cuts(bdt_values, variables_dict, BDT, MANUAL, mode):
                passed_events += t_data.event_weight
                all_vars = variables + spectators
                clone.Fill()
                for name, var in all_vars:
                    for v in var:
                        if v == -999:
                            break
                        histo_dict[name].Fill(v, t_data.event_weight)

    clone.Write()
    f_clone.Close()

    f_bdt = ROOT.TFile("plots/h_bdt_%s.root" % tree_name, "RECREATE")
    for h in h_bdts:
        h_bdts[h].Write()
    f_bdt.Close()

    for h in histograms:
        if MANUAL:
            folder = "_cuts"
        elif BDT:
            folder = "_bdt"
        else:
            folder = ""
        if mode == "numu":
            folder = "_numu"
        elif mode == "nc":
            folder = "_nc"
        f = ROOT.TFile("plots%s/%s_%s.root" % (folder, h.GetName(), tree_name),
                       "RECREATE")
        h.Write()
        f.Close()

    return passed_events

def fill_histos(chain, histo_dict, h_bdt_types, option="", mode="nue"):
    ROOT.TMVA.Tools.Instance()
    reader = ROOT.TMVA.Reader(":".join([
        "!V",
        "!Silent",
        "Color"]))

    bdt_ntuple = ROOT.TNtuple(
        "bdt_ntuple", "bdt_ntuple", "single"+":".join(bdt_types) + ":category")

    variables_dict = load_variables(chain)
    load_bdt(reader)

    passed_events = 0
    entries = int(chain.GetEntries() / 1)
    f = open("selected_events/%s_passed.txt" % option, "w")
    clone = chain.GetTree().CloneTree(0)
    for i in tqdm(range(entries)):
        chain.GetEntry(i)
        category = int(chain.category)

        # cuts_response = reader.EvaluateMVA("Cuts method", rectangular_cut)
        # likelihood_response = reader.EvaluateMVA("Likelihood method")

        bdts = {}
        for bdt_name in bdt_types:
            bdts[bdt_name] = reader.EvaluateMVA("BDT%s" % bdt_name)

        numu = mode == "numu"

        if pre_cuts(variables_dict, numu):
            for bdt_name in bdt_types:
                h_bdt_types[bdt_name][category].Fill(bdts[bdt_name], chain.event_weight)
        else:
            continue

        bdt_ntuple.Fill(array("f", list(bdts.values())+[chain.category]))

        if apply_cuts(bdts, variables_dict, BDT, MANUAL, mode):
            if not chain.true_nu_is_fidvol and category != 0 and category != 6 and category != 1 and category != 7:
                category = 5
            clone.Fill()

            print(int(chain.run), int(chain.subrun), int(chain.event), int(category), chain.reco_energy, file=f)
            if option == "nue":
                if category == 2:
                    passed_events += chain.event_weight
            else:
                passed_events += chain.event_weight

            all_vars = variables + spectators
            var_dict = dict(all_vars)
            for name, var in all_vars:
                for i_v, v in enumerate(var):
                    if v == -999:
                        break
                    if "track" in name and name != "no_tracks":
                        pdg_code = int(var_dict["track_pdg"][i_v])
                        if option != "bnbext" and abs(pdg_code) == 2147483648:
                            pdg_code = int(var_dict["shower_pdg"][0])
                    elif "shower" in name:
                        pdg_code = int(var_dict["shower_pdg"][i_v])
                        if option != "bnbext" and abs(pdg_code) == 2147483648:
                            pdg_code = int(var_dict["track_pdg"][0])

                    if option != "bnbext" and abs(pdg_code) == 2147483648:
                        pdg_code = 99999

                    if abs(pdg_code) in h_pdgs[name]:
                        h_pdgs[name][abs(pdg_code)].Fill(v, chain.event_weight)
                    else:
                        h_pdgs[name][99999].Fill(v, chain.event_weight)

                    histo_dict[name][category].Fill(v, chain.event_weight)

    f.close()
    f = ROOT.TFile("plots/bdt_ntuple_%s.root" % option,
                    "RECREATE")
    clone.Write()
    bdt_ntuple.Write()
    f.Close()

    return passed_events


if __name__ == "__main__":
    h_pdgs = {}
    variables_dict = dict(variables + spectators)

    for n, b in variables_dict.items():
        pdg_dict = {}
        for pdg in pdgs:
            pdg_dict[pdgs[pdg]] = ROOT.TH1F("h_%s_%s" % (
                n, pdg), pdg+labels[n], binning[n][0], binning[n][1], binning[n][2])
        h_pdgs[n] = pdg_dict

    h_int = {}
    for name_int in interactions:
        h_int[interactions[name_int]] = ROOT.TH1F(
            name_int, ";E_{corr};N. Entries / 0.05 GeV", len(bins) - 1, bins)

    mc_chain = ROOT.TChain("mc_tree")
    nue_chain = ROOT.TChain("nue_tree")
    bnbext_chain = ROOT.TChain("bnbext_tree")

    mc_chain.Add("root_files/mc_file.root")
    mc_chain.Add("root_files/dirt.root")

    nue_chain.Add("root_files/nue_file.root")
    bnbext_chain.Add("root_files/bnbext_file.root")

    kinds = []
    categories = ["intime", "cosmic", "nu_e", "nu_mu", "nc", "dirt", "data",
                "mixed", "cc0pi"]
    stacked_histos = []


    for i, n in enumerate(variables_dict.keys()):
        histos = []

        h_stack = ROOT.THStack("h_" + n, labels[n])

        for c in categories:
            if n != "reco_energy" and n != "nu_E":
                h = ROOT.TH1F("h_%s_%s" % (n, c), labels[n],
                            binning[n][0], binning[n][1], binning[n][2])
            else:
                h = ROOT.TH1F("h_%s_%s" % (n, c), labels[n], len(bins) - 1, bins)
            histos.append(h)

        h_stack.Add(histos[0])
        h_stack.Add(histos[3])
        h_stack.Add(histos[4])
        h_stack.Add(histos[5])
        h_stack.Add(histos[7])
        h_stack.Add(histos[1])
        h_stack.Add(histos[6])
        h_stack.Add(histos[8])
        h_stack.Add(histos[2])
        stacked_histos.append(h_stack)
        kinds.append(histos)


    histo_dict = dict(zip(variables_dict.keys(), kinds))
    h_bdt_types = {}
    h_bdt_stacks = {}
    for bdt_type in bdt_types:
        h_bdt_stack = ROOT.THStack("h_bdt_%s" % bdt_type, ";BDT response; N. Entries / 0.04")
        h_bdts = []

        for i, c in enumerate(categories):
            h = ROOT.TH1F("h_bdt_%s_%s" % (bdt_type, c), "BDT response; N. Entries / 0.04", 50, -1, 1)
            h.SetLineColor(1)
            h.SetFillColor(colors[i])
            if i != 0:
                h.SetLineWidth(0)
            h_bdts.append(h)

        h_bdt_stack.Add(h_bdts[0])
        h_bdt_stack.Add(h_bdts[3])
        h_bdt_stack.Add(h_bdts[4])
        h_bdt_stack.Add(h_bdts[5])
        h_bdt_stack.Add(h_bdts[7])
        h_bdt_stack.Add(h_bdts[1])
        h_bdt_stack.Add(h_bdts[6])
        h_bdt_stack.Add(h_bdts[8])
        h_bdt_stack.Add(h_bdts[2])
        h_bdt_stacks[bdt_type] = h_bdt_stack
        h_bdt_types[bdt_type] = h_bdts

    print("Data events", fill_histos_data("bnb", mode))
    print("nu_e events", fill_histos(nue_chain, histo_dict, h_bdt_types, "nue", mode))
    print("LEE events", fill_histos_data("lee", mode))
    print("BNB + cosmic events", fill_histos(mc_chain, histo_dict, h_bdt_types, "mc", mode))
    print("EXT events", fill_histos(bnbext_chain, histo_dict, h_bdt_types, "bnbext", mode))

    f_bdt = ROOT.TFile("plots/h_bdt.root", "RECREATE")
    for bdt_type in bdt_types:
        h_bdt_stacks[bdt_type].Write()
    f_bdt.Close()

    for i, histos in enumerate(kinds):
        for j, h in enumerate(histos):
            h.SetLineColor(1)
            if j != 0:
                h.SetLineWidth(0)
            h.SetFillColor(colors[j])

    for h in stacked_histos:
        if MANUAL:
            folder = "_cuts"
        elif BDT:
            folder = "_bdt"
        else:
            folder = ""
        if mode == "numu":
            folder = "_numu"
        elif mode == "nc":
            folder = "_nc"
        f = ROOT.TFile("plots%s/%s_mc.root" % (folder, h.GetName()), "RECREATE")
        h.Write()
        f.Close()

    pickle.dump(h_pdgs, open("plots%s/pdg_plots.p" % folder, "wb"))
