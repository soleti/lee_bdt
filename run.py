#!/usr/local/bin/python3

import ROOT
from bdt_common import binning, labels, variables, spectators, fill_histos_data
from bdt_common import bins, colors, interactions, load_variables
from bdt_common import pre_cuts, fix_binning
from bdt_common import bdt_types, apply_cuts, load_bdt, printProgressBar
from bdt_common import pdgs
from array import array
import pickle

ROOT.gStyle.SetOptStat(0)

DRAW_SYS = True
DRAW_PLOTS = True

h_pdgs = {}
vars = dict(variables + spectators)

for n, b in vars.items():
    pdg_dict = {}
    for pdg in pdgs:
        pdg_dict[pdgs[pdg]] = ROOT.TH1F("h_%s_%s" % (n, pdg), pdg+labels[n], binning[n][0], binning[n][1], binning[n][2])
    h_pdgs[n] = pdg_dict

h_int = {}
for name in interactions:
    h_int[interactions[name]] = ROOT.TH1F(name, ";E_{corr};N. Entries / 0.05 GeV", len(bins) - 1, bins)

h_angle_energy_bkg = ROOT.TH2F("h_angle_energy_bkg",
                               ";Shower opening angle [#circ];Shower energy [GeV]",
                               45, 0, 45, 20, 0.05, 0.5)
h_angle_energy_sig = ROOT.TH2F("h_angle_energy_sig",
                               ";Shower opening angle [#circ];Shower energy [GeV]",
                               45, 0, 45, 20, 0.05, 0.5)

def fill_histos(chain, histo_dict, h_bdt_types, option=""):
    ROOT.TMVA.Tools.Instance()
    reader = ROOT.TMVA.Reader(":".join([
        "!V",
        "!Silent",
        "Color"]))

    bdt_ntuple = ROOT.TNtuple(
        "bdt_ntuple", "bdt_ntuple", "single"+":".join(bdt_types) + ":category")

    variables_dict = load_variables(chain)
    load_bdt(reader)

    h_energy_sig = ROOT.TH1F("h_energy_sig_%s" % option,
                             ";E_{#nu_{e}};N. Entries / 0.05 GeV",
                             len(bins) - 1,
                             bins)
    h_energy_bkg = ROOT.TH1F("h_energy_bkg_%s" % option,
                             ";E_{#nu_{e}};N. Entries / 0.05 GeV",
                             len(bins) - 1,
                             bins)
    h_reco_sig = ROOT.TH1F("h_reco_sig_%s" % option,
                           ";Reco. energy [GeV];N. Entries / 0.05 GeV", len(bins) - 1, bins)
    h_reco_bkg = ROOT.TH1F("h_reco_bkg_%s" % option,
                           ";Reco. energy [GeV];N. Entries / 0.05 GeV", len(bins) - 1, bins)


    passed_events = 0
    entries = int(chain.GetEntries() / 1)
    for i in range(entries):

        # printProgressBar(i, entries, prefix="Progress:",
        #                  suffix="Complete", length=20)

        chain.GetEntry(i)
        category = int(chain.category)

        # cuts_response = reader.EvaluateMVA("Cuts method", rectangular_cut)
        # likelihood_response = reader.EvaluateMVA("Likelihood method")

        bdts = {}
        for bdt_name in bdt_types:
            bdts[bdt_name] = reader.EvaluateMVA("BDT%s" % bdt_name)

        if pre_cuts(variables_dict):
            for bdt_name in bdt_types:
                h_bdt_types[bdt_name][category].Fill(bdts[bdt_name], chain.event_weight)
        else:
            continue

        bdt_ntuple.Fill(array("f", list(bdts.values())+[chain.category]))

        if apply_cuts(bdts, variables_dict):
            if not chain.true_nu_is_fidvol and category != 0 and category != 6 and category != 1 and category != 7:
                category = 5
            # if category == 4:
            #     print(int(chain.run), int(chain.subrun), int(chain.event))
            if option == "nue":
                if category == 2 and chain.reco_energy < 3:
                    passed_events += chain.event_weight
            else:
                passed_events += chain.event_weight
            corr = 1
            # if category == 2 and not manual:
            #     corr = 1.2

            if chain.true_nu_is_fidvol:
                h_energy_sig.Fill(chain.nu_E, chain.event_weight * corr)

            h_reco_sig.Fill(chain.reco_energy, chain.event_weight * corr)

            if (option == "mc" or option == "bnbext" or (option == "nue" and category != 2)):
                if chain.true_nu_is_fidvol:
                    h_energy_bkg.Fill(chain.nu_E, chain.event_weight * corr)
                h_reco_bkg.Fill(chain.reco_energy, chain.event_weight * corr)

            all = variables + spectators
            var_dict = dict(all)
            for name, var in all:
                for i, v in enumerate(var):
                    if v > -999:
                        if "track" in name and name != "no_tracks":
                            pdg_code = int(var_dict["track_pdg"][i])
                            if option != "bnbext" and abs(pdg_code) == 2147483648:
                                pdg_code = int(var_dict["shower_pdg"][0])
                        elif "shower" in name:
                            pdg_code = int(var_dict["shower_pdg"][i])
                            if option != "bnbext" and abs(pdg_code) == 2147483648:
                                pdg_code = int(var_dict["track_pdg"][0])

                        if option != "bnbext" and abs(pdg_code) == 2147483648:
                            pdg_code = 99999

                        if abs(pdg_code) in h_pdgs[name]:
                            h_pdgs[name][abs(pdg_code)].Fill(v, chain.event_weight)
                        else:
                            h_pdgs[name][99999].Fill(v, chain.event_weight)

                        histo_dict[name][category].Fill(v, chain.event_weight)


    f = ROOT.TFile("plots/bdt_ntuple_%s.root" % option,
                    "RECREATE")
    bdt_ntuple.Write()
    f.Close()

    return passed_events

mc_chain = ROOT.TChain("mc_tree")
nue_chain = ROOT.TChain("nue_tree")
bnbext_chain = ROOT.TChain("bnbext_tree")

mc_chain.Add("root_files/mc_file.root")
nue_chain.Add("root_files/nue_file.root")
bnbext_chain.Add("root_files/bnbext_file.root")

kinds = []
categories = ["intime", "cosmic", "nu_e", "nu_mu", "nc", "dirt", "data",
              "mixed", "cc0pi"]
stacked_histos = []

variables_dict = dict(variables + spectators)

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

print("nu_e events", fill_histos(nue_chain, histo_dict, h_bdt_types, "nue"))
print("LEE events", fill_histos_data("lee"))
print("BNB + cosmic events", fill_histos(mc_chain, histo_dict, h_bdt_types, "mc"))
print("Data events", fill_histos_data("bnb"))
print("EXT events", fill_histos(bnbext_chain, histo_dict, h_bdt_types, "bnbext"))

h_interactions = ROOT.THStack(
    "h_interactions", ";E_{corr};N. Entries / 0.05 GeV")
l_interactions = ROOT.TLegend(0.09, 0.7747, 0.904, 0.8968, "", "brNDC")

l_interactions.SetNColumns(3)
for interaction in h_int:
    histo = h_int[interaction]
    if histo.Integral() > 0:
        histo.SetLineColor(1)
        h_interactions.Add(histo)

h_true_e = ROOT.THStack(
    "h_true_e", ";E_{corr} [GeV]; N. Entries / 0.05 GeV")

for j in range(h_interactions.GetNhists()):
    h_clone = h_interactions.GetHists()[j].Clone()
    h_clone.Scale(1)
    fix_binning(h_clone)

    h_fixed = ROOT.TH1F("h_fixed%i" % j, "", len(bins) - 1, 0, len(bins) - 1)

    for i in range(1, h_clone.GetNbinsX() + 1):
        h_fixed.SetBinContent(i, h_clone.GetBinContent(i))
        h_fixed.SetBinError(i, h_clone.GetBinError(i))
        h_fixed.GetXaxis().SetBinLabel(i, "")

    h_fixed.SetLineColor(ROOT.kBlack)
    l_interactions.AddEntry(h_fixed, "%s: %.1f events" % (
        h_clone.GetName(), h_interactions.GetHists()[j].Integral()), "f")

    h_true_e.Add(h_fixed)

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
    f = ROOT.TFile("plots/%s_mc.root" % h.GetName(), "RECREATE")
    h.Write()
    f.Close()

pickle.dump(h_pdgs, open("plots/pdg_plots.p", "wb"))
