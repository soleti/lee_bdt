#!/usr/local/bin/python3

import ROOT

from bdt_common import bdt_cut, binning, labels, variables, spectators, fill_histos_data
from bdt_common import sigmaCalc, manual_cuts, bdt, manual, bins, colors

print("BDT: ", bdt, "Manual: ", manual)


def fill_histos(chain, histo_dict, h_bdts):
    events = 0

    for i in range(chain.GetEntries()):
        chain.GetEntry(i)
        category = int(chain.category)
        h_bdts[category].Fill(chain.BDT, chain.event_weight)

        if bdt:
            apply_bdt = chain.BDT > bdt_cut
        else:
            apply_bdt = True

        if manual:
            apply_manual = manual_cuts(chain)
        else:
            apply_manual = True

        if apply_bdt and apply_manual:
            if chain.category == 2 and bins[0] < chain.reco_energy < bins[-1]:
                events += ttrain.event_weight

            for name, var in variables:
                histo_dict[name][category].Fill(var[0], chain.event_weight)

            for name, var in spectators:
                histo_dict[name][category].Fill(var[0], chain.event_weight)

            # if chain.category != 2 and 0 < chain.reco_energy < 2:
            #     print("nu_mu",
            #           int(chain.run), int(chain.subrun), int(chain.event),
            #           int(chain.category),
            #           inv_interactions[int(chain.interaction_type)])

    return events


tout = ROOT.TChain("dataset/TestTree")
ttrain = ROOT.TChain("dataset/TrainTree")
ttrain.Add("test.root")
tout.Add("test.root")

kinds = []
categories = ["intime", "cosmic", "nu_e", "nu_mu", "nc", "dirt", "data",
              "mixed"]
stacked_histos = []

for name, var in variables:
    tout.SetBranchAddress(name, var)
    ttrain.SetBranchAddress(name, var)

for name, var in spectators:
    tout.SetBranchAddress(name, var)
    ttrain.SetBranchAddress(name, var)

variables_dict = dict(variables + spectators)

for i, n in enumerate(variables_dict.keys()):
    histos = []

    h_stack = ROOT.THStack("h_" + n, labels[n])

    for c in categories:
        if n != "reco_energy":
            h = ROOT.TH1F("h_%s_%s" % (n, c), labels[n],
                          binning[n][0], binning[n][1], binning[n][2])
        else:
            h = ROOT.TH1F("h_%s_%s" % (n, c), labels[n], len(bins) - 1, bins)
        histos.append(h)

    h_stack.Add(histos[2])
    h_stack.Add(histos[3])
    h_stack.Add(histos[4])
    h_stack.Add(histos[5])
    h_stack.Add(histos[7])
    h_stack.Add(histos[1])
    h_stack.Add(histos[0])
    h_stack.Add(histos[6])

    stacked_histos.append(h_stack)
    kinds.append(histos)

histo_dict = dict(zip(variables_dict.keys(), kinds))

h_bdt_stack = ROOT.THStack("h_bdt", ";BDT response; N. Entries / 0.05")
h_bdts = []

for i, c in enumerate(categories):
    h = ROOT.TH1F("h_bdt_%s" % c, "BDT response; N. Entries / 0.05", 40, -1, 1)
    h.SetLineColor(1)
    h.SetFillColor(colors[i])
    h.SetLineWidth(2)
    h_bdts.append(h)

h_bdt_stack.Add(h_bdts[2])
h_bdt_stack.Add(h_bdts[3])
h_bdt_stack.Add(h_bdts[4])
h_bdt_stack.Add(h_bdts[5])
h_bdt_stack.Add(h_bdts[7])
h_bdt_stack.Add(h_bdts[1])
h_bdt_stack.Add(h_bdts[0])
h_bdt_stack.Add(h_bdts[6])

h_bkg = ROOT.TH1F("h_bkg", "", 16, 0.2, 1)
h_sig = ROOT.TH1F("h_sig", "", 16, 0.2, 1)

best_sigma = 0
chosen_cut = 0

measure_sigma = False

if measure_sigma:
    for i in range(200):

        cut = 0.2 + i * 0.004

        h_bkg.Reset()
        h_sig.Reset()

        for i in range(tout.GetEntries()):
            tout.GetEntry(i)

            if tout.BDT > cut:

                if tout.category == 2:
                    h_sig.Fill(tout.reco_energy, tout.event_weight)
                else:
                    h_bkg.Fill(tout.reco_energy, tout.event_weight)

        for i in range(ttrain.GetEntries()):
            ttrain.GetEntry(i)

            if ttrain.BDT > cut:
                if ttrain.category == 2:
                    h_sig.Fill(ttrain.reco_energy, ttrain.event_weight)
                else:
                    h_bkg.Fill(ttrain.reco_energy, ttrain.event_weight)

        f_cosmic_in_time = ROOT.TFile("bdt_cosmic/cut%.3f.root" % cut)
        h_cosmic_in_time = f_cosmic_in_time.Get("h_reco")
        h_bkg.Add(h_cosmic_in_time)

        if sigmaCalc(h_sig, h_bkg) > best_sigma:
            best_sigma = sigmaCalc(h_sig, h_bkg)
            h_chosen_bkg = h_bkg.Clone()
            chosen_cut = cut
        f_cosmic_in_time.Close()
    print(best_sigma, chosen_cut)


print(fill_histos(ttrain, histo_dict, h_bdts))
print(fill_histos(tout, histo_dict, h_bdts))

f_bdt = ROOT.TFile("plots/h_bdt.root", "RECREATE")
h_bdt_stack.Write()
f_bdt.Close()

for i, histos in enumerate(kinds):
    for j, h in enumerate(histos):
        h.SetLineColor(1)
        h.SetLineWidth(2)
        h.SetFillColor(colors[j])

for h in stacked_histos:
    f = ROOT.TFile("plots/%s_mc.root" % h.GetName(), "RECREATE")
    h.Write()
    f.Close()

fill_histos_data("bnb", bdt, manual)
fill_histos_data("bnbext", bdt, manual)
fill_histos_data("lee", bdt, manual)
