#!/usr/local/bin/python3

import ROOT
from tqdm import tqdm
from array import array
import sys
from glob import glob

sample = sys.argv[1]
events_file = sys.argv[2]
sbnfit_chain = ROOT.TChain("robertoana/sbnfit")
files = glob("data_files/%s/*.root" % sample)

for f in files:
    sbnfit_chain.Add(f)

print(sample)
filtered_file = ROOT.TFile("root_files/filtered_events_%s.root" % sample, "RECREATE")
filtered_trees = []
branches = []
entries = int(sbnfit_chain.GetEntries() / 1)

categories = ["bnbext", "cosmic", "nu_e_cc0pinp", "nu_mu",
              "nc", "outfv", "bnb", "contaminated", "nu_e_cc", "other", "lee"]

reco_energy = array("d", [0])

for i in range(11):
    filtered_tree = sbnfit_chain.GetTree().CloneTree(0)
    branch = filtered_tree.Branch("reco_energy", reco_energy, "reco_energy/D")
    branches.append(branch)
    filtered_tree.SetName("%s" % categories[i])
    filtered_trees.append(filtered_tree)


selected_events = open("selected_events/%s.txt" % events_file).readlines()
events = [" ".join(s.split()[:3]) for s in selected_events]
cats = [s.split()[3] for s in selected_events]
energies = [s.split()[4] for s in selected_events]

for i in tqdm(range(entries)):
    sbnfit_chain.GetEntry(i)
    event = "%i %i %i" % (sbnfit_chain.run, sbnfit_chain.subrun, sbnfit_chain.event)
    if event in events:
        cat = int(cats[events.index(event)])
        reco_energy[0] = float(energies[events.index(event)])
        branches[cat].Fill()
        filtered_trees[cat].Fill()

filtered_file.cd()
for i in range(11):
    filtered_trees[i].Write()

