#!/usr/bin/env python3.4

import math
import ROOT
from glob import glob

ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetPalette(87)
ROOT.gStyle.SetNumberContours(99)

nue_cosmic = glob("nue_files_fidvol/*/*.root")

chain_nue = ROOT.TChain("robertoana/pandoratree")
chain_nue_pot = ROOT.TChain("robertoana/pot")

for f in nue_cosmic:
    chain_nue.Add(f)
    chain_nue_pot.Add(f)

total_pot = 6.6e20
total_nue_pot = 0
for i in range(chain_nue_pot.GetEntries()):
    chain_nue_pot.GetEntry(i)
    total_nue_pot += chain_nue_pot.pot
print("Total POT v_e", total_nue_pot)

e_energy = ROOT.TEfficiency("e_energy",";#nu_{e} energy [GeV];Efficiency",20,0,2)

is_fiducial = 0
eNp = 0
yesshower_notrack = 0
noshower_yestrack = 0
noshower_notrack = 0
passed = 0
not_passed = 0
flash_passed = 0

entries = chain_nue.GetEntries()
for i in range(entries):
    chain_nue.GetEntry(i)

    protons = sum(1 for i,pdg in enumerate(chain_nue.nu_daughters_pdg) if abs(pdg) == 2212)
    electrons = sum(1 for i,pdg in enumerate(chain_nue.nu_daughters_pdg) if abs(pdg) == 11)
    photons = sum(1 for i in chain_nue.nu_daughters_pdg if i == 22)
    pions = sum(1 for i in chain_nue.nu_daughters_pdg if abs(i) == 211 or abs(i) == 111)


    if electrons > 0 and photons == 0 and pions == 0 and protons > 0:
        eNp+=1

        if chain_nue.true_nu_is_fiducial:
            is_fiducial += 1

            if chain_nue.flash_passed:
                flash_passed += 1

            p = False
            p_track = False
            p_shower = False

            if protons == chain_nue.n_tracks and chain_nue.n_tracks == chain_nue.nu_matched_tracks: p_track = True

            if electrons == chain_nue.n_showers and chain_nue.n_showers == chain_nue.nu_matched_showers: p_shower = True

            if p_track and p_shower:
                p = True

            proton_energy = sum([chain_nue.nu_daughters_E[i] for i,pdg in enumerate(chain_nue.nu_daughters_pdg) if pdg == 2212])
            electron_energy = max([chain_nue.nu_daughters_E[i] for i,pdg in enumerate(chain_nue.nu_daughters_pdg) if abs(pdg) == 11])

            if chain_nue.passed:
                passed+=1
            else:
                not_passed+=1
                if chain_nue.shower_passed > 0 and chain_nue.track_passed <= 0:
                    yesshower_notrack += 1
                if chain_nue.track_passed > 0 and chain_nue.shower_passed <= 0:
                    noshower_yestrack += 1
                if chain_nue.flash_passed and chain_nue.track_passed <= 0 and chain_nue.shower_passed <= 0:
                    noshower_notrack += 1

            e_energy.Fill(chain_nue.passed, chain_nue.nu_E)

print("Entries", entries)
print("1eNp", eNp)
print("1eNp + Is fiducial", is_fiducial)
print("1eNp + Is fiducial + Flash passed", flash_passed)
print("Passed", passed)

eff = passed/is_fiducial
eff_err = math.sqrt((eff*(1-eff))/eNp)

print("Efficiency: ({0:.1f} +- {1:.1f}) %".format(eff*100, eff_err*100))

f_energy = ROOT.TFile("f_energy.root","RECREATE")
e_energy.Write()
f_energy.Close()

pt = ROOT.TPaveText(0.1,0.91,0.45,0.97, "ndc")
pt.AddText("MicroBooNE Preliminary")
pt.SetFillColor(0)
pt.SetBorderSize(0)
pt.SetShadowColor(0)

c_energy = ROOT.TCanvas("c_energy")
e_energy.Draw("apl")
e_energy.SetMarkerStyle(20)
e_energy.SetLineColor(ROOT.kRed+1)
e_energy.SetLineWidth(2)
c_energy.Update()
e_energy.GetPaintedGraph().GetXaxis().SetRangeUser(0.1,2)
e_energy.GetPaintedGraph().GetYaxis().SetRangeUser(0,1.1)
pt.Draw()
c_energy.SaveAs("plots/energy.pdf")
c_energy.Draw()
input()
