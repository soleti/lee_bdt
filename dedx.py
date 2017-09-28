#!/usr/bin/env python

import math
import ROOT
import landau
from glob import glob
from numpy import array


def langau(var, par):
    x = array([var[0]])
    A = par[0]
    mpv = par[1]
    landau_sigma = par[2]
    gauss_sigma = par[3]
    return landau.gauss_landau(x, mpv, landau_sigma, gauss_sigma, scale=A)[0]


def langau_lin(var, par):
    x = array([var[0]])
    A = par[0]
    mpv = par[1]
    landau_sigma = par[2]
    gauss_sigma = par[3]
    c = par[4]
    return landau.gauss_landau(x, mpv, landau_sigma, gauss_sigma,
                               scale=A)[0] + c


ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptFit(0)

h_dedx_electron = ROOT.TH1F("h_dedx_electron",
                            "Electrons;dE/dx [MeV/cm];Area normalized",
                            50, 0, 8)
h_dedx_photon = ROOT.TH1F("h_dedx_photon",
                          "Photons;dE/dx [MeV/cm];Area normalized",
                          50, 0, 8)

h_dedx_hits_electron = ROOT.TH1F("h_dedx_hits_electron",
                                 "Electron hits;dE/dx [MeV/cm];Area normalized",
                                 50, 0, 8)

h_dedx_hits_photon = ROOT.TH1F("h_dedx_hits_photon",
                               "Photon hits;dE/dx [MeV/cm];Area normalized",
                               50, 0, 8)

h_dedx_hits_energy = ROOT.TH2F("h_dedx_hits_energy",
                               "Electron hits;dE/dx [MeV/cm];#pi/#gamma energy [GeV]",
                               50, 0, 8, 50, 0, 1)

h_p_dedx = ROOT.TH2F("h_dedx_p", ";p [GeV/c];dE/dx [MeV/cm]",
                     100, 0, 2, 100, 1, 8)

nue_cosmic = glob("mc_nue_nofidvol/*/Pandora*.root")
chain = ROOT.TChain("robertoana/pandoratree")


for f in nue_cosmic:
    chain.Add(f)

with open("mc_passed.txt","r") as f:
    lines = f.readlines()

events = [" ".join(line.split()[:3]) for line in lines]
weights = [float(line.split()[3]) for line in lines]

entries = chain.GetEntries()
for i in range(entries):
    chain.GetEntry(i)
    gamma = 0
    e = 0
    pions = 0
    gamma_energy = 0
    pions_energy = 0
    run_subrun_event = "{} {} {}".format(int(chain.run), int(chain.subrun), int(chain.event))
    if 1:
        weight = 1#weights[events.index(run_subrun_event)]

        for i in range(len(chain.nu_daughters_E)):
            if abs(chain.nu_daughters_pdg[i]) == 11:
                e += 1
            if abs(chain.nu_daughters_pdg[i]) == 22:
                gamma += 1
                gamma_energy += chain.nu_daughters_E[i]
            if abs(chain.nu_daughters_pdg[i]) == 111:
                pions += 1
                pions_energy += chain.nu_daughters_E[i]

        reco_gamma = 0
        reco_electron = 0

        if chain.passed:
            neutrino_vertex = [chain.vx, chain.vy, chain.vz]
            for ish in range(chain.n_showers):
                pdg = chain.matched_showers[ish]
                dedx = chain.shower_dEdx[ish][2]
                print(chain.matched_showers_process[ish])
                if abs(pdg) == 22:

                    for ihit in range(len(chain.dEdx_hits[ish])):
                        h_dedx_hits_photon.Fill(chain.dEdx_hits[ish][ihit], weight)

                    h_dedx_photon.Fill(dedx, weight)

                if abs(pdg) == 11:

                    for ihit in range(len(chain.dEdx_hits[ish])):
                        h_dedx_hits_electron.Fill(chain.dEdx_hits[ish][ihit], weight)

                    h_dedx_electron.Fill(dedx, weight)

        #
        # for i in range(chain.n_showers):
        #     shower_vertex = [chain.shower_start_x[i],
        #                      chain.shower_start_y[i],
        #                      chain.shower_start_z[i]]
        #
        #     distance = math.sqrt(sum([
        #         (s - n)**2 for s, n in zip(shower_vertex, neutrino_vertex)]))
        #
        #     pdg = chain.matched_showers[i]
        #     # print(chain.matched_showers_process[i])
        #     dedx = chain.shower_dEdx[i][2]
        #     energy = chain.shower_energy[i]
            # if abs(pdg) == 22:
            #     h_p_dedx.Fill(energy, dedx)
            #     h_dedx_electron.Fill(dedx)
            #     for ish in range(chain.n_showers):
            #         for ihit in range(len(chain.dEdx_hits[ish])):
            #             h_dedx_hits.Fill(chain.dEdx_hits[ish][ihit])
            #
            #             h_dedx_hits_energy.Fill(chain.dEdx_hits[ish][ihit],
            #                                     pions_energy + gamma_energy)


c_dedx = ROOT.TCanvas("c_dedx")
hs = ROOT.THStack("hs", ";dE/dx [MeV/cm];Area normalized")
h_dedx_photon.Scale(h_dedx_electron.Integral() / h_dedx_photon.Integral())
h_dedx_photon.SetLineColor(ROOT.kRed + 1)
h_dedx_electron.SetLineColor(ROOT.kBlue + 1)
hs.Add(h_dedx_electron)
hs.Add(h_dedx_photon)
hs.Draw("hist nostack")
c_dedx.Update()

# c_p_dedx = ROOT.TCanvas("c_p_dedx")
# h_p_dedx.Draw("colz")
# c_p_dedx.Update()

c_dedx_hits = ROOT.TCanvas("c_dedx_hits")
h_dedx_hits_electron.SetMarkerStyle(20)
h_dedx_hits_electron.SetLineColor(ROOT.kBlack)
h_dedx_hits_electron.Draw("ep")
h_dedx_hits_electron.Sumw2()
h_dedx_hits_electron.Scale(1 / h_dedx_hits_electron.Integral())
f_dedx_hits_electron = ROOT.TFile("f_dedx_hits_electron.root", "RECREATE")
h_dedx_hits_electron.Write()
f_dedx_hits_electron.Close()

h_dedx_hits_photon.SetMarkerStyle(25)
h_dedx_hits_photon.SetLineColor(ROOT.kBlack)
h_dedx_hits_photon.Draw("ep same")
h_dedx_hits_photon.Sumw2()
h_dedx_hits_photon.Scale(1 / h_dedx_hits_photon.Integral())
f_dedx_hits_photon = ROOT.TFile("f_dedx_hits_photon.root", "RECREATE")
h_dedx_hits_photon.Write()
f_dedx_hits_photon.Close()

f_langau = ROOT.TF1("f_langau", langau_lin, 0, 8, 5)

f_langau.SetNpx(10000)
f_langau.SetParNames("A", "#mu", "Landau #sigma", "Gaussian #sigma", "c")
f_langau.SetParLimits(0, 0, 1)
f_langau.SetParLimits(2, 0, 1)
f_langau.SetParLimits(3, 0, 1)
f_langau.SetParLimits(4, 0, 0.1)
f_langau.SetParameters(1000, 1.7, 0.1, 0.2, 0.1)

# h_dedx_hits_electron.Fit(f_langau, "R")

l_langau = ROOT.TLegend(0.37, 0.65, 0.84, 0.85)
l_langau.AddEntry(f_langau, "Fitted Landau", "l")
l_langau.AddEntry(f_langau, "#mu = %.2f MeV/cm"
                  % f_langau.GetParameter(1), "")
l_langau.AddEntry(f_langau, "Landau #sigma = %.2f MeV/cm"
                  % f_langau.GetParameter(2), "")
l_langau.AddEntry(f_langau, "Gaussian #sigma = %.2f MeV/cm"
                  % f_langau.GetParameter(3), "")
# l_langau.Draw()

c_dedx_hits.Update()
#
# c_dedx_hits_energy = ROOT.TCanvas("c_dedx_hits_energy")
# h_dedx_hits_energy.Draw("colz")
# c_dedx_hits_energy.Update()
raw_input()
