#!/usr/bin/env python

import ROOT
import landau
from glob import glob
from numpy import array
import sys
import math


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


def choose_shower(root_chain):
    most_energetic_shower = 0
    shower_id = 0
    for ish in range(root_chain.n_showers):
        if root_chain.shower_energy[ish] < 3:
            if root_chain.shower_energy[ish] > most_energetic_shower:
                most_energetic_shower = root_chain.shower_energy[ish]
                shower_id = ish
    return shower_id


ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptFit(0)

h_dedx_electron = ROOT.TH1F("h_dedx_electron",
                            "Electrons;dE/dx [MeV/cm];Area normalized",
                            50, 0, 8)
h_dedx_photon = ROOT.TH1F("h_dedx_photon",
                          "Photons;dE/dx [MeV/cm];Area normalized",
                          50, 0, 8)

h_dedx_data = ROOT.TH1F("h_dedx_data", "Data;dE/dx [MeV/cm];Area normalized",
                        50, 0, 8)

h_dedx_hits_electron = ROOT.TH1F("h_dedx_hits_electron",
                                 "Electron hits;dE/dx [MeV/cm];a.u.",
                                 50, 0, 8)

h_dedx_hits_photon = ROOT.TH1F("h_dedx_hits_photon",
                               "Photon hits;dE/dx [MeV/cm];a.u.",
                               50, 0, 8)

h_dedx_hits_other = ROOT.TH1F("h_dedx_hits_other",
                              "Other hits;dE/dx [MeV/cm];a.u.", 50, 0, 8)

h_dedx_hits_data = ROOT.TH1F("h_dedx_hits_data",
                             "Data hits;dE/dx [MeV/cm];a.u.",
                             50, 0, 8)

h_dedx_hits_energy = ROOT.TH2F("h_dedx_hits_energy",
                               "Electron hits;dE/dx [MeV/cm];Energy [GeV]",
                               50, 0, 8, 50, 0, 1)

h_p_dedx = ROOT.TH2F("h_dedx_p", ";p [GeV/c];dE/dx [MeV/cm]",
                     100, 0, 2, 100, 1, 8)


if (len(sys.argv) > 1):
    sample = sys.argv[1]
else:
    sample = "mc"

weight = 1

if sample == "mc":
    print("Monte Carlo")
    weigth = 1.66
    files = glob("mc_bnb_mcc83/*/Pandora*.root")
elif sample == "databnb":
    print("Data BNB")
    weight = 13.87
    files = glob("data_bnb_mcc83/*/Pandora*.root")
elif sample == "dataext":
    print("Data EXT")
    weight = 18.02
    files = glob("data_ext_mcc83/*/Pandora*.root")

chain = ROOT.TChain("robertoana/pandoratree")

for f in files:
    chain.Add(f)


entries = chain.GetEntries()
for i in range(entries):
    chain.GetEntry(i)
    gamma = 0
    e = 0
    pions = 0
    gamma_energy = 0
    pions_energy = 0

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

    neutrino_vertex = [chain.vx, chain.vy, chain.vz]
    for ish in range(chain.n_showers):
        pdg = chain.matched_showers[ish]
        dedx = chain.shower_dEdx[ish][2]
        most_proton_track = -1
        for itrk in range(chain.n_tracks):
            if chain.predict_p[itrk] > most_proton_track:
                most_proton_track = chain.predict_p[itrk]
                track_id = itrk

        shower_id = choose_shower(chain)

        track_vertex = [
            chain.track_start_x[track_id],
            chain.track_start_y[track_id],
            chain.track_start_z[track_id]]

        track_end = [
            chain.track_end_x[track_id],
            chain.track_end_y[track_id],
            chain.track_end_z[track_id]]

        shower_vertex = [
            chain.shower_start_x[shower_id],
            chain.shower_start_y[shower_id],
            chain.shower_start_z[shower_id]]

        neutrino_vertex = [chain.vx, chain.vy, chain.vz]

        shower_vertex_d = math.sqrt(
            sum([(s - n)**2
                 for s, n in zip(shower_vertex, neutrino_vertex)]))
        track_vertex_d = math.sqrt(
            sum([(t - n)**2
                 for t, n in zip(track_vertex, neutrino_vertex)]))
        track_end_d = math.sqrt(
            sum([(t - n)**2 for t, n in zip(track_end, neutrino_vertex)]))

        if sample == "mc":
            if abs(pdg) == 22:

                for ihit in range(len(chain.dEdx_hits[ish])):
                    hit_dedx = chain.dEdx_hits[ish][ihit]
                    h_dedx_hits_photon.Fill(hit_dedx, weight)
                h_dedx_photon.Fill(dedx, weight)

            elif abs(pdg) == 11:

                for ihit in range(len(chain.dEdx_hits[ish])):
                    hit_dedx = chain.dEdx_hits[ish][ihit]
                    h_dedx_hits_electron.Fill(hit_dedx, weight)

                h_dedx_electron.Fill(dedx, weight)

            else:
                for ihit in range(len(chain.dEdx_hits[ish])):
                    hit_dedx = chain.dEdx_hits[ish][ihit]
                    h_dedx_hits_other.Fill(hit_dedx, weight)
        else:
            for ihit in range(len(chain.dEdx_hits[ish])):
                hit_dedx = chain.dEdx_hits[ish][ihit]
                h_dedx_hits_data.Fill(hit_dedx, weight)

            h_dedx_data.Fill(dedx, weight)


if sample == "mc":
    c_dedx = ROOT.TCanvas("c_dedx")
    hs = ROOT.THStack("hs", ";dE/dx [MeV/cm];a.u.")
    h_dedx_photon.Scale(1 / h_dedx_photon.Integral())
    h_dedx_electron.Scale(1 / h_dedx_electron.Integral())
    l_dedx = ROOT.TLegend(0.37, 0.65, 0.84, 0.85)
    l_dedx.AddEntry(h_dedx_photon, "Photon", "f")
    l_dedx.AddEntry(h_dedx_electron, "Electron", "f")
    h_dedx_photon.SetLineColor(ROOT.kRed + 1)
    h_dedx_electron.SetLineColor(ROOT.kBlue + 1)
    hs.Add(h_dedx_electron)
    hs.Add(h_dedx_photon)
    hs.Draw("hist nostack")
    l_dedx.Draw()
    c_dedx.Update()

    c_dedx_hits = ROOT.TCanvas("c_dedx_hits")
    h_dedx_hits_electron.SetLineColor(ROOT.kBlack)
    h_dedx_hits_electron.Draw("ep")
    h_dedx_hits_electron.Sumw2()
    f_dedx_hits_electron = ROOT.TFile("f_dedx_hits_electron.root", "RECREATE")
    h_dedx_hits_electron.Write()
    f_dedx_hits_electron.Close()
    h_dedx_hits_electron.SetFillColor(ROOT.kAzure + 1)

    h_dedx_hits_photon.SetLineColor(ROOT.kBlack)
    h_dedx_hits_photon.Draw("ep same")
    h_dedx_hits_photon.Sumw2()
    f_dedx_hits_photon = ROOT.TFile("f_dedx_hits_photon.root", "RECREATE")
    h_dedx_hits_photon.Write()
    h_dedx_hits_photon.SetFillColor(ROOT.kRed - 4)
    f_dedx_hits_photon.Close()

    h_dedx_hits_other.SetLineColor(ROOT.kBlack)
    h_dedx_hits_other.Draw("ep same")
    h_dedx_hits_other.Sumw2()
    f_dedx_hits_other = ROOT.TFile("f_dedx_hits_other.root", "RECREATE")
    h_dedx_hits_other.Write()
    h_dedx_hits_other.SetFillColor(ROOT.kGray)

    f_dedx_hits_other.Close()

    h_dedx_stack = ROOT.THStack("h_dedx_stack", "")
    h_dedx_stack.Add(h_dedx_hits_electron)
    h_dedx_stack.Add(h_dedx_hits_photon)
    h_dedx_stack.Add(h_dedx_hits_other)
    f_langau = ROOT.TF1("f_langau", langau_lin, 0.5, 3.5, 5)

    f_langau.SetNpx(10000)
    f_langau.SetParNames("A", "#mu", "Landau #sigma", "Gaussian #sigma", "c")
    f_langau.SetParLimits(0, 0, 1)
    f_langau.SetParLimits(2, 0, 1)
    f_langau.SetParLimits(3, 0, 1)
    f_langau.SetParLimits(4, 0, 0.1)
    f_langau.SetParameters(1000, 1.7, 0.1, 0.2, 0.05)

    #h_dedx_hits_electron.Fit(f_langau, "R")

    l_langau = ROOT.TLegend(0.37, 0.65, 0.84, 0.85)
    l_langau.AddEntry(f_langau, "Fitted Landau", "l")
    l_langau.AddEntry(f_langau, "#mu = %.2f MeV/cm"
                      % f_langau.GetParameter(1), "")
    l_langau.AddEntry(f_langau, "Landau #sigma = %.2f MeV/cm"
                      % f_langau.GetParameter(2), "")
    l_langau.AddEntry(f_langau, "Gaussian #sigma = %.2f MeV/cm"
                      % f_langau.GetParameter(3), "")
    l_langau.Draw()
    h_dedx_stack.Draw("hist")
    c_dedx_hits.Update()

else:
    c_dedx_hits_data = ROOT.TCanvas("c_dedx_hits_data")
    if sample == "dataext":
        f_dedx_hits_data = ROOT.TFile("f_dedx_hits_dataext.root", "RECREATE")
        h_dedx_hits_data.Write()
    elif sample == "databnb":
        f_dedx_hits_data = ROOT.TFile("f_dedx_hits_databnb.root", "RECREATE")
        h_dedx_hits_data.Write()

    f_dedx_hits_data.Close()
    h_dedx_hits_data.Draw()
    c_dedx_hits_data.Update()

raw_input()
