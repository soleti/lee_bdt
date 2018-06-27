#!/usr/local/bin/python3

import ROOT
from glob import glob
from bdt_common import is_fiducial, gauss_exp, printProgressBar, choose_shower, bins, pre_cuts, manual_cuts
from array import array
from root_numpy import hist2array
from proton_energy import length2energy
import math
import pickle
import numpy as np


ROOT.gStyle.SetOptTitle(1)
ROOT.gStyle.SetTitleFillColor(ROOT.kWhite)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetPalette(87)
ROOT.gStyle.SetNumberContours(99)
ROOT.gStyle.SetStatW(0.16)
ROOT.gStyle.SetOptFit(0)
PROTON_MASS = 0.938

nue_cosmic = glob("mc_nue_cali/*.root")
c = ROOT.TChain("robertoana/pandoratree")

for f in nue_cosmic:
    c.Add(f)

entries = c.GetEntries()
h_reco_nu = []
h_true_nu = []

h_res_nu = []
h_res_nu_total = ROOT.TH1F("h_res_nu_total",
                           ";(E_{dep} - E_{reco})/E_{dep}[GeV]; N. Entries / 0.0125",
                           80, -0.5, 0.5)

h_nu_reco_true = ROOT.TH2F("h_nu_reco_true",
                           ";E_{#nu} [GeV];E_{deposited} [GeV]",
                           50, 0, 3,
                           50, 0, 3)

h_reco_spectrum_corr = ROOT.TH1F("h_reco_spectrum_corr", "Reco. energy [GeV]", len(bins) - 1, bins)
h_true_spectrum_corr = ROOT.TH1F("h_true_spectrum_corr", "True energy [GeV]", len(bins) - 1, bins)

h_shower_proton = ROOT.TH1F("h_shower_proton",";Shower angle [#circ]; N. Entries / 4#circ",45,0,180)
h_shower_electron = ROOT.TH1F(
    "h_shower_electron", ";Shower angle [#circ]; N. Entries / 4#circ", 45, 0, 180)

h_shower_proton_dedx = ROOT.TH1F(
    "h_shower_proton_dedx", ";Shower dE/dx [MeV/cm]; N. Entries / 0.25 MeV/cm", 40, 0, 10)
h_shower_electron_dedx = ROOT.TH1F(
    "h_shower_electron_dedx", ";Shower dE/dx [MeV/cm]; N. Entries / 0.25 MeV/cm", 40, 0, 10)

h_shower = ROOT.TH1F("h_shower", ";Angle [#circ]; N. Entries / 4#circ", 45, 0, 180)
h_shower_dedx = ROOT.TH1F("h_shower_dedx", ";Shower dE/dx [MeV/cm]; N. Entries / 0.25 MeV/cm", 40, 0, 10)

h_track_proton = ROOT.TH1F(
    "h_track_proton", ";Track angle [#circ]; N. Entries / 4#circ", 45, 0, 180)
h_track_electron = ROOT.TH1F(
    "h_track_electron", ";Track angle [#circ]; N. Entries / 4#circ", 45, 0, 180)

h_track_proton_dedx = ROOT.TH1F(
    "h_track_proton_dedx", ";Track dE/dx [MeV/cm]; N. Entries / 0.25 MeV/cm", 40, 0, 6)
h_track_electron_dedx = ROOT.TH1F(
    "h_track_electron_dedx", ";Track dE/dx [MeV/cm]; N. Entries / 0.25 MeV/cm", 40, 0, 6)

h_track = ROOT.TH1F(
    "h_track", ";Track angle [#circ]; N. Entries / 4#circ", 45, 0, 180)
h_track_dedx = ROOT.TH1F(
    "h_track_dedx", ";Track dE/dx [MeV/cm]; N. Entries / 0.25 MeV/cm", 40, 0, 6)

h_angle = ROOT.TH1F("h_angle", "", 180, 0, 180)
h_angle_proton = ROOT.TH1F("h_angle_proton", "", 180, 0, 180)

for i in range(len(bins) - 1):
    h_reco_nu.append(ROOT.TH1F("h_reco_nu%i" % i,
                               "%.2f GeV < E_{dep} < %.2f GeV;E_{reco} [GeV];N. Entries / 0.05 GeV"
                               % (bins[i] , bins[i + 1]),
                               60, 0, 3))
    h_true_nu.append(ROOT.TH1F("h_true_nu%i" % i,
                               "%.2f GeV < E_{dep} < %.2f GeV;E_{reco} [GeV];N. Entries / 0.05 GeV"
                               % (bins[i], bins[i + 1]),
                               60, 0, 3))
    h_res_nu.append(ROOT.TH1F("h_res_nu%i" % i,
                              "%.2f GeV < E_{dep} < %.2f GeV;(E_{reco} - E_{dep})/E_{dep} [GeV]; N. Entries / 0.04"
                              % (bins[i], bins[i + 1]),
                              40, -0.8, 0.8))

print("Entries", entries)

PROTON_THRESHOLD = 0.040
ELECTRON_THRESHOLD = 0.020

f_fit = ROOT.TFile("plots/f_fit.root")
pol = f_fit.Get("f_pol2")
f_fit.Close()

p_a = pol.GetParameter(0)
p_b = pol.GetParameter(1)
p_c = pol.GetParameter(2)


for i in range(int(entries / 5)):
    printProgressBar(i, entries, prefix='Progress:', suffix='Complete', length=20)
    c.GetEntry(i)
    p = 0
    e = 0
    photons = 0
    pions = 0
    proton_energy = 0
    n_tracks = 0
    electron_energy = 0
    for i_pdg, energy in enumerate(c.nu_daughters_E):

        p_start = [
            c.nu_daughters_vx[i_pdg],
            c.nu_daughters_vy[i_pdg],
            c.nu_daughters_vz[i_pdg],
        ]

        p_end = [
            c.nu_daughters_endx[i_pdg],
            c.nu_daughters_endy[i_pdg],
            c.nu_daughters_endz[i_pdg],
        ]

        if abs(c.nu_daughters_pdg[i_pdg]) == 2212:
            if not (is_fiducial(p_start) and is_fiducial(p_end)):
                p = 0
                break

            if energy - PROTON_MASS > PROTON_THRESHOLD:
                p += 1
                proton_energy += energy - PROTON_MASS

        if abs(c.nu_daughters_pdg[i_pdg]) == 11:
            if not (is_fiducial(p_start)):
                e = 0
                break
            if energy > ELECTRON_THRESHOLD:
                e += 1
                electron_energy += energy

        if c.nu_daughters_pdg[i_pdg] == 22:
            photons += 1

        if c.nu_daughters_pdg[i_pdg] == 111 or abs(c.nu_daughters_pdg[i_pdg]) ==  211:
            pions += 1

    eNp = e == 1 and photons == 0 and pions == 0 and p >= 1
    nu_energy = c.nu_E #electron_energy + proton_energy

    true_neutrino_vertex = [c.true_vx_sce,
                            c.true_vy_sce,
                            c.true_vz_sce]

    if not c.passed or not eNp or not is_fiducial(true_neutrino_vertex) or not c.category == 2:
        continue


    hit_index = 2
    track_like_shower_id = -1
    shower_id = choose_shower(c, hit_index)

    no_tracks = False
    max_pca = 0
    min_pca = 1
    for ish, pca in enumerate(c.shower_pca):
        if pca > max_pca:
            max_pca = pca
            if ish != shower_id:
                track_like_shower_id = ish
        if pca < min_pca:
            min_pca = pca


    if c.n_tracks == 0 and c.n_showers > 1:
        no_tracks = True


    reco_nu_energy = 0
    corr_energy = 0
    bin = h_true_spectrum_corr.FindBin(nu_energy) - 1

    vec_shower = ROOT.TVector3(
        c.shower_dir_x[shower_id],
        c.shower_dir_y[shower_id],
        c.shower_dir_z[shower_id])

    for i_sh in range(c.n_showers):
        shower_vertex = [c.shower_start_x[i_sh],
                        c.shower_start_y[i_sh],
                        c.shower_start_z[i_sh]]

        v_sh = ROOT.TVector3(
            c.shower_dir_x[i_sh],
            c.shower_dir_y[i_sh],
            c.shower_dir_z[i_sh])

        cos = v_sh.Dot(vec_shower) / (v_sh.Mag() * vec_shower.Mag())
        shower_angle = math.degrees(math.acos(min(cos, 1)))

        h_shower.Fill(shower_angle)
        h_shower_dedx.Fill(c.shower_dEdx[i_sh][2])
        if c.matched_showers[i_sh] == 2212:
            h_shower_proton_dedx.Fill(c.shower_dEdx[i_sh][2])
            h_shower_proton.Fill(shower_angle)
        elif abs(c.matched_showers[i_sh]) == 11:
            h_shower_electron.Fill(shower_angle)
            h_shower_electron_dedx.Fill(c.shower_dEdx[i_sh][2])

        if (i_sh != shower_id and shower_angle > 15) or (no_tracks and i_sh == track_like_shower_id):
            length_e = length2energy(c.shower_length[i_sh])
            reco_nu_energy += length_e
            corr_energy += length_e
        else:
            reco_nu_energy += c.shower_energy[i_sh][hit_index]
            corr_energy += c.shower_energy[i_sh][hit_index] / 0.78 + 0.02

    for i_tr in range(c.n_tracks):

        v_tr = ROOT.TVector3(
            c.track_dir_x[i_tr],
            c.track_dir_y[i_tr],
            c.track_dir_z[i_tr])

        cos = v_tr.Dot(vec_shower) / (v_tr.Mag() * vec_shower.Mag())
        track_angle = math.degrees(math.acos(min(cos, 1)))
        h_track_dedx.Fill(c.track_dEdx[i_tr][2])

        if c.matched_tracks[i_tr] == 2212:
            h_track_proton.Fill(track_angle)
            h_track_proton_dedx.Fill(c.track_dEdx[i_tr][2])
        elif abs(c.matched_tracks[i_tr]) == 11:
            h_track_electron.Fill(track_angle)
            h_track_electron_dedx.Fill(c.track_dEdx[i_tr][2])

        h_track.Fill(track_angle)

        # if 20 < track_angle < 175 or c.track_dEdx[i_tr][2] > 2.3 or c.track_dEdx[i_tr][2] < 0.9:
        length_e = length2energy(c.track_len[i_tr])
        reco_nu_energy += length_e
        corr_energy += length_e 
        # else:
        #     reco_nu_energy += c.track_energy_hits[i_tr][hit_index] / 0.69
        #     corr_energy += c.track_energy_hits[i_tr][hit_index] / 0.69
    if 0 < reco_nu_energy < 3 and 0 < nu_energy < 3 and c.category == 2:
        p_c2 = p_c - corr_energy
        delta = p_b * p_b - 4 * p_a * p_c2
        corr_energy2 = 0

        if delta > 0:
            corr_energy2 = (- p_b + math.sqrt(delta)) / (2 * p_a)
        else:
            corr_energy2 = corr_energy

        corr_energy2 = corr_energy

        if h_reco_spectrum_corr.FindBin(corr_energy2) == h_true_spectrum_corr.FindBin(nu_energy):
            h_reco_spectrum_corr.Fill(reco_nu_energy)

        h_true_spectrum_corr.Fill(reco_nu_energy)
        h_reco_nu[bin].Fill(corr_energy2)
        h_true_nu[bin].Fill(nu_energy)
        h_nu_reco_true.Fill(nu_energy, corr_energy2)
        h_res_nu[bin].Fill((corr_energy2 - nu_energy) / nu_energy)
        h_res_nu_total.Fill((corr_energy2 - nu_energy) / nu_energy)


if __name__ == "__main__":
    c_angle = ROOT.TCanvas("c_angle")
    e_angle = ROOT.TEfficiency(h_angle_proton, h_angle)
    e_angle.Draw("ap")
    c_angle.Update()
    # ****************************************************
    # Reco vs true plot
    # ****************************************************

    c_nu_reco_true = ROOT.TCanvas("c_nu_reco_true")
    a_nu_reco_true = array("f", [])

    a_bins = array("f", [])
    
    x_errs_low = array("f", [])
    x_errs_high = array("f", [])
    y_errs_low = array("f", [])
    y_errs_high = array("f", [])

    for i, bin in enumerate(h_reco_nu):
        a_nu_reco_true.append(bin.GetXaxis().GetBinCenter(bin.GetMaximumBin()))
        y_errs_low.append(bin.GetStdDev() / 2)
        y_errs_high.append(bin.GetStdDev() / 2)
        x_value = h_true_nu[i].GetMean()#Xaxis().GetBinCenter(h_true_nu[i].GetMaximumBin())
        a_bins.append(x_value)
        bin_index = h_true_spectrum_corr.FindBin(x_value)
        x_err_h = h_true_spectrum_corr.GetXaxis().GetBinUpEdge(bin_index) - x_value
        x_err_l = x_value - h_true_spectrum_corr.GetXaxis().GetBinLowEdge(bin_index)
        x_errs_high.append(x_err_h)
        x_errs_low.append(x_err_l)

    g_nu_reco_true = ROOT.TGraphAsymmErrors(
        len(bins) - 1, a_bins, a_nu_reco_true, x_errs_low, x_errs_high, y_errs_low, y_errs_high)

    h_nu_reco_true.SetMinimum(-0.001)
    h_nu_reco_true.Draw("col")
    h_nu_reco_true.GetYaxis().SetTitleOffset(1.1)
    h_nu_reco_true.GetXaxis().SetTitleOffset(1.1)
        
    g_nu_reco_true.SetMarkerStyle(20)
    g_nu_reco_true.Draw("p same")

    f_line = ROOT.TF1("f_pol2", "pol1", 0, 1)
    f_line.SetParNames("a", "b")
    f_line.SetParameters(0, 0)
    f_line.SetParLimits(0, -1, 0)
    f_line.SetParameters(1, 0.84)
    l_p_true_reco = ROOT.TLegend(0.11, 0.913, 0.900, 0.968)
    l_p_true_reco.SetNColumns(2)
    l_p_true_reco.AddEntry(g_nu_reco_true, "Most probable values", "lep")
    g_nu_reco_true.Fit(f_line)
    # f_fit = ROOT.TFile("plots/f_fit.root", "RECREATE")
    # f_line.Write()
    # f_fit.Close()
    l_p_true_reco.AddEntry(f_line, "E_{deposited} = %.2f E_{#nu} + %.2f" %
                        (f_line.GetParameter(1), f_line.GetParameter(0)), "l")
    l_p_true_reco.Draw()
    c_nu_reco_true.SetLeftMargin(0.12)
    c_nu_reco_true.SetBottomMargin(0.13)

    c_nu_reco_true.Update()
    c_nu_reco_true.SaveAs("plots/h_nu_slope.pdf")

    # ****************************************************
    # Reco, true per interval 
    # ****************************************************
    ROOT.gStyle.SetTitleSize(0.6)
    c_nu_energy = ROOT.TCanvas("c_nu_energy", "", 1000, 700)
    c_nu_energy.Divide(5, 2)

    for i in range(len(bins) - 1):
        c_nu_energy.cd(i + 1)
        h_true_nu[i].Draw()
        h_reco_nu[i].Draw("same")
        if i == 0:
            l_true_reco = ROOT.TLegend(0.54, 0.76, 0.85, 0.83)
            l_true_reco.AddEntry(h_reco_nu[i], "Reco. energy", "l")
            l_true_reco.AddEntry(h_true_nu[i], "True energy", "l")
            l_true_reco.Draw()
        h_true_nu[i].SetLineColor(ROOT.kBlack)
        h_reco_nu[i].GetXaxis().SetTitleSize(0.05)
        h_reco_nu[i].GetXaxis().SetTitleOffset(1.05)
        h_reco_nu[i].GetYaxis().SetTitleOffset(1.3)
        ROOT.gPad.SetLeftMargin(0.15)
        ROOT.gPad.SetBottomMargin(0.18)

    c_nu_energy.Update()
    c_nu_energy.SaveAs("plots/h_nu_energy.pdf")

    # ROOT.gStyle.SetOptFit(1)


    # ****************************************************
    # nu energy resolution per interval 
    # ****************************************************
    c_nu_res = ROOT.TCanvas("c_nu_res", "", 1200, 700)
    c_nu_res.Divide(5, 2)
    a_res = array("f", [])
    a_res_err_h = array("f", [])
    a_res_err_l = array("f", [])
        
    for i in range(len(bins) - 1):
        c_nu_res.cd(i + 1)

        h_res_nu[i].Draw()
        h_res_nu[i].GetXaxis().SetTitleSize(0.05)
        h_res_nu[i].GetXaxis().SetTitleOffset(1.05)
        h_res_nu[i].GetYaxis().SetTitleOffset(1.3)

        f_gausexp = ROOT.TF1("f_gausexp", gauss_exp, -1, 1, 4)

        f_gausexp.SetParLimits(1, -0.1, 0.1)
        f_gausexp.SetParLimits(2, 0.001, 0.6)
        f_gausexp.SetParLimits(3, 0.001, 1)

        f_gausexp.SetParameters(h_res_nu[i].GetMaximum(),
                                0,
                                h_res_nu[i].GetRMS(),
                                0.4)
        f_gausexp.SetParNames("A", "#mu", "#sigma", "k")

        if i == 9:
            f_gausexp.SetParLimits(2, 6.49288e-02, 8.09288e-02)
            # h_res_nu[i].Fit(f_gausexp,
            #                         "RQ",
            #                         "",
            #                         -0.35,
            #                         0.16)
        elif i == 8:
            f_gausexp.SetParLimits(2, 6.69288e-02, 8.09288e-02)
            # h_res_nu[i].Fit(f_gausexp,
            #                 "RQ",
            #                 "",
            #                 -1.6 * h_res_nu[i].GetRMS(),
            #                       1.1 * h_res_nu[i].GetRMS())
        else:
            print("fit")
            # h_res_nu[i].Fit(f_gausexp,
            #                 "RQ",
            #                 "",
            #                 -1.6 * max(0.1, h_res_nu[i].GetRMS()),
            #                 1.1)

        a_res.append(f_gausexp.GetParameter(2) * 100)
        a_res_err_h.append(f_gausexp.GetParError(2) * 100 / 2)
        a_res_err_l.append(f_gausexp.GetParError(2) * 100 / 2)

        ROOT.gPad.SetLeftMargin(0.13)
        ROOT.gPad.SetBottomMargin(0.18)
    c_nu_res.Update()
    c_nu_res.SaveAs("plots/h_nu_res.pdf")

    # ****************************************************
    # nu energy resolution vs true energy
    # ****************************************************
    ROOT.gStyle.SetOptTitle(0)

    c_e_res_e = ROOT.TCanvas("c_nu_res_e")
    a_bins.pop(0)
    a_res.pop(0)
    x_errs_low.pop(0)
    x_errs_high.pop(0)
    a_res_err_l.pop(0)
    a_res_err_h.pop(0)
    g_e_res_e = ROOT.TGraphAsymmErrors(len(bins) - 2, a_bins, a_res, x_errs_low, x_errs_high, a_res_err_l, a_res_err_h)
    g_e_res_e.Draw("ap")
    g_e_res_e.GetYaxis().SetTitle("#sigma [%]")
    g_e_res_e.GetXaxis().SetTitle("#nu_{e} [GeV]")
    g_e_res_e.GetXaxis().SetTitleSize(0.05)
    g_e_res_e.GetXaxis().SetTitleOffset(1)
    g_e_res_e.GetYaxis().SetTitleOffset(0.9)
    c_e_res_e.SetBottomMargin(0.13)
    g_e_res_e.SetMarkerStyle(20)
    f_res = ROOT.TF1("f_res", "sqrt(([0]/sqrt(x))**2+([1]/x)**2+[2]**2)", 0, 2)
    f_res.SetParNames("a", "b", "c")
    # f_res.SetParLimits(0, 0.1, 1)
    f_res.SetParLimits(1, 0.1, 0.98)
    f_res.SetParLimits(2, 0.1, 1.94)

    g_e_res_e.Fit(f_res)
    l_res = ROOT.TLegend(0.37, 0.68, 0.80, 0.86)
    l_res.AddEntry(f_res, "(#frac{%.2f}{#sqrt{E / GeV}} #oplus #frac{%.2f}{E / GeV} #oplus %.2f) %%" %
                   (f_res.GetParameter(0),
                    f_res.GetParameter(1),
                    f_res.GetParameter(2)), "l")
    l_res.Draw()
    c_e_res_e.Update()
    f_e_res_e = ROOT.TFile("plots/g_nu_res_e.root", "RECREATE")
    g_e_res_e.Write()
    f_e_res_e.Close()
    c_e_res_e.SaveAs("plots/h_nu_res_e.pdf")

    c_res_total = ROOT.TCanvas("c_res_total")
    h_res_nu_total.Draw("hist")
    f_gausexp.SetParLimits(1, -0.2, 0.2)
    f_gausexp.SetParLimits(2, 0.01, 0.15)
    f_gausexp.SetParLimits(3, 0, 1)
    # h_res_nu_total.Fit(f_gausexp,
    #                        "",
    #                        "",
    #                        -0.2,
    #                        0.3)
    f_gausexp.Draw("same")
    c_res_total.SetBottomMargin(0.13)
    h_res_nu_total.GetXaxis().SetTitleOffset(1.1)
    c_res_total.SetLeftMargin(0.13)

    c_res_total.Update()
    c_res_total.SaveAs("plots/h_nu_res_total.pdf")


    c_eff = ROOT.TCanvas("c_eff")
    e_energy = ROOT.TEfficiency(h_reco_spectrum_corr, h_true_spectrum_corr)
    e_energy.SetFillStyle(3002)
    e_energy.SetFillColor(1)
    e_energy.Draw("a2")
    e_energy.Draw("p same")
    c_eff.Update()

    e_energy.GetPaintedGraph().GetXaxis().SetTitle("Reco. energy [GeV]")
    e_energy.GetPaintedGraph().GetYaxis().SetTitle("Bin purity")
    e_energy.GetPaintedGraph().GetYaxis().SetRangeUser(0.001, 1)
    c_eff.Update()

    c_shower = ROOT.TCanvas("c_shower")
    h_shower_proton.Divide(h_shower)
    h_shower_electron.Divide(h_shower)

    h_shower_proton.Draw("hist")
    h_shower_proton.GetYaxis().SetRangeUser(0.001, 1)

    h_shower_electron.SetLineColor(ROOT.kRed + 1)
    h_shower_electron.Draw("hist same")

    c_shower.Update()

    c_track = ROOT.TCanvas("c_track")
    h_track_proton.Divide(h_track)
    h_track_electron.Divide(h_track)

    h_track_proton.Draw("hist")
    h_track_proton.GetYaxis().SetRangeUser(0.001,1)

    h_track_electron.SetLineColor(ROOT.kRed + 1)
    h_track_electron.Draw("hist same")

    c_track.Update()


    c_track_dedx = ROOT.TCanvas("c_track_dedx")
    h_track_proton_dedx.Divide(h_track_dedx)
    h_track_electron_dedx.Divide(h_track_dedx)

    h_track_proton_dedx.Draw("hist")
    h_track_proton_dedx.GetYaxis().SetRangeUser(0.001, 1)

    h_track_electron_dedx.SetLineColor(ROOT.kRed + 1)
    h_track_electron_dedx.Draw("hist same")

    c_track_dedx.Update()

    c_shower_dedx = ROOT.TCanvas("c_shower_dedx")
    h_shower_proton_dedx.Divide(h_shower_dedx)
    h_shower_electron_dedx.Divide(h_shower_dedx)

    h_shower_proton_dedx.Draw("hist")
    h_shower_proton_dedx.GetYaxis().SetRangeUser(0.001, 1)
        
    h_shower_electron_dedx.SetLineColor(ROOT.kRed + 1)
    h_shower_electron_dedx.Draw("hist same")

    c_shower_dedx.Update()


input()
