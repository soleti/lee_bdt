#!/usr/local/bin/python3

import ROOT
from glob import glob
from bdt_common import is_fiducial, gauss_exp, printProgressBar
from array import array
from proton_energy import length2energy
ROOT.gStyle.SetOptTitle(1)
ROOT.gStyle.SetTitleFillColor(ROOT.kWhite)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetPalette(87)
ROOT.gStyle.SetNumberContours(99)
ROOT.gStyle.SetStatW(0.16)
ROOT.gStyle.SetOptFit(0)
PROTON_MASS = 0.938

nue_cosmic = glob("mc_nue_ubxsec/*.root")
c = ROOT.TChain("robertoana/pandoratree")

for f in nue_cosmic:
    c.Add(f)

entries = c.GetEntries()
n_bins = 10
h_reco_nu = []
h_true_nu = []

h_res_nu = []
h_res_nu_total = ROOT.TH1F("h_res_nu_total",
                                 ";(#nu_{E_{reco}} - #nu_{E})/#nu_{E}[GeV]; N. Entries / 0.0125",
                                 80, -0.5, 0.5)

h_nu_reco_true = ROOT.TH2F("h_nu_reco_true",
                               ";#nu_{E} [GeV];#nu_{E_{reco}}",
                               50, 0, 1,
                               50, 0, 1)

for i in range(n_bins):
    h_reco_nu.append(ROOT.TH1F("h_reco_nu%i" % i,
                                   "%.2f GeV < #nu_{E} < %.2f GeV;Reco. #nu_{E} [GeV];N. Entries / 0.02 GeV"
                                   % (0.1 * i , 0.1 * (i + 1)),
                                   60, 0, 1.2))
    h_true_nu.append(ROOT.TH1F("h_true_nu%i" % i,
                                   "%.2f GeV < #nu_{E} < %.2f GeV;Reco. #nu_{E_{reco}} [GeV];N. Entries / 0.02 GeV"
                                   % (0.1 * i, 0.1 * (i + 1)),
                                   60, 0, 1.2))
    h_res_nu.append(ROOT.TH1F("h_res_nu%i" % i,
                                  "%.2f GeV < #nu_{E} < %.2f GeV;(#nu_{E_{reco}} - #nu_{E})/#nu_{E} [GeV]; N. Entries / 0.05"
                                  % (0.1 * i, 0.1 * (i + 1)),
                                  50, -0.8, 0.8))

print("Entries", entries)

PROTON_THRESHOLD = 0.040
ELECTRON_THRESHOLD = 0.020

OFFSET = 0.8
BIAS = 0.017

for i in range(int(entries / 1)):
    printProgressBar(i, entries, prefix='Progress:', suffix='Complete', length=20)
    c.GetEntry(i)
    
    p = 0
    e = 0
    photons = 0
    pions = 0
    proton_energy = 0
    n_tracks = 0
    nu_energy = c.nu_E
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

        if c.nu_daughters_pdg[i_pdg] == 22:
            photons += 1

        if c.nu_daughters_pdg[i_pdg] == 111:
            pions += 1

    eNp = e == 1 and photons == 0 and pions == 0 and p == 1

    true_neutrino_vertex = [c.true_vx_sce,
                            c.true_vy_sce,
                            c.true_vz_sce]

    if not c.passed or not eNp or not is_fiducial(true_neutrino_vertex):
        continue

    reco_nu_energy = 0
    bin = int(nu_energy / (1 / n_bins))

    for i_sh in range(c.n_showers):
        reco_nu_energy += c.shower_energy[i_sh][2]

    for i_tr in range(c.n_tracks):
        reco_nu_energy += length2energy(c.track_len[i_tr])

    if bin < n_bins and reco_nu_energy:
        h_reco_nu[bin].Fill(reco_nu_energy)
        h_true_nu[bin].Fill(nu_energy)
        h_nu_reco_true.Fill(nu_energy, reco_nu_energy)
        h_res_nu[bin].Fill((reco_nu_energy / OFFSET + BIAS - nu_energy) / nu_energy)
        h_res_nu_total.Fill((reco_nu_energy / OFFSET + BIAS - nu_energy) / nu_energy)


if __name__ == "__main__":

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
        a_nu_reco_true.append(bin.GetMaximumBin() * 0.02 - 0.01)
        y_errs_low.append(bin.GetStdDev() / 2)
        y_errs_high.append(bin.GetStdDev() / 2)
        x_value = h_true_nu[i].FindBin(h_true_nu[i].GetMean()) * 0.02 - 0.01
        a_bins.append(x_value)
        x_err_h = (i + 1) * 0.1 - x_value
        x_err_l = x_value - i * 0.1
        x_errs_high.append(x_err_h)
        x_errs_low.append(x_err_l)

    g_nu_reco_true = ROOT.TGraphAsymmErrors(
        10, a_bins, a_nu_reco_true, x_errs_low, x_errs_high, y_errs_low, y_errs_high)

    h_nu_reco_true.Draw("colz")
    h_nu_reco_true.GetYaxis().SetTitleOffset(1.1)
    h_nu_reco_true.GetXaxis().SetTitleOffset(1.1)
        
    g_nu_reco_true.SetMarkerStyle(20)
    g_nu_reco_true.Draw("p same")

    f_line = ROOT.TF1("f_line", "[0]*x+[1]", 0, 1)
    f_line.SetParNames("m", "q")
    f_line.SetParameters(1, 0)
    l_p_true_reco = ROOT.TLegend(0.11, 0.913, 0.900, 0.968)
    l_p_true_reco.SetNColumns(2)
    l_p_true_reco.AddEntry(g_nu_reco_true, "Most probable values", "lep")
    g_nu_reco_true.Fit(f_line)
    l_p_true_reco.AddEntry(f_line, "#nu_{E_{reco}} = %.2f #nu_{E} + %.2f GeV" %
                        (f_line.GetParameter(0), f_line.GetParameter(1)), "l")
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

    for i in range(n_bins):
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
        
    for i in range(n_bins):
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
            h_res_nu[i].Fit(f_gausexp,
                                    "RQ",
                                    "",
                                    -0.35,
                                    0.16)
        elif i == 8:
            f_gausexp.SetParLimits(2, 6.69288e-02, 8.09288e-02)
            h_res_nu[i].Fit(f_gausexp,
                            "RQ",
                            "",
                            -1.6 * h_res_nu[i].GetRMS(),
                                  1.1 * h_res_nu[i].GetRMS())
        else:
            print(h_res_nu[i].GetRMS())
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
    g_e_res_e = ROOT.TGraphAsymmErrors(9, a_bins, a_res, x_errs_low, x_errs_high, a_res_err_l, a_res_err_h)
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
    h_res_nu_total.Fit(f_gausexp,
                           "",
                           "",
                           -0.2,
                           0.3)
    f_gausexp.Draw("same")
    c_res_total.SetBottomMargin(0.13)
    h_res_nu_total.GetXaxis().SetTitleOffset(1.1)
    c_res_total.SetLeftMargin(0.13)

    c_res_total.Update()
    c_res_total.SaveAs("plots/h_nu_res_total.pdf")
input()
