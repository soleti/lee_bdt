#!/usr/local/bin/python3

import ROOT
from glob import glob
from bdt_common import is_fiducial, gauss_exp, printProgressBar
from array import array

ROOT.gStyle.SetOptTitle(1)
ROOT.gStyle.SetTitleFillColor(ROOT.kWhite)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetPalette(87)
ROOT.gStyle.SetNumberContours(99)
ROOT.gStyle.SetStatW(0.16)
ROOT.gStyle.SetOptFit(0)
ELECTRON_MASS = 0.00052

nue_cosmic = glob("data_files/mc_nue_pid/*.root")
c = ROOT.TChain("robertoana/pandoratree")

for f in nue_cosmic:
    c.Add(f)

entries = c.GetEntries()
n_bins = 10
h_reco_electron = []
h_true_electron = []

h_res_electron = []
h_res_electron_total = ROOT.TH1F("h_res_electron_total",
                                 ";(e_{E_{reco}} - e_{E})/e_{E}[GeV]; N. Entries / 0.0125",
                                 80, -0.5, 0.5)

h_electron_reco_true = ROOT.TH2F("h_electron_reco_true",
                                 ";E^{e} [GeV];E^{e}_{reco} [GeV]",
                               50, 0, 1,
                               50, 0, 1)

for i in range(n_bins):
    h_reco_electron.append(ROOT.TH1F("h_reco_electron%i" % i,
                                   "%.2f GeV < e_{E} < %.2f GeV;Reco. e_{E} [GeV];N. Entries / 0.02 GeV"
                                   % (0.1 * i , 0.1 * (i + 1)),
                                   60, 0, 1.2))
    h_true_electron.append(ROOT.TH1F("h_true_electron%i" % i,
                                   "%.2f GeV < e_{E} < %.2f GeV;Reco. e_{E_{reco}} [GeV];N. Entries / 0.02 GeV"
                                   % (0.1 * i, 0.1 * (i + 1)),
                                   60, 0, 1.2))
    h_res_electron.append(ROOT.TH1F("h_res_electron%i" % i,
                                  "%.2f GeV < e_{E} < %.2f GeV;(e_{E_{reco}} - e_{E})/e_{E} [GeV]; N. Entries / 0.05"
                                  % (0.1 * i, 0.1 * (i + 1)),
                                  50, -0.8, 0.8))

print("Entries", entries)

ELECTRON_THRESHOLD = 0.020
OFFSET = 0.8
BIAS = 0.017

for i in range(int(entries / 1)):
    printProgressBar(i, entries, prefix='Progress:', suffix='Complete', length=20)
    c.GetEntry(i)
    
    e = 0
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

        if abs(c.nu_daughters_pdg[i_pdg]) == 11:
            if not (is_fiducial(p_start) and is_fiducial(p_end) and p_start[2] < 700 and p_start[2] > 300 and p_end[2] < 700 and p_end[2] > 300):
                e = 0
                break
            if energy > ELECTRON_THRESHOLD:
                electron_energy = energy
                e += 1

    true_neutrino_vertex = [c.true_vx_sce,
                            c.true_vy_sce,
                            c.true_vz_sce]

    if not c.passed or not e or not is_fiducial(true_neutrino_vertex):
        continue

    reco_electron_energy = 0
    bin = int(electron_energy / (1 / n_bins))

    # if c.n_showers > 1:
    #     continue

    for i_sh in range(c.n_showers):
        if c.matched_showers[i_sh] == 11:
            reco_electron_energy += c.shower_energy[i_sh][2] * c.shower_energy_cali[i_sh][2]

    if bin < n_bins and reco_electron_energy:
        h_reco_electron[bin].Fill(reco_electron_energy)
        h_true_electron[bin].Fill(electron_energy)
        h_electron_reco_true.Fill(electron_energy, reco_electron_energy)
        h_res_electron[bin].Fill((reco_electron_energy / OFFSET + BIAS - electron_energy) / electron_energy)
        h_res_electron_total.Fill((reco_electron_energy / OFFSET + BIAS - electron_energy) / electron_energy)


if __name__ == "__main__":

    # ****************************************************
    # Reco vs true plot
    # ****************************************************

    c_electron_reco_true = ROOT.TCanvas("c_electron_reco_true")
    a_electron_reco_true = array("f", [])

    a_bins = array("f", [])
    
    x_errs_low = array("f", [])
    x_errs_high = array("f", [])
    y_errs_low = array("f", [])
    y_errs_high = array("f", [])

    for i, bin in enumerate(h_reco_electron):
        # bin.GetMaximumBin() * 0.02 - 0.01
        a_electron_reco_true.append(bin.GetMaximumBin() * 0.02 - 0.01)
        y_errs_low.append(bin.GetStdDev() / 2)
        y_errs_high.append(bin.GetStdDev() / 2)
        x_value = h_true_electron[i].FindBin(h_true_electron[i].GetMean()) * 0.02 - 0.01
        a_bins.append(x_value)
        x_err_h = (i + 1) * 0.1 - x_value
        x_err_l = x_value - i * 0.1
        x_errs_high.append(x_err_h)
        x_errs_low.append(x_err_l)

    g_electron_reco_true = ROOT.TGraphAsymmErrors(
        10, a_bins, a_electron_reco_true, x_errs_low, x_errs_high, y_errs_low, y_errs_high)

    h_electron_reco_true.Draw("col")
    h_electron_reco_true.SetMinimum(-0.001)
    h_electron_reco_true.GetYaxis().SetTitleOffset(1.1)
    h_electron_reco_true.GetXaxis().SetTitleOffset(1.1)
        
    g_electron_reco_true.SetMarkerStyle(20)
    g_electron_reco_true.Draw("p same")

    f_line = ROOT.TF1("f_line", "[0]*x+[1]", 0, 1)
    f_line.SetParNames("m", "q")
    f_line.SetParameters(0.784, -0.017)
    l_p_true_reco = ROOT.TLegend(0.11, 0.913, 0.900, 0.968)
    l_p_true_reco.SetNColumns(2)
    l_p_true_reco.AddEntry(g_electron_reco_true, "Most probable values", "lep")
    g_electron_reco_true.Fit(f_line)
    l_p_true_reco.AddEntry(f_line, "E_{reco}^{e} = %.2f E^{e} - %.2f GeV" % (f_line.GetParameter(0), -f_line.GetParameter(1)), "l")
    f_line.SetLineColor(ROOT.kRed + 1)
    f_line.Draw("same")
    l_p_true_reco.Draw()
    c_electron_reco_true.SetLeftMargin(0.12)
    c_electron_reco_true.SetBottomMargin(0.13)

    c_electron_reco_true.Update()
    c_electron_reco_true.SaveAs("plots/h_electron_slope.pdf")

    # ****************************************************
    # Reco, true per interval 
    # ****************************************************
    ROOT.gStyle.SetTitleSize(0.6)
    c_electron_energy = ROOT.TCanvas("c_electron_energy", "", 1000, 700)
    c_electron_energy.Divide(5, 2)

    for i in range(n_bins):
        c_electron_energy.cd(i + 1)
        h_true_electron[i].Draw()
        h_reco_electron[i].Draw("same")
        if i == 0:
            l_true_reco = ROOT.TLegend(0.54, 0.76, 0.85, 0.83)
            l_true_reco.AddEntry(h_reco_electron[i], "Reco. energy", "l")
            l_true_reco.AddEntry(h_true_electron[i], "True energy", "l")
            l_true_reco.Draw()
        h_true_electron[i].SetLineColor(ROOT.kBlack)
        h_reco_electron[i].GetXaxis().SetTitleSize(0.05)
        h_reco_electron[i].GetXaxis().SetTitleOffset(1.05)
        h_reco_electron[i].GetYaxis().SetTitleOffset(1.3)
        ROOT.gPad.SetLeftMargin(0.15)
        ROOT.gPad.SetBottomMargin(0.18)

    c_electron_energy.Update()
    c_electron_energy.SaveAs("plots/h_electron_energy.pdf")

    # ROOT.gStyle.SetOptFit(1)


    # ****************************************************
    # electron energy resolution per interval 
    # ****************************************************
    c_electron_res = ROOT.TCanvas("c_electron_res", "", 1200, 700)
    c_electron_res.Divide(5, 2)
    a_res = array("f", [])
    a_res_err_h = array("f", [])
    a_res_err_l = array("f", [])
        
    for i in range(n_bins):
        c_electron_res.cd(i + 1)

        h_res_electron[i].Draw()
        h_res_electron[i].GetXaxis().SetTitleSize(0.05)
        h_res_electron[i].GetXaxis().SetTitleOffset(1.05)
        h_res_electron[i].GetYaxis().SetTitleOffset(1.3)

        f_gausexp = ROOT.TF1("f_gausexp", gauss_exp, -1, 1, 4)

        f_gausexp.SetParLimits(1, -0.1, 0.1)
        f_gausexp.SetParLimits(2, 0.001, 0.6)
        f_gausexp.SetParLimits(3, 0.001, 1)

        f_gausexp.SetParameters(h_res_electron[i].GetMaximum(),
                                0,
                                h_res_electron[i].GetRMS(),
                                0.4)
        f_gausexp.SetParNames("A", "#mu", "#sigma", "k")

        if i == 9:
            f_gausexp.SetParLimits(2, 6.49288e-02, 8.09288e-02)
            h_res_electron[i].Fit(f_gausexp,
                                    "RQ",
                                    "",
                                    -0.35,
                                    0.16)
        elif i == 8:
            f_gausexp.SetParLimits(2, 6.69288e-02, 8.09288e-02)
            h_res_electron[i].Fit(f_gausexp,
                                  "RQ",
                                  "",
                                  -1.6 * h_res_electron[i].GetRMS(),
                                  1.1 * h_res_electron[i].GetRMS())
        else:
            h_res_electron[i].Fit(f_gausexp,
                                  "RQ",
                                  "",
                                  -1.6 * h_res_electron[i].GetRMS(),
                                  1.1 * h_res_electron[i].GetRMS())

        a_res.append(f_gausexp.GetParameter(2) * 100)
        a_res_err_h.append(f_gausexp.GetParError(2) * 100 / 2)
        a_res_err_l.append(f_gausexp.GetParError(2) * 100 / 2)

        ROOT.gPad.SetLeftMargin(0.13)
        ROOT.gPad.SetBottomMargin(0.18)
    c_electron_res.Update()
    c_electron_res.SaveAs("plots/h_electron_res.pdf")

    # ****************************************************
    # electron energy resolution vs true energy
    # ****************************************************
    ROOT.gStyle.SetOptTitle(0)

    c_e_res_e = ROOT.TCanvas("c_e_res_e")
    a_bins.pop(0)
    a_res.pop(0)
    x_errs_low.pop(0)
    x_errs_high.pop(0)
    a_res_err_l.pop(0)
    a_res_err_h.pop(0)
    g_e_res_e = ROOT.TGraphAsymmErrors(9, a_bins, a_res, x_errs_low, x_errs_high, a_res_err_l, a_res_err_h)
    g_e_res_e.Draw("ap")
    g_e_res_e.GetYaxis().SetTitle("#sigma [%]")
    g_e_res_e.GetXaxis().SetTitle("E_{e} [GeV]")
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
    f_e_res_e = ROOT.TFile("plots/g_e_res_e.root", "RECREATE")
    g_e_res_e.Write()
    f_e_res_e.Close()
    c_e_res_e.SaveAs("plots/h_electron_res_e.pdf")

    c_res_total = ROOT.TCanvas("c_res_total")
    h_res_electron_total.Draw("hist")
    f_gausexp.SetParLimits(1, -0.2, 0.2)
    f_gausexp.SetParLimits(2, 0.01, 0.15)
    f_gausexp.SetParLimits(3, 0, 1)
    h_res_electron_total.Fit(f_gausexp,
                           "",
                           "",
                           -0.2,
                           0.3)
    f_gausexp.Draw("same")
    c_res_total.SetBottomMargin(0.13)
    h_res_electron_total.GetXaxis().SetTitleOffset(1.1)
    c_res_total.SetLeftMargin(0.13)

    c_res_total.Update()
    c_res_total.SaveAs("plots/h_electron_res_total.pdf")
input()
