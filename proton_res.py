#!/usr/local/bin/python3

import ROOT
from glob import glob
from bdt_common import is_fiducial, gauss_exp, printProgressBar
from proton_energy import length2energy
from array import array

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
n_bins = 6
h_reco_proton = []
h_true_proton = []

h_res_proton = []
h_res_proton_total = ROOT.TH1F("h_res_proton_total",
                               ";(p_{E_{reco}} - p_{E_{k}})/p_{E_{k}}[GeV]; N. Entries / 0.0125",
                               80, -0.5, 0.5)

h_proton_reco_true = ROOT.TH2F("h_proton_reco_true",
                               ";p_{E_{k}} [GeV];p_{E_{reco}}",
                               50, 0, 0.6,
                               50, 0, 0.6)

for i in range(n_bins):
    h_reco_proton.append(ROOT.TH1F("h_reco_proton%i" % i,
                                   "%.2f GeV < p_{E_{k}} < %.2f GeV;Reco. p_{E_{k}} [GeV];N. Entries / 0.02 GeV"
                                   % (0.1 * i , 0.1 * (i + 1)),
                                   40, 0, 0.8))
    h_true_proton.append(ROOT.TH1F("h_true_proton%i" % i,
                                   "%.2f GeV < p_{E_{k}} < %.2f GeV;Reco. p_{E_{k}} [GeV];N. Entries / 0.02 GeV"
                                   % (0.1 * i, 0.1 * (i + 1)),
                                   40, 0, 0.8))
    h_res_proton.append(ROOT.TH1F("h_res_proton%i" % i,
                                  "%.2f GeV < p_{E_{k}} < %.2f GeV;(p_{E_{reco}} - p_{E_{k}})/p_{E_{k}}[GeV]; N. Entries / 0.0125"
                                  % (0.1 * i, 0.1 * (i + 1)),
                                  80, -0.5, 0.5))

print("Entries", entries)

PROTON_THRESHOLD = 0.040
ELECTRON_THRESHOLD = 0.020
OFFSET = 0.96
single_tracks_events = 0
more_tracks = 0
for i in range(int(entries / 10)):
    printProgressBar(i, entries, prefix='Progress:', suffix='Complete', length=20)
    c.GetEntry(i)
    
    p = 0
    e = 0
    photons = 0
    pions = 0
    proton_energy = 0
    n_tracks = 0
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
            if not (is_fiducial(p_start) and is_fiducial(p_end) and p_start[2] < 700 and p_start[2] > 300 and p_end[2] < 700 and p_end[2] > 300):
                p = 0 
                break

            if energy - PROTON_MASS > PROTON_THRESHOLD:
                p += 1
                proton_energy += energy - PROTON_MASS

        if abs(c.nu_daughters_pdg[i_pdg]) == 11:
            if energy > ELECTRON_THRESHOLD:
                e += 1

        if c.nu_daughters_pdg[i_pdg] == 22:
            photons += 1

        if c.nu_daughters_pdg[i_pdg] == 111:
            pions += 1

    eNp = e == 1 and photons == 0 and pions == 0 and p == 1

    if not c.passed or not eNp:
        continue

    true_neutrino_vertex = [c.true_vx_sce,
                            c.true_vy_sce,
                            c.true_vz_sce]

    if not is_fiducial(true_neutrino_vertex):
        continue

    reco_proton_energy = 0
    bin = int(proton_energy / (0.6 / n_bins))
    
    for i_tr in range(c.n_tracks):
        if c.matched_tracks[i_tr] == 2212:
            reco_proton_energy += length2energy(c.track_len[i_tr])
            n_tracks += 1
    if bin < n_bins and reco_proton_energy:
        h_reco_proton[bin].Fill(reco_proton_energy)
        h_true_proton[bin].Fill(proton_energy)
        h_proton_reco_true.Fill(proton_energy, reco_proton_energy)
        h_res_proton[bin].Fill((reco_proton_energy / OFFSET - proton_energy) / proton_energy)
        h_res_proton_total.Fill(
            (reco_proton_energy / OFFSET - proton_energy) / proton_energy)

    if n_tracks == 1:
        single_tracks_events +=1
    if n_tracks > 1:
        more_tracks +=1

print(single_tracks_events, more_tracks)

if __name__ == "__main__":

    # ****************************************************
    # Reco vs true plot
    # ****************************************************

    c_proton_reco_true = ROOT.TCanvas("c_proton_reco_true")
    a_proton_reco_true = array("f", [])

    a_bins = array("f", [])
    
    x_errs_low = array("f", [])
    x_errs_high = array("f", [])
    y_errs_low = array("f", [])
    y_errs_high = array("f", [])

    for i, bin in enumerate(h_reco_proton):
        a_proton_reco_true.append(bin.GetMaximumBin() * 0.02 - 0.01)
        y_errs_low.append(bin.GetStdDev() / 2)
        y_errs_high.append(bin.GetStdDev() / 2)

        # x_value = h_true_proton[i].FindBin(h_true_proton[i].GetMean()) * 0.02 - 0.01
        x_value = h_true_proton[i].GetMaximumBin() * 0.02 - 0.01

        a_bins.append(x_value)
        x_err_h = (i + 1) * 0.1 - x_value
        x_err_l = x_value - i * 0.1
        x_errs_high.append(x_err_h)
        x_errs_low.append(x_err_l)

    g_proton_reco_true = ROOT.TGraphAsymmErrors(
        6, a_bins, a_proton_reco_true, x_errs_low, x_errs_high, y_errs_low, y_errs_high)

    h_proton_reco_true.Draw("colz")
    h_proton_reco_true.SetMinimum(-0.001)
    h_proton_reco_true.GetYaxis().SetTitleOffset(1.1)
    g_proton_reco_true.SetMarkerStyle(20)
    g_proton_reco_true.Draw("p same")

    f_line = ROOT.TF1("f_line", "[0]*x+[1]", 0, 1)
    f_line.SetParNames("m", "q")
    f_line.SetParameters(1, 0)
    l_p_true_reco = ROOT.TLegend(0.11, 0.913, 0.900, 0.968)
    l_p_true_reco.SetNColumns(2)
    l_p_true_reco.AddEntry(g_proton_reco_true, "Most probable value", "lep")
    g_proton_reco_true.Fit(f_line, "R", "0", 0, 0.5)
    l_p_true_reco.AddEntry(f_line, "E_{k}^{reco} = %.2f E_{k}^{true} + %.2f GeV" %
                        (0.99, 0), "l")
    l_p_true_reco.Draw()
    f_line.Draw("same")
    c_proton_reco_true.SetLeftMargin(0.12)
    c_proton_reco_true.SetBottomMargin(0.13)

    c_proton_reco_true.Update()
    # c_proton_reco_true.SaveAs("plots/h_proton_slope.pdf")

    # ****************************************************
    # Reco, true per interval 
    # ****************************************************
    ROOT.gStyle.SetTitleSize(0.6)
    c_proton_energy = ROOT.TCanvas("c_proton_energy", "", 1000, 700)
    c_proton_energy.Divide(3, 2)

    for i in range(n_bins):
        c_proton_energy.cd(i + 1)
        h_true_proton[i].Draw()
        h_reco_proton[i].Draw("same")
        if i == 0:
            l_true_reco = ROOT.TLegend(0.54, 0.76, 0.85, 0.83)
            l_true_reco.AddEntry(h_reco_proton[i], "Reco. energy", "l")
            l_true_reco.AddEntry(h_true_proton[i], "True energy", "l")
            l_true_reco.Draw()
        h_true_proton[i].SetLineColor(ROOT.kBlack)
        h_reco_proton[i].GetXaxis().SetTitleSize(0.05)
        h_reco_proton[i].GetXaxis().SetTitleOffset(1.05)
        h_reco_proton[i].GetYaxis().SetTitleOffset(1.3)
        ROOT.gPad.SetLeftMargin(0.13)
        ROOT.gPad.SetBottomMargin(0.18)

    c_proton_energy.Update()
    # c_proton_energy.SaveAs("plots/h_proton_energy.pdf")

    # ROOT.gStyle.SetOptFit(1)


    # ****************************************************
    # Proton energy resolution per interval 
    # ****************************************************
    c_proton_res = ROOT.TCanvas("c_proton_res", "", 1200, 700)
    c_proton_res.Divide(3, 2)
    a_res = array("f", [])
    a_res_err_h = array("f", [])
    a_res_err_l = array("f", [])
        
    for i in range(n_bins):
        c_proton_res.cd(i + 1)

        h_res_proton[i].Draw()
        h_res_proton[i].GetXaxis().SetTitleSize(0.05)
        h_res_proton[i].GetXaxis().SetTitleOffset(1.05)
        h_res_proton[i].GetYaxis().SetTitleOffset(1.3)

        f_gausexp = ROOT.TF1("f_gausexp", gauss_exp, -1, 1, 4)

        f_gausexp.SetParLimits(1, -0.1, 0.1)
        f_gausexp.SetParLimits(2, 0.001, 0.6)
        f_gausexp.SetParLimits(3, 0.001, 1)

        f_gausexp.SetParameters(h_res_proton[i].GetMaximum(),
                                0,
                                h_res_proton[i].GetRMS(),
                                0.4)
        f_gausexp.SetParNames("A", "#mu", "#sigma", "k")

        if i >= 3:
            f_gausexp.SetParameters(h_res_proton[i].GetMaximum(),
                                    0,
                                    0.018,
                                    0.39)

            h_res_proton[i].Fit(f_gausexp,
                                "RQ",
                                "",
                                -0.05,
                                0.08)
        elif i == 0:
            h_res_proton[i].Fit(f_gausexp,
                                "RQ",
                                "",
                                -0.3,
                                0.2)
        elif i == 2:
            f_gausexp.SetParameters(h_res_proton[i].GetMaximum(),
                                    0,
                                    0.018,
                                    0.39)
            h_res_proton[i].Fit(f_gausexp,
                                "RQ",
                                "",
                                -0.1,
                                0.12)
        else:
            h_res_proton[i].Fit(f_gausexp,
                                "RQ",
                                "",
                                -h_res_proton[i].GetRMS(),
                                h_res_proton[i].GetRMS())

        a_res.append(f_gausexp.GetParameter(2) * 100)
        a_res_err_h.append(f_gausexp.GetParError(2) * 100 / 2)
        a_res_err_l.append(f_gausexp.GetParError(2) * 100 / 2)

        ROOT.gPad.SetLeftMargin(0.13)
        ROOT.gPad.SetBottomMargin(0.18)
    c_proton_res.Update()
    # c_proton_res.SaveAs("plots/h_proton_res.pdf")

    # ****************************************************
    # Proton energy resolution vs true energy
    # ****************************************************
    ROOT.gStyle.SetOptTitle(0)

    c_p_res_e = ROOT.TCanvas("c_p_res_e")
    g_p_res_e = ROOT.TGraphAsymmErrors(6, a_bins, a_res, x_errs_low, x_errs_high, a_res_err_l, a_res_err_h)
    g_p_res_e.Draw("ap")
    g_p_res_e.GetYaxis().SetTitle("#sigma [%]")
    g_p_res_e.GetXaxis().SetTitle("E_{p_{k}} [GeV]")
    g_p_res_e.GetXaxis().SetTitleSize(0.05)
    g_p_res_e.GetXaxis().SetTitleOffset(1)
    g_p_res_e.GetYaxis().SetTitleOffset(0.9)
    c_p_res_e.SetBottomMargin(0.13)
    g_p_res_e.SetMarkerStyle(20)
    f_res = ROOT.TF1("f_res", "sqrt(([0]/sqrt(x))**2+([1]/x)**2+[2]**2)", 0, 2)
    f_res.SetParNames("a", "b", "c")
    f_res.SetParLimits(0, 0.1, 1)
    f_res.SetParLimits(1, 0.1, 1)
    f_res.SetParLimits(2, 0.1, 2)

    g_p_res_e.Fit(f_res)
    l_res = ROOT.TLegend(0.37, 0.68, 0.80, 0.86)
    l_res.AddEntry(f_res, "(#frac{%.2f}{#sqrt{E / GeV}} #oplus #frac{%.2f}{E / GeV} #oplus %.2f) %%" %
                   (f_res.GetParameter(0),
                    f_res.GetParameter(1),
                    f_res.GetParameter(2)), "l")
    l_res.Draw()
    f_p_res_e = ROOT.TFile("plots/g_p_res_e.root", "RECREATE")
    g_p_res_e.Write()
    f_p_res_e.Close()
    c_p_res_e.Update()
    c_p_res_e.SaveAs("plots/h_p_res_e.pdf")

    c_res_total = ROOT.TCanvas("c_res_total")
    h_res_proton_total.Draw("hist")
    h_res_proton_total.Fit(f_gausexp,
                           "",
                           "",
                           -0.3,
                           0.12)
    f_gausexp.SetNpx(1000)
    c_res_total.SetBottomMargin(0.13)
    h_res_proton_total.GetXaxis().SetTitleOffset(1.1)
    c_res_total.SetLeftMargin(0.13)
    f_gausexp.Draw("same")
    c_res_total.Update()

input()
