#!/usr/local/bin/python3

from glob import glob
from array import array
from tqdm import tqdm_notebook

import ROOT
from bdt_common import is_fiducial, gauss_exp
from settings import ELECTRON_MASS, ELECTRON_THRESHOLD
from settings import PROTON_MASS, PROTON_THRESHOLD
from settings import variables, spectators
from settings import fill_kin_branches, pre_cuts
from proton_energy import length2energy


ROOT.gStyle.SetOptTitle(1)
ROOT.gStyle.SetTitleFillColor(ROOT.kWhite)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetNumberContours(99)
ROOT.gStyle.SetStatW(0.16)
ROOT.gStyle.SetOptFit(0)


def calibration_graph_proton(h_reco, h_true):
    a_reco_true = array("f", [])

    a_bins = array("f", [])

    x_errs_low = array("f", [])
    x_errs_high = array("f", [])
    y_errs_low = array("f", [])
    y_errs_high = array("f", [])

    f_gaussexp = ROOT.TF1("f_gausexp", gauss_exp, 0, 2, 4)
    f_gaussexp.SetParNames("n", "mu", "sigma", "k")

    f_gaussexp.SetParLimits(3, 0.001, 1)
    f_gaussexp.SetParameter(3, 0.5)

    x_width = (h_true[-1].GetXaxis().GetXmax() - h_true[-1].GetXaxis().GetXmin())  / len(h_reco)

    for i, i_bin in enumerate(h_reco):
        f_gaussexp.SetParameter(2, h_reco[i].GetStdDev())
        f_gaussexp.SetParLimits(2, h_reco[i].GetStdDev()-0.1, h_reco[0].GetStdDev()+0.1)
        true_center = h_true[i].GetBinCenter(h_true[i].GetMaximumBin())
        f_gaussexp.SetParLimits(1, true_center-0.05, true_center+0.01)
        if i == 7:
            f_gaussexp.SetParLimits(1, true_center-0.05, 0.42)
        f_gaussexp.SetParameter(1, true_center)
        maximum_x = i_bin.GetMean()
        f_gaussexp.SetParameter(1, h_true[i].GetMean())
        i_bin.Fit(f_gaussexp, "Q")
        a_reco_true.append(f_gaussexp.GetParameter(1))
        print(i, f_gaussexp.GetParameter(1))
        hm_1 = f_gaussexp.GetX(f_gaussexp.GetMaximum()/2, 0, maximum_x)
        hm_2 = f_gaussexp.GetX(f_gaussexp.GetMaximum()/2, maximum_x, 2)
        y_errs_low.append(maximum_x - hm_1)
        y_errs_high.append(hm_2 - maximum_x)
        x_value = true_center
        a_bins.append(x_value)
        x_err_h = (i + 1) * x_width - x_value + PROTON_THRESHOLD
        x_err_l = x_value - i * x_width - PROTON_THRESHOLD
        x_errs_high.append(x_err_h)
        x_errs_low.append(x_err_l)

    g_reco_true = ROOT.TGraphAsymmErrors(len(h_reco),
                                         a_bins,
                                         a_reco_true,
                                         x_errs_low, x_errs_high,
                                         y_errs_low, y_errs_high)

    return g_reco_true, a_bins, x_errs_low, x_errs_high


def calibration_graph_neutrino(h_reco, h_true):
    a_reco_true = array("f", [])

    a_bins = array("f", [])

    x_errs_low = array("f", [])
    x_errs_high = array("f", [])
    y_errs_low = array("f", [])
    y_errs_high = array("f", [])

    f_gaussexp = ROOT.TF1("f_gausexp", gauss_exp, 0, 2, 4)
    f_gaussexp.SetParNames("n", "mu", "sigma", "k")
    f_gaussexp.SetParameter(2, h_reco[0].GetStdDev())
    f_gaussexp.SetParameter(3, 0.5)

    f_gaussexp.SetParLimits(2, 0.01, 2)
    f_gaussexp.SetParLimits(3, 0.001, 1)

    x_width = (h_true[-1].GetXaxis().GetXmax() -
               h_true[-1].GetXaxis().GetXmin()) / len(h_reco)
    for i, i_bin in enumerate(h_reco):
        if i == 2:
            f_gaussexp.SetParLimits(1, 0, 0.55)
        else:
            f_gaussexp.SetParLimits(1, 0, 2)

        if i > 0:
            i_bin.Fit(f_gaussexp, "Q")
        maximum = f_gaussexp.GetMaximum()
        maximum_x = f_gaussexp.GetX(maximum)
        hm_1 = f_gaussexp.GetX(f_gaussexp.GetMaximum()/2, 0, maximum_x)
        hm_2 = f_gaussexp.GetX(f_gaussexp.GetMaximum()/2,
                               maximum_x, h_true[-1].GetXaxis().GetXmax())

        a_reco_true.append(maximum_x)
        y_errs_low.append(maximum_x - hm_1)
        y_errs_high.append(hm_2 - maximum_x)
        x_value = h_true[i].GetMean()
        if i < 1:
            x_value = -1
        a_bins.append(x_value)
        x_err_h = (i + 1) * x_width - x_value
        x_err_l = x_value - i * x_width
        x_errs_high.append(x_err_h)
        x_errs_low.append(x_err_l)

    g_reco_true = ROOT.TGraphAsymmErrors(len(h_reco),
                                         a_bins,
                                         a_reco_true,
                                         x_errs_low, x_errs_high,
                                         y_errs_low, y_errs_high)

    return g_reco_true, a_bins, x_errs_low, x_errs_high


def calibration_graph_deposited(h_reco, h_true):
    a_reco_true = array("f", [])

    a_bins = array("f", [])

    x_errs_low = array("f", [])
    x_errs_high = array("f", [])
    y_errs_low = array("f", [])
    y_errs_high = array("f", [])

    f_gaussexp = ROOT.TF1("f_gausexp", gauss_exp, 0, 2, 4)
    f_gaussexp.SetParNames("n", "mu", "sigma", "k")
    f_gaussexp.SetParameter(2, h_reco[0].GetStdDev())
    f_gaussexp.SetParameter(3, 0.5)

    f_gaussexp.SetParLimits(2, 0.01, 2)
    f_gaussexp.SetParLimits(3, 0.001, 1)

    x_width = (h_true[-1].GetXaxis().GetXmax() -
               h_true[-1].GetXaxis().GetXmin()) / len(h_reco)
    for i, i_bin in enumerate(h_reco):
        if i == 2:
            f_gaussexp.SetParLimits(1, 0, 0.55)
        else:
            f_gaussexp.SetParLimits(1, 0, 2)

        if i > 0:
            i_bin.Fit(f_gaussexp, "Q")
        maximum = f_gaussexp.GetMaximum()
        maximum_x = f_gaussexp.GetX(maximum)
        hm_1 = f_gaussexp.GetX(f_gaussexp.GetMaximum()/2, 0, maximum_x)
        hm_2 = f_gaussexp.GetX(f_gaussexp.GetMaximum()/2,
                               maximum_x, h_true[-1].GetXaxis().GetXmax())

        a_reco_true.append(maximum_x)
        y_errs_low.append(maximum_x - hm_1)
        y_errs_high.append(hm_2 - maximum_x)
        x_value = h_true[i].GetMean()
        if i < 1:
            x_value = -1
        a_bins.append(x_value)
        x_err_h = (i + 1) * x_width - x_value
        x_err_l = x_value - i * x_width
        x_errs_high.append(x_err_h)
        x_errs_low.append(x_err_l)

    g_reco_true = ROOT.TGraphAsymmErrors(len(h_reco),
                                         a_bins,
                                         a_reco_true,
                                         x_errs_low, x_errs_high,
                                         y_errs_low, y_errs_high)

    return g_reco_true, a_bins, x_errs_low, x_errs_high


def calibration_graph(h_reco, h_true):
    a_reco_true = array("f", [])

    a_bins = array("f", [])

    x_errs_low = array("f", [])
    x_errs_high = array("f", [])
    y_errs_low = array("f", [])
    y_errs_high = array("f", [])

    f_gaussexp = ROOT.TF1("f_gausexp", gauss_exp, 0, 2, 4)
    f_gaussexp.SetParNames("n", "mu", "sigma", "k")
    f_gaussexp.SetParameter(2, h_reco[0].GetStdDev())
    f_gaussexp.SetParameter(3, 0.5)

    f_gaussexp.SetParLimits(1, 0, 1.5)
    f_gaussexp.SetParLimits(2, 0.01, 2)
    f_gaussexp.SetParLimits(3, 0.001, 1)

    x_width = (h_true[-1].GetXaxis().GetXmax() -
               h_true[-1].GetXaxis().GetXmin()) / len(h_reco)

    for i, i_bin in enumerate(h_reco):
        i_bin.Fit(f_gaussexp, "Q")
        maximum = f_gaussexp.GetMaximum()
        maximum_x = f_gaussexp.GetX(maximum)
        hm_1 = f_gaussexp.GetX(f_gaussexp.GetMaximum()/2, 0, maximum_x)
        hm_2 = f_gaussexp.GetX(f_gaussexp.GetMaximum()/2, maximum_x, 2)
        a_reco_true.append(maximum_x)
        y_errs_low.append(maximum_x - hm_1)
        y_errs_high.append(hm_2 - maximum_x)
        x_value = h_true[i].GetMean()
        a_bins.append(x_value)
        x_err_h = (i + 1) * x_width - x_value + ELECTRON_THRESHOLD
        x_err_l = x_value - i * x_width - ELECTRON_THRESHOLD
        x_errs_high.append(x_err_h)
        x_errs_low.append(x_err_l)

    g_reco_true = ROOT.TGraphAsymmErrors(len(h_reco),
                                         a_bins,
                                         a_reco_true,
                                         x_errs_low, x_errs_high,
                                         y_errs_low, y_errs_high)

    return g_reco_true, a_bins, x_errs_low, x_errs_high


def select_proton_event(c, var_dict):
    p = 0
    proton_energy = 0
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
            if (energy - PROTON_MASS) > PROTON_THRESHOLD and is_fiducial(p_start) and is_fiducial(p_end):
                proton_energy = energy - PROTON_MASS
                p += 1

    true_neutrino_vertex = [c.true_vx_sce,
                            c.true_vy_sce,
                            c.true_vz_sce]

    if not c.passed or not p == 1 or not is_fiducial(true_neutrino_vertex):
        return False

    if not c.category == 2:
        return False

    # If there are no tracks we require at least two showers
    showers_2_tracks_0 = True
    if c.n_tracks == 0 and c.n_showers == 1:
        showers_2_tracks_0 = False

    if not showers_2_tracks_0:
        return False

    fill_kin_branches(c, 1, var_dict, "nue", False)

    if not pre_cuts(var_dict):
        return False

    return proton_energy


def select_electron_event(c, var_dict):
    e = 0
    electron_energy = 0
    for i_pdg, energy in enumerate(c.nu_daughters_E):

        p_start = [
            c.nu_daughters_vx[i_pdg],
            c.nu_daughters_vy[i_pdg],
            c.nu_daughters_vz[i_pdg],
        ]

        if abs(c.nu_daughters_pdg[i_pdg]) == 11:
            if not (is_fiducial(p_start)):
                e = 0
                break
            if energy - ELECTRON_MASS > ELECTRON_THRESHOLD:
                electron_energy = energy
                e += 1

    true_neutrino_vertex = [c.true_vx_sce,
                            c.true_vy_sce,
                            c.true_vz_sce]

    if not c.passed or not e or not is_fiducial(true_neutrino_vertex):
        return False

    if not c.category == 2:
        return False

    # If there are no tracks we require at least two showers
    showers_2_tracks_0 = True
    if c.n_tracks == 0 and c.n_showers == 1:
        showers_2_tracks_0 = False

    if not showers_2_tracks_0:
        return False

    fill_kin_branches(c, 1, var_dict, "nue", False)

    if not pre_cuts(var_dict):
        return False

    return electron_energy


def electron_resolution(file_path, bias, offset, n_bins=10, max_energy=2, scale=1):
    files = glob("%s/*.root" % file_path)
    c = ROOT.TChain("robertoana/pandoratree")

    for f in files:
        c.Add(f)

    entries = int(c.GetEntries() / scale)
    h_res_electron = []

    h_binning = ROOT.TH1F("h_binning", "", n_bins,
                            ELECTRON_THRESHOLD, max_energy)
    for i in range(n_bins):
        h_res_electron.append(ROOT.TH1F("h_res_electron%i" % i,
                                        "%.0f MeV < E^{e} < %.0f MeV;(E_{corr}^{e} - E^{e})/E^{e}; N. Entries / 0.05"
                                        % (h_binning.GetBinLowEdge(i + 1) * 1000, h_binning.GetBinLowEdge(i + 2) * 1000),
                                        50, -1, 1))

    var_dict = dict(variables + spectators)

    for i in tqdm_notebook(range(entries)):
        c.GetEntry(i)

        electron_energy = select_electron_event(c, var_dict)
        if not electron_energy:
            continue

        reco_electron_energy = 0
        i_bin = int(electron_energy / (max_energy / n_bins))

        for i_sh in range(c.n_showers):
            if c.matched_showers[i_sh] == 11:
                reco_electron_energy += c.shower_energy[i_sh][2] * \
                    c.shower_energy_cali[i_sh][2]

        if i_bin < n_bins and reco_electron_energy:
            corr_energy = (reco_electron_energy - offset) / bias
            h_res_electron[i_bin].Fill((corr_energy - electron_energy) / electron_energy)

    return h_res_electron


def electron_calibration(file_path, n_bins=10, max_energy=2, scale=1):
    files = glob("%s/*.root" % file_path)
    c = ROOT.TChain("robertoana/pandoratree")

    for f in files:
        c.Add(f)

    entries = int(c.GetEntries() / scale)
    h_reco_true = ROOT.TH2F("h_reco_true",
                            ";E^{e} [GeV];E^{e}_{reco} [GeV]",
                            50, ELECTRON_THRESHOLD, max_energy,
                            50, ELECTRON_THRESHOLD, max_energy)

    h_reco = []
    h_true = []
    h_binning = ROOT.TH1F("h_binning", "", n_bins,
                          ELECTRON_THRESHOLD, max_energy)

    for i in range(n_bins):
        h_reco.append(ROOT.TH1F("h_reco_electron%i" % i,
                                "%.0f MeV < E^{e} < %.0f MeV;Energy[GeV];N. Entries / 0.02 GeV"
                                % (h_binning.GetBinLowEdge(i + 1) * 1000, h_binning.GetBinLowEdge(i + 2) * 1000),
                                60, ELECTRON_THRESHOLD, max_energy))

        h_true.append(ROOT.TH1F("h_true_electron%i" % i,
                                "%.0f MeV < E^{e} < %.0f MeV;Energy[GeV];N. Entries / 0.02 GeV"
                                % (h_binning.GetBinLowEdge(i + 1) * 1000, h_binning.GetBinLowEdge(i + 2) * 1000),
                                60, ELECTRON_THRESHOLD, max_energy))

    var_dict = dict(variables + spectators)

    for i in tqdm_notebook(range(entries)):
        c.GetEntry(i)

        electron_energy = select_electron_event(c, var_dict)
        if not electron_energy:
            continue

        reco_electron_energy = 0
        i_bin = h_binning.FindBin(electron_energy)

        for i_sh in range(c.n_showers):
            if c.matched_showers[i_sh] == 11:
                reco_electron_energy += c.shower_energy[i_sh][2] * \
                    c.shower_energy_cali[i_sh][2]

        if 0 < i_bin <= n_bins and reco_electron_energy:
            h_reco[i_bin-1].Fill(reco_electron_energy)
            h_true[i_bin-1].Fill(electron_energy)
            h_reco_true.Fill(electron_energy, reco_electron_energy)

    return h_reco_true, h_true, h_reco


def proton_calibration(file_path, n_bins=10, max_energy=1, scale=1):
    files = glob("%s/*.root" % file_path)
    c = ROOT.TChain("robertoana/pandoratree")

    for f in files:
        c.Add(f)

    entries = int(c.GetEntries() / scale)
    h_reco_true = ROOT.TH2F("h_reco_true",
                            ";E^{p} [GeV];E^{p}_{reco} [GeV]",
                            50, PROTON_THRESHOLD, max_energy,
                            50, PROTON_THRESHOLD, max_energy)

    h_reco = []
    h_true = []

    h_binning = ROOT.TH1F("h_binning", "", n_bins,
                          PROTON_THRESHOLD, max_energy)
    for i in range(n_bins):
        h_reco.append(ROOT.TH1F("h_reco_electron%i" % i,
                                "%.0f MeV < E^{p} < %.0f MeV;Energy[GeV];N. Entries / 0.02 GeV"
                                % (h_binning.GetBinLowEdge(i + 1) * 1000, h_binning.GetBinLowEdge(i + 2) * 1000),
                                30, PROTON_THRESHOLD, max_energy))

        h_true.append(ROOT.TH1F("h_true_electron%i" % i,
                                "%.0f MeV < E^{p} < %.0f MeV;Energy[GeV];N. Entries / 0.02 GeV"
                                % (h_binning.GetBinLowEdge(i + 1) * 1000, h_binning.GetBinLowEdge(i + 2) * 1000),
                                30, PROTON_THRESHOLD, max_energy))

    var_dict = dict(variables + spectators)

    for i in tqdm_notebook(range(entries)):
        c.GetEntry(i)

        proton_energy = select_proton_event(c, var_dict)
        if not proton_energy:
            continue

        reco_proton_energy = 0
        i_bin = h_binning.FindBin(proton_energy)

        for i_tr in range(c.n_tracks):
            if c.matched_tracks[i_tr] == 2212:
                reco_proton_energy += length2energy(c.track_len[i_tr])

        if 0 < i_bin <= n_bins and reco_proton_energy:
            h_reco[i_bin-1].Fill(reco_proton_energy)
            h_true[i_bin-1].Fill(proton_energy)
            h_reco_true.Fill(proton_energy, reco_proton_energy)

    return h_reco_true, h_true, h_reco


def check_cc0pinp(chain_nue):
    """Checks if the event is a contained electron neutrino
    CC0π-Np interaction

    Args:
        chain_nue: ROOT TChain of the events

    Returns:
        True if the event is a contained CC0π-Np interaction
    """
    protons = 0
    electrons = 0
    photons = 0
    pions = 0
    electron_energy = 0
    proton_energy = 0

    for i, energy in enumerate(chain_nue.nu_daughters_E):
        if chain_nue.nu_daughters_pdg[i] == 2212:
            p_vertex = [chain_nue.nu_daughters_vx[i],
                        chain_nue.nu_daughters_vy[i],
                        chain_nue.nu_daughters_vz[i]]

            p_end = [chain_nue.nu_daughters_endx[i],
                     chain_nue.nu_daughters_endy[i],
                     chain_nue.nu_daughters_endz[i]]

            if energy - PROTON_MASS > PROTON_THRESHOLD:
                protons += 1
                proton_energy += energy - PROTON_MASS
            if not is_fiducial(p_vertex) or not is_fiducial(p_end):
                protons = 0
                break

        if chain_nue.nu_daughters_pdg[i] == 11:
            e_vertex = [chain_nue.nu_daughters_endx[i],
                        chain_nue.nu_daughters_endy[i],
                        chain_nue.nu_daughters_endz[i]]

            if energy - ELECTRON_MASS > ELECTRON_THRESHOLD and is_fiducial(e_vertex):
                electron_energy += energy
                electrons += 1

        if chain_nue.nu_daughters_pdg[i] == 22:
            # if energy > 0.035:
            photons += 1

        if chain_nue.nu_daughters_pdg[i] == 111 or chain_nue.nu_daughters_pdg[i] == 211:
            # if energy > 0.06:
            pions += 1

    if electrons == 1 and pions == 0 and protons >= 1 and photons == 0 and chain_nue.ccnc == 0:
        return electron_energy + proton_energy

    return False


def neutrino_calibration(file_path, n_bins=10, max_energy=2, scale=1,
                         proton_calib=(1, 0), electron_calib=(0.77, 0.02), neutrino_calib=(0.9, 0)):
    files = glob("%s/*.root" % file_path)
    c = ROOT.TChain("robertoana/pandoratree")

    for f in files:
        c.Add(f)

    entries = int(c.GetEntries() / scale)
    h_reco_true = ROOT.TH2F("h_reco_true",
                            ";E^{#nu_{e}} [GeV];E^{#nu_{e}}_{reco} [GeV]",
                            200, 0, max_energy*5,
                            200, 0, max_energy*5)
    h_reco_true_corr = ROOT.TH2F("h_reco_true_corr",
                                 ";E^{#nu_{e}} [GeV];E^{#nu_{e}}_{reco} [GeV]",
                                 200, 0, max_energy*5,
                                 200, 0, max_energy*5)

    h_reco = []
    h_true = []
    h_reco_corr = []

    h_binning = ROOT.TH1F("h_binning", "", n_bins,
                          0, max_energy)
    for i in range(n_bins):
        h_reco.append(ROOT.TH1F("h_reco_neutrino%i" % i,
                                "%.0f MeV < E^{#nu_{e}} < %.0f MeV;Energy[GeV];N. Entries / 0.02 GeV"
                                % (h_binning.GetBinLowEdge(i + 1) * 1000, h_binning.GetBinLowEdge(i + 2) * 1000),
                                30, 0, max_energy))
        h_reco_corr.append(ROOT.TH1F("h_reco_corr_neutrino%i" % i,
                                     "%.0f MeV < E^{#nu_{e}} < %.0f MeV;Energy[GeV];N. Entries / 0.02 GeV"
                                     % (h_binning.GetBinLowEdge(i + 1) * 1000, h_binning.GetBinLowEdge(i + 2) * 1000),
                                     30, 0, max_energy))
        h_true.append(ROOT.TH1F("h_true_neutrinon%i" % i,
                                "%.0f MeV < E^{#nu_{e}} < %.0f MeV;Energy[GeV];N. Entries / 0.02 GeV"
                                % (h_binning.GetBinLowEdge(i + 1) * 1000, h_binning.GetBinLowEdge(i + 2) * 1000),
                                30, 0, max_energy))

    var_dict = dict(variables + spectators)

    for i in tqdm_notebook(range(entries)):
        c.GetEntry(i)

        eNp_interaction = check_cc0pinp(c)

        if not eNp_interaction:
            continue

        fill_kin_branches(c, 1, var_dict, "nue", False)
        if not pre_cuts(var_dict):
            continue

        neutrino_energy = c.nu_E
        reco_neutrino_energy = 0
        i_bin = h_binning.FindBin(neutrino_energy)

        for i_tr in range(c.n_tracks):
            if c.matched_tracks[i_tr] == 2212:
                reco_neutrino_energy += (length2energy(c.track_len[i_tr]) - proton_calib[1]) / proton_calib[0]

        for i_sh in range(c.n_showers):
            if c.matched_showers[i_sh] == 11:
                reco_neutrino_energy += (c.shower_energy[i_sh][2] * c.shower_energy_cali[i_sh][2] - electron_calib[1]) / electron_calib[0]

        if reco_neutrino_energy > 0.02:
            reco_corr = (reco_neutrino_energy - neutrino_calib[1]) / neutrino_calib[0]
            if 0 < i_bin <= n_bins:
                h_reco[i_bin-1].Fill(reco_neutrino_energy)
                h_true[i_bin-1].Fill(neutrino_energy)
                h_reco_corr[i_bin-1].Fill(reco_corr)

            h_reco_true_corr.Fill(neutrino_energy, reco_corr)
            h_reco_true.Fill(neutrino_energy, reco_neutrino_energy)

    return h_reco_true, h_true, h_reco, h_reco_true_corr, h_reco_corr


def deposited_calibration(file_path, n_bins=10, max_energy=2, scale=1,
                         proton_calib=(1, 0), electron_calib=(0.77, 0.02), calib=(0.97, 0)):
    files = glob("%s/*.root" % file_path)
    c = ROOT.TChain("robertoana/pandoratree")

    for f in files:
        c.Add(f)

    entries = int(c.GetEntries() / scale)
    h_reco_true = ROOT.TH2F("h_reco_true",
                            ";E_{k} [GeV];E_{corr} [GeV]",
                            300, 0, 3,
                            300, 0, 3)
    h_reco_true_corr = ROOT.TH2F("h_reco_true_corr",
                                 ";E_{k} [GeV];E_{corr} [GeV]",
                                 300, 0, 3,
                                 300, 0, 3)

    h_reco = []
    h_true = []
    h_reco_corr = []

    h_binning = ROOT.TH1F("h_binning", "", n_bins,
                          0, max_energy)
    for i in range(n_bins):
        h_reco.append(ROOT.TH1F("h_reco_neutrino%i" % i,
                                "%.0f MeV < E_{deposited} < %.0f MeV;Energy[GeV];N. Entries / 0.02 GeV"
                                % (h_binning.GetBinLowEdge(i + 1) * 1000, h_binning.GetBinLowEdge(i + 2) * 1000),
                                30, 0, max_energy))
        h_reco_corr.append(ROOT.TH1F("h_reco_corr_neutrino%i" % i,
                                     "%.0f MeV < E_{deposited} < %.0f MeV;Energy[GeV];N. Entries / 0.02 GeV"
                                     % (h_binning.GetBinLowEdge(i + 1) * 1000, h_binning.GetBinLowEdge(i + 2) * 1000),
                                     30, 0, max_energy))
        h_true.append(ROOT.TH1F("h_true_neutrinon%i" % i,
                                "%.0f MeV < E_{deposited} < %.0f MeV;Energy[GeV];N. Entries / 0.02 GeV"
                                % (h_binning.GetBinLowEdge(i + 1) * 1000, h_binning.GetBinLowEdge(i + 2) * 1000),
                                30, 0, max_energy))

    var_dict = dict(variables + spectators)

    for i in tqdm_notebook(range(entries)):
        c.GetEntry(i)

        eNp_interaction = check_cc0pinp(c)

        if not eNp_interaction:
            continue

        fill_kin_branches(c, 1, var_dict, "nue", False)
        if not pre_cuts(var_dict):
            continue

        neutrino_energy = eNp_interaction
        reco_neutrino_energy = 0
        i_bin = h_binning.FindBin(neutrino_energy)

        for i_tr in range(c.n_tracks):
            if c.matched_tracks[i_tr] == 2212:
                reco_neutrino_energy += (length2energy(
                    c.track_len[i_tr]) - proton_calib[1]) / proton_calib[0]

        for i_sh in range(c.n_showers):
            if c.matched_showers[i_sh] == 11:
                reco_neutrino_energy += (
                    c.shower_energy[i_sh][2] * c.shower_energy_cali[i_sh][2] - electron_calib[1]) / electron_calib[0]

        if reco_neutrino_energy > 0.02:
            reco_corr = (reco_neutrino_energy -
                         calib[1]) / calib[0]
            if 0 < i_bin <= n_bins:
                h_reco[i_bin-1].Fill(reco_neutrino_energy)
                h_true[i_bin-1].Fill(neutrino_energy)
                h_reco_corr[i_bin-1].Fill(reco_corr)

            h_reco_true_corr.Fill(neutrino_energy, reco_corr)
            h_reco_true.Fill(neutrino_energy, reco_neutrino_energy)

    return h_reco_true, h_true, h_reco, h_reco_true_corr, h_reco_corr
