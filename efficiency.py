#!/usr/local/bin/python3

import math
import statistics
from glob import glob

import ROOT
from tqdm import tqdm_notebook
from settings import binning, labels
from settings import is_fiducial, variables, spectators
from settings import pre_cuts, load_bdt, apply_cuts, fill_kin_branches
from settings import N_UNI, ELECTRON_THRESHOLD, PROTON_THRESHOLD
from settings import PROTON_MASS, ELECTRON_MASS

ELECTRON_THRESHOLD = 0.020
PROTON_THRESHOLD = 0.040

ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptFit(0)
ROOT.gStyle.SetPalette(87)
ROOT.gStyle.SetNumberContours(99)


class Efficiency:
    """Class describing several selection efficiencies"""

    def __init__(self, tot, selected={}, var="", pot=0, total_sys_err=0, sys_err=0, sys_2d=0):
        """Initialize the Efficiency object

        Args:
            tot: TH1F object of the denominator
            selected: dictionary of the different numerators
            var: Efficiency variable
            pot: Number of protons-on-target
            total_sys_err: Total systematic error
            sys_err: Array of systematic errors (one per bin)
            sys_2d: 2D histogram of the systematic errors
        """
        self.tot = tot
        self.pot = pot
        self.selected = selected
        self.efficiency = {}
        self.efficiency_err = {}
        self.variable = var
        self.sys_err = sys_err
        self.sys_2d = sys_2d
        self.total_sys_err = total_sys_err
        for sel in selected:
            eff_value = selected[sel].Integral()/tot.Integral()
            self.efficiency[sel] = eff_value
            self.efficiency_err[sel] = 1./tot.Integral() * math.sqrt(selected[sel].Integral()*(1-eff_value))

    def histo_divided(self, num):
        h_num = self.selected[num].Clone()
        h_num.Divide(self.tot)

        return h_num

    def tefficiency(self, num):
        """Creates a ROOT TEfficiency object

        Args:
            num: name of the selected numerator

        Returns:
            ROOT TEfficiency object for the num selection
        """
        clone_tot = self.tot.Clone()
        clone_tot.GetYaxis().SetTitleOffset(0.9)
        clone_tot.SetTitle(";%s;%s" % (clone_tot.GetXaxis().GetTitle(), "Efficiency"))
        teff = ROOT.TEfficiency(self.selected[num], clone_tot)
        return teff

    def tefficiency_err(self, num):
        """Creates a TH1F with bin content corresponding to the efficiency and
        bin error correspondin to the quadratic sum of systematic error and
        statistical error

        Args:
            num: name of the selected numerator

        Returns:
            ROOT TH1F histogram
        """
        teff = self.tefficiency(num)
        teff_sys = self.selected[num].Clone()
        teff_sys.Divide(self.tot)

        for i_bin in range(1, self.selected[num].GetNbinsX()+1):
            stat_error = teff.GetEfficiencyErrorLow(i_bin) + teff.GetEfficiencyErrorUp(i_bin)
            teff_sys.SetBinError(i_bin, math.sqrt(stat_error**2 + self.sys_err[i_bin - 1]**2))

        return teff_sys


    def draw(self, systematics=True):
        """Draws the efficiency in a ROOT Canvas with statistical
        and systematic errors

        Args:
            systematics: draw the systematic error or not

        Returns:
            The canvas, the TEfficiency objects and the legend
        """
        c_eff_sys = ROOT.TCanvas("c_eff_sys_%s" % self.tot.GetName(),
                                 self.tot.GetName(),
                                 640, 480)
        eff_stat = self.tefficiency("passed")
        eff_sys = self.tefficiency_err("passed")
        eff_stat.SetLineWidth(2)
        eff_stat.Draw("ap")
        eff_stat.SetTitle("Efficiency")
        eff_sys.SetLineColor(ROOT.kRed + 1)
        eff_stat.SetLineColor(ROOT.kRed + 1)

        if systematics:
            self.sys_2d.Draw("colz same")
            eff_sys.Draw("e same")
            eff_stat.Draw("p same")

        l_eff = ROOT.TLegend(0.097, 0.854, 0.9, 0.9)

        if systematics:
            l_eff.AddEntry(eff_stat,
                           "Selection efficiency: (%.1f #pm %.1f (stat) #pm %.1f (sys)) %%"
                           % (self.efficiency["passed"] * 100, self.efficiency_err["passed"] * 100, self.total_sys_err * 100),
                           "le")
        else:
            l_eff.AddEntry(eff_stat,
                           "Selection efficiency: (%.1f #pm %.1f) %%"
                           % (self.efficiency["passed"] * 100, self.efficiency_err["passed"] * 100),
                           "le")

        l_eff.Draw()
        c_eff_sys.SetTopMargin(0.18)

        if systematics:
            return c_eff_sys, eff_stat, eff_sys, self.sys_2d, l_eff
        else:
            return c_eff_sys, eff_stat, l_eff

def draw_true(histo):
    """Draws the histogram in a ROOT Canvas with a legend
    showing the number of events

    Args:
        histo: Histogram to be drawn

    Returns:
        The canvas, the histogram and the legend
    """
    tot = histo.Clone()
    tot.SetLineColor(ROOT.TColor.GetColor("#62b570"))
    tot.SetLineWidth(3)
    leg = ROOT.TLegend(0.097, 0.854, 0.9, 0.91)
    leg.SetNColumns(2)
    leg.AddEntry(tot,
                 "Generated #nu_{e} CC0#pi-Np: %.1f events" % tot.Integral(),
                 "f")
    canvas = ROOT.TCanvas("c_tot", "True %s" % histo.GetName(), 640, 480)
    tot.Draw("hist")
    tot.GetYaxis().SetTitleOffset(0.9)
    canvas.SetTopMargin(0.16)

    return canvas, tot, leg


def stacked(eff, h_stack, colors):
    """Draws a stacked histogram of every selection in the Efficiency
    object

    Args:
        eff: Efficiency object
        h_stack: THStack object to be drawn
        colors: Hexadecimal colors of the THStack histograms

    Returns:
        The THStack and the legend
    """
    h_cats = []
    leg = ROOT.TLegend(0.1, 0.79, 0.9, 0.91)

    for cat in eff.selected:
        h_cat = eff.selected[cat].Clone()
        eff_value = eff.efficiency[cat]
        h_cat.Divide(eff.tot)
        if colors:
            h_cat.SetFillColor(ROOT.TColor.GetColor(colors[cat]))
        h_cat.SetLineColor(ROOT.TColor.GetColor("#555555"))
        h_cat.SetLineWidth(1)
        if cat == "passed":
            h_cat.SetLineWidth(3)
            h_cat.SetLineColor(1)
        h_cats.append(h_cat)
        leg.AddEntry(h_cat,
                     "%s: %.1f%%" % (h_cat.GetName(), eff_value*100),
                     "f")
        h_stack.Add(h_cat)

    return h_cats, leg


def check_reco_fidvol(chain_nue):
    """Checks if the neutrino vertex, the reconstructed tracks
    and the starting point of reconstructed showers are contained
    in the fiducial volume

    Args:
        chain_nue: ROOT TChain of the events

    Returns:
        True if the event is fully contained
    """
    track_fidvol = True

    for i_tr in range(chain_nue.n_tracks):
        track_start = [
            chain_nue.track_start_x[i_tr],
            chain_nue.track_start_y[i_tr],
            chain_nue.track_start_z[i_tr]
        ]
        track_end = [
            chain_nue.track_end_x[i_tr],
            chain_nue.track_end_y[i_tr],
            chain_nue.track_end_z[i_tr]
        ]

        track_fidvol = track_fidvol and is_fiducial(
            track_start) and is_fiducial(track_end)
        if not track_fidvol:
            break

    shower_fidvol = True

    for i_sh in range(chain_nue.n_showers):
        shower_start = [
            chain_nue.shower_start_x[i_sh],
            chain_nue.shower_start_y[i_sh],
            chain_nue.shower_start_z[i_sh]
        ]

        shower_fidvol = shower_fidvol and is_fiducial(shower_start)
        if not shower_fidvol:
            break

    neutrino_vertex = [chain_nue.vx, chain_nue.vy, chain_nue.vz]
    neutrino_fidvol = is_fiducial(neutrino_vertex)

    return track_fidvol and shower_fidvol and neutrino_fidvol


def check_1e1p(chain_nue):
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
        vertex = [chain_nue.nu_daughters_start_v[i][0],
                  chain_nue.nu_daughters_start_v[i][1],
                  chain_nue.nu_daughters_start_v[i][2]]

        end = [chain_nue.nu_daughters_end_v[i][0],
                chain_nue.nu_daughters_end_v[i][1],
                chain_nue.nu_daughters_end_v[i][2]]

        if chain_nue.nu_daughters_pdg[i] == 2212:
            if energy - PROTON_MASS > PROTON_THRESHOLD and is_fiducial(vertex) and is_fiducial(end):
                protons += 1
                proton_energy = max(energy - PROTON_MASS, proton_energy)

        if chain_nue.nu_daughters_pdg[i] == 11:
            if energy - ELECTRON_MASS > ELECTRON_THRESHOLD and is_fiducial(vertex):
                electron_energy += energy - ELECTRON_MASS
                electrons += 1

        if chain_nue.nu_daughters_pdg[i] == 22:
            if energy > ELECTRON_THRESHOLD:
                photons += 1

        if chain_nue.nu_daughters_pdg[i] == 111:
            if energy - 0.135 > ELECTRON_THRESHOLD and is_fiducial(vertex):
                pions += 1

        if abs(chain_nue.nu_daughters_pdg[i]) == 211:
            if energy - 0.140 > PROTON_THRESHOLD and is_fiducial(vertex) and is_fiducial(end):
                pions += 1


    is_1e1p = electrons == 1 and chain_nue.ccnc == 0
    if is_1e1p:
        return (proton_energy, electron_energy)

    return False

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
        vertex = [chain_nue.nu_daughters_start_v[i][0],
                  chain_nue.nu_daughters_start_v[i][1],
                  chain_nue.nu_daughters_start_v[i][2]]

        end = [chain_nue.nu_daughters_end_v[i][0],
               chain_nue.nu_daughters_end_v[i][1],
               chain_nue.nu_daughters_end_v[i][2]]

        if chain_nue.nu_daughters_pdg[i] == 2212:
            if energy - PROTON_MASS > PROTON_THRESHOLD and is_fiducial(vertex) and is_fiducial(end):
                protons += 1
                proton_energy = max(energy - PROTON_MASS, proton_energy)

        if chain_nue.nu_daughters_pdg[i] == 11:
            if energy - ELECTRON_MASS > ELECTRON_THRESHOLD and is_fiducial(vertex):
                electron_energy += energy - ELECTRON_MASS
                electrons += 1

        if chain_nue.nu_daughters_pdg[i] == 22:
            if energy > ELECTRON_THRESHOLD and is_fiducial(vertex):
                photons += 1

        if chain_nue.nu_daughters_pdg[i] == 111:
            if energy - 0.135 > ELECTRON_THRESHOLD and is_fiducial(vertex):
                pions += 1

        if abs(chain_nue.nu_daughters_pdg[i]) == 211:
            if energy - 0.140 > PROTON_THRESHOLD and is_fiducial(vertex) and is_fiducial(end):
                pions += 1

    is_CC0piNp = electrons == 1 and pions == 0 and protons >= 1 and photons == 0
    if is_CC0piNp:
        return (proton_energy, electron_energy)

    return False


def lepton_angles(chain_filter):
    """Return the true theta and phi angles of the electron

    Args:
        chain_filter: ROOT TChain of the events

    Returns:
        Theta and phi angles of the electron (in degrees)
    """
    theta, phi = -999, -999

    for i, pdg in enumerate(chain_filter.nu_daughters_pdg):
        if pdg == 11:
            px = chain_filter.nu_daughters_p[i][0]
            py = chain_filter.nu_daughters_p[i][1]
            pz = chain_filter.nu_daughters_p[i][2]
            theta = math.degrees(math.acos(pz/math.sqrt(px**2+py**2+pz**2)))
            phi = math.degrees(math.atan2(py, px))
            break

    return theta, phi


def eff_total_sys_error(h_tot_sys, h_passed_sys):
    ratios = []
    for uni in range(N_UNI):
        ratio = h_passed_sys[uni].Integral()/h_tot_sys[uni].Integral()
        ratios.append(ratio)

    return statistics.stdev(ratios)

def eff_2d(h_tot_sys, h_passed_sys, variable):
    h_2d = ROOT.TH2F("h_2d",
                     ";%s;%s" %
                     (h_tot_sys[0].GetXaxis().GetTitle(), "Efficiency"),
                     binning[variable][0], binning[variable][1], binning[variable][2],
                     150, 0, 1)
    for uni in range(N_UNI):
        eff = h_passed_sys[uni].Clone()
        eff.Divide(h_tot_sys[uni])
        for i_bin in range(1, eff.GetNbinsX() + 1):
            value = eff.GetBinContent(i_bin)
            center = eff.GetBinCenter(i_bin)
            h_2d.Fill(center, value)

    return h_2d

def eff_sys_error(h_tot_sys, h_passed_sys):
    sys_error = []

    for i_bin in range(1, h_tot_sys[0].GetNbinsX()+1):
        err = []
        for uni in range(N_UNI):
            if h_tot_sys[uni].GetBinContent(i_bin):
                num = h_passed_sys[uni].GetBinContent(i_bin)
                den = h_tot_sys[uni].GetBinContent(i_bin)
                err.append(num/den)

        if len(err) > 1:
            sys_error.append(statistics.stdev(err))
        else:
            sys_error.append(0)

    return sys_error


def bkg_efficiency(files_path, scale=1):

    tot_bkg = {}
    passed_bkg = {}
    for i in range(10):
        tot_bkg[i] = ROOT.TH1F("h_tot_bkg_%i" % i,
                               labels["energy"],
                               binning["energy"][0],
                               binning["energy"][1],
                               binning["energy"][2])

        passed_bkg[i] = ROOT.TH1F("h_passed_bkg_%i" % i,
                                  labels["energy"],
                                  binning["energy"][0],
                                  binning["energy"][1],
                                  binning["energy"][2])

    bnb_cosmic = glob(files_path+"/*/*.root")
    nue_cosmic = glob("data_files/mc_nue_sbnfit/*.root")
    dirt = glob("data_files/dirt/*.root")
    samples = [bnb_cosmic, nue_cosmic, dirt]

    for s in tqdm_notebook(samples):
        chain_nue = ROOT.TChain("robertoana/pandoratree")
        chain_filter = ROOT.TChain("nueFilter/filtertree")
        chain_pot = ROOT.TChain("nueFilter/pot")
        for filename in s:
            chain_filter.Add(filename)
            chain_nue.Add(filename)
            chain_pot.Add(filename)

        # for filename in nue_cosmic:
        #     chain_filter.Add(filename)
        #     chain_nue.Add(filename)
        #     chain_pot.Add(filename)

        entries = int(chain_filter.GetEntries() / scale)
        pot_entries = int(chain_pot.GetEntries() / scale)

        total_pot = 0
        for i in range(pot_entries):
            chain_pot.GetEntry(i)
            total_pot += chain_pot.pot

        total_pot /= 1.028
        weight = 4.341e+19 / total_pot
        if s == dirt:
            weight *= 0.5
        var_dict = dict(variables + spectators)

        i_nue = -1
        for i_evt in tqdm_notebook(range(entries)):
            chain_filter.GetEntry(i_evt)
            true_neutrino_vertex = [chain_filter.true_vx,
                                    chain_filter.true_vy,
                                    chain_filter.true_vz]

            true_category = 0
            if not is_fiducial(true_neutrino_vertex):
                true_category = 5

            if chain_filter.passed:
                i_nue += 1

            if is_fiducial(true_neutrino_vertex):
                if chain_filter.ccnc == 1:
                    true_category = 4
                elif chain_filter.ccnc == 0 and chain_filter.nu_pdg == 14:
                    true_category = 3


            energy = chain_filter.nu_energy
            tot_bkg[true_category].Fill(energy, weight)

            if chain_filter.passed:
                chain_nue.GetEntry(i_nue)
                true_neutrino_vertex_sce = [chain_nue.true_vx_sce,
                                            chain_nue.true_vy_sce,
                                            chain_nue.true_vz_sce]

                if chain_nue.numu_passed == 1:
                    continue

                showers_2_tracks_0 = True
                if chain_nue.n_tracks == 0 and chain_nue.n_showers == 1:
                    showers_2_tracks_0 = False

                if not showers_2_tracks_0:
                    continue

                if not check_reco_fidvol(chain_nue):
                    continue

                fill_kin_branches(chain_nue, chain_nue.bnbweight, var_dict, "bnb", False)

                if not pre_cuts(var_dict):
                    continue

                category = var_dict["category"][0]

                if not is_fiducial(true_neutrino_vertex_sce) and category != 0 and category != 6 and category != 1 and category != 7:
                    category = 5

                passed_bkg[category].Fill(chain_nue.nu_E, chain_nue.bnbweight * weight)

    effs = {}
    for t in tot_bkg:
        if t == 2 or t == 0:
            continue
        if tot_bkg[t].Integral():
            effs[t] = Efficiency(tot_bkg[t], {"passed": passed_bkg[t]}, "energy", total_pot)

    return effs


def bkg_efficiency_nue(files_path, scale=1):
    bnb_cosmic = glob(files_path+"/*.root")

    chain_nue = ROOT.TChain("robertoana/pandoratree")
    chain_filter = ROOT.TChain("nueFilter/filtertree")
    chain_pot = ROOT.TChain("nueFilter/pot")
    for filename in bnb_cosmic:
        chain_filter.Add(filename)
        chain_nue.Add(filename)
        chain_pot.Add(filename)

    entries = int(chain_filter.GetEntries() / scale)
    pot_entries = int(chain_pot.GetEntries() / scale)

    total_pot = 0
    for i in range(pot_entries):
        chain_pot.GetEntry(i)
        total_pot += chain_pot.pot

    total_pot /= 1.028

    tot_bkg = {}
    passed_bkg = {}
    for i in range(10):
        tot_bkg[i] = ROOT.TH1F("h_tot_bkg_%i" % i,
                               labels["energy"],
                               binning["energy"][0],
                               binning["energy"][1],
                               binning["energy"][2])

        passed_bkg[i] = ROOT.TH1F("h_passed_bkg_%i" % i,
                                  labels["energy"],
                                  binning["energy"][0],
                                  binning["energy"][1],
                                  binning["energy"][2])

    i_nue = -1
    var_dict = dict(variables + spectators)
    for i_evt in tqdm_notebook(range(entries)):
        chain_filter.GetEntry(i_evt)
        true_neutrino_vertex = [chain_filter.true_vx,
                                chain_filter.true_vy,
                                chain_filter.true_vz]

        true_category = 0
        if not is_fiducial(true_neutrino_vertex):
            true_category = 5

        if chain_filter.passed:
            i_nue += 1


        if is_fiducial(true_neutrino_vertex) and not check_cc0pinp(chain_filter) and chain_filter.ccnc == 0:
            true_category = 8
            # print(list(chain_filter.nu_daughters_pdg), list(chain_filter.nu_daughters_E))

        energy = chain_filter.nu_energy
        tot_bkg[true_category].Fill(energy, chain_filter.bnbweight)

        if chain_filter.passed:
            chain_nue.GetEntry(i_nue)

            true_neutrino_vertex_sce = [chain_nue.true_vx_sce,
                                        chain_nue.true_vy_sce,
                                        chain_nue.true_vz_sce]

            if chain_nue.numu_passed == 1:
                continue

            if chain_nue.category not in (2, 7, 8):
                continue


            if chain_nue.ccnc == 1:
                continue

            showers_2_tracks_0 = True
            if chain_nue.n_tracks == 0 and chain_nue.n_showers == 1:
                showers_2_tracks_0 = False

            if not showers_2_tracks_0:
                continue

            if not check_reco_fidvol(chain_nue):
                continue

            if check_cc0pinp(chain_filter):
                continue

            fill_kin_branches(chain_nue, chain_filter.bnbweight, var_dict, "nue_cc", False)
            if not pre_cuts(var_dict):
                continue

            reco_category = var_dict["category"][0]

            if not is_fiducial(true_neutrino_vertex_sce):
                reco_category = 5

            passed_bkg[reco_category].Fill(energy, chain_filter.bnbweight)

    effs = {}
    for t in tot_bkg:
        if t == 2 or t == 0:
            continue
        if tot_bkg[t].Integral():
            effs[t] = Efficiency(
                tot_bkg[t], {"passed": passed_bkg[t]}, "energy", total_pot)

    return effs


def efficiency(files_path, eff_variables=[], systematics=False, scale=1, is_1e1p=False):
    nue_cosmic = glob(files_path+"/*10.root")
    chain_nue = ROOT.TChain("robertoana/pandoratree")
    chain_filter = ROOT.TChain("nueFilter/filtertree")
    chain_pot = ROOT.TChain("nueFilter/pot")

    for filename in nue_cosmic:
        chain_filter.Add(filename)
        chain_nue.Add(filename)
        chain_pot.Add(filename)

    entries = int(chain_filter.GetEntries() / scale)
    pot_entries = int(chain_pot.GetEntries() / scale)

    total_pot = 0
    for i in range(pot_entries):
        chain_pot.GetEntry(i)
        total_pot += chain_pot.pot

    total_pot /= 1.028
    categories = ["passed", "quality cuts", "CC #nu_{#mu} selected", "not contained",
                  "cosmic selected", "1 shower", "no showers", "no flash", "no data products"]

    h_tot = {}
    h_tot_sys = {}

    for v in eff_variables:
        h_tot[v] = ROOT.TH1F("h_tot_%s" % v,
                             labels[v],
                             binning[v][0],
                             binning[v][1],
                             binning[v][2])
        h_tot_uni = []
        for u in range(N_UNI):
            h_tot_uni.append(ROOT.TH1F("h_tot_%s_%i" % (v, u),
                                       labels[v],
                                       binning[v][0],
                                       binning[v][1],
                                       binning[v][2]))
        h_tot_sys[v] = h_tot_uni

    h_selected = {}
    h_passed_sys = {}
    for v in eff_variables:
        h_sel = {}
        for c in categories:
            h_sel[c] = ROOT.TH1F("%s" % c,
                                 labels[v],
                                 binning[v][0],
                                 binning[v][1],
                                 binning[v][2])
        h_selected[v] = h_sel

        h_passed_uni = []
        for u in range(N_UNI):
            h_passed_uni.append(ROOT.TH1F("h_tot_%s_%i" % (v, u),
                                          labels[v],
                                          binning[v][0],
                                          binning[v][1],
                                          binning[v][2]))
        h_passed_sys[v] = h_passed_uni

    i_nue = -1
    eff_vars = {}
    var_dict = dict(variables + spectators)
    ROOT.TMVA.Tools.Instance()
    reader = ROOT.TMVA.Reader(":".join([
        "!V",
        "Silent",
        "Color"]))
    load_bdt(reader)

    for i_evt in tqdm_notebook(range(entries)):
        chain_filter.GetEntry(i_evt)

        if is_1e1p:
            eNp_interaction = check_1e1p(chain_filter)
        else:
            eNp_interaction = check_cc0pinp(chain_filter)

        true_neutrino_vertex = [chain_filter.true_vx_sce,
                                chain_filter.true_vy_sce,
                                chain_filter.true_vz_sce]

        if chain_filter.passed:
            i_nue += 1

        if not (is_fiducial(true_neutrino_vertex) and eNp_interaction):
            continue
        proton_energy, electron_energy = eNp_interaction[0], eNp_interaction[1]
        eff_vars["proton_energy"] = proton_energy
        eff_vars["electron_energy"] = electron_energy
        eff_vars["energy"] = chain_filter.nu_energy
        eff_vars["vx"] = chain_filter.true_vx_sce
        eff_vars["vy"] = chain_filter.true_vy_sce
        eff_vars["vz"] = chain_filter.true_vz_sce

        eff_vars["theta"], eff_vars["phi"] = lepton_angles(chain_filter)

        if systematics:
            flux_weights = [1]*len(chain_filter.flux_weights[0])
            for flux in chain_filter.flux_weights:
                flux_weights = [a*b for a, b in zip(flux_weights, flux)]
            genie_weights = chain_filter.genie_weights[0]
        weight = chain_filter.bnbweight

        for v in eff_variables:
            h_tot[v].Fill(eff_vars[v], weight)

            for u in range(N_UNI):
                h_tot_sys[v][u].Fill(eff_vars[v], weight * flux_weights[u] * genie_weights[u])

            if not chain_filter.passed:
                if chain_filter.selection_result == 2:
                    h_selected[v]["no showers"].Fill(eff_vars[v], weight)
                if chain_filter.selection_result == 4:
                    h_selected[v]["no flash"].Fill(eff_vars[v], weight)
                if chain_filter.selection_result == 5:
                    h_selected[v]["no data products"].Fill(eff_vars[v], weight)
                continue

            chain_nue.GetEntry(i_nue)
            contaminated = chain_nue.cosmic_fraction < 0.5 and chain_nue.category == 7
            selected = (chain_nue.category == 2 or contaminated) and chain_nue.passed == 1
            if not selected:
                h_selected[v]["cosmic selected"].Fill(eff_vars[v], weight)
                continue

            # If there are no tracks we require at least two showers
            showers_2_tracks_0 = True
            if chain_nue.n_tracks == 0 and chain_nue.n_showers == 1:
                showers_2_tracks_0 = False

            if not showers_2_tracks_0:
                h_selected[v]["1 shower"].Fill(eff_vars[v], weight)
                continue

            if chain_nue.numu_passed == 1:
                h_selected[v]["CC #nu_{#mu} selected"].Fill(eff_vars[v], weight)
                continue

            if not check_reco_fidvol(chain_nue):
                h_selected[v]["not contained"].Fill(eff_vars[v], weight)
                continue

            fill_kin_branches(chain_nue, 1, var_dict, "nue", False)

            if not pre_cuts(var_dict):
                h_selected[v]["quality cuts"].Fill(eff_vars[v], weight)
                continue

            if v == "proton_energy":
                there_is_proton = False
                for i_sh in range(chain_nue.n_showers):
                    if chain_nue.matched_showers[i_sh] == 2212:
                        there_is_proton = True
                        break
                for i_tr in range(chain_nue.n_tracks):
                    if there_is_proton:
                        break
                    if chain_nue.matched_tracks[i_tr] == 2212:
                        there_is_proton = True
                        break
                if there_is_proton:
                    h_selected[v]["passed"].Fill(eff_vars[v], weight)
            elif v == "electron_energy":
                there_is_electron = False
                for i_sh in range(chain_nue.n_showers):
                    if chain_nue.matched_showers[i_sh] == 11 and electron_energy > 0.02:
                        there_is_electron = True
                        break
                for i_tr in range(chain_nue.n_tracks):
                    if there_is_electron:
                        break
                    if chain_nue.matched_tracks[i_tr] == 11 and electron_energy > 0.02:
                        there_is_electron = True
                        break
                if there_is_electron:
                    h_selected[v]["passed"].Fill(eff_vars[v], weight)
            else:
                h_selected[v]["passed"].Fill(eff_vars[v], weight)

            for u in range(N_UNI):
                h_passed_sys[v][u].Fill(eff_vars[v], weight * flux_weights[u] * genie_weights[u])

            # bdt_values = {}
            # for bdt_name in BDT_TYPES:
            #     bdt_values[bdt_name] = reader.EvaluateMVA("BDT%s" % bdt_name)

        #     if apply_cuts(bdt_values, var_dict, bdt=False, manual=True):
        #         h_selected["cuts"].Fill(eff_vars[v], weight)
        #     if apply_cuts(bdt_values, var_dict, bdt=True, manual=False):
        #         h_selected["bdt"].Fill(eff_vars[v], weight)

    eff_objects = {}
    for v in eff_variables:
        h_2d = eff_2d(h_tot_sys[v], h_passed_sys[v], v)
        sys_err = eff_sys_error(h_tot_sys[v], h_passed_sys[v])
        total_sys_err = eff_total_sys_error(h_tot_sys[v], h_passed_sys[v])
        eff_objects[v] = Efficiency(h_tot[v], h_selected[v], v, total_pot, total_sys_err, sys_err, h_2d)

    return eff_objects

