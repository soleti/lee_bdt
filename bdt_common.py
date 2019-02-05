from array import array
import numpy as np
import math
import ROOT
from scipy.stats import poisson
import sys
import collections
import os
from tqdm import tqdm

np.set_printoptions(formatter={'float': '{: 0.3f}'.format})

ELECTRON_MASS = 0.51e-3
PROTON_MASS = 0.938
PROTON_THRESHOLD = 0.040
ELECTRON_THRESHOLD = 0.020

MAX_N_TRACKS = 15
MAX_N_SHOWERS = 15
N_UNI = 100
BDT, MANUAL = False, False
# Number to be obtained from Zarko's POT counting tool
total_data_bnb_pot = 4.341e+19


bdt_nc_cut = 0.165 # thesis
bdt_numu_cut = 0.105 # thesis
bdt_cut_cosmic = 0.09 # thesis

# NEW CUTS
bdt_nc_cut = 0.15
bdt_numu_cut = 0.15
bdt_cut_cosmic = 0.2
# NEW MAX CUTS
bdt_nc_cut = 0.19
bdt_numu_cut = 0.17
bdt_cut_cosmic = 0.22

bdt_nc_cut = 0.18
bdt_numu_cut = 0.17
bdt_cut_cosmic = 0.21

# COLLAB MEETING MAX CUTS
# bdt_nc_cut = 0.175
# bdt_numu_cut = 0.19
# bdt_cut_cosmic = 0.116

bdt_cut_neutrino = 0.16
rectangular_cut = 0.03
bdt_types = ["", "cosmic", "numu", "nc", "neutrino"]

miniboone_bins = np.array([0, 0.200,  0.300,  0.375,  0.475,  0.550,  0.675,  0.800,  0.950,  1.100,  1.300,  1.500,  3.000])
miniboone_bins2 = np.array([0, 0.200,  0.300,  0.375,  0.475,  0.550,  0.675,  0.800,  0.950,  1.100,  1.300,  1.500,  1.8])

energy_bins = np.array([0, 0.2, 0.4, 0.6, 0.9, 1.25, 1.9, 3])
energy_bins2 = np.array([0, 0.2, 0.4, 0.6, 0.9, 1.25, 1.9, 2.5])

# bins = array("f", [0, 0.200, 0.350, 0.5, 0.65, 0.8,
#                    0.950, 1.200, 1.500, 3.000])

# bins = np.array([0, 0.200, 0.350, 0.5, 0.65, 0.850,
#                  1.050, 1.300, 1.550, 3.000])

# bins2 = np.array([0, 0.200, 0.350, 0.5, 0.65, 0.850,
#                   1.050, 1.300, 1.550, 1.700])

bins = energy_bins
bins2 = energy_bins2

# bins = np.array([0, 0.2, 0.4, 0.6, 0.9, 1.25, 1.9, 3])

# bins2 = np.array([0, 0.2, 0.4, 0.6, 0.9, 1.25, 1.9, 2.5])

svm_loaded = False

pdg_colors = {
    2147483648: "#ffffff",
    22: '#e6194b',
    11: '#f58231',
    2112: '#ffe119',
    13: '#bfef45',
    211: '#3cb44b',
    2212: '#42d4f4',
    999999: '#4363db',
    99999: '#a9a9a9'
}

pdgs = [
    ("Data off-beam", 2147483648),
    ("#gamma", 22),
    ("e^{#pm}", 11),
    ("n", 2112),
    ("#mu^{#pm}", 13),
    ("#pi^{#pm}", 211),
    ("p", 2212),
    ("Overlay cosmic", 999999),
    ("Other", 99999),
]


pdgs = collections.OrderedDict(pdgs)

inv_pdgs = {v: k for k, v in pdgs.items()}


def load_svm():
    if not svm_loaded:
        reader_svm = ROOT.TMVA.Reader(":".join([
            "!V",
            "!Silent",
            "Color"]))

        bdt_numu = array("f", [0])
        bdt_nc = array("f", [0])
        bdt_cosmic = array("f", [0])
        reader_svm.AddVariable("cosmic", bdt_cosmic)
        reader_svm.AddVariable("numu", bdt_numu)
        reader_svm.AddVariable("nc", bdt_nc)

        reader_svm.BookMVA("SVM", "dataset/weights/TMVAClassification_SVM.weights.xml")
        svm_loaded = True


def fixed_width_histo_2d(h):
    h_clone = h.Clone()
    h_fixed = ROOT.TH2F(h.GetName()+"_fixed",
                        labels["reco_energy"],
                        len(bins) - 1, bins2,
                        len(bins) - 1, bins2)

    for i in range(1, h_clone.GetNbinsX() + 1):
        for j in range(1, h_clone.GetNbinsY() + 1):
            h_fixed.SetBinContent(i, j, h_clone.GetBinContent(i, j))
            h_fixed.SetBinError(i, j, h_clone.GetBinError(i, j))
    h_fixed.GetYaxis().SetTitle(h_fixed.GetXaxis().GetTitle())

    return h_fixed

def fixed_width_histo(h):
    h_clone = h.Clone()
    h_fixed = ROOT.TH1F(h.GetName()+"_fixed", labels["reco_energy"], len(bins) - 1, bins2)

    for i in range(1, h_clone.GetNbinsX() + 1):
        h_fixed.SetBinContent(i, h_clone.GetBinContent(i))
        h_fixed.SetBinError(i, h_clone.GetBinError(i))

    return h_fixed


def is_1eNp(c):
    p = 0
    e = 0
    photons = 0
    pions = 0

    for i_pdg, energy in enumerate(c.nu_daughters_E):
        p_start = [
            c.nu_daughters_vx[i_pdg],
            c.nu_daughters_vy[i_pdg],
            c.nu_daughters_vz[i_pdg]
        ]
        p_end = [
            c.nu_daughters_endx[i_pdg],
            c.nu_daughters_endy[i_pdg],
            c.nu_daughters_endz[i_pdg]
        ]
        if abs(c.nu_daughters_pdg[i_pdg]) == 2212:
            if not (is_fiducial(p_start) and is_fiducial(p_end)):
                p = 0
                break

            if energy - PROTON_MASS > PROTON_THRESHOLD:
                p += 1

        if abs(c.nu_daughters_pdg[i_pdg]) == 11:
            if not is_fiducial(p_start):
                e = 0
                break
            if energy > ELECTRON_THRESHOLD:
                e += 1

        if c.nu_daughters_pdg[i_pdg] == 22:
            photons += 1

        if c.nu_daughters_pdg[i_pdg] == 111 or abs(c.nu_daughters_pdg[i_pdg]) == 211:
            pions += 1

    eNp = e == 1 and photons == 0 and pions == 0 and p >= 1
    true_neutrino_vertex = [c.true_vx_sce,
                            c.true_vy_sce,
                            c.true_vz_sce]

    return eNp and is_fiducial(true_neutrino_vertex)

def choose_shower(root_chain, plane):
    most_energetic_shower = 0
    shower_id = 0
    for ish in range(root_chain.n_showers):
        if root_chain.shower_energy[ish][plane] < 3:
            if root_chain.shower_energy[ish][plane] > most_energetic_shower:
                most_energetic_shower = root_chain.shower_energy[ish][plane]
                shower_id = ish
    return shower_id


def distance(p1, p2):
    return math.sqrt(
        sum([(t - n)**2 for t, n in zip(p1, p2)]))

def gauss_exp(var, par):
    """
    n:par[0]
    mu:par[1]
    sigma:par[2]
    k:par[3]
    """
    n = par[0]
    mu = par[1]
    sigma = max(0.001, par[2])
    k = par[3]
    x = var[0]

    if (x - mu) / sigma >= -k:
        return n * math.exp(-0.5 * ((x - mu) / sigma)**2)
    else:
        return n * math.exp(k**2 / 2 + k * ((x - mu) / sigma))

def fix_binning(histogram, width=0.05):
    for bin_i in range(1, histogram.GetNbinsX() + 1):
        bin_width = histogram.GetBinWidth(bin_i)
        histogram.SetBinError(bin_i, histogram.GetBinError(
            bin_i) / (bin_width / width))
        histogram.SetBinContent(bin_i, histogram.GetBinContent(
            bin_i) / (bin_width / width))


def printProgressBar(iteration, total, prefix='', suffix='', decimals=1, length=100, fill='█'):
    """
    Call in a loop to create terminal progress bar
    @params:
        iteration   - Required  : current iteration (Int)
        total       - Required  : total iterations (Int)
        prefix      - Optional  : prefix string (Str)
        suffix      - Optional  : suffix string (Str)
        decimals    - Optional  : positive number of decimals in percent complete (Int)
        length      - Optional  : character length of bar (Int)
        fill        - Optional  : bar fill character (Str)
    """
    percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
    filledLength = int(length * iteration // total)
    bar = fill * filledLength + '-' * (length - filledLength - 1)
    print('\r%s |%s| %s%% %s' % (prefix, bar, percent, suffix), end='\r')
    # Print New Line on Complete
    if iteration == total:
        print()


def bin_size(name):
    return abs(binning[name][2] - binning[name][1]) / binning[name][0]


def FC_histo(h_bkg, bins=array("f", [])):
    if len(bins) > 0:
        h_limit = ROOT.TH1F("h_limit", "", len(bins) - 1, bins)
    else:
        h_limit = ROOT.TH1F("h_limit", "",
                            h_bkg.GetNbinsX(),
                            h_bkg.GetXaxis().GetXmin(),
                            h_bkg.GetXaxis().GetXmax())

    for i in range(1, h_limit.GetNbinsX() + 1):
        h_limit.SetBinContent(i, h_bkg.GetBinContent(i)+upper_limit_FC(h_bkg.GetBinContent(i)))
    return h_limit


def upper_limit_FC(bkg, cf=0.9999994):
    FC = ROOT.TFeldmanCousins(cf)
    FC.SetMuMax(200)

    sigma = math.sqrt(bkg)
    nmin = max(0, int(bkg - 5. * sigma))   # Use 0 if nmin<0
    nmax = max(20, int(bkg + 5. * sigma)+1) # Use at least 20 for low means
    po = poisson(bkg)
    UL = 0.
    for i in range(nmin, nmax):
        pmf = po.pmf(i)
        ul = FC.CalculateUpperLimit(i, bkg)
        #print "i=%i, Po(i)=%f, U(i,b)=%f" % (i, pmf, ul)
        UL += po.pmf(i) * ul
    return UL

colors = [ROOT.TColor.GetColor("#000000"), ROOT.TColor.GetColor("#e7623d"),
          ROOT.TColor.GetColor("#62b570"), ROOT.TColor.GetColor("#8779b1"),
          ROOT.TColor.GetColor("#609ec6"), ROOT.TColor.GetColor("#c46d27"),
          ROOT.TColor.GetColor("#ffffff"), ROOT.TColor.GetColor("#e7623d"),
          ROOT.TColor.GetColor("#1e7a2e")]

def pre_cuts(var_dict, numu=False):
    tr_id = int(var_dict["track_id"][0])
    sh_id = int(var_dict["shower_id"][0])

    if numu:
        numu = int(var_dict["numu_score"][0]) == 1
    else:
        numu = int(var_dict["numu_score"][0]) == 0
    hits = var_dict["track_hits"][0] > 6 and \
           var_dict["shower_hits"][0] > 5 and \
           var_dict["total_hits_y"][0] > 0 and \
           var_dict["total_hits_u"][0] > 0 and \
           var_dict["total_hits_v"][0] > 0
    shower_track_energy = var_dict["total_shower_energy_cali"][0] > 0.01 and \
                          var_dict["total_track_energy_length"][0] > 0 and \
                          var_dict["shower_energy"][0] > 0.01
    track_pid = var_dict["track_pidchipr"][tr_id] < 80
    ene = True #
    category = var_dict["category"][0]

    if category == 6 and 0.6 < var_dict["reco_energy"][0] < 1.25:
        ene = ene and var_dict["n_showers"][0] < 3
    if not var_dict["true_nu_is_fidvol"][0] and category != 0 and category != 6 and category != 1 and category != 7:
        category = 5
    shower_dedx = var_dict["shower_dedx"][sh_id] < 3.5
    for i_sh in range(int(var_dict["n_showers"][0])):
        shower_dedx = shower_dedx and var_dict["shower_dedx"][i_sh] < 4.5
    return numu and hits and shower_track_energy and var_dict["reco_energy"][0] < 3# and shower_dedx and ene

def is_active(point):
    ok_y = y_start < point[1] < y_end
    ok_x = x_start < point[0] < x_end
    ok_z = z_start < point[2] < z_end
    return ok_y and ok_x and ok_z

def is_fiducial(point):
    ok_y = point[1] > y_start + 15 and point[1] < y_end - 15
    ok_x = point[0] > x_start + 10 and point[0] < x_end - 10
    ok_z = point[2] > z_start + 10 and point[2] < z_end - 40
    return ok_y and ok_x and ok_z


def load_bdt(reader):
    for name, var in variables:
        reader.AddVariable(name, var)

    for name, var in spectators:
        reader.AddSpectator(name, var)

    for type in bdt_types:
        reader.BookMVA("BDT%s" % type,
                       "dataset/weights/TMVAClassification_BDT%s.weights.xml" % type)


def load_variables(chain):
    for name, var in variables:
        chain.SetBranchAddress(name, var)

    for name, var in spectators:
        chain.SetBranchAddress(name, var)

    return dict(variables + spectators)


def bdt_cut(bdt_dict, bdt_values=[bdt_nc_cut, bdt_numu_cut, bdt_cut_cosmic]):
    # SVM_response = reader_svm.EvaluateMVA("SVM")
    # return SVM_response > 0.85
    return bdt_dict["nc"] > bdt_values[0] and bdt_dict["numu"] > bdt_values[1] and bdt_dict["cosmic"] > bdt_values[2]

def apply_cuts(bdt_dict, var_dict, bdt=BDT, manual=MANUAL, mode="nue"):
    if bdt:
        # if int(var_dict["category"][0]) == 6:
        #     apply_bdt = bdt_cut(bdt_dict, [bdt_nc_cut, bdt_numu_cut, 0.109])
        # else:
        apply_bdt = bdt_cut(bdt_dict)
    else:
        apply_bdt = True

    if manual:
        apply_manual = manual_cuts(var_dict, mode)
    else:
        apply_manual = True

    return apply_bdt and apply_manual


def find_interaction(dictionary, interaction):
    for name, id_int in dictionary.items():
        if id_int == interaction:
            return name


def manual_cuts(var_dict, mode="nue"):
    tr_id = int(var_dict["track_id"][0])
    sh_id = int(var_dict["shower_id"][0])
    if var_dict["no_tracks"][0] == 1:
        tr_id = 0
        shower_distance = var_dict["shower_distance"][sh_id] < 2
        track_distance = var_dict["track_distance"][tr_id] < 2
    ratio = var_dict["hits_ratio"][0] > 0.55
    ratio_photon = var_dict["hits_ratio"][0] > 0.7
    track_res = var_dict["track_res_std"][tr_id] < 2
    shower_dedx_photon = True
    for i_sh in range(int(var_dict["n_showers"][0])):
        shower_distance = var_dict["shower_distance"][i_sh] < 3.9
        if var_dict["category"][0] == 6 and 0.5 < var_dict["reco_energy"][0] < 0.85:
            shower_dedx_photon = shower_dedx_photon and 2.91/3.85e-5 < var_dict["shower_dqdx"][i_sh] < 5/3.85e-5
        else:
            shower_dedx_photon = shower_dedx_photon and 3.1/3.85e-5 < var_dict["shower_dqdx"][i_sh] < 5/3.85e-5
    track_distance = var_dict["track_distance"][tr_id] < 5
    shower_open_angle = 1 < var_dict["shower_open_angle"][sh_id] < 20
    track_length = var_dict["track_length"][tr_id] < 80
    track_length_numu = var_dict["track_length"][tr_id] > 20
    total_hits_y = var_dict["total_hits_y"][0] > 100
    total_hits_y_numu = var_dict["total_hits_y"][0] > 150
    shower_energy = var_dict["shower_energy"][sh_id] > 0.05
    shower_hits = var_dict["shower_hits"][0] > 100
    total_shower_energy = var_dict["total_shower_energy_cali"][0] > 0.1
    total_shower_energy_numu = var_dict["total_shower_energy_cali"][0] > 0.05
    shower_dedx = 1/3.85e-5 < var_dict["shower_dqdx"][sh_id] < 3.2/3.85e-5
    for i_sh in range(int(var_dict["n_showers"][0])):
        shower_dedx = shower_dedx and var_dict["shower_dqdx"][i_sh] < 3.5/3.85e-5
    track_shower_angle = -0.95 < var_dict["track_shower_angle"][tr_id]
    track_hits = var_dict["track_hits"][0] > 100
    track_proton_chi2 = True
    track_mu_chi2 = False
    for i_tr in range(int(var_dict["n_tracks"][0])):
        track_proton_chi2 = track_proton_chi2 and 0 < var_dict["track_pidchipr"][tr_id] < 80
        if var_dict["track_pidchipr"][i_tr] > 40 and var_dict["track_pidchipr"][i_tr] < 220:
            track_mu_chi2 = True
    track_proton_chi2_numu = track_mu_chi2
    shower_pi_chi2 = not (0 < var_dict["shower_pidchipi"][sh_id] < 12)
    pt_ratio = var_dict["pt_ratio"][0] < 0.85

    # total_hits_y = var_dict["total_hits_y"][0] > 85 # remove for thesis
    # shower_open_angle = 1 < var_dict["shower_open_angle"][sh_id] < 30

    # collabmeeting_cuts = [track_res, ratio, shower_distance, pt_ratio, shower_open_angle,
    #                       track_length, track_distance, total_hits_y, track_shower_angle,
    #                       shower_energy, shower_dedx, track_proton_chi2, shower_pi_chi2]


    # total hits, shower energy, ratio, dedx, track distance, shower distance, proton chi2, track shower angle, track length

    collabmeeting_cuts = [track_res, ratio, shower_distance, shower_open_angle,
                          track_length, track_distance, total_hits_y, track_shower_angle,
                          shower_energy, shower_dedx, track_proton_chi2, shower_pi_chi2]

    # collabmeeting_cuts = [total_hits_y, track_res, shower_energy, ratio, shower_dedx, shower_open_angle,
    #                       track_distance, shower_distance, track_proton_chi2, shower_pi_chi2, track_shower_angle, track_length]
    passed_collab = len(collabmeeting_cuts) == sum(collabmeeting_cuts)

    photon_cuts = [ratio, track_length, track_distance, total_hits_y, track_shower_angle,
                   shower_energy, shower_dedx_photon, track_proton_chi2]

    passed_photon = len(photon_cuts) == sum(photon_cuts)

    numu_cuts = [shower_distance, track_length_numu, shower_energy,
                 track_shower_angle, shower_dedx, track_proton_chi2_numu]

    passed_numu = len(numu_cuts) == sum(numu_cuts)

    if mode == "numu":
        return passed_numu
    elif mode == "nc":
        return passed_photon

    return passed_collab


def sigma_calc_matrix(h_signal, h_background, scale_factor=1, systematics=False, mode="selection", error_scale=1):
    #it is just, Δχ2 = (number of events signal in Energy bins in a 1D matrix)
    #(2D Matrix - (statistical uncertainty)^2 in a the diagonal of the matrix)^-1
    #(number of events signal in Energy bins in a 1D matrix)^Transpose

    sig_array = h_signal * scale_factor
    bkg_array = (h_background) * scale_factor
    if mode == "selection":
        folder = ""
    elif mode == "bdt":
        folder = "_bdt"
    elif mode == "cuts":
        folder = "_cuts"
    nbins = len(sig_array)

    emtx = np.zeros((nbins, nbins))
    for x in range(nbins):
        emtx[x][x] = bkg_array[x]

    if systematics:
        fname_flux = "plots%s/sys/h_cov_reco_energy_flux.root" % folder
        fname_genie = "plots%s/sys/h_cov_reco_energy_genie.root" % folder
        fname_det = "plots/sys/h_reco_energy_det_sys.root"

        if os.path.isfile(fname_genie):
            f_genie = ROOT.TFile(fname_genie)
            h_genie = f_genie.Get("h_cov_reco_energy")
            fill_cov_matrix(emtx, h_genie, scale_factor * error_scale)
            f_genie.Close()
        else:
            print("GENIE covariance matrix %s not available" % fname_genie)

        if os.path.isfile(fname_flux):
            f_flux = ROOT.TFile(fname_flux)
            h_flux = f_flux.Get("h_cov_reco_energy")
            fill_cov_matrix(emtx, h_flux, scale_factor * error_scale)
            f_flux.Close()
        else:
            print("Flux covariance matrix %s not available" % fname_flux)

        if os.path.isfile(fname_det):
            f_det = ROOT.TFile(fname_det)
            h_det = f_det.Get("h_frac_reco_energy")
            fill_cov_matrix_det(emtx, h_det, bkg_array, scale_factor * error_scale)
            f_det.Close()
        else:
            print("Det. sys. covariance matricx %s not available" % fname_det)

    emtxinv = np.linalg.inv(emtx)
    chisq = float(sig_array.dot(emtxinv).dot(sig_array.T))

    #print "Sqrt of that (==sigma?) is ",np.sqrt(chisq)
    return np.sqrt(chisq)

def sigma_calc_likelihood_sys(sig, bkg, delta_b):
    s = sig
    b = bkg
    db = delta_b

    ln_1 = math.log((s+b)*(b+db**2)/(b**2+(s+b)*db**2))
    ln_2 = math.log(1+db**2*s/(b*b+b*db**2))

    sigma = math.sqrt(2*((s+b)*ln_1-b**2/db**2*ln_2))

    return sigma

def sigma_calc_likelihood(h_signal, h_background, scale_factor=1):
    sig_array = h_signal * scale_factor
    bkg_array = (h_background) * scale_factor

    sum_sigma = 0
    for sig, bkg, in zip(sig_array, bkg_array):
        n0 = sig+bkg
        sigma2 = 2 * n0 * math.log(1+sig/bkg) - 2*sig
        sum_sigma += sigma2
    return math.sqrt(sum_sigma)

def fill_cov_matrix(matrix, histo, scale_factor=1):
    for x_bin in range(1, histo.GetNbinsX()+1):
        for y_bin in range(1, histo.GetNbinsY()+1):
            matrix[x_bin - 1][y_bin - 1] += histo.GetBinContent(x_bin, y_bin)*(scale_factor**2)

def fill_cov_matrix_det(matrix, histo, bkg_array, scale_factor=1):
    for x_bin in range(1, histo.GetNbinsX()+1):
        for y_bin in range(1, histo.GetNbinsY()+1):
            matrix[x_bin - 1][y_bin - 1] += histo.GetBinContent(x_bin, y_bin)*bkg_array[x_bin-1]*bkg_array[y_bin-1]/2

def save_histo_sbnfit(histo, name, scale_factor=1):
    h_tosave = ROOT.TH1D(name, "", len(bins) - 1, bins)
    for i_bin in range(1, histo.GetNbinsX()+1):
        h_tosave.SetBinContent(i_bin, histo.GetBinContent(i_bin)*scale_factor)
    return h_tosave



total_pot = total_data_bnb_pot

description2 = ["Cosmic in-time",
                "Cosmic",
                "#nu_{e} CC0#pi-Np",
                "Beam Intrinsic #nu_{#mu}",
                "Beam Intrinsic NC",
                "Outside fid. vol.",
                "Data",
                "Cosmic contaminated",
                "#nu_{e} CC",
                "Low-energy excess"]

description = ["Data beam-off",
               "Beam Intrinsic #nu_{#mu}",
               "Beam Intrinsic NC",
               "Outside fid. vol.",
               "Cosmic contaminated",
               "Cosmic",
               "Data beam-on",
               "#nu_{e} CC",
               "#nu_{e} CC0#pi-Np",
               "Low-energy excess"]

interactions = {
    "kUnknownInteraction": -1,
    "QE": 0,
    "Resonant": 1,
    "DIS": 2,
    "Coherent": 3,
    "MEC": 4,
    "kElectronScattering": 5,
    "kIMDAnnihilation": 6,
    "kInverseBetaDecay": 7,
    "kGlashowResonance": 8,
    "kAMNuGamma": 9,
    "kCohElastic": 10,
    "kDiffractive": 11,
    "kEM": 12,
    "kWeakMix": 13,
    "kNuanceOffset": 1000,
    "kCCQE": 1000 + 1,
    "kNCQE": 1000 + 2,
    "kResCCNuProtonPiPlus": 1000 + 3,
    "kResCCNuNeutronPi0": 1000 + 4,
    "kResCCNuNeutronPiPlus": 1000 + 5,
    "kResNCNuProtonPi0": 1000 + 6,
    "kResNCNuProtonPiPlus": 1000 + 7,
    "kResNCNuNeutronPi0": 1000 + 8,
    "kResNCNuNeutronPiMinus": 1000 + 9,
    "kResCCNuBarNeutronPiMinus": 1000 + 10,
    "kResCCNuBarProtonPi0": 1000 + 11,
    "kResCCNuBarProtonPiMinus": 1000 + 12,
    "kResNCNuBarProtonPi0": 1000 + 13,
    "kResNCNuBarProtonPiPlus": 1000 + 14,
    "kResNCNuBarNeutronPi0": 1000 + 15,
    "kResNCNuBarNeutronPiMinus": 1000 + 16,
    "kResCCNuDeltaPlusPiPlus": 1000 + 17,
    "kResCCNuDelta2PlusPiMinus": 1000 + 21,
    "kResCCNuBarDelta0PiMinus": 1000 + 28,
    "kResCCNuBarDeltaMinusPiPlus": 1000 + 32,
    "kResCCNuProtonRhoPlus": 1000 + 39,
    "kResCCNuNeutronRhoPlus": 1000 + 41,
    "kResCCNuBarNeutronRhoMinus": 1000 + 46,
    "kResCCNuBarNeutronRho0": 1000 + 48,
    "kResCCNuSigmaPluskaonPlus": 1000 + 53,
    "kResCCNuSigmaPluskaon0": 1000 + 55,
    "kResCCNuBarSigmaMinuskaon0": 1000 + 60,
    "kResCCNuBarSigma0kaon0": 1000 + 62,
    "kResCCNuProtonEta": 1000 + 67,
    "kResCCNuBarNeutronEta": 1000 + 70,
    "kResCCNukaonPlusLambda0": 1000 + 73,
    "kResCCNuBarkaon0Lambda0": 1000 + 76,
    "kResCCNuProtonPiPlusPiMinus": 1000 + 79,
    "kResCCNuProtonPi0Pi0": 1000 + 80,
    "kResCCNuBarNeutronPiPlusPiMinus": 1000 + 85,
    "kResCCNuBarNeutronPi0Pi0": 1000 + 86,
    "kResCCNuBarProtonPi0Pi0": 1000 + 90,
    "kCCDIS": 1000 + 91,
    "kNCDIS": 1000 + 92,
    "kUnUsed1": 1000 + 93,
    "kUnUsed2": 1000 + 94,
    "k1eNpHyperon": 1000 + 95,
    "kNCCOH": 1000 + 96,
    "kCCCOH": 1000 + 97,
    "kNuElectronElastic": 1000 + 98,
    "kInverseMuDecay": 1000 + 99
}

inv_interactions = {v: k for k, v in interactions.items()}

x_start = 0
x_end = 256.35
y_start = -116.5
y_end = 116.5
z_start = 0
z_end = 1036.8
pi0_mass = array("f", [0])

genie_weights = array("f", N_UNI * [1.])
flux_weights = array("f", N_UNI * [1.])

shower_angle = array("f", MAX_N_SHOWERS * [-sys.float_info.max])
track_angle = array("f", [0])
shower_pdg = array("f", MAX_N_SHOWERS * [-sys.float_info.max])
track_pdg = array("f", MAX_N_TRACKS * [-sys.float_info.max])
track_energy_length = array("f", MAX_N_TRACKS * [-sys.float_info.max])
theta_ratio = array("f", [0])
track_length = array("f", MAX_N_TRACKS * [-sys.float_info.max])
track_theta = array("f", MAX_N_TRACKS * [-sys.float_info.max])
track_phi = array("f", MAX_N_TRACKS * [-sys.float_info.max])
track_start_y = array("f", MAX_N_TRACKS * [-sys.float_info.max])
track_end_y = array("f", MAX_N_TRACKS * [-sys.float_info.max])
track_start_x = array("f", MAX_N_TRACKS * [-sys.float_info.max])
track_end_x = array("f", MAX_N_TRACKS * [-sys.float_info.max])
track_start_z = array("f", MAX_N_TRACKS * [-sys.float_info.max])
track_end_z = array("f", MAX_N_TRACKS * [-sys.float_info.max])
track_res_mean = array("f", MAX_N_TRACKS * [-sys.float_info.max])
track_res_std = array("f", MAX_N_TRACKS * [-sys.float_info.max])
track_pca = array("f", MAX_N_TRACKS * [-sys.float_info.max])
track_pidchipr = array("f", MAX_N_TRACKS * [-sys.float_info.max])
track_pidchimu = array("f", MAX_N_TRACKS * [-sys.float_info.max])

track_likelihood = array("f", MAX_N_TRACKS * [-sys.float_info.max])
track_mip_likelihood = array("f", MAX_N_TRACKS * [-sys.float_info.max])
track_p_likelihood = array("f", MAX_N_TRACKS * [-sys.float_info.max])
track_pidchi = array("f", MAX_N_TRACKS * [-sys.float_info.max])
track_pida = array("f", MAX_N_TRACKS * [-sys.float_info.max])
track_shower_angle = array("f", MAX_N_TRACKS * [-sys.float_info.max])
track_distance = array("f", MAX_N_TRACKS * [-sys.float_info.max])
track_dqdx = array("f", MAX_N_TRACKS * [-sys.float_info.max])

shower_theta = array("f", MAX_N_SHOWERS * [-sys.float_info.max])
shower_phi = array("f", MAX_N_SHOWERS * [-sys.float_info.max])
shower_score = array("f", MAX_N_SHOWERS * [-sys.float_info.max])
shower_start_y = array("f", MAX_N_SHOWERS * [-sys.float_info.max])
shower_start_x = array("f", MAX_N_SHOWERS * [-sys.float_info.max])
shower_start_z = array("f", MAX_N_SHOWERS * [-sys.float_info.max])
shower_energy = array("f", MAX_N_SHOWERS * [-sys.float_info.max])
shower_res_mean = array("f", MAX_N_SHOWERS * [-sys.float_info.max])
shower_res_std = array("f", MAX_N_SHOWERS * [-sys.float_info.max])
shower_dedx = array("f", MAX_N_SHOWERS * [-sys.float_info.max])
shower_dqdx = array("f", MAX_N_SHOWERS * [-sys.float_info.max])
shower_dedx_cali = array("f", MAX_N_SHOWERS * [-sys.float_info.max])
shower_dedx_u = array("f", MAX_N_SHOWERS * [-sys.float_info.max])
shower_dedx_v = array("f", MAX_N_SHOWERS * [-sys.float_info.max])
shower_pca = array("f", MAX_N_SHOWERS * [-sys.float_info.max])
shower_open_angle = array("f", MAX_N_SHOWERS * [-sys.float_info.max])
shower_distance = array("f", MAX_N_SHOWERS * [-sys.float_info.max])
shower_pidchimu = array("f", MAX_N_SHOWERS * [-sys.float_info.max])
shower_pidchipr = array("f", MAX_N_SHOWERS * [-sys.float_info.max])
shower_pidchipi = array("f", MAX_N_SHOWERS * [-sys.float_info.max])


b_shower_id = array("f", [0])
b_track_id = array("f", [0])

pt = array("f", [0])
p = array("f", [0])
pt_ratio = array("f", [0])
n_tracks = array("f", [0])
n_showers = array("f", [0])
n_tracks_before = array("f", [0])
n_showers_before = array("f", [0])
reco_energy = array("f", [0])
event_weight = array("f", [0])
category = array("f", [0])
event = array("f", [0])
run = array("f", [0])
subrun = array("f", [0])
interaction_type = array("f", [0])
is_signal = array("f", [0])

numu_score = array("f", [0])
total_shower_energy = array("f", [0])
total_shower_energy_cali = array("f", [0])
total_track_energy = array("f", [0])
shower_hits = array("f", [0])
hits_ratio = array("f", [0])
track_hits = array("f", [0])
total_track_energy_length = array("f", [0])
n_objects = array("f", [0])
total_hits = array("f", [0])
total_hits_u = array("f", [0])
total_hits_v = array("f", [0])
total_hits_y = array("f", [0])
no_tracks = array("f", [0])
shower_length = array("f", [0])
nu_E = array("f", [0])
E_dep = array("f", [0])
true_nu_is_fidvol = array("f", [0])

shower_dedx_bdt = array("f", [0])

spectators = [
    ("track_pdg", track_pdg),
    ("shower_pdg", shower_pdg),
    ("shower_dedx_bdt", shower_dedx_bdt),
    ("true_nu_is_fidvol", true_nu_is_fidvol),
    ("total_hits_u", total_hits_u),
    ("total_hits_v", total_hits_v),
    ("track_mip_likelihood", track_mip_likelihood),
    ("track_p_likelihood", track_p_likelihood),
    ("genie_weights", genie_weights),
    ("flux_weights", flux_weights),
    # ("track_pidchi", track_pidchi),
    # ("track_pida", track_pida),
    ("nu_E", nu_E),
    ("E_dep", nu_E),
    ("category", category),
    ("event_weight", event_weight),
    ("event", event),
    ("run", run),
    ("subrun", subrun),
    ("interaction_type", interaction_type),
    ("is_signal", is_signal),
    ("reco_energy", reco_energy),
    ("total_track_energy", total_track_energy),
    ("total_hits_y", total_hits_y),
    ("track_start_x", track_start_x),
    ("track_end_x", track_end_x),
    ("track_start_z", track_start_z),
    ("track_end_z", track_end_z),
    ("shower_start_z", shower_start_z),
    ("track_end_y", track_end_y),
    ("shower_hits", shower_hits),
    ("track_hits", track_hits),
    ("shower_start_x", shower_start_x),
    # ("shower_length", shower_length),
    ("n_objects", n_objects),
    ("pt", pt),
    ("p", p),

    ("track_pca", track_pca),
    ("shower_pca", shower_pca),
    ("track_res_std", track_res_std),
    ("shower_res_mean", shower_res_mean),
    ("shower_res_std", shower_res_std),
    ("total_hits", total_hits),
    ("total_shower_energy", total_shower_energy),
    ("track_energy_length", track_energy_length),
    ("total_track_energy_length", total_track_energy_length),
    ("numu_score", numu_score),
    ("shower_dedx", shower_dedx),
    ("shower_dedx_cali", shower_dedx_cali),
    ("shower_dedx_u", shower_dedx_u),
    ("shower_dedx_v", shower_dedx_v),
    # ("track_angle", track_angle),
    ("shower_id", b_shower_id),
    ("track_id", b_track_id),
    ("n_showers", n_showers),
    ("n_tracks", n_tracks),
    ("n_showers_before", n_showers_before),
    ("n_tracks_before", n_tracks_before),
    ("total_shower_energy_cali", total_shower_energy_cali),
    ("shower_energy", shower_energy),
    ("no_tracks", no_tracks)
    # ("track_pidchimu", track_pidchimu),
]

variables = [
    ("pt_ratio", pt_ratio),
    ("theta_ratio", theta_ratio),
    ("shower_pidchimu", shower_pidchimu),
    ("shower_pidchipr", shower_pidchipr),
    ("shower_pidchipi", shower_pidchipi),
    ("shower_angle", shower_angle),
    ("track_dqdx", track_dqdx),
    ("shower_dqdx", shower_dqdx),
    ("track_pidchipr", track_pidchipr),
    ("track_start_y", track_start_y),
    ("shower_start_y", shower_start_y),
    ("track_length", track_length),
    ("shower_open_angle", shower_open_angle),
    ("hits_ratio", hits_ratio),
    ("track_distance", track_distance),
    ("shower_distance", shower_distance),
    ("track_res_mean", track_res_mean),
    ("track_theta", track_theta),
    ("track_shower_angle", track_shower_angle),
    ("shower_theta", shower_theta),
    ("shower_phi", shower_phi),
    ("track_phi", track_phi),
    ("track_likelihood", track_likelihood),
]

binning = {
    "track_pdg": [100, 0, 100],
    "shower_pdg": [100, 0, 100],
    "shower_id": [15, 0, 15],
    "track_id": [15, 0, 15],
    "genie_weights": [50, 0, 1],
    "flux_weights": [50, 0, 1],

    "shower_dedx_bdt": [50, -1, 1],
    "shower_angle": [20, 0, 180],
    "track_angle": [20, 0, 180],
    "pi0_mass": [20, 0, 0.4],
    "shower_track_d": [30, 0, 15],
    "true_nu_is_fidvol": [2, 0, 2],
    "total_hits_u": [1, 0, 1],
    "total_hits_v": [1, 0, 1],
    "track_pidchipr": [30, 0, 300],
    "track_pidchimu": [40, 0, 80],

    "shower_pidchipr": [40, 0, 300],
    "shower_pidchipi": [20, 0, 60],
    "shower_pidchimu": [20, 0, 60],

    "track_likelihood": [40, -10, 10],
    "track_mip_likelihood": [40, 0, 1.1],
    "track_p_likelihood": [40, 0, 1.1],

    "track_pidchi": [20, 0, 40],
    "track_pida": [20, 0, 25],
    "track_res_mean": [20, -2, 2],
    "track_res_std": [20, 0, 5],
    "shower_res_mean": [20, -2, 2],
    "shower_res_std": [20, 0, 5],
    "no_tracks": [2, 0, 2],
    "n_objects": [8, 2, 10],
    "n_tracks": [5, 1, 6],
    "n_showers": [5, 1, 6],
    "n_tracks_before": [5, 1, 6],
    "n_showers_before": [5, 1, 6],
    "track_theta": [18, 0, 180],
    "track_phi": [18, -180, 180],
    "shower_theta": [9, 0, 180],
    "shower_phi": [9, -180, 180],
    "shower_distance": [20, 0, 15],
    "track_distance": [20, 0, 15],
    "track_shower_angle": [20, -1, 1],
    "track_start_y": [20, y_start, y_end],
    "track_start_z": [20, z_start, z_end],
    "track_start_x": [20, x_start, x_end],
    "track_end_y": [20, y_start, y_end],
    "track_end_z": [20, z_start, z_end],
    "track_end_x": [20, x_start, x_end],
    "shower_start_y": [10, y_start, y_end],
    "shower_start_z": [10, z_start, z_end],
    "shower_start_x": [10, x_start, x_end],
    "shower_end_y": [20, y_start, y_end],
    "shower_end_z": [20, z_start, z_end],
    "shower_end_x": [20, x_start, x_end],
    "track_length": [20, 0, 300],
    "shower_length": [20, 0, 150],
    "proton_score": [10, 0, 1],
    "shower_energy": [20, 0, 0.5],
    "shower_energy_y": [20, 0, 0.5],
    "total_hits": [30, 0, 300],
    "total_hits_u": [20, 0, 1000],
    "total_hits_v": [20, 0, 1000],
    "total_hits_y": [20, 0, 1000],
    "pt": [20, 0, 1],
    "pt_ratio": [20, 0, 1],

    "p": [20, 0, 1],

    "reco_energy": [10, 0, 2],
    "shower_open_angle": [23, 0, 46],
    "shower_dedx": [10, 0.3, 6],
    "shower_dedx_cali": [10, 0, 6],
    "shower_dqdx": [18, 0, 155000],

    "shower_dedx_u": [20, 0.3, 6],
    "shower_dedx_v": [20, 0.3, 6],

    "shower_dedx_merged": [20, 0, 5],
    "numu_score": [20, 0, 20],
    "category": [11, 0, 11],
    "event_weight": [20, 0, 100],
    "event": [20, 0, 10000],
    "run": [20, 0, 10],
    "subrun": [20, 0, 1000],
    "interaction_type": [5, 0, 5],
    "is_signal": [1, -1.5, 1.5],
    "shower_pca": [20, 0.9, 1],
    "track_pca": [20, 0.99, 1],
    "total_shower_energy": [20, 0, 0.5],
    "total_shower_energy_cali": [20, 0, 0.5],

    "total_shower_energy_y": [20, 0, 0.5],
    "total_track_energy": [20, 0, 1],
    "shower_hits": [20, 0, 200],
    "hits_ratio": [20, 0, 1],
    "theta_ratio": [40, 0, 5],
    "energy": [20, 0, 3],
    "vx": [30, x_start, x_end],
    "vy": [30, y_start, y_end],
    "vz": [30, z_start, z_end],
    "theta": [18, 0, 180],
    "phi": [18, -180, 180],
    "nu_E": [40, 0, 4],
    "E_dep": [30, 0, 3],
    "track_hits": [20, 0, 400],
    "track_dqdx": [20, 0, 1200],
    "total_track_energy_length": [12, 0, 1.2],
    "track_energy_length": [12, 0, 1.2]

}

if MANUAL or BDT:
    binning["shower_theta"] = [9, 0, 180]
    binning["shower_phi"] = [9, -180, 180]
    binning["shower_start_y"] = [10, y_start, y_end]
    binning["shower_start_z"] = [10, z_start, z_end]
    binning["shower_start_x"] = [10, x_start, x_end]


labels = {
    "shower_id": "shower_id",
    "track_id": "track_id",
    "track_pdg": "track pdg",
    "shower_pdg": "shower pdg",
    "shower_dedx_bdt": ";dE/dx BDT;",
    "shower_track_d": ";shower track d;",
    "shower_angle": ";Secondary shower/main shower angle [#circ];",
    "track_angle": ";track_angle;",
    "pi0_mass": ";pi0 mass;",
    "true_nu_is_fidvol": ";true_nu_is_fidvol;",
    "total_hits_u": ";total hits u;",
    "total_hits_v": "; total hits u;",
    "no_tracks": ";No tracks",
    "flux_weights": "; flux;",
    "genie_weights": "; genie;",

    "no_tracks": ";No tracks",
    "track_pidchipr": ";Track proton PID #chi^{2};N. Entries / %.1f" % bin_size("track_pidchipr"),
    "track_pidchimu": ";Track #mu PID #chi^{2};N. Entries / %.1f" % bin_size("track_pidchimu"),

    "shower_pidchipr": ";Shower proton PID #chi^{2};N. Entries / %.1f" % bin_size("shower_pidchipr"),
    "shower_pidchipi": ";Shower #pi PID #chi^{2};N. Entries / %.1f" % bin_size("shower_pidchipi"),
    "shower_pidchimu": ";Shower #mu PID #chi^{2};N. Entries / %.1f" % bin_size("shower_pidchimu"),

    "track_likelihood": ";log(L_{MIP}/L_{p});N. Entries / %.1f" % bin_size("track_likelihood"),
    "track_mip_likelihood": ";L_{MIP};N. Entries / %.1f" % bin_size("track_mip_likelihood"),
    "track_p_likelihood": ";L_{p};N. Entries / %.1f" % bin_size("track_p_likelihood"),

    "track_pidchi": ";Track min. PID #chi^{2};N. Entries / %.1f" % bin_size("track_pidchi"),
    "track_pida": ";Track PIDa; N. Entries / %.1f" % bin_size("track_pida"),
    "track_res_mean": ";Track res. #mu [cm]; N. Entries / %.2f" % bin_size("track_res_mean"),
    "track_res_std": ";Track res. #sigma [cm]; N. Entries / %.2f" % bin_size("track_res_std"),
    "shower_res_mean": ";Shower res. #mu [cm]; N. Entries / %.2f" % bin_size("shower_res_mean"),
    "shower_res_std": ";Shower res. #sigma [cm]; N. Entries / %.2f" % bin_size("shower_res_std"),
    "nu_E": ";#nu_{e} energy [GeV];N. Entries / %.1f" % bin_size("nu_E"),
    "energy": ";#nu_{e} energy [GeV];N. Entries / %.2f GeV" % bin_size("energy"),
    "E_dep": ";E_{deposited} [GeV];N. Entries / %.1f" % bin_size("nu_E"),
    "n_objects": ";# objects;N.Entries / %i" % bin_size("n_objects"),
    "n_tracks": ";# tracks;N.Entries / %i" % bin_size("n_tracks"),
    "n_showers": ";# showers;N.Entries / %i" % bin_size("n_showers"),
    "n_tracks_before": ";# tracks;N.Entries / %i" % bin_size("n_tracks"),
    "n_showers_before": ";# showers;N.Entries / %i" % bin_size("n_showers"),
    "track_theta": ";Track #theta [#circ];N. Entries / %.1f#circ" % bin_size("track_theta"),
    "track_phi": ";Track #phi [#circ];N. Entries / %.1f#circ" % bin_size("track_phi"),
    "shower_theta": ";Shower #theta [#circ];N. Entries / %.1f#circ" % bin_size("shower_theta"),
    "shower_phi": ";Shower #phi [#circ];N. Entries / %.1f#circ" % bin_size("shower_phi"),
    "shower_distance": ";Shower distance [cm];N. Entries / %.1f cm" % bin_size("shower_distance"),
    "track_distance": ";Track distance [cm];N. Entries / %.1f cm" % bin_size("track_distance"),
    "track_shower_angle": ";cos#alpha [#circ];N. Entries / %.2f" % bin_size("track_shower_angle"),
    "track_start_y": ";Track start y [cm];N. Entries / %.1f cm" % bin_size("track_start_y"),
    "track_start_z": ";Track start z [cm];N. Entries / %.1f cm" % bin_size("track_start_z"),
    "track_start_x": ";Track start x [cm];;N. Entries / %.1f cm" % bin_size("track_start_x"),
    "track_end_y": ";Track end y [cm];N. Entries / %.1f cm" % bin_size("track_end_y"),
    "track_end_z": ";Track end z [cm];N. Entries / %.1f cm" % bin_size("track_end_z"),
    "track_end_x": ";Track end x [cm];N. Entries / %.1f cm" % bin_size("track_end_x"),
    "shower_start_y": ";Shower start y [cm];N. Entries / %.1f cm" % bin_size("shower_start_y"),
    "shower_start_z": ";Shower start z [cm];N. Entries / %.1f cm" % bin_size("shower_start_z"),
    "shower_start_x": ";Shower start x [cm];N. Entries / %.1f cm" % bin_size("shower_start_x"),
    "shower_end_y": ";Shower end y [cm];N. Entries / %.1f cm" % bin_size("shower_end_y"),
    "shower_end_z": ";Shower end z [cm];N. Entries / %.1f cm" % bin_size("shower_end_z"),
    "shower_end_x": ";Shower end x [cm];N. Entries / %.1f cm" % bin_size("shower_end_x"),
    "track_length": ";Track length [cm];N. Entries / %.1f cm" % bin_size("track_length"),
    "shower_length": ";Shower length [cm];N. Entries / %.1f cm" % bin_size("shower_length"),
    "proton_score": ";Proton score; N. Entries / %.1f" % bin_size("proton_score"),
    "shower_energy": ";Shower energy [GeV]; N. Entries / %.2f GeV" % bin_size("shower_energy"),
    "shower_energy_y": ";Shower energy (Y plane) [GeV]; N. Entries / %.2f GeV" % bin_size("shower_energy_y"),
    "pt": ";p_{t} [GeV/c];N. Entries / %.2f GeV/c" % bin_size("pt"),
    "p": ";p [GeV/c];N. Entries / %.2f GeV/c" % bin_size("p"),
    "pt_ratio": ";p_{t}/p [GeV/c];N. Entries / %.2f GeV/c" % bin_size("pt_ratio"),
    "vx": ";x [cm];",
    "vy": ";y [cm];",
    "vz": ";z [cm];",
    "reco_energy": ";E_{deposited} [GeV]; N. Entries / 0.05 GeV",
    "shower_open_angle": ";Shower open angle [#circ]; N. Entries / %.1f#circ" % bin_size("shower_open_angle"),
    "shower_dedx": ";Shower dE/dx (not calibrated) [MeV/cm]; N. Entries / %.1f MeV/cm" % bin_size("shower_dedx"),
    "shower_dedx_cali": ";Shower dE/dx [MeV/cm]; N. Entries / %.1f MeV/cm" % bin_size("shower_dedx_cali"),
    "shower_dqdx": ";Shower dQ/dx; N. Entries / %.1f " % bin_size("shower_dqdx"),

    "shower_dedx_u": ";Shower dE/dx (U plane) [MeV/cm]; N. Entries / %.1f MeV/cm" % bin_size("shower_dedx_u"),
    "shower_dedx_v": ";Shower dE/dx (V plane) [MeV/cm]; N. Entries / %.1f MeV/cm" % bin_size("shower_dedx_v"),
    "shower_dedx_merged": ";Shower dE/dx (merged hits) [MeV/cm]; N. Entries / %.1f MeV/cm" % bin_size("shower_dedx_merged"),
    "numu_score": ";#nu_{#mu} selection score; N. Entries / %i" % bin_size("numu_score"),
    "category": ";category",
    "event_weight": ";event_weight",
    "event": ";event",
    "run": ";run",
    "subrun": ";subrun",
    "interaction_type": ";;N. Entries",
    "is_signal": ";Selected events; N. Entries",
    "shower_pca": ";Shower PCA;N. Entries / %.3f" % bin_size("shower_pca"),
    "track_pca": ";Track PCA;N. Entries / %.3f" % bin_size("track_pca"),
    "total_shower_energy": ";Total shower E (not calibrated) [GeV]; N. Entries / %.2f GeV" % bin_size("total_shower_energy"),
    "total_shower_energy_cali": ";Total shower E [GeV]; N. Entries / %.2f GeV" % bin_size("total_shower_energy_cali"),
    "total_track_energy": ";Total track E [GeV]; N. Entries / %.2f GeV" % bin_size("total_track_energy"),
    "track_energy_length": ";Track E (stopping power) [GeV]; N. Entries / %.2f GeV" % bin_size("track_energy_length"),
    "shower_hits": ";Shower hits; N. Entries / %i" % bin_size("shower_hits"),
    "hits_ratio": ";Shower hits/total hits; N. Entries / %i" % bin_size("hits_ratio"),
    "theta_ratio": ";#theta_{e}/#theta_{p}; N. Entries / %i" % bin_size("theta_ratio"),
    "theta": ";Lepton #theta [#circ];N. Entries / %.1f#circ" % bin_size("track_theta"),
    "phi": ";Lepton #phi [#circ];N. Entries / %.1f#circ" % bin_size("track_theta"),
    "total_hits": ";Total hits; N. Entries / %i" % bin_size("shower_hits"),
    "total_hits_u": ";Total hits (U plane); N. Entries / %i" % bin_size("shower_hits"),
    "total_hits_v": ";Total hits (V plane); N. Entries / %i" % bin_size("shower_hits"),
    "total_hits_y": ";Total hits (Y plane); N. Entries / %i" % bin_size("shower_hits"),
    "track_hits": ";Track hits; N. Entries / %i" % bin_size("track_hits"),
    "track_dqdx": ";Track dQ/dx [a.u.]; N. Entries / %.2f a.u." % bin_size("track_dqdx"),
    "total_track_energy_length": ";Total track E (stopping power) [GeV]; N. Entries / %.2f GeV" % bin_size("total_track_energy_length")
}
