from array import array
import numpy as np
import math
import ROOT
from root_numpy import hist2array
from scipy.stats import poisson

ELECTRON_MASS = 0.51e-3 
PROTON_MASS = 0.938
PROTON_THRESHOLD = 0.040
ELECTRON_THRESHOLD = 0.020

bdt, manual = False, False
# Number to be obtained from Zarko's POT counting tool
total_data_bnb_pot_mcc83 = 4.903e+19
total_data_bnb_pot_mcc86 = 4.841e+19
total_data_bnb_pot_nogoodruns = 4.793e+19
total_data_bnb_goodruns_remap = 4.407e+19

total_data_bnb_pot_goodruns = 4.356e+19
total_data_bnb_pot_numu = 1.627e20
total_data_bnb_pot = total_data_bnb_goodruns_remap

# bdt_cut = 0.1470

bdt_cut = -0
rectangular_cut = 0.03

bins = array("f", [0, 0.200, 0.300, 0.375, 0.475, 0.550,
                   0.675, 0.800, 0.950, 1.100, 1.300, 1.500, 3.000])
bins = array("f", [0, 0.200, 0.300, 0.375, 0.475, 0.550,
                   0.675, 0.800, 0.950, 1.100, 1.300, 1.500, 3.000])

bins2 = array("f", [0, 0.200, 0.300, 0.375, 0.475, 0.550,
                   0.675, 0.800, 0.950, 1.100, 1.300, 1.500, 1.700])

# bins = array("f", [0, 0.2, 0.6, 1, 1.700, 3.000])

# bins = array("f", [0, 0.200, 0.650, 1.100, 1.500, 3.000])

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
            if not (is_fiducial(p_start)):
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
    sigma = par[2]
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

def printProgressBar (iteration, total, prefix = '', suffix = '', decimals = 1, length = 100, fill = '█'):
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
    print('\r%s |%s| %s%% %s' % (prefix, bar, percent, suffix), end = '\r')
    # Print New Line on Complete
    if iteration == total: 
        print()


def bin_size(name):
    return abs(binning[name][2] - binning[name][1]) / binning[name][0]


def FC_histo(h_bkg, bins=array("f", [])):
    if len(bins) > 0:
        h_limit = ROOT.TH1F("h_limit", "", len(bins) - 1, bins)
    else:
        h_limit = ROOT.TH1F("h_limit", "", h_bkg.GetNbinsX(), h_bkg.GetXaxis().GetXmin(), h_bkg.GetXaxis().GetXmax())

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

# bins = array("f", [0.2, 0.25, 0.3, 0.35,
#                    0.4, 0.45, 0.5, 0.6, 0.8, 1])

# bins = array("f", [0.2, 0.4, 0.5, 0.65, 0.8, 1])

# bins = array("f", [0, 0.200, 0.300, 0.375, 0.475, 0.550, 0.675, 0.800, 0.950, 1.100,  1.300, 1.500, 3.000])

def pre_cuts(chain):
    # dedx = chain.dedx < 3.5 
    # shower_distance = chain.shower_distance < 5
    # track_distance = chain.track_distance < 5
    # reco_energy = bins[0] < chain.reco_energy < bins[-1]
    # hits = chain.track_hits > 5 and chain.shower_hits > 40
    # track_start = [chain.track_start_x, chain.track_start_y, chain.track_start_z]
    # track_end = [chain.track_end_x, chain.track_end_y, chain.track_end_z]
    # shower_start = [chain.shower_start_x, chain.shower_start_y, chain.shower_start_z]
    # fiducial = is_fiducial(track_start) and is_fiducial(track_end) and is_fiducial(shower_start)
    numu = -13 <= chain.numu_score < 170 or manual
    track_pca = chain.track_pca > 0.9
    n_showers = chain.n_showers == 2
    # aa = 0.3 < chain.total_shower_energy/0.78 + 0.02 + chain.total_track_energy_length < 0.375
    shower_track_energy = chain.shower_energy > 0.01 and chain.total_shower_energy > 0.01 and chain.shower_energy_y > 0.01 and chain.total_track_energy_length > 0 
    reco_energy = bins[0] < chain.reco_energy < bins[-1]
    # reco_energy = 0.55 < chain.reco_energy < 0.675
    hits = chain.track_hits > 5 and chain.shower_hits_y > 5 and chain.total_hits_y > 0 and chain.total_hits_u > 0 and chain.total_hits_v > 0
    shower_angle = chain.track_shower_angle > -0.98 
    pre_bdt = 1 < chain.dedx < 3 and chain.shower_distance < 5 and chain.shower_track_d < 5 and chain.track_distance < 5 and chain.shower_distance < 5
    # reco_energy = 0.55 < chain.reco_energy < 0.675
    # dqdx_length = not (chain.track_dedx < 50000 and chain.track_length > 25)

    corrected_energy = 0.200 < ((chain.total_shower_energy_cali + 1.36881e-02) /
                                7.69908e-01) + chain.total_track_energy_length < 0.8
    return hits and numu and reco_energy and shower_track_energy# and chain.track_shower_angle > -0.95# and (not 100 < chain.track_start_z < 150)


def is_fiducial(point):
    ok_y = point[1] > y_start + 20 and point[1] < y_end - 20
    ok_x = point[0] > x_start + 10 and point[0] < x_end - 10
    ok_z = point[2] > z_start + 10 and point[2] < z_end - 50
    return ok_y and ok_x and ok_z


def fill_histos_data(tree_name, bdt, manual):
    f_data = ROOT.TFile("%s_file.root" % tree_name)
    t_data = f_data.Get("%s_tree" % tree_name)

    ROOT.TMVA.Tools.Instance()
    reader = ROOT.TMVA.Reader(":".join([
        "!V",
        "!Silent",
        "Color"]))

    for name, var in variables:
        t_data.SetBranchAddress(name, var)

    for name, var in spectators:
        t_data.SetBranchAddress(name, var)

    for name, var in variables:
        reader.AddVariable(name, var)

    for name, var in spectators:
        reader.AddSpectator(name, var)

    reader.BookMVA("BDT method",
                   "dataset/weights/TMVAClassification_BDT.weights.xml")
    # reader.BookMVA("Likelihood method",
    #                "dataset/weights/TMVAClassification_Likelihood.weights.xml")
    # reader.BookMVA("Cuts method",
    #                "dataset/weights/TMVAClassification_Cuts.weights.xml")

    variables_dict = dict(variables + spectators)

    histograms = []

    for i, n in enumerate(variables_dict.keys()):
        if n != "reco_energy":
            h = ROOT.TH1F("h_%s" % n, labels[n],
                          binning[n][0], binning[n][1], binning[n][2])
        else:
            h = ROOT.TH1F("h_%s" % n, labels[n], len(bins) - 1, bins)
        histograms.append(h)

    histo_dict = dict(zip(variables_dict.keys(), histograms))

    h_bdt = ROOT.TH1F("h_bdt_%s" % tree_name,
                      "BDT response; N. Entries / 0.05", 20, -1, 1)
    passed_events = 0

    for i in range(t_data.GetEntries()):
        t_data.GetEntry(i)
        BDT_response = reader.EvaluateMVA("BDT method")
        # likelihood_response = reader.EvaluateMVA("Likelihood method")
        # cuts_response = reader.EvaluateMVA("Cuts method", rectangular_cut)

        if pre_cuts(t_data):
            h_bdt.Fill(BDT_response, t_data.event_weight)

            if bdt:
                apply_bdt = BDT_response > bdt_cut
            else:
                apply_bdt = True

            if manual:
                apply_manual = manual_cuts(t_data)
            else:
                apply_manual = True

            if apply_bdt and apply_manual:
                # print(t_data.run, t_data.subrun, t_data.event, t_data.total_shower_energy / 0.78 +
                #       0.02 + t_data.total_track_energy_length)
                # if t_data.n_objects >= 5 and t_data.category == 6:
                passed_events += t_data.event_weight
                for name, var in variables:
                    histo_dict[name].Fill(var[0], t_data.event_weight)
                for name, var in spectators:
                    if name == "reco_energy":
                        corrected_energy_cali = ((t_data.total_shower_energy_cali * 0.987 + 1.36881e-02) /
                                                  7.69908e-01) + t_data.total_track_energy_length
                        corrected_energy = ((t_data.total_shower_energy + 1.23986e-02) /
                                            7.87131e-01) + t_data.total_track_energy_length

                        corrected_energy_hits = ((t_data.total_shower_energy_cali * 0.987 + 1.36881e-02) /
                                                 7.69908e-01) + (t_data.total_track_energy + 3.57033e-02) / 7.70870e-01
                                    
                        histo_dict[name].Fill(corrected_energy_cali, t_data.event_weight)
                        # histo_dict[name][category].Fill(chain.reco_energy, chain.event_weight)
                    else:
                        histo_dict[name].Fill(var[0], t_data.event_weight)

    f_bdt = ROOT.TFile("plots/h_bdt_%s.root" % tree_name, "RECREATE")
    h_bdt.Write()
    f_bdt.Close()

    for h in histograms:
        f = ROOT.TFile("plots/%s_%s.root" % (h.GetName(), tree_name),
                       "RECREATE")
        h.Write()
        f.Close()

    return passed_events


def find_interaction(dictionary, interaction):
    for name, id_int in dictionary.items():
        if id_int == interaction:
            return name


def manual_cuts(chain):
    track_shower_angle = -0.9 < chain.track_shower_angle 
    shower_distance = chain.shower_distance < 5
    shower_open_angle = 1 < chain.shower_open_angle < 19
    dedx = 1 < chain.dedx_merged < 3.2
    total_hits_y = chain.total_hits_y > 100
    shower_energy = chain.total_shower_energy_cali > 0.06
    track_length = chain.track_length < 80
    track_distance = chain.track_distance < 5
    track_res = chain.track_res_std < 2
    shower_angle = not (5 < chain.shower_angle < 45)
    dqdx = chain.dqdx_bdt_max > 0.1# and (chain.dqdx_bdt > 0.1 or chain.dqdx_bdt < -0.2)
    corrected_energy =  0.300 < (chain.total_shower_energy_cali + 0.02)/ \
        0.78 + chain.total_track_energy_length < 0.425

    numu_exclude = chain.numu_score != 17 and chain.numu_score != 2

    ratio = chain.hits_ratio > 0.55
    dedx_bdt  = chain.dedx_bdt > -10

    if chain.no_tracks == 1:
        dedx_bdt = chain.dedx_bdt > -0.1
        shower_distance = chain.shower_distance < 2
        track_distance = chain.track_distance < 2


    collabmeeting_cuts = [dedx, total_hits_y, track_shower_angle, numu_exclude, shower_energy, ratio, dqdx, track_res, dedx_bdt,
                          track_distance, shower_distance, shower_open_angle, track_length]

    # collabmeeting_cuts = [total_hits_y, numu_exclude, shower_energy, ratio, dqdx, track_res, dedx_bdt,
    #                       track_shower_angle, shower_open_angle, track_length]

    # collabmeeting_cuts = [dedx_bdt, track_res, ratio, shower_distance, numu_exclude, shower_open_angle,
    #                     track_length, track_distance, dqdx, total_hits_y,
    #                     track_shower_angle, shower_energy, dedx]

    dedx_inverted = 3.2 < chain.dedx < 6
    inverted_cuts = [dedx_inverted, shower_open_angle, total_hits_y, track_length, track_distance,
                     track_shower_angle, shower_energy, dqdx]

                     
    dedx_harder = 1.5 < chain.dedx < 2.6
    bdt_harder = chain.dqdx_bdt > 0.1
    energy_harder = chain.shower_energy > 0.1
    harder_cuts = collabmeeting_cuts + [energy_harder, dedx_harder, bdt_harder, chain.no_tracks == 0]

    track_length2 = chain.track_length > 20
    
    numu_cuts = [dedx, chain.dedx < 3, shower_open_angle, total_hits_y, track_length2,
    not numu_exclude, chain.numu_score != 2,  track_shower_angle, shower_energy, not dqdx]

    passed_numu = len(numu_cuts) == sum(numu_cuts)

    passed_collab = len(collabmeeting_cuts) == sum(collabmeeting_cuts)
    passed_inverted = len(inverted_cuts) == sum(inverted_cuts)
    passed_harder = len(harder_cuts) == sum(harder_cuts)
    return passed_collab

def sigma_calc_matrix(h_signal, h_background, scale_factor=1, sys=0):
    #it is just, Δχ2 = (number of events signal in Energy bins in a 1D matrix)
    #(2D Matrix - (statistical uncertainty)^2 in a the diagonal of the matrix)^-1
    #(number of events signal in Energy bins in a 1D matrix)^Transpose
    tries = 50
    sig_array = h_signal * scale_factor
    bkg_array = h_background * scale_factor

    nbins = len(sig_array)

    emtx = np.zeros((nbins, nbins))
    sig_err = np.zeros((tries, nbins))
    emtx_err = np.zeros((tries, nbins, nbins))

    for x in range(nbins):
        emtx[x][x] = bkg_array[x] + (sys*bkg_array[x])**2
        for i in range(tries):
            r = poisson.rvs(bkg_array[x] / 0.02442477422674749)
            s = poisson.rvs(sig_array[x] / 0.02158221715003847)
            sig_err[i][x] = s
            emtx_err[i][x][x] = max(1, r + (sys*r)**2)

    emtxinv = np.linalg.inv(emtx)
    sigmas = []
    for i in range(tries):
        try:
            emtxinv_err = np.linalg.inv(emtx_err[i])
            chisq_err = float(sig_err[i].dot(emtxinv_err).dot(sig_err[i].T))
            sigmas.append(math.sqrt(chisq_err))
        except np.linalg.linalg.LinAlgError:
            print("Singular matrix")
            continue

    sigma_err = np.std(sigmas)
    chisq = float(sig_array.dot(emtxinv).dot(sig_array.T))

    #print "Sqrt of that (==sigma?) is ",np.sqrt(chisq)
    return np.sqrt(chisq), sigma_err / np.sqrt(tries)


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
    "kQE": 0,
    "kRes": 1,
    "kDIS": 2,
    "kCoh": 3,
    "kCohElastic": 4,
    "kElectronScattering": 5,
    "kIMDAnnihilation": 6,
    "kInverseBetaDecay": 7,
    "kGlashowResonance": 8,
    "kAMNuGamma": 9,
    "kMEC": 10,
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

shower_angle = array("f", [0])
track_angle = array("f", [0])
shower_pdg = array("f", [0])
track_pdg = array("f", [0])

track_length = array("f", [0])
track_theta = array("f", [0])
track_phi = array("f", [0])
shower_theta = array("f", [0])
shower_phi = array("f", [0])
shower_energy = array("f", [0])
shower_energy_y = array("f", [0])

pt = array("f", [0])
n_tracks = array("f", [0])
n_showers = array("f", [0])
track_shower_angle = array("f", [0])
track_id = array("f", [0])
shower_start_x = array("f", [0])
track_start_x = array("f", [0])
shower_distance = array("f", [0])
shower_true_distance = array("f", [0])

proton_score = array("f", [0])
track_distance = array("f", [0])
track_start_y = array("f", [0])
track_end_y = array("f", [0])
track_start_x = array("f", [0])
track_end_x = array("f", [0])
track_start_z = array("f", [0])
track_end_z = array("f", [0])
shower_start_y = array("f", [0])
shower_end_y = array("f", [0])
shower_start_x = array("f", [0])
shower_end_x = array("f", [0])
shower_start_z = array("f", [0])
shower_end_z = array("f", [0])
reco_energy = array("f", [0])
event_weight = array("f", [0])
category = array("f", [0])
event = array("f", [0])
run = array("f", [0])
subrun = array("f", [0])
interaction_type = array("f", [0])
is_signal = array("f", [0])
shower_open_angle = array("f", [0])
dedx = array("f", [0])
dedx_cali = array("f", [0])
dedx_u = array("f", [0])
dedx_v = array("f", [0])

numu_score = array("f", [0])
shower_pca = array("f", [0])
track_pca = array("f", [0])
total_shower_energy = array("f", [0])
total_shower_energy_cali = array("f", [0])

total_shower_energy_y = array("f", [0])

total_track_energy = array("f", [0])
track_energy = array("f", [0])
shower_hits = array("f", [0])
hits_ratio = array("f", [0])

shower_hits_y = array("f", [0])
track_hits = array("f", [0])
track_dedx = array("f", [0])
track_energy_length = array("f", [0])
total_track_energy_length = array("f", [0])
n_objects = array("f", [0])
dedx_merged = array("f", [0])
total_hits = array("f", [0])
total_hits_u = array("f", [0])
total_hits_v = array("f", [0])
total_hits_y = array("f", [0])
no_tracks = array("f", [0])
shower_length = array("f", [0])
nu_E = array("f", [0])
E_dep = array("f", [0])
track_pidchipr = array("f", [0])
track_pidchi = array("f", [0])
track_pida = array("f", [0])
track_res_mean = array("f", [0])
track_res_std = array("f", [0])
shower_res_mean = array("f", [0])
shower_res_std = array("f", [0])
true_nu_is_fidvol = array("f", [0])
dqdx_bdt = array("f", [0])
dqdx_pion = array("f", [0])

dqdx_bdt_max = array("f", [0])

dedx_bdt = array("f", [0])
shower_track_d = array("f", [0])

spectators = [
    # ("track_pdg", track_pdg),
    # ("shower_pdg", shower_pdg),
    ("dedx_bdt", dedx_bdt),
    ("true_nu_is_fidvol", true_nu_is_fidvol),
    ("total_hits_u", total_hits_u),
    ("total_hits_v", total_hits_v),
    ("track_pidchipr", track_pidchipr),
    ("track_pidchi", track_pidchi),

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
    # ("total_shower_energy_y", total_shower_energy_y),
    ("total_hits_y", total_hits_y),
    ("track_energy_length", track_energy_length),
    ("n_tracks", n_tracks),
    ("track_energy", track_energy),
    ("track_start_x", track_start_x),
    ("track_end_x", track_end_x),
    ("track_start_z", track_start_z),
    ("track_end_z", track_end_z),
    ("shower_start_z", shower_start_z),
    ("track_start_y", track_start_y),
    ("track_end_y", track_end_y),
    ("shower_start_y", shower_start_y),
    ("shower_phi", shower_phi),
    ("track_phi", track_phi),
    ("shower_hits", shower_hits),    
    ("shower_hits_y", shower_hits_y),
    ("track_hits", track_hits),
    ("shower_energy_y", shower_energy_y),
    ("track_dedx", track_dedx),
    # ("proton_score", proton_score),
    ("shower_start_x", shower_start_x),
    # ("shower_length", shower_length),
    ("no_tracks", no_tracks),
    ("dedx_merged", dedx_merged),
    ("n_objects", n_objects),
    ("pt", pt),
    ("shower_theta", shower_theta),
    ("track_theta", track_theta),
    ("track_pca", track_pca),
    ("shower_pca", shower_pca),
    ("track_res_mean", track_res_mean),
    ("track_res_std", track_res_std),
    ("shower_res_mean", shower_res_mean),
    ("shower_res_std", shower_res_std),
    ("total_hits", total_hits),
    ("total_shower_energy", total_shower_energy),
    ("total_shower_energy_cali", total_shower_energy_cali),

    ("total_track_energy_length", total_track_energy_length),    
    ("shower_energy", shower_energy),
    # ("dqdx_pion", dqdx_pion),
    # ("shower_true_distance", shower_true_distance),
    ("numu_score", numu_score),
    ("dedx_cali", dedx_cali),
    ("dedx_u", dedx_u),
    ("dedx_v", dedx_v),
    # ("track_angle", track_angle),

]

variables = [ 
    ("n_showers", n_showers),
    ("shower_track_d", shower_track_d),
    ("shower_open_angle", shower_open_angle),
    ("track_distance", track_distance),
    ("shower_distance", shower_distance),
    ("dedx", dedx),
    ("track_length", track_length),
    ("track_shower_angle", track_shower_angle),
    ("dqdx_bdt", dqdx_bdt),
    ("dqdx_bdt_max", dqdx_bdt_max),
    ("hits_ratio", hits_ratio),
    ("shower_angle", shower_angle),
]

binning = {
    "track_pdg": [100, 0, 100],
    "shower_pdg": [100, 0, 100],

    "dqdx_bdt": [20, -0.5, 0.5],
    "dqdx_pion": [50, -1, 1],
    "dqdx_bdt_max": [50, -1, 1],
    "dedx_bdt": [50, -1, 1],
    "shower_angle": [20, 0, 180],
    "track_angle": [20, 0, 180],
    "pi0_mass": [20, 0, 0.4],
    "shower_track_d": [30, 0, 15],
    "true_nu_is_fidvol": [2, 0, 2],
    "total_hits_u": [1,0,1],
    "total_hits_v": [1,0,1],
    "track_pidchipr": [20, 0, 40],
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
    "track_theta": [18, 0, 180],
    "track_phi": [18, -180, 180],
    "shower_theta": [18, 0, 180],
    "shower_phi": [18, -180, 180],
    "shower_distance": [20, 0, 15],
    "shower_true_distance": [20, 0, 15],
    "track_distance": [20, 0, 15],
    "track_shower_angle": [20, -1, 1],
    "track_start_y": [20, y_start, y_end],
    "track_start_z": [20, z_start, z_end],
    "track_start_x": [20, x_start, x_end],
    "track_end_y": [20, y_start, y_end],
    "track_end_z": [20, z_start, z_end],
    "track_end_x": [20, x_start, x_end],
    "shower_start_y": [20, y_start, y_end],
    "shower_start_z": [20, z_start, z_end],
    "shower_start_x": [20, x_start, x_end],
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
    "reco_energy": [40, 0, 2],
    "shower_open_angle": [46, 0, 46],
    "dedx": [20, 0, 5],
    "dedx_cali": [20, 0, 5],
    "dedx_u": [20, 0, 5],
    "dedx_v": [20, 0, 5],

    "dedx_merged": [20, 0, 5],
    "numu_score": [20, 0, 20],
    "category": [11, 0, 11],
    "event_weight": [20, 0, 100],
    "event": [20, 0, 10000],
    "run": [20, 0, 10],
    "subrun": [20, 0, 1000],
    "interaction_type": [100, 1000, 1100],
    "is_signal": [2, 0, 1],
    "shower_pca": [20, 0.9, 1],
    "track_pca": [20, 0.99, 1],
    "total_shower_energy": [20, 0, 0.5],
    "total_shower_energy_cali": [20, 0, 0.5],

    "total_shower_energy_y": [20, 0, 0.5],
    "total_track_energy": [20, 0, 1],
    "shower_hits": [20, 0, 200],
    "hits_ratio": [20, 0, 1],

    "shower_hits_y": [20, 0, 500],
    "nu_E": [40, 0, 4],
    "E_dep": [30, 0, 3],

    "track_hits": [20, 0, 400],
    "track_energy": [20, 0, 1],
    "track_dedx": [40, 0, 700],
    "track_energy_length": [20, 0, 1],
    "total_track_energy_length": [20, 0, 2]
}

if manual:
    binning["shower_theta"] = [9, 0, 180]
    binning["shower_phi"] = [9, -180, 180]
    binning["shower_start_y"] = [10, y_start, y_end]
    binning["shower_start_z"] = [10, z_start, z_end]
    binning["shower_start_x"] = [10, x_start, x_end]


labels = {
    "track_pdg": "track pdg",
    "shower_pdg": "shower pdg",

    "dqdx_bdt": ";Proton BDT;N. Entries / %.2f" % bin_size("dqdx_bdt"),
    "dqdx_pion": ";dQ/dx BDT;",
    
    "dqdx_bdt_max": ";dQ/dx BDT;",

    "dedx_bdt": ";dE/dx BDT;",
    "shower_track_d": ";shower track d;",

    "shower_angle": ";shower_angle;",
    "track_angle": ";track_angle;",

    "pi0_mass": ";pi0 mass;",

    "true_nu_is_fidvol": ";true_nu_is_fidvol;",
    "total_hits_u": ";total hits u;",
    "total_hits_v": "; total hits u;",
    "no_tracks": ";No tracks",
    "track_pidchipr": ";Track proton PID #chi^{2};N. Entries / %.1f" % bin_size("track_pidchipr"),
    "track_pidchi": ";Track min. PID #chi^{2};N. Entries / %.1f" % bin_size("track_pidchi"),

    "track_pida": ";Track PIDa; N. Entries / %.1f" % bin_size("track_pida"),
    "track_res_mean": ";Track res. #mu [cm]; N. Entries / %.2f" % bin_size("track_res_mean"),
    "track_res_std": ";Track res. #sigma [cm]; N. Entries / %.2f" % bin_size("track_res_std"),
    "shower_res_mean": ";Shower res. #mu [cm]; N. Entries / %.2f" % bin_size("shower_res_mean"),
    "shower_res_std": ";Shower res. #sigma [cm]; N. Entries / %.2f" % bin_size("shower_res_std"),
    "nu_E": ";True E #nu;N. Entries / %.1f" % bin_size("nu_E"),
    "E_dep": ";True E #nu;N. Entries / %.1f" % bin_size("nu_E"),

    "n_objects": ";# objects;N.Entries / %i" % bin_size("n_objects"),
    "n_tracks": ";# tracks;N.Entries / %i" % bin_size("n_tracks"),
    "n_showers": ";# showers;N.Entries / %i" % bin_size("n_showers"),
    "track_theta": ";Track #theta [#circ];N. Entries / %.1f#circ" % bin_size("track_theta"),
    "track_phi": ";Track #phi [#circ];N. Entries / %.1f#circ" % bin_size("track_phi"),
    "shower_theta": ";Shower #theta [#circ];N. Entries / %.1f#circ" % bin_size("shower_theta"),
    "shower_phi": ";Shower #phi [#circ];N. Entries / %.1f#circ" % bin_size("shower_phi"),
    "shower_distance": ";Track distance [cm];N. Entries / %.1f cm" % bin_size("shower_distance"),
    "shower_true_distance": ";Shower-true #nu distance [cm];N. Entries / %.1f cm" % bin_size("shower_true_distance"),

    "track_distance": ";Shower distance [cm];N. Entries / %.1f cm" % bin_size("track_distance"),
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
    "reco_energy": ";E_{deposited} [GeV]; N. Entries / 0.05 GeV",
    "shower_open_angle": ";Shower open angle [#circ]; N. Entries / %.1f#circ" % bin_size("shower_open_angle"),
    "dedx": ";Shower dE/dx (not calibrated) [MeV/cm]; N. Entries / %.1f MeV/cm" % bin_size("dedx"),
    "dedx_cali": ";Shower dE/dx [MeV/cm]; N. Entries / %.1f MeV/cm" % bin_size("dedx_cali"),
    "dedx_u": ";Shower dE/dx (U plane) [MeV/cm]; N. Entries / %.1f MeV/cm" % bin_size("dedx_u"),
    "dedx_v": ";Shower dE/dx (V plane) [MeV/cm]; N. Entries / %.1f MeV/cm" % bin_size("dedx_v"),
    "dedx_merged": ";Shower dE/dx (merged hits) [MeV/cm]; N. Entries / %.1f MeV/cm" % bin_size("dedx_merged"),
    "numu_score": ";#nu_{#mu} selection score; N. Entries / %i" % bin_size("numu_score"),
    "category": ";category",
    "event_weight": ";event_weight",
    "event": ";event",
    "run": ";run",
    "subrun": ";subrun",
    "interaction_type": ";interaction_type",
    "is_signal": ";is_signal",
    "shower_pca": ";Shower PCA;N. Entries / %.3f" % bin_size("shower_pca"),
    "track_pca": ";Track PCA;N. Entries / %.3f" % bin_size("track_pca"),
    "total_shower_energy": ";Total shower E (not calibrated) [GeV]; N. Entries / %.2f GeV" % bin_size("total_shower_energy"),
    "total_shower_energy_cali": ";Total shower E [GeV]; N. Entries / %.2f GeV" % bin_size("total_shower_energy_cali"),

    "total_track_energy": ";Total track E [GeV]; N. Entries / %.2f GeV" % bin_size("total_track_energy"),
    "shower_hits": ";Shower hits; N. Entries / %i" % bin_size("shower_hits"),
    "hits_ratio": ";Shower hits/total hits; N. Entries / %i" % bin_size("hits_ratio"),

    "total_hits": ";Total hits; N. Entries / %i" % bin_size("shower_hits"),
    "total_hits_u": ";Total hits (U plane); N. Entries / %i" % bin_size("shower_hits"),
    "total_hits_v": ";Total hits (V plane); N. Entries / %i" % bin_size("shower_hits"),
    "total_hits_y": ";Total hits (Y plane); N. Entries / %i" % bin_size("shower_hits"),

    "shower_hits_y": ";Shower hits (Y plane); N. Entries / %i" % bin_size("shower_hits_y"),
    "track_hits": ";Track hits; N. Entries / %i" % bin_size("track_hits"),
    "track_energy": ";Track E [GeV]; N. Entries / %.2f GeV" % bin_size("track_energy"),
    "track_dedx": ";Track dQ/dx [a.u.]; N. Entries / %.2f a.u." % bin_size("track_dedx"),
    "track_energy_length": ";Track E (stopping power) [GeV]; N. Entries / %.2f GeV" % bin_size("track_energy_length"),
    "total_track_energy_length": ";Total track E (stopping power) [GeV]; N. Entries / %.2f GeV" % bin_size("total_track_energy_length")

}
