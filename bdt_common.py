from array import array
import math
import ROOT
from scipy.stats import poisson

ELECTRON_MASS = 0.51e-3
PROTON_MASS = 0.938


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
    bar = fill * filledLength + '-' * (length - filledLength)
    print('\r%s |%s| %s%% %s' % (prefix, bar, percent, suffix), end = '\r')
    # Print New Line on Complete
    if iteration == total: 
        print('\r%s |%s| %s%% %s' % (prefix, bar, percent, suffix), end='\r')
        print()

def bin_size(name):
    return binning[name][2] / (binning[name][0] - binning[name][1])

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
    nmin = max(0,  int(bkg - 5. * sigma))   # Use 0 if nmin<0
    nmax = max(20, int(bkg + 5. * sigma)+1) # Use at least 20 for low means
    po = poisson(bkg)
    UL = 0.
    for i in range(nmin, nmax):
        pmf = po.pmf(i)
        ul = FC.CalculateUpperLimit(i, bkg)
        #print "i=%i, Po(i)=%f, U(i,b)=%f" % (i, pmf, ul)
        UL += po.pmf(i) * ul
    return UL

colors = [ROOT.TColor.GetColor("#efac3a"), ROOT.TColor.GetColor("#e7623d"),
          ROOT.TColor.GetColor("#62b570"), ROOT.TColor.GetColor("#8779b1"),
          ROOT.TColor.GetColor("#609ec6"), ROOT.TColor.GetColor("#c46d27"),
          ROOT.TColor.GetColor("#ffffff"), ROOT.TColor.GetColor("#e7623d"),
          ROOT.TColor.GetColor("#1e7a2e")]

bdt, manual = False, False
# Number to be obtained from Zarko's POT counting tool
total_data_bnb_pot_mcc83 = 4.903e+19
total_data_bnb_pot_mcc86 = 4.857e+19
total_data_bnb_pot = total_data_bnb_pot_mcc86
# bdt_cut = 0.1470

bdt_cut = 0.1131
rectangular_cut = 0.565

SIGNAL_INTERVAL = [0, 6]

bins = array("f", [0.2, 0.25, 0.3, 0.35,
                   0.4, 0.45, 0.5, 0.6, 0.8, 1])

# bins = array("f", [0.2, 0.4, 0.5, 0.65, 0.8, 1])

bins = array("f", [0.200,  0.300,  0.375,  0.475,  0.550,  0.675, 0.800,  0.950,  1.100,   1.300,  1.500,  3.000])

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
    numu = 12 < chain.numu_score
    track_pca = chain.track_pca > 0.9
    shower_track_energy = chain.shower_energy > 0.01 and chain.total_shower_energy > 0.01 and chain.shower_energy_y > 0.01
    reco_energy = bins[0] < chain.reco_energy < bins[-1]
    hits = chain.track_hits > 5 and chain.shower_hits_y > 0 and chain.total_hits_u > 0 and chain.total_hits_v > 0 and chain.total_hits_y > 0
    shower_angle = chain.track_shower_angle > -0.98
    return reco_energy and hits and shower_angle  and shower_track_energy and chain.total_track_energy_length > 0



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
                      "BDT response; N. Entries / 0.05", 40, -1, 1)
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
                passed_events += t_data.event_weight
                for name, var in variables:
                    histo_dict[name].Fill(var[0], t_data.event_weight)
                for name, var in spectators:
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
    shower_energy = chain.shower_energy > 0.1
    track_energy = chain.total_track_energy > 0.015
    dedx = 1.4 < chain.dedx < 2.7
    shower_distance = chain.shower_distance < 3.5
    track_distance = chain.track_distance < 4
    track_length = chain.track_length < 100
    proton_score = chain.proton_score > 0.95
    open_angle = 1 < chain.shower_open_angle < 15
    shower_theta = chain.shower_theta < 90
    n_tracks = chain.n_tracks < 3
    track_shower_angle = -0.9 < chain.track_shower_angle < 0.9
    track_theta = chain.track_theta < 130
    shower_nhits = chain.shower_hits > 70
    shower_pca = chain.shower_pca < 0.995
    pt = chain.pt < 0.7
    track_nhits = chain.track_hits < 300
    # dedx = 0.6 < chain.dedx < 3.07
    # proton_score = chain.proton_score > 0.57
    # shower_distance = chain.shower_distance < 2.23
    # track_distance = chain.track_distance < 9.55
    # open_angle = 1 < chain.shower_open_angle < 23.7
    # shower_theta = 15 < chain.shower_theta < 99

    cuts = [shower_energy, dedx, proton_score, open_angle,
            shower_theta, shower_distance, track_distance,
            shower_pca, track_shower_angle, pt, shower_nhits,
            track_nhits]

    dedx = 0.8576 < chain.dedx < 5.907
    proton_score = 0.1426 < chain.proton_score < 1
    shower_distance = 0.01 < chain.shower_distance < 4
    track_distance = 3.0483806859534085e-03 < chain.track_distance < 8.9818492209332348e+00
    shower_open_angle = 1.1999870880917283e+00 < chain.shower_open_angle < 3.3774471946669600e+01
    shower_theta = 1.7068450029571980e+00 < chain.shower_theta < 1.7098164676405759e+02
    track_phi = -1.7474272785368544e+02 < chain.track_phi < 180
    shower_phi = -180 < chain.shower_phi < 170
    track_theta = 1.1442699092190036e+01 < chain.track_theta < 1.6352903537803476e+02
    track_shower_angle = -1 < chain.track_shower_angle < 9.4899730606947252e-01

    optimized_cuts = [dedx, proton_score, shower_distance, track_distance, 
                      shower_open_angle, shower_theta, track_phi, shower_phi, track_theta, track_shower_angle]

    passed = len(cuts) == sum(cuts)
    passed_optimized = len(optimized_cuts) == sum(optimized_cuts)
    no_tracks = chain.no_tracks < 1
    shower_hits = chain.shower_hits > 40
    track_hits = chain.track_hits > 5
    track_shower_angle = -0.9 < chain.track_shower_angle 
    track_distance = chain.track_distance < 3
    shower_pca = 0.835 < chain.shower_pca < 0.995
    n_tracks = chain.n_tracks < 4
    n_showers = chain.n_showers < 4
    shower_distance = chain.shower_distance < 2.8
    track_energy = 0.01 < chain.total_track_energy
    track_dedx = 1 < chain.track_dedx 
    shower_open_angle = 1 < chain.shower_open_angle < 25
    shower_energy = chain.total_shower_energy > 0.16
    dedx = 1 < chain.dedx 
    pt = chain.pt < 0.7
    total_hits = chain.total_hits > 500
    total_hits_y = chain.total_hits_y > 170
    shower_start_z = chain.shower_start_z < 1000
    proton_score = chain.proton_score > 0.1 or chain.proton_score < 0
    track_theta = chain.track_theta < 155
    shower_theta = chain.shower_theta < 90
    shower_energy2 = chain.shower_energy > 0.15
    category = chain.category != 2 and chain.category != 8
    track_end_y = chain.track_end_y < 110
    dedx_y = 1.2 < chain.dedx_y < 3.2
    track_length = chain.track_length < 120
    numu_cut = 11 < chain.numu_score < 16 or chain.numu_score < 10
    track_distance = chain.track_distance < 3
    # collabmeeting_cuts = [dedx, total_hits, shower_theta, track_length, track_dedx, shower_energy, shower_energy2, dedx_y, numu_cut, shower_distance, track_shower_angle]
    collabmeeting_cuts = [dedx_y, numu_cut, shower_distance]

    shower_energy = chain.shower_energy > 0.2
    shower_energy2 = chain.total_shower_energy > 0.2

    dedx = 1.2 < chain.dedx < 5.5
    proton_score = chain.proton_score > 0.9
    openangle = 1 < chain.shower_open_angle < 20
    theta = chain.shower_theta < 90
    shower_distance = chain.shower_distance < 4
    track_distance = chain.track_distance < 4

    oscmeeting_cuts = [n_showers, shower_pca, shower_energy2, shower_energy, dedx, proton_score, openangle, theta, shower_distance, track_distance]
    
    theta = chain.shower_theta < 90
    dedx = 0.9 < chain.dedx < 3
    showerangle = chain.track_shower_angle > -0.93 
    bdt_cuts = [dedx, shower_theta, showerangle]
    bdt_cuts = [dedx]

    passed_bdt = len(bdt_cuts) == sum(bdt_cuts)
    passed_collab = len(collabmeeting_cuts) == sum(collabmeeting_cuts)
    passed_osc = len(oscmeeting_cuts) == sum(oscmeeting_cuts)

    return passed_collab

def sigmaCalc(h_signal, h_background, sys_err=0):
    for i in range(h_signal.GetNbinsX()):
        print(h_signal.GetBinContent(i))
        print(h_background.GetBinContent(i))
    chi2 = sum(
        [h_signal.GetBinContent(i)**2 / (h_background.GetBinContent(i))
         for i in range(1, h_signal.GetNbinsX() + 1)
         if h_signal.GetBinContent(i) > 0])

    return math.sqrt(chi2)


total_pot = total_data_bnb_pot

description = ["#nu_{e} CC0#pi-Np",
               "#nu_{e} CC",
               "Beam Intrinsic #nu_{#mu}",
               "Beam Intrinsic NC",
               "Dirt",
               "Cosmic contaminated",
               "Cosmic",
               "Cosmic in-time",
               "Data",
               "Low-energy excess"]

interactions = {
    "kQE": 0,
    "kRes": 1,
    "kDIS": 2,
    "kCoh": 3,
    "kNuanceOffset": 1000,
    "k1eNp": 1000 + 1,
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
numu_score = array("f", [0])
dedx_hits = array("f", [0])
shower_pca = array("f", [0])
track_pca = array("f", [0])
total_shower_energy = array("f", [0])
total_shower_energy_y = array("f", [0])

total_track_energy = array("f", [0])
total_track_energy_dedx = array("f", [0])
track_energy = array("f", [0])
track_energy_dedx = array("f", [0])
shower_hits = array("f", [0])
shower_hits_y = array("f", [0])
track_hits = array("f", [0])
track_dedx = array("f", [0])
track_energy_length = array("f", [0])
total_track_energy_length = array("f", [0])
n_objects = array("f", [0])
dedx_y = array("f", [0])
total_hits = array("f", [0])
total_hits_u = array("f", [0])
total_hits_v = array("f", [0])
total_hits_y = array("f", [0])
no_tracks = array("f", [0])
shower_length = array("f", [0])
spectators = [
    ("category", category),
    ("event_weight", event_weight),
    ("event", event),
    ("run", run),
    ("subrun", subrun),
    ("interaction_type", interaction_type),
    ("is_signal", is_signal),
    ("dedx_hits", dedx_hits),
    ("reco_energy", reco_energy),
    ("total_track_energy", total_track_energy),
    ("total_track_energy_dedx", total_track_energy),
    ("track_energy_dedx", track_energy_dedx),
    ("total_shower_energy", total_shower_energy),
    ("total_shower_energy_y", total_shower_energy_y),
    ("total_hits", total_hits),
    ("total_hits_u", total_hits_u),
    ("total_hits_v", total_hits_v),
    ("total_hits_y", total_hits_y),

    ("track_energy_length", track_energy_length),
    ("total_track_energy_length", total_track_energy_length),
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
    ("shower_energy", shower_energy),
    ("shower_energy_y", shower_energy_y),
    ("pt", pt),
    ("shower_theta", shower_theta),
    ("track_theta", track_theta),
    ("n_showers", n_showers),
    ("track_dedx", track_dedx),
    ("proton_score", proton_score),
    ("shower_start_x", shower_start_x),
    ("dedx", dedx),
    ("shower_length", shower_length),
    ("no_tracks", no_tracks)

]

variables = [ 
    ("dedx_y", dedx_y),
    ("n_objects", n_objects),
    ("numu_score", numu_score),
    ("shower_open_angle", shower_open_angle),
    ("track_length", track_length),
    ("track_distance", track_distance),
    ("shower_distance", shower_distance),
    ("shower_pca", shower_pca),
    ("track_pca", track_pca),
    ("track_shower_angle", track_shower_angle)
]

binning = {
    "no_tracks": [2, 0, 2],
    "n_objects": [8, 2, 10],
    "n_tracks": [5, 1, 6],
    "n_showers": [5, 1, 6],
    "track_theta": [36, 0, 180],
    "track_phi": [36, -180, 180],
    "shower_theta": [18, 0, 180],
    "shower_phi": [36, -180, 180],
    "shower_distance": [25, 0, 15],
    "track_distance": [25, 0, 15],
    "track_shower_angle": [40, -1, 1],
    "track_start_y": [30, y_start, y_end],
    "track_start_z": [30, z_start, z_end],
    "track_start_x": [30, x_start, x_end],
    "track_end_y": [30, y_start, y_end],
    "track_end_z": [30, z_start, z_end],
    "track_end_x": [30, x_start, x_end],
    "shower_start_y": [30, y_start, y_end],
    "shower_start_z": [30, z_start, z_end],
    "shower_start_x": [30, x_start, x_end],
    "shower_end_y": [30, y_start, y_end],
    "shower_end_z": [30, z_start, z_end],
    "shower_end_x": [30, x_start, x_end],
    "track_length": [50, 0, 200],
    "shower_length": [50, 0, 200],
    "proton_score": [10, 0, 1],
    "shower_energy": [100, 0, 0.5],
    "shower_energy_y": [100, 0, 0.5],
    "total_hits": [50, 0, 1000],
    "total_hits_u": [50, 0, 500],
    "total_hits_v": [50, 0, 500],
    "total_hits_y": [50, 0, 500],
    "pt": [20, 0, 1],
    "reco_energy": [40, 0, 2],
    "shower_open_angle": [46, 0, 46],
    "dedx": [30, 0, 6],
    "dedx_y": [30, 0, 6],
    "numu_score": [20, 0, 20],
    "category": [7, 0, 7],
    "event_weight": [20, 0, 100],
    "event": [20, 0, 10000],
    "run": [20, 0, 10],
    "subrun": [20, 0, 1000],
    "interaction_type": [100, 1000, 1100],
    "is_signal": [2, 0, 1],
    "dedx_hits": [50, 0, 1],
    "shower_pca": [50, 0.8, 1],
    "track_pca": [50, 0.99, 1],
    "total_shower_energy": [50, 0, 0.5],
    "total_shower_energy_y": [50, 0, 0.5],
    "total_track_energy": [40, 0, 1],
    "total_track_energy_dedx": [10, 0, 1],
    "shower_hits": [50, 0, 200],
    "shower_hits_y": [50, 0, 500],

    "track_hits": [50, 0, 400],
    "track_energy": [50, 0, 1],
    "track_energy_dedx": [50, 0, 1],
    "track_dedx": [40, 0.2, 6],
    "track_energy_length": [50, 0, 1],
    "total_track_energy_length": [30, 0, 1]
}

labels = {
    "no_tracks": ";No tracks",

    "n_objects": ";# objects;N.Entries / %i" % bin_size("n_objects"),
    "n_tracks": ";# tracks;N.Entries / %i" % bin_size("n_tracks"),
    "n_showers": ";# showers;N.Entries / %i" % bin_size("n_showers"),
    "track_theta": ";Track #theta [#circ];N. Entries / %.1f#circ" % bin_size("track_theta"),
    "track_phi": ";Track #phi [#circ];N. Entries / %.1f#circ" % bin_size("track_phi"),
    "shower_theta": ";Shower #theta [#circ];N. Entries / %.1f#circ" % bin_size("shower_theta"),
    "shower_phi": ";Shower #phi [#circ];N. Entries / %.1f#circ" % bin_size("shower_phi"),
    "shower_distance": ";Shower distance [cm];N. Entries / %.1f cm" % bin_size("shower_distance"),
    "track_distance": ";Track distance [cm];N. Entries / %.1f cm" % bin_size("track_distance"),
    "track_shower_angle": ";cos#theta [#circ];N. Entries / %.2f" % bin_size("track_shower_angle"),
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
    "reco_energy": ";Reco. energy [GeV]; N. Entries / 0.05 GeV",
    "shower_open_angle": ";Shower open angle [#circ]; N. Entries / %.1f#circ" % bin_size("shower_open_angle"),
    "dedx": ";Shower dE/dx [MeV/cm]; N. Entries / %.1f MeV/cm" % bin_size("dedx"),
    "dedx_y": ";Shower dE/dx (Y plane) [MeV/cm]; N. Entries / %.1f MeV/cm" % bin_size("dedx_y"),
    "numu_score": ";#nu_{#mu} selection score; N. Entries / %i" % bin_size("numu_score"),
    "category": ";category",
    "event_weight": ";event_weight",
    "event": ";event",
    "run": ";run",
    "subrun": ";subrun",
    "interaction_type": ";interaction_type",
    "is_signal": ";is_signal",
    "dedx_hits": ";dedx_hits",
    "shower_pca": ";Shower PCA;N. Entries / %.3f" % bin_size("shower_pca"),
    "track_pca": ";Track PCA;N. Entries / %.3f" % bin_size("track_pca"),
    "total_shower_energy": ";Total shower E [GeV]; N. Entries / %.2f GeV" % bin_size("total_shower_energy"),
    "total_shower_energy_y": ";Total shower E (Y plane) [GeV]; N. Entries / %.2f GeV" % bin_size("total_shower_energy_y"),
    "total_track_energy": ";Total track E [GeV]; N. Entries / %.2f GeV" % bin_size("total_track_energy"),
    "total_track_energy_dedx": ";Total track E (dE/dx) [GeV]; N. Entries / %.2f GeV" % bin_size("total_track_energy_dedx"),
    "shower_hits": ";Shower hits; N. Entries / %i" % bin_size("shower_hits"),
    "total_hits": ";Total hits; N. Entries / %i" % bin_size("shower_hits"),
    "total_hits_u": ";Total hits (U plane); N. Entries / %i" % bin_size("shower_hits"),
    "total_hits_v": ";Total hits (V plane); N. Entries / %i" % bin_size("shower_hits"),
    "total_hits_y": ";Total hits (Y plane); N. Entries / %i" % bin_size("shower_hits"),

    "shower_hits_y": ";Shower hits (Y plane); N. Entries / %i" % bin_size("shower_hits_y"),
    "track_hits": ";Track hits; N. Entries / %i" % bin_size("track_hits"),
    "track_energy": ";Track E [GeV]; N. Entries / %.2f GeV" % bin_size("track_energy"),
    "track_energy_dedx": ";Track E (dE/dx) [GeV]; N. Entries / %.2f GeV" % bin_size("track_energy_dedx"),
    "track_dedx": ";Track dE/dx [MeV/cm]; N. Entries / %.2f MeV/cm" % bin_size("track_dedx"),
    "track_energy_length": ";Track E (stopping power) [GeV]; N. Entries / %.2f GeV" % bin_size("track_energy_length"),
    "total_track_energy_length": ";Total track E (stopping power) [GeV]; N. Entries / %.2f GeV" % bin_size("total_track_energy_length")

}
