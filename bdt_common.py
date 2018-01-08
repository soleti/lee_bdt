from array import array
import math
import ROOT

bdt, manual = False, True
total_data_bnb_pot = 4.758e19


def fill_histos(tree_name, bdt, manual):
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

    variables_dict = dict(variables + spectators)

    histograms = []

    for i, n in enumerate(variables_dict.keys()):
        if n != "reco_energy":
            h = ROOT.TH1F("h_%s" % n, labels[n],
                          binning[n][0], binning[n][1], binning[n][2])
        else:
            bins = array("f", [0.2, 0.25, 0.3, 0.35, 0.4,
                               0.45, 0.5, 0.6, 0.8, 1])
            h = ROOT.TH1F("h_%s" % n, labels[n], len(bins) - 1, bins)
        histograms.append(h)

    histo_dict = dict(zip(variables_dict.keys(), histograms))

    h_bdt = ROOT.TH1F("h_bdt_%s" % tree_name,
                      "BDT response; N. Entries / 0.05", 40, -1, 1)
    passed_events = 0

    for i in range(t_data.GetEntries()):
        t_data.GetEntry(i)
        BDT_response = reader.EvaluateMVA("BDT method")
        h_bdt.Fill(BDT_response, t_data.event_weight)

        if bdt:
            apply_bdt = BDT_response > bdt_cut
        else:
            apply_bdt = True

        if manual:
            apply_manual = manual_cuts(t_data)
        else:
            apply_manual = True

        if apply_bdt and apply_manual and 0.1 < t_data.reco_energy < 1:
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
    shower_energy = chain.shower_energy > 0.2
    dedx = 1.4 < chain.dedx < 3
    shower_distance = chain.shower_distance < 2.5
    track_distance = chain.track_distance < 2.5
    proton_score = chain.proton_score > 0.98
    open_angle = 2 < chain.shower_open_angle < 15
    shower_theta = chain.shower_theta < 90
    n_tracks = chain.n_tracks < 3
    shower_pca = chain.shower_pca < 0.985
    track_shower_angle = -0.9 < chain.track_shower_angle < 0.9
    track_theta = chain.track_theta < 130

    dedx = 0.6 < chain.dedx < 3.07
    proton_score = chain.proton_score > 0.57
    shower_distance = chain.shower_distance < 2.23
    track_distance = chain.track_distance < 9.55
    open_angle = 1 < chain.shower_open_angle < 23.7
    shower_theta = 15 < chain.shower_theta < 99


    cuts = [dedx, proton_score, open_angle,
            shower_theta, shower_distance, track_distance]

    # optimized_cuts = [dedx, proton_score, shower_distance, track_distance,
    #                   shower_open_angle, pt, shower_theta, track_length]

    passed = len(cuts) == sum(cuts)
    # passed_optimized = len(optimized_cuts) == sum(optimized_cuts)
    return passed


def sigmaCalc(h_signal, h_background, sys_err=0):
    chi2 = sum(
        [h_signal.GetBinContent(i)**2 /
         (h_background.GetBinContent(i) +
          math.pow(sys_err * h_background.GetBinContent(i), 2))
         for i in range(1, h_signal.GetNbinsX() - 1)
         if h_background.GetBinContent(i) > 0])

    return math.sqrt(chi2)


total_pot = total_data_bnb_pot

description = ["Beam Intrinsic #nu_{e}",
               "Cosmic in-time",
               "Cosmic",
               "Cosmic contaminated",
               "Beam Intrinsic #nu_{#mu}",
               "Beam Intrinsic NC",
               "Dirt",
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

bdt_cut = 0.47
track_length = array("f", [0])
track_theta = array("f", [0])
track_phi = array("f", [0])
shower_theta = array("f", [0])
shower_phi = array("f", [0])
shower_energy = array("f", [0])
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
total_track_energy = array("f", [0])

spectators = [
    ("category", category),
    ("event_weight", event_weight),
    ("event", event),
    ("run", run),
    ("subrun", subrun),
    ("interaction_type", interaction_type),
    ("is_signal", is_signal),
    ("dedx_hits", dedx_hits),
    ("shower_energy", shower_energy),
    ("reco_energy", reco_energy),

    ("pt", pt),
    ("track_shower_angle", track_shower_angle),

    ("track_theta", track_theta),
    ("track_start_y", track_start_y),
    ("track_end_y", track_end_y),
    ("track_start_x", track_start_x),
    ("track_end_x", track_end_x),
    ("track_start_z", track_start_z),
    ("track_end_z", track_end_z),
    ("shower_start_y", shower_start_y),
    ("shower_start_x", shower_start_x),
    ("shower_start_z", shower_start_z),
    ("track_length", track_length),
    ("track_phi", track_phi),
    ("shower_phi", shower_phi),
    ("n_tracks", n_tracks),
    ("n_showers", n_showers),
    ("shower_pca", shower_pca),
    ("track_pca", track_pca),
    ("total_track_energy", total_track_energy),
    ("numu_score", numu_score)
]

variables = [
    ("dedx", dedx),
    ("proton_score", proton_score),
    ("shower_distance", shower_distance),
    ("track_distance", track_distance),
    ("shower_open_angle", shower_open_angle),
    ("shower_theta", shower_theta),
]

labels = {
    "n_tracks": ";# tracks;N.Entries / 1",
    "n_showers": ";# showers;N.Entries / 1",
    "track_theta": ";Track #theta [#circ];N. Entries / 20#circ",
    "track_phi": ";Track #phi [#circ];N. Entries / 40#circ",
    "shower_theta": ";Shower #theta [#circ];N. Entries / 20#circ",
    "shower_phi": ";Shower #phi [#circ];N. Entries / 40#circ",
    "shower_distance": ";Shower distance [cm];N. Entries / 1 cm",
    "track_distance": ";Track distance [cm];N. Entries / 1 cm",
    "track_shower_angle": ";cos#theta [#circ];N. Entries / 0.2",
    "track_start_y": ";Track start y [cm];",
    "track_start_z": ";Track start z [cm];",
    "track_start_x": ";Track start x [cm];",
    "track_end_y": ";Track end y [cm];",
    "track_end_z": ";Track end z [cm];",
    "track_end_x": ";Track end x [cm];",
    "shower_start_y": ";Shower start y [cm]",
    "shower_start_z": ";Shower start z [cm]",
    "shower_start_x": ";Shower start x [cm]",
    "shower_end_y": ";Shower end y [cm]",
    "shower_end_z": ";Shower end z [cm]",
    "shower_end_x": ";Shower end x [cm]",
    "track_length": ";Track length [cm];N. Entries / 2 cm",
    "proton_score": ";Proton score; N. Entries / 0.025",
    "shower_energy": ";Shower energy [GeV]; N. Entries / 0.1 GeV",
    "pt": ";p_{t} [GeV/c];N. Entries / 0.1 GeV/c",
    "reco_energy": ";Reco. energy [GeV]; N. Entries / 0.05 GeV",
    "shower_open_angle": ";Shower open angle [#circ]; N. Entries / 2#circ",
    "dedx": ";dE/dx [MeV/cm]; N. Entries / 0.3 MeV/cm",
    "numu_score": ";#nu_{#mu} selection score; N. Entries / 0.01",
    "category": ";category",
    "event_weight": ";event_weight",
    "event": ";event",
    "run": ";run",
    "subrun": ";subrun",
    "interaction_type": ";interaction_type",
    "is_signal": ";is_signal",
    "dedx_hits": ";dedx_hits",
    "shower_pca": ";Shower PCA;N. Entries / 0.025",
    "track_pca": ";Track PCA;N. Entries / 0.025",
    "total_shower_energy": ";Total shower E [GeV]; N. Entries / 0.025 GeV",
    "total_track_energy": ";Total track E [GeV]; N. Entries / 0.025 GeV"

}

binning = {
    "n_tracks": [5, 1, 6],
    "n_showers": [5, 1, 6],
    "track_theta": [36, 0, 180],
    "track_phi": [36, -180, 180],
    "shower_theta": [36, 0, 180],
    "shower_phi": [36, -180, 180],
    "shower_distance": [50, 0, 10],
    "track_distance": [50, 0, 10],
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
    "track_length": [40, 0, 200],
    "proton_score": [40, 0.9, 1],
    "shower_energy": [40, 0, 1],
    "pt": [20, 0, 2],
    "reco_energy": [40, 0, 2],
    "shower_open_angle": [46, 0, 46],
    "dedx": [38, 0.3, 6],
    "numu_score": [20, 0, 1],
    "category": [7, 0, 7],
    "event_weight": [20, 0, 100],
    "event": [20, 0, 10000],
    "run": [20, 0, 10],
    "subrun": [20, 0, 1000],
    "interaction_type": [100, 1000, 1100],
    "is_signal": [2, 0, 1],
    "dedx_hits": [200, 0, 1],
    "shower_pca": [50, 0.9, 1],
    "track_pca": [50, 0.9, 1],
    "total_shower_energy": [40, 0, 0.2],
    "total_track_energy": [40, 0, 0.2]

}
