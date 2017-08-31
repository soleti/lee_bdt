from array import array

def find_interaction(dictionary,interaction):
    for name,id_int in dictionary.items():
        if id_int == interaction:
            return name

total_pot = 5e19

description = ["Other", "Cosmic", "Beam Intrinsic #nu_{e}", "Beam Intrinsic #nu_{#mu}", "Beam Intrinsic NC", "Dirt", "Cosmic contaminated"]

interactions = {
    "kQE":0,
    "kRes":1,
    "kDIS":2,
    "kCoh":3,
    "kNuanceOffset":1000,
    "k1eNp":1000+1,
    "kNCQE":1000+2,
    "kResCCNuProtonPiPlus":1000+3,
    "kResCCNuNeutronPi0":1000+4,
    "kResCCNuNeutronPiPlus":1000+5,
    "kResNCNuProtonPi0":1000+6,
    "kResNCNuProtonPiPlus":1000+7,
    "kResNCNuNeutronPi0":1000+8,
    "kResNCNuNeutronPiMinus":1000+9,
    "kResCCNuBarNeutronPiMinus":1000+10,
    "kResCCNuBarProtonPi0":1000+11,
    "kResCCNuBarProtonPiMinus":1000+12,
    "kResNCNuBarProtonPi0":1000+13,
    "kResNCNuBarProtonPiPlus":1000+14,
    "kResNCNuBarNeutronPi0":1000+15,
    "kResNCNuBarNeutronPiMinus":1000+16,
    "kResCCNuDeltaPlusPiPlus":1000+17,
    "kResCCNuDelta2PlusPiMinus":1000+21,
    "kResCCNuBarDelta0PiMinus":1000+28,
    "kResCCNuBarDeltaMinusPiPlus":1000+32,
    "kResCCNuProtonRhoPlus":1000+39,
    "kResCCNuNeutronRhoPlus":1000+41,
    "kResCCNuBarNeutronRhoMinus":1000+46,
    "kResCCNuBarNeutronRho0":1000+48,
    "kResCCNuSigmaPluskaonPlus":1000+53,
    "kResCCNuSigmaPluskaon0":1000+55,
    "kResCCNuBarSigmaMinuskaon0":1000+60,
    "kResCCNuBarSigma0kaon0":1000+62,
    "kResCCNuProtonEta":1000+67,
    "kResCCNuBarNeutronEta":1000+70,
    "kResCCNukaonPlusLambda0":1000+73,
    "kResCCNuBarkaon0Lambda0":1000+76,
    "kResCCNuProtonPiPlusPiMinus":1000+79,
    "kResCCNuProtonPi0Pi0":1000+80,
    "kResCCNuBarNeutronPiPlusPiMinus":1000+85,
    "kResCCNuBarNeutronPi0Pi0":1000+86,
    "kResCCNuBarProtonPi0Pi0":1000+90,
    "kCCDIS":1000+91,
    "kNCDIS":1000+92,
    "kUnUsed1":1000+93,
    "kUnUsed2":1000+94,
    "k1eNpHyperon":1000+95,
    "kNCCOH":1000+96,
    "kCCCOH":1000+97,
    "kNuElectronElastic":1000+98,
    "kInverseMuDecay":1000+99
}


x_start = 0
x_end = 256.35
y_start = -116.5
y_end = 116.5
z_start = 0
z_end = 1036.8

bdt_cut = -10.05

track_length = array("f", [ 0 ] )
track_theta = array("f", [ 0 ] )
track_phi = array("f", [ 0 ] )
shower_theta = array("f", [ 0 ] )
shower_phi = array("f", [ 0 ] )
shower_energy = array("f", [ 0 ] )
pt = array("f", [ 0 ] )
n_tracks = array("f", [ 0 ] )
n_showers = array("f", [ 0 ] )
track_shower_angle = array("f", [ 0 ] )
track_id = array("f", [ 0 ] )
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
reco_energy = array("f", [ 0 ] )
event_weight = array("f", [ 0 ] )
category = array("f", [ 0 ] )
event = array("f", [0])
run = array("f", [0])
subrun = array("f", [0])
interaction_type = array("f", [0])
is_signal = array("f", [ 0 ] )
shower_open_angle = array("f", [0])
dedx = array("f", [0])

spectators = [
    ("category", category),
    ("event_weight", event_weight),
    ("event", event),
    ("run", run),
    ("subrun", subrun),
    ("interaction_type", interaction_type),
    ("is_signal", is_signal)

]

variables = [
    ("track_length",track_length),
    ("track_theta",track_theta),
    ("track_phi",track_phi),
    ("shower_energy",shower_energy),
    ("shower_theta",shower_theta),
    ("shower_phi",shower_phi),
    ("pt",pt),
    ("n_tracks",n_tracks),
    ("n_showers",n_showers),
    ("track_shower_angle",track_shower_angle),
    ("proton_score",proton_score),
    ("shower_distance",shower_distance),
    ("track_distance",track_distance),
    ("track_start_y",track_start_y),
    ("track_end_y",track_end_y),
    ("track_start_x",track_start_x),
    ("track_end_x",track_end_x),
    ("track_start_z",track_start_z),
    ("track_end_z",track_end_z),
    ("shower_start_y",shower_start_y),
    ("shower_start_x",shower_start_x),
    ("shower_start_z",shower_start_z),
    ("reco_energy", reco_energy),
    ("shower_open_angle", shower_open_angle),
    ("dedx",dedx)
]
labels = {
    "n_tracks":";# tracks;N.Entries / 1",
    "n_showers":";# showers;N.Entries / 1",
    "track_theta":";Track #theta [#circ];N. Entries / 20#circ",
    "track_phi":";Track #phi [#circ];N. Entries / 40#circ",
    "shower_theta":";Shower #theta [#circ];N. Entries / 20#circ",
    "shower_phi":";Shower #phi [#circ];N. Entries / 40#circ",
    "shower_distance":";Shower distance [cm];N. Entries / 1 cm",
    "track_distance":";Track distance [cm];N. Entries / 1 cm",
    "track_shower_angle":";cos#theta [#circ];N. Entries / 0.2",
    "track_start_y":";Track start y [cm];",
    "track_start_z":";Track start z [cm];",
    "track_start_x":";Track start x [cm];",
    "track_end_y":";Track end y [cm];",
    "track_end_z":";Track end z [cm];",
    "track_end_x":";Track end x [cm];",
    "shower_start_y":";Shower start y [cm]",
    "shower_start_z":";Shower start z [cm]",
    "shower_start_x":";Shower start x [cm]",
    "shower_end_y":";Shower end y [cm]",
    "shower_end_z":";Shower end z [cm]",
    "shower_end_x":";Shower end x [cm]",
    "track_length":";Track length [cm];N. Entries / 2 cm",
    "proton_score":";Proton score; N. Entries / 0.1",
    "shower_energy":";Shower energy [GeV]; N. Entries / 0.1 GeV",
    "pt":";p_{t} [GeV/c];N. Entries / 0.1 GeV/c",
    "reco_energy":";Reco. energy [GeV]; N. Entries / 0.1 GeV",
    "shower_open_angle":";Shower open angle [#circ]; N. Entries / 2#circ",
    "dedx":";dE/dx [MeV/cm]; N. Entries / 0.2 MeV/cm"
}

binning = {
    "n_tracks":[5,1,6],
    "n_showers":[5,1,6],
    "track_theta":[9,0,180],
    "track_phi":[9,-180,180],
    "shower_theta":[9,0,180],
    "shower_phi":[9,-180,180],
    "shower_distance":[10,0,10],
    "track_distance":[10,0,10],
    "track_shower_angle":[10,-1,1],
    "track_start_y":[10,y_start,y_end],
    "track_start_z":[10,z_start,z_end],
    "track_start_x":[10,x_start,x_end],
    "track_end_y":[10,y_start,y_end],
    "track_end_z":[10,z_start,z_end],
    "track_end_x":[10,x_start,x_end],
    "shower_start_y":[10,y_start,y_end],
    "shower_start_z":[10,z_start,z_end],
    "shower_start_x":[10,x_start,x_end],
    "shower_end_y":[10,y_start,y_end],
    "shower_end_z":[10,z_start,z_end],
    "shower_end_x":[10,x_start,x_end],
    "track_length":[10,0,100],
    "proton_score":[10,0,1],
    "shower_energy":[20,0,2],
    "pt":[20,0,2],
    "reco_energy":[19,0.1,2],
    "shower_open_angle":[45,0,90],
    "dedx":[40,0,8]
}
