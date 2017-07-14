#!/usr/bin/env python3.4

import ROOT
from array import array
f_data = ROOT.TFile("data_file.root")
t_data = f_data.Get("data_tree")

ROOT.TMVA.Tools.Instance()
reader = ROOT.TMVA.Reader(":".join([
    "!V",
    "!Silent",
    "Color",]))

track_length = array("f", [ 0 ] )
track_theta = array("f", [ 0 ] )
track_phi = array("f", [ 0 ] )
track_z = array("f", [ 0 ] )
shower_theta = array("f", [ 0 ] )
shower_phi = array("f", [ 0 ] )
shower_energy = array("f", [ 0 ] )
shower_z = array("f", [ 0 ] )
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

reco_energy = array("f", [ 0 ] )
event_weight = array("f", [ 0 ] )
category = array("f", [ 0 ] )
event = array("f", [0])
run = array("f", [0])
subrun = array("f", [0])
interaction_type = array("f", [0])

t_data.SetBranchAddress("track_length", track_length)
t_data.SetBranchAddress("track_theta",track_theta)
t_data.SetBranchAddress("track_phi",track_phi)
t_data.SetBranchAddress("shower_energy",shower_energy)
t_data.SetBranchAddress("shower_theta",shower_theta)
t_data.SetBranchAddress("shower_phi",shower_phi)
t_data.SetBranchAddress("pt",pt)
t_data.SetBranchAddress("n_tracks",n_tracks)
t_data.SetBranchAddress("n_showers",n_showers)
t_data.SetBranchAddress("track_shower_angle",track_shower_angle)
t_data.SetBranchAddress("proton_score",proton_score)
t_data.SetBranchAddress("shower_distance",shower_distance)
t_data.SetBranchAddress("track_distance",track_distance)
t_data.SetBranchAddress("track_start_y",track_start_y)
t_data.SetBranchAddress("track_end_y",track_end_y)
t_data.SetBranchAddress("track_start_x",track_start_x)
t_data.SetBranchAddress("track_end_x",track_end_x)
t_data.SetBranchAddress("track_start_z",track_start_z)
t_data.SetBranchAddress("track_end_z",track_end_z)

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
    ("track_end_z",track_end_z)
]

for name, var in variables:
    reader.AddVariable(name, var)

spectators = [
    ("reco_energy", reco_energy),
    ("category", category),
    ("event_weight", event_weight),
    ("event", event),
    ("run", run),
    ("subrun", subrun),
    ("interaction_type", interaction_type)
]

for name, var in spectators:
    reader.AddSpectator(name, var)

reader.BookMVA("BDT method","dataset/weights/TMVAClassification_BDT.weights.xml")

for i in range(t_data.GetEntries()):
    t_data.GetEntry(i)
    BDT_response = reader.EvaluateMVA("BDT method")
    print(BDT_response)
