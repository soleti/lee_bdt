#!/usr/bin/env python3.4

from ROOT import TChain, TH1F, TH2F, gStyle, TCanvas, THStack, TLegend, gPad, TCut
from ROOT import TLine, TVector3, TPaveText, TTree, TFile
from ROOT import kRed, kGreen, kBlue, kOrange, kGray, kWhite
from array import array
from glob import glob
import math


bnb_cosmic = glob("nu_files_fidvol/*/*.root")
nue_cosmic = glob("nue_files_fidvol/*/*.root")
data_bnb = glob("data_files_bnb_mcc81/*/*.root")
data_bnbext = glob("data_files_bnbext_mcc81/*/*.root")

chain = TChain("robertoana/pandoratree")
chain_pot = TChain("robertoana/pot")
chain_numu = TChain("UBXSec/tree")

chain_nue = TChain("robertoana/pandoratree")
chain_nue_pot = TChain("robertoana/pot")

chain_data_bnb = TChain("robertoana/pandoratree")
chain_data_bnb_pot = TChain("robertoana/pot")
chain_data_bnb_numu = TChain("UBXSec/tree")

chain_data_bnbext = TChain("robertoana/pandoratree")
chain_data_bnbext_pot = TChain("robertoana/pot")
chain_data_bnbext_numu = TChain("UBXSec/tree")


gStyle.SetOptStat(0)
gStyle.SetPalette(87)
gStyle.SetNumberContours(99)

for f in bnb_cosmic:
    chain.Add(f)
    chain_pot.Add(f)
    chain_numu.Add(f)

for f in data_bnb:
    chain_data_bnb.Add(f)
    chain_data_bnb_pot.Add(f)
    chain_data_bnb_numu.Add(f)

for f in data_bnbext:
    chain_data_bnbext.Add(f)
    chain_data_bnbext_pot.Add(f)
    chain_data_bnbext_numu.Add(f)

for f in nue_cosmic:
    chain_nue.Add(f)
    chain_nue_pot.Add(f)

print(chain_data_bnb.GetEntries(),chain_data_bnbext.GetEntries())


total_bnb_pot = 0
for i in range(chain_pot.GetEntries()):
    chain_pot.GetEntry(i)
    total_bnb_pot += chain_pot.pot

print("Total POT BNB", total_bnb_pot)

total_nue_pot = 0
for i in range(chain_nue_pot.GetEntries()):
    chain_nue_pot.GetEntry(i)
    total_nue_pot += chain_nue_pot.pot
print("Total POT v_e", total_nue_pot)

total_data_bnb_pot = 0
for i in range(chain_data_bnb_pot.GetEntries()):
    chain_data_bnb_pot.GetEntry(i)
    total_data_bnb_pot += chain_data_bnb_pot.pot
total_data_bnb_pot *= 1e12
print("Total data BNB POT", total_data_bnb_pot)

total_data_bnbext_pot = 0
for i in range(chain_data_bnbext_pot.GetEntries()):
    chain_data_bnbext_pot.GetEntry(i)
    total_data_bnbext_pot += chain_data_bnbext_pot.pot
total_data_bnbext_pot *= 1e12
print("Total data BNB EXT POT", total_data_bnbext_pot)

total_pot = 6.6e20

numu_selected_events = {}
for i in range(chain_numu.GetEntries()):
    chain_numu.GetEntry(i)
    numu_selected_events[chain_numu.event] = chain_numu.muon_is_reco

numu_selected_events_data = {}
for i in range(chain_data_bnb_numu.GetEntries()):
    chain_data_bnb_numu.GetEntry(i)
    numu_selected_events_data[chain_data_bnb_numu.event] = chain_data_bnb_numu.muon_is_reco

numu_selected_events_data_ext = {}
for i in range(chain_data_bnbext_numu.GetEntries()):
    chain_data_bnbext_numu.GetEntry(i)
    numu_selected_events_data_ext[chain_data_bnbext_numu.event] = chain_data_bnbext_numu.muon_is_reco

h_n_candidates = TH1F("h_n_candidates",";# candidates;N. Entries / 1",5,0,5)
h_true_reco_e = TH2F("h_true_reco_e",";True energy [GeV];Reco. energy [GeV]",18,0.2,2,18,0.2,2)
h_diff = TH1F("h_diff",";(True energy - Reco. energy)/True energy;N.Entries / 0.1 GeV", 30,-1,2)

h_vx_diff = TH1F("h_vx_diff",";#Deltax [cm];N.Entries / 0.5 cm",80,-20,20)
h_vy_diff = TH1F("h_vy_diff",";#Deltay [cm];N.Entries / 0.5 cm",80,-20,20)
h_vz_diff = TH1F("h_vz_diff",";#Deltaz [cm];N.Entries / 0.5 cm",80,-20,20)

h_vx_diff_cut = TH1F("h_vx_diff_cut",";#Deltax [cm];N.Entries / 0.5 cm",80,-20,20)
h_vy_diff_cut = TH1F("h_vy_diff_cut",";#Deltay [cm];N.Entries / 0.5 cm",80,-20,20)
h_vz_diff_cut = TH1F("h_vz_diff_cut",";#Deltaz [cm];N.Entries / 0.5 cm",80,-20,20)

colors = [kGray+2, kRed - 3, kGreen - 2, kBlue - 5, kBlue - 9, kOrange+3, kWhite, kWhite]
description = ["Other", "Cosmic", "Beam Intrinsic #nu_{e}", "Beam Intrinsic #nu_{#mu}", "Beam Intrinsic NC", "Dirt", "Data BNB - BNB EXT"]

l_energy = TLegend(0.455,0.53,0.705,0.85)
l_plots = TLegend(0.48,0.55,0.84,0.84)


def fill_kin_branches(root_chain, weight, variables):

    if root_chain.category > 6: print(root_chain.event)

    longest_track = 0
    longest_track_id = 0
    most_proton_track_id = 0
    most_proton_track = 0
    most_proton_track_length = 0
    for itrk in range(root_chain.n_tracks):
        if root_chain.predict_p[itrk] > most_proton_track:
            most_proton_track = root_chain.predict_p[itrk]
            most_proton_track_id = itrk
            most_proton_track_length = root_chain.track_len[itrk]
        if root_chain.track_len[itrk] > longest_track:
            longest_track = root_chain.track_len[itrk]
            longest_track_id = itrk

    most_energetic_shower = 0
    most_energetic_shower_id = 0
    for ish in range(root_chain.n_showers):
        if root_chain.shower_energy[ish] < 3:
            if root_chain.shower_energy[ish] > most_energetic_shower:
                most_energetic_shower = root_chain.shower_energy[ish]
                most_energetic_shower_id = ish


    v_track = TVector3(root_chain.track_dir_x[longest_track_id],root_chain.track_dir_y[longest_track_id],root_chain.track_dir_z[longest_track_id])
    v_shower = TVector3(root_chain.shower_dir_x[most_energetic_shower_id],root_chain.shower_dir_y[most_energetic_shower_id],root_chain.shower_dir_z[most_energetic_shower_id])
    costheta_shower_track = v_track.Dot(v_shower)/(v_track.Mag()*v_shower.Mag())

    signal = 0
    if root_chain.category == 2: signal = 1

    track_vertex = [root_chain.track_start_x[longest_track_id],root_chain.track_start_y[longest_track_id],root_chain.track_start_z[longest_track_id]]
    shower_vertex = [root_chain.shower_start_x[most_energetic_shower_id],root_chain.shower_start_y[most_energetic_shower_id],root_chain.shower_start_z[most_energetic_shower_id]]
    neutrino_vertex = [root_chain.vx,root_chain.vy,root_chain.vz]
    true_neutrino_vertex = [root_chain.true_vx_sce,root_chain.true_vy_sce,root_chain.true_vz_sce]

    shower_vertex_d = math.sqrt(sum([(s-n)**2 for s,n in zip(shower_vertex,neutrino_vertex)]))
    track_vertex_d = math.sqrt(sum([(t-n)**2 for t,n in zip(track_vertex,neutrino_vertex)]))

    variables["is_signal"][0] = signal
    variables["track_length"][0] = most_proton_track_length
    variables["track_phi"][0] = math.degrees(root_chain.track_phi[most_proton_track_id])
    variables["track_theta"][0] = math.degrees(root_chain.track_theta[most_proton_track_id])
    variables["shower_energy"][0] = most_energetic_shower
    variables["shower_theta"][0] = math.degrees(root_chain.shower_theta[most_energetic_shower_id])
    variables["shower_phi"][0] = math.degrees(root_chain.shower_phi[most_energetic_shower_id])
    variables["shower_distance"][0] = shower_vertex_d
    variables["track_distance"][0] = track_vertex_d

    variables["track_start_x"][0] = root_chain.track_start_x[most_proton_track_id]
    variables["track_end_x"][0] = root_chain.track_end_x[most_proton_track_id]

    variables["track_start_y"][0] = root_chain.track_start_y[most_proton_track_id]
    variables["track_end_y"][0] = root_chain.track_end_y[most_proton_track_id]

    variables["track_start_z"][0] = root_chain.track_start_z[most_proton_track_id]
    variables["track_end_z"][0] = root_chain.track_end_z[most_proton_track_id]

    variables["shower_start_x"][0] = root_chain.shower_start_x[most_energetic_shower_id]
    #variables["shower_end_x"][0] = root_chain.shower_end_x[most_energetic_shower_id]

    variables["shower_start_y"][0] = root_chain.shower_start_y[most_energetic_shower_id]
    #variables["shower_end_y"][0] = root_chain.shower_end_y[most_energetic_shower_id]

    variables["shower_start_z"][0] = root_chain.shower_start_z[most_energetic_shower_id]
    #variables["shower_end_z"][0] = root_chain.shower_end_z[most_energetic_shower_id]

    variables["reco_energy"][0] = root_chain.E
    variables["category"][0] = root_chain.category
    variables["event_weight"][0] = weight
    variables["pt"][0] = pt_plot(root_chain)

    variables["n_tracks"][0] = root_chain.n_tracks
    variables["n_showers"][0] = root_chain.n_showers
    variables["track_shower_angle"][0] = costheta_shower_track

    variables["event"][0] = root_chain.event
    variables["run"][0] = root_chain.run
    variables["subrun"][0] = root_chain.subrun
    variables["proton_score"][0] = root_chain.predict_p[most_proton_track_id]
    variables["interaction_type"][0] = root_chain.interaction_type


def costheta_plot(root_chain, weight=1):
    electrons = sum(1 for i,pdg in enumerate(root_chain.nu_daughters_pdg) if abs(pdg) == 11)

    if electrons == 1:
        e_dir = []

        for i, pdg in enumerate(root_chain.nu_daughters_pdg):
            if abs(pdg) == 11:
                p = math.sqrt(root_chain.nu_daughters_px[i]**2+root_chain.nu_daughters_py[i]**2+root_chain.nu_daughters_pz[i]**2)
                e_dir.append([root_chain.nu_daughters_px[i]/p, root_chain.nu_daughters_py[i]/p, root_chain.nu_daughters_pz[i]/p])

        for ish in range(root_chain.n_showers):
            v_reco = TVector3(root_chain.shower_dir_x[ish], root_chain.shower_dir_y[ish], root_chain.shower_dir_z[ish])
            v_true = TVector3(e_dir[0][0], e_dir[0][1], e_dir[0][2])
            costheta = v_reco.Dot(v_true)/(v_reco.Mag()*v_true.Mag())
            h_costhetas[root_chain.category].Fill(costheta, weight)


def pt_plot(root_chain):
    p_showers = []
    for ish in range(root_chain.n_showers):
        if root_chain.shower_energy[ish] > 0:
            p_vector = TVector3(root_chain.shower_dir_x[ish],root_chain.shower_dir_y[ish],root_chain.shower_dir_z[ish])
            if (root_chain.shower_energy[ish] < 10):
                p_vector.SetMag(math.sqrt((root_chain.shower_energy[ish]+0.052)**2-0.052**2))
            p_showers.append(p_vector)

    p_tracks = []
    for itr in range(root_chain.n_tracks):
        if root_chain.track_energy[itr] > 0:
            p_vector = TVector3(root_chain.track_dir_x[itr], root_chain.track_dir_y[itr], root_chain.track_dir_z[itr])
            p_vector.SetMag(math.sqrt((root_chain.track_energy[itr]+0.938)**2 - 0.938**2))
            p_tracks.append(p_vector)

    p_track_sum = TVector3()
    if len(p_tracks) > 0:
        p_track_sum = p_tracks[0]
        for i in p_tracks[1:]:
            p_track_sum += i

    p_shower_sum = TVector3()
    if len(p_showers) > 0:
        p_shower_sum = p_showers[0]
        for i in p_showers[1:]:
            p_shower_sum += i

    pt = (p_track_sum+p_shower_sum).Perp()
    return pt


is_signal = array("f", [ 0 ] )
reco_energy = array("f", [ 0 ] )
track_length = array("f", [ 0 ] )
track_theta = array("f", [ 0 ] )
track_phi = array("f", [ 0 ] )

shower_theta = array("f", [ 0 ] )
shower_phi = array("f", [ 0 ] )
shower_energy = array("f", [ 0 ] )

event_weight = array("f", [ 0 ] )
category = array("f", [ 0 ] )
pt = array("f", [ 0 ] )
n_tracks = array("f", [ 0 ] )
n_showers = array("f", [ 0 ] )
track_shower_angle = array("f", [ 0 ] )
event = array("f", [0])
run = array("f", [0])
subrun = array("f", [0])
shower_distance = array("f", [0])
interaction_type = array("f", [0])
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
variables = {"reco_energy":reco_energy,
    "track_length":track_length,
    "track_theta":track_theta,
    "track_phi":track_phi,
    "shower_theta":shower_theta,
    "shower_phi":shower_phi,
    "shower_energy":shower_energy,
    "event_weight":event_weight,
    "category":category,
    "is_signal":is_signal,
    "pt":pt,
    "n_tracks":n_tracks,
    "n_showers":n_showers,
    "track_shower_angle":track_shower_angle,
    "event":event,
    "shower_distance":shower_distance,
    "run":run,
    "subrun":subrun,
    "interaction_type":interaction_type,
    "proton_score":proton_score,
    "track_distance":track_distance,
    "track_start_y":track_start_y,
    "track_end_y":track_end_y,
    "track_start_x":track_start_x,
    "track_end_x":track_end_x,
    "track_start_z":track_start_z,
    "track_end_z":track_end_z,
    "shower_start_y":shower_start_y,
    "shower_start_x":shower_start_x,
    "shower_start_z":shower_start_z,
}

kin_tree = TTree("kin_tree","kin_tree")

bnb_tree = TTree("bnb_tree","bnb_tree")
bnbext_tree = TTree("bnbext_tree","bnbext_tree")

for n,b in variables.items():
    kin_tree.Branch(n,b,n+"/f")
    bnb_tree.Branch(n,b,n+"/f")
    bnbext_tree.Branch(n,b,n+"/f")

# NU_E INTRINSIC + COSMIC SAMPLE
for i in range(chain_nue.GetEntries()):
    chain_nue.GetEntry(i)

    if chain_nue.passed:
        h_vx_diff.Fill(chain_nue.vx-chain_nue.true_vx)
        h_vy_diff.Fill(chain_nue.vy-chain_nue.true_vy)
        h_vz_diff.Fill(chain_nue.vz-chain_nue.true_vz)
        h_vx_diff_cut.Fill(chain_nue.vx-chain_nue.true_vx_sce)
        h_vy_diff_cut.Fill(chain_nue.vy-chain_nue.true_vy_sce)
        h_vz_diff_cut.Fill(chain_nue.vz-chain_nue.true_vz_sce)
        h_n_candidates.Fill(chain_nue.n_candidates)

        if chain_nue.category != 1:
            h_diff.Fill((chain_nue.nu_E-chain_nue.E)/chain_nue.nu_E)
            h_true_reco_e.Fill(chain_nue.nu_E,chain_nue.E)

        proton = True
        for itrk in range(chain_nue.n_tracks):
            if chain_nue.predict_p[itrk] < -0.75:
                proton = False
                break

        tracks_contained = min([chain_nue.track_is_fiducial[i] for i in range(chain_nue.n_tracks)])
        showers_contained = min([chain_nue.shower_is_fiducial[i] for i in range(chain_nue.n_showers)])
        if proton:
            fill_kin_branches(chain_nue, total_pot/total_nue_pot, variables)
            kin_tree.Fill()


# BNB + COSMIC SAMPLE
for i in range(chain.GetEntries()):
    chain.GetEntry(i)
    numu_passed = False

    if chain.event in numu_selected_events:
        numu_passed = numu_selected_events[chain.event]

    if chain.passed and not numu_passed and abs(chain.nu_pdg) != 12:
        proton = True
        for itrk in range(chain.n_tracks):
            if chain.predict_p[itrk] < -0.75:
                proton = False
                break

        tracks_contained = min([chain.track_is_fiducial[i] for i in range(chain.n_tracks)])
        showers_contained = min([chain.shower_is_fiducial[i] for i in range(chain.n_showers)])

        if proton:
            fill_kin_branches(chain, total_pot/total_bnb_pot, variables)
            kin_tree.Fill()
            h_n_candidates.Fill(chain.n_candidates)

for i in range(chain_data_bnb.GetEntries()):
    chain_data_bnb.GetEntry(i)
    numu_passed = False

    if chain_data_bnb.event in numu_selected_events_data:
        numu_passed = numu_selected_events_data[chain_data_bnb.event]

    tracks_contained = 0
    showers_contained = 0
    if chain_data_bnb.n_tracks > 0:
        tracks_contained = min([chain_data_bnb.track_is_fiducial[i] for i in range(chain_data_bnb.n_tracks)])
        showers_contained = min([chain_data_bnb.shower_is_fiducial[i] for i in range(chain_data_bnb.n_showers)])

    if chain_data_bnb.passed and not numu_passed:
        fill_kin_branches(chain_data_bnb,total_pot/total_data_bnb_pot,variables)
        bnb_tree.Fill()

for i in range(chain_data_bnbext.GetEntries()):
    chain_data_bnbext.GetEntry(i)
    numu_passed = False

    if chain_data_bnbext.event in numu_selected_events_data_ext:
        numu_passed = numu_selected_events_data_ext[chain_data_bnbext.event]

    tracks_contained = 0
    showers_contained = 0
    if chain_data_bnbext.n_tracks > 0:
        tracks_contained = min([chain_data_bnbext.track_is_fiducial[i] for i in range(chain_data_bnbext.n_tracks)])
        showers_contained = min([chain_data_bnbext.shower_is_fiducial[i] for i in range(chain_data_bnbext.n_showers)])

    if chain_data_bnbext.passed and not numu_passed:
        fill_kin_branches(chain_data_bnbext,1.2300 * (382718/chain_data_bnbext.GetEntries()) * (chain_data_bnb.GetEntries()/547616)*total_pot/5e19,variables)
        bnbext_tree.Fill()


kin_file = TFile("kin_file.root", "RECREATE")
kin_tree.Write()
kin_file.Close()

bnb_file = TFile("bnb_file.root", "RECREATE")
bnb_tree.Write()
bnb_file.Close()

bnbext_file = TFile("bnbext_file.root", "RECREATE")
bnbext_tree.Write()
bnbext_file.Close()

pt2 = TPaveText(0.1,0.91,0.45,0.97, "ndc")
pt2.AddText("MicroBooNE Preliminary")
pt2.SetFillColor(0)
pt2.SetBorderSize(0)
pt2.SetShadowColor(0)


c_reco_true = TCanvas("c_reco_true")
h_true_reco_e.Draw("colz")
line = TLine(0.2,0.2,2,2)
line.SetLineStyle(2)
line.SetLineWidth(3)
line.SetLineColor(kRed+1)
line.Draw()
c_reco_true.Update()
pt2.Draw()
c_reco_true.SaveAs("plots/reco_true.pdf")
c_reco_true.Draw()

input()
