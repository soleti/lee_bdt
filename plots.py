#!/usr/bin/env python3.4

from ROOT import TChain, TH1F, TH2F, gStyle, TCanvas, THStack, TLegend, gPad, TEfficiency, TCut
from ROOT import TLine, TVector3, TPaveText, TTree, TFile
from ROOT import kRed, kGreen, kBlue, kOrange, kGray
from array import array
from glob import glob
import math

bnb_cosmic = glob("nu_files/*/*.root")
nue_cosmic = glob("nue_files/*/*.root")

chain = TChain("robertoana/pandoratree")
chain_pot = TChain("robertoana/pot")
chain_numu = TChain("UBXSec/tree")

chain_nue = TChain("robertoana/pandoratree")
chain_nue_pot = TChain("robertoana/pot")

gStyle.SetOptStat(0)
gStyle.SetPalette(87)
gStyle.SetNumberContours(99)


# In[2]:


for f in bnb_cosmic:
    chain.Add(f)
    chain_pot.Add(f)
    chain_numu.Add(f)
    chain_numu_proton.Add(f)

for f in nue_cosmic:
    chain_nue.Add(f)
    chain_nue_pot.Add(f)
    chain_nue_proton.Add(f)

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

total_pot = 6.6e20


# In[3]:


numu_selected_events = {}
for i in range(chain_numu.GetEntries()):
    chain_numu.GetEntry(i)
    numu_selected_events[chain_numu.event] = chain_numu.muon_is_reco


# In[4]:


h_true_reco_e = TH2F("h_true_reco_e",";True energy [GeV];Reco. energy [GeV]",18,0.2,2,18,0.2,2)
h_diff = TH1F("h_diff",";(True energy - Reco. energy)/True energy;N.Entries / 0.1 GeV", 30,-1,2)

e_energy = TEfficiency("e_energy",";#nu_{e} energy [GeV];Efficiency",20,0,2)
e_energy_perfect = TEfficiency("e_energy_perfect",";#nu_{e} energy [GeV];Efficiency",20,0,2)
e_energy_track = TEfficiency("e_energy_track",";#nu_{e} energy [GeV];Efficiency",20,0,2)
e_energy_shower = TEfficiency("e_energy_shower",";#nu_{e} energy [GeV];Efficiency",20,0,2)

h_vx_diff = TH1F("h_vx_diff",";#Deltax [cm];N.Entries / 0.5 cm",80,-20,20)
h_vy_diff = TH1F("h_vy_diff",";#Deltay [cm];N.Entries / 0.5 cm",80,-20,20)
h_vz_diff = TH1F("h_vz_diff",";#Deltaz [cm];N.Entries / 0.5 cm",80,-20,20)

h_vx_diff_cut = TH1F("h_vx_diff_cut",";#Deltax [cm];N.Entries / 0.5 cm",80,-20,20)
h_vy_diff_cut = TH1F("h_vy_diff_cut",";#Deltay [cm];N.Entries / 0.5 cm",80,-20,20)
h_vz_diff_cut = TH1F("h_vz_diff_cut",";#Deltaz [cm];N.Entries / 0.5 cm",80,-20,20)

h_energy = THStack("h_energy",";Reco. energy [GeV];N.Entries / 0.1 GeV")

h_cosmic = TH1F("h_cosmic",";Reco. energy [GeV];N.Entries / 0.1 GeV",20,0,2)
h_dirt = TH1F("h_dirt",";Reco. energy [GeV];N.Entries / 0.1 GeV",20,0,2)
h_nu_e = TH1F("h_nu_e",";Reco. energy [GeV];N.Entries / 0.1 GeV",20,0,2)
h_nu_mu = TH1F("h_nu_mu",";Reco. energy [GeV];N.Entries / 0.1 GeV",20,0,2)
h_other = TH1F("h_other",";Reco. energy [GeV];N.Entries / 0.1 GeV",20,0,2)
h_nc = TH1F("h_nc",";Reco. energy [GeV];N.Entries / 0.1 GeV",20,0,2)

h_energies = [h_other,h_cosmic,h_nu_e,h_nu_mu,h_nc,h_dirt]

h_n_tracks_other = TH1F("h_n_tracks_other",";# tracks;N.Entries / 1",5,1,6)
h_n_tracks_cosmic = TH1F("h_n_tracks_cosmic",";# tracks;N.Entries / 1",5,1,6)
h_n_tracks_nu_e = TH1F("h_n_tracks_nu_e",";# tracks;N.Entries / 1",5,1,6)
h_n_tracks_nu_mu = TH1F("h_n_tracks_nu_mu",";# tracks;N.Entries / 1",5,1,6)
h_n_tracks_nc = TH1F("h_n_tracks_nc",";# tracks;N.Entries / 1",5,1,6)
h_n_tracks_dirt = TH1F("h_n_tracks_dirt",";# tracks;N.Entries / 1",5,1,6)

h_n_trackss = [h_n_tracks_other,h_n_tracks_cosmic,h_n_tracks_nu_e,h_n_tracks_nu_mu,h_n_tracks_nc,h_n_tracks_dirt]

h_n_showers_other = TH1F("h_n_showers_other",";# showers;N.Entries / 1",5,1,6)
h_n_showers_cosmic = TH1F("h_n_showers_cosmic",";# showers;N.Entries / 1",5,1,6)
h_n_showers_nu_e = TH1F("h_n_showers_nu_e",";# showers;N.Entries / 1",5,1,6)
h_n_showers_nu_mu = TH1F("h_n_showers_nu_mu",";# showers;N.Entries / 1",5,1,6)
h_n_showers_nc = TH1F("h_n_showers_nc",";# showers;N.Entries / 1",5,1,6)
h_n_showers_dirt = TH1F("h_n_showers_dirt",";# showers;N.Entries / 1",5,1,6)

h_n_showerss = [h_n_showers_other,h_n_showers_cosmic,h_n_showers_nu_e,h_n_showers_nu_mu,h_n_showers_nc,h_n_showers_dirt]

h_vy_other = TH1F("h_vy_other",";y [cm];N. Entries / 13 cm", 10,-130,130)
h_vy_cosmic = TH1F("h_vy_cosmic",";y [cm];N. Entries / 13 cm", 10,-130,130)
h_vy_nu_e = TH1F("h_vy_nu_e",";y [cm];N. Entries / 13 cm", 10,-130,130)
h_vy_nu_mu = TH1F("h_vy_nu_mu",";y [cm];N. Entries / 13 cm", 10,-130,130)
h_vy_nc = TH1F("h_vy_nc",";y [cm];N. Entries / 13 cm", 10,-130,130)
h_vy_dirt = TH1F("h_vy_dirt",";y [cm];N. Entries / 13 cm", 10,-130,130)

h_vys = [h_vy_other,h_vy_cosmic,h_vy_nu_e,h_vy_nu_mu,h_vy_nc,h_vy_dirt]

h_trkz_other = TH1F("h_trkz_other",";cos#theta;N. Entries / 0.1", 20,-1,1)
h_trkz_cosmic = TH1F("h_trkz_cosmic",";cos#theta;N. Entries / 0.1", 20,-1,1)
h_trkz_nu_e = TH1F("h_trkz_nu_e",";cos#theta;N. Entries / 0.1", 20,-1,1)
h_trkz_nu_mu = TH1F("h_trkz_nu_mu",";cos#theta;N. Entries / 0.1", 20,-1,1)
h_trkz_nc = TH1F("h_trkz_nc",";cos#theta;N. Entries / 0.1", 20,-1,1)
h_trkz_dirt = TH1F("h_trkz_dirt",";cos#theta;N. Entries / 0.1", 20,-1,1)

h_trkzs = [h_trkz_other,h_trkz_cosmic,h_trkz_nu_e,h_trkz_nu_mu,h_trkz_nc,h_trkz_dirt]

h_shwrz_other = TH1F("h_shwrz_other",";cos#theta;N. Entries / 0.1", 20,-1,1)
h_shwrz_cosmic = TH1F("h_shwrz_cosmic",";cos#theta;N. Entries / 0.1", 20,-1,1)
h_shwrz_nu_e = TH1F("h_shwrz_nu_e",";cos#theta;N. Entries / 0.1", 20,-1,1)
h_shwrz_nu_mu = TH1F("h_shwrz_nu_mu",";cos#theta;N. Entries / 0.1", 20,-1,1)
h_shwrz_nc = TH1F("h_shwrz_nc",";cos#theta;N. Entries / 0.1", 20,-1,1)
h_shwrz_dirt = TH1F("h_shwrz_dirt",";cos#theta;N. Entries / 0.1", 20,-1,1)

h_shwrzs = [h_shwrz_other,h_shwrz_cosmic,h_shwrz_nu_e,h_shwrz_nu_mu,h_shwrz_nc,h_shwrz_dirt]

h_trklen_other = TH1F("h_trklen_other",";L [cm];N. Entries / 5 cm", 24,0,120)
h_trklen_cosmic = TH1F("h_trklen_cosmic",";L [cm];N. Entries / 5 cm", 24,0,120)
h_trklen_nu_e = TH1F("h_trklen_nu_e",";L [cm];N. Entries / 5 cm", 24,0,120)
h_trklen_nu_mu = TH1F("h_trklen_nu_mu",";L [cm];N. Entries / 5 cm", 24,0,120)
h_trklen_nc = TH1F("h_trklen_nc",";L [cm];N. Entries / 5 cm", 24,0,120)
h_trklen_dirt = TH1F("h_trklen_dirt",";L [cm];N. Entries / 5 cm", 24,0,120)

h_trklens = [h_trklen_other,h_trklen_cosmic,h_trklen_nu_e,h_trklen_nu_mu,h_trklen_nc,h_trklen_dirt]

h_shower_e_other = TH1F("h_shower_e_other",";Energy [GeV];N. Entries / 0.1 GeV", 20,0,1)
h_shower_e_cosmic = TH1F("h_shower_e_cosmic",";Energy [GeV];N. Entries / 0.1 GeV", 20,0,1)
h_shower_e_nu_e = TH1F("h_shower_e_nu_e",";Energy [GeV];N. Entries / 0.1 GeV", 20,0,1)
h_shower_e_nu_mu = TH1F("h_shower_e_nu_mu",";Energy [GeV];N. Entries / 0.1 GeV", 20,0,1)
h_shower_e_nc = TH1F("h_shower_e_nc",";Energy [GeV];N. Entries / 0.1 GeV", 20,0,1)
h_shower_e_dirt = TH1F("h_shower_e_dirt",";Energy [GeV];N. Entries / 0.1 GeV", 20,0,1)

h_shower_es = [h_shower_e_other,h_shower_e_cosmic,h_shower_e_nu_e,h_shower_e_nu_mu,h_shower_e_nc,h_shower_e_dirt]

h_theta_other = TH1F("h_theta_other",";#theta [#circ];N. Entries / 5#circ", 18,0,180)
h_theta_cosmic = TH1F("h_theta_cosmic",";#theta [#circ];N. Entries / 5#circ", 18,0,180)
h_theta_nu_e = TH1F("h_theta_nu_e",";#theta [#circ];N. Entries / 5#circ", 18,0,180)
h_theta_nu_mu = TH1F("h_theta_nu_mu",";#theta [#circ];N. Entries / 5#circ", 18,0,180)
h_theta_nc = TH1F("h_theta_nc",";#theta [#circ];N. Entries / 5#circ", 18,0,180)
h_theta_dirt = TH1F("h_theta_dirt",";#theta [#circ];N. Entries / 5#circ", 18,0,180)

h_thetas = [h_theta_other,h_theta_cosmic,h_theta_nu_e,h_theta_nu_mu,h_theta_nc,h_theta_dirt]

h_phi_other = TH1F("h_phi_other",";#phi [#circ];N. Entries / 10#circ", 18,-180,180)
h_phi_cosmic = TH1F("h_phi_cosmic",";#phi [#circ];N. Entries / 10#circ", 18,-180,180)
h_phi_nu_e = TH1F("h_phi_nu_e",";#phi [#circ];N. Entries / 10#circ", 18,-180,180)
h_phi_nu_mu = TH1F("h_phi_nu_mu",";#phi [#circ];N. Entries / 10#circ", 18,-180,180)
h_phi_nc = TH1F("h_phi_nc",";#phi [#circ];N. Entries / 10#circ", 18,-180,180)
h_phi_dirt = TH1F("h_phi_dirt",";#phi [#circ];N. Entries / 10#circ", 18,-180,180)

h_phis = [h_phi_other,h_phi_cosmic,h_phi_nu_e,h_phi_nu_mu,h_phi_nc,h_phi_dirt]

h_track_theta_other = TH1F("h_track_theta_other",";#theta [#circ];N. Entries / 5#circ", 18,0,180)
h_track_theta_cosmic = TH1F("h_track_theta_cosmic",";#theta [#circ];N. Entries / 5#circ", 18,0,180)
h_track_theta_nu_e = TH1F("h_track_theta_nu_e",";#theta [#circ];N. Entries / 5#circ", 18,0,180)
h_track_theta_nu_mu = TH1F("h_track_theta_nu_mu",";#theta [#circ];N. Entries / 5#circ", 18,0,180)
h_track_theta_nc = TH1F("h_track_theta_nc",";#theta [#circ];N. Entries / 5#circ", 18,0,180)
h_track_theta_dirt = TH1F("h_track_theta_dirt",";#theta [#circ];N. Entries / 5#circ", 18,0,180)

h_track_thetas = [h_track_theta_other,h_track_theta_cosmic,h_track_theta_nu_e,h_track_theta_nu_mu,h_track_theta_nc,h_track_theta_dirt]

h_track_phi_other = TH1F("h_track_phi_other",";#phi [#circ];N. Entries / 10#circ", 18,-180,180)
h_track_phi_cosmic = TH1F("h_track_phi_cosmic",";#phi [#circ];N. Entries / 10#circ", 18,-180,180)
h_track_phi_nu_e = TH1F("h_track_phi_nu_e",";#phi [#circ];N. Entries / 10#circ", 18,-180,180)
h_track_phi_nu_mu = TH1F("h_track_phi_nu_mu",";#phi [#circ];N. Entries / 10#circ", 18,-180,180)
h_track_phi_nc = TH1F("h_track_phi_nc",";#phi [#circ];N. Entries / 10#circ", 18,-180,180)
h_track_phi_dirt = TH1F("h_track_phi_dirt",";#phi [#circ];N. Entries / 10#circ", 18,-180,180)

h_track_phis = [h_track_phi_other,h_track_phi_cosmic,h_track_phi_nu_e,h_track_phi_nu_mu,h_track_phi_nc,h_track_phi_dirt]


h_pt_other = TH1F("h_pt_other",";p_{t} [GeV/c^{2}];N. Entries / 0.1 GeV/c^{2}", 10,0,1)
h_pt_cosmic = TH1F("h_pt_cosmic",";p_{t} [GeV/c^{2}];N. Entries / 0.1 GeV/c^{2}", 10,0,1)
h_pt_nu_e = TH1F("h_pt_nu_e",";p_{t} [GeV/c^{2}];N. Entries / 0.1 GeV/c^{2}", 10,0,1)
h_pt_nu_mu = TH1F("h_pt_nu_mu",";p_{t} [GeV/c^{2}];N. Entries / 0.1 GeV/c^{2}", 10,0,1)
h_pt_nc = TH1F("h_pt_nc",";p_{t} [GeV/c^{2}];N. Entries / 0.1 GeV/c^{2}", 10,0,1)
h_pt_dirt = TH1F("h_pt_dirt",";p_{t} [GeV/c^{2}];N. Entries / 0.1 GeV/c^{2}", 10,0,1)

h_pts = [h_pt_other,h_pt_cosmic,h_pt_nu_e,h_pt_nu_mu,h_pt_nc,h_pt_dirt]

h_distance_other = TH1F("h_distance_other",";Distance [cm];N. Entries / 10 cm", 80,0,800)
h_distance_cosmic = TH1F("h_distance_cosmic",";Distance [cm];N. Entries / 10 cm", 80,0,800)
h_distance_nu_e = TH1F("h_distance_nu_e",";Distance [cm];N. Entries / 10 cm", 80,0,800)
h_distance_nu_mu = TH1F("h_distance_nu_mu",";Distance [cm];N. Entries / 10 cm", 80,0,800)
h_distance_nc = TH1F("h_distance_nc",";Distance [cm];N. Entries / 10 cm", 80,0,800)
h_distance_dirt = TH1F("h_distance_dirt",";Distance [cm];N. Entries / 10 cm", 80,0,800)

h_distances = [h_distance_other,h_distance_cosmic,h_distance_nu_e,h_distance_nu_mu,h_distance_nc,h_distance_dirt]

h_costheta_other = TH1F("h_costheta_other",";cos#theta [#circ];N. Entries / 0.2", 10,-1,1)
h_costheta_cosmic = TH1F("h_costheta_cosmic",";cos#theta [#circ];N. Entries / 0.2", 10,-1,1)
h_costheta_nu_e = TH1F("h_costheta_nu_e",";cos#theta [#circ];N. Entries / 0.2", 10,-1,1)
h_costheta_nu_mu = TH1F("h_costheta_nu_mu",";cos#theta [#circ];N. Entries / 0.2", 10,-1,1)
h_costheta_nc = TH1F("h_costheta_nc",";cos#theta [#circ];N. Entries / 0.2", 10,-1,1)
h_costheta_dirt = TH1F("h_costheta_dirt",";cos#theta [#circ];N. Entries / 0.2", 10,-1,1)

h_costhetas = [h_costheta_other,h_costheta_cosmic,h_costheta_nu_e,h_costheta_nu_mu,h_costheta_nc,h_costheta_dirt]

h_shower_track_angle_other = TH1F("h_shower_track_angle_other",";cos#theta [#circ];N. Entries / 0.2", 10,-1,1)
h_shower_track_angle_cosmic = TH1F("h_shower_track_angle_cosmic",";cos#theta [#circ];N. Entries / 0.2", 10,-1,1)
h_shower_track_angle_nu_e = TH1F("h_shower_track_angle_nu_e",";cos#theta [#circ];N. Entries / 0.2", 10,-1,1)
h_shower_track_angle_nu_mu = TH1F("h_shower_track_angle_nu_mu",";cos#theta [#circ];N. Entries / 0.2", 10,-1,1)
h_shower_track_angle_nc = TH1F("h_shower_track_angle_nc",";cos#theta [#circ];N. Entries / 0.2", 10,-1,1)
h_shower_track_angle_dirt = TH1F("h_shower_track_angle_dirt",";cos#theta [#circ];N. Entries / 0.2", 10,-1,1)

h_shower_track_angles = [h_shower_track_angle_other,h_shower_track_angle_cosmic,h_shower_track_angle_nu_e,h_shower_track_angle_nu_mu,h_shower_track_angle_nc,h_shower_track_angle_dirt]

kinematic_plots = [h_distances,h_shower_track_angles,h_costhetas,h_pts,h_phis,h_thetas,h_shower_es,h_trklens,h_shwrzs,h_trkzs,h_vys,h_n_trackss,h_n_showerss,h_track_thetas,h_track_phis]

colors = [kGray+2, kRed - 3, kGreen - 2, kBlue - 5, kBlue - 9, kOrange+3]
description = ["Other", "Cosmic", "Beam Intrinsic #nu_{e}", "Beam Intrinsic #nu_{#mu}", "Beam Intrinsic NC", "Dirt"]

l_energy = TLegend(0.48,0.55,0.84,0.84)
l_plots = TLegend(0.48,0.55,0.84,0.84)

passed = 0
total = 0
perfect = 0
evt_shower = 0
evt_track = 0
perfect_tracks = 0
perfect_showers = 0
flash_passed = 0
no_photons = 0


# In[5]:

def fill_kin_branches(root_chain, weight, variables):
    longest_track = 0
    longest_track_id = 0

    for itrk in range(root_chain.n_tracks):
        if root_chain.track_len[itrk] > longest_track:
            longest_track = root_chain.track_len[itrk]
            longest_track_id = itrk

    most_energetic_shower = 0
    most_energetic_shower_id = 0
    for ish in range(root_chain.n_showers):
        if root_chain.shower_energy[ish] > most_energetic_shower:
            most_energetic_shower = root_chain.shower_energy[ish]
            most_energetic_shower_id = ish


    v_track = TVector3(root_chain.track_dir_x[longest_track_id],root_chain.track_dir_y[longest_track_id],root_chain.track_dir_z[longest_track_id])
    v_shower = TVector3(root_chain.shower_dir_x[most_energetic_shower_id],root_chain.shower_dir_y[most_energetic_shower_id],root_chain.shower_dir_z[most_energetic_shower_id])
    costheta_shower_track = v_track.Dot(v_shower)/(v_track.Mag()*v_shower.Mag())

    signal = 0
    if root_chain.category == 2: signal = 1

    shower_vertex = [root_chain.shower_start_x[most_energetic_shower_id],root_chain.shower_start_y[most_energetic_shower_id],root_chain.shower_start_z[most_energetic_shower_id]]
    neutrino_vertex = [root_chain.vx,root_chain.vy,root_chain.vz]
    shower_vertex_d = math.sqrt(sum([(s-n)**2 for s,n in zip(shower_vertex,neutrino_vertex)]))

    variables["is_signal"][0] = signal
    variables["track_length"][0] = longest_track
    variables["track_phi"][0] = root_chain.track_phi[longest_track_id]
    variables["track_theta"][0] = root_chain.track_theta[longest_track_id]
    variables["track_z"][0] = root_chain.track_dir_z[longest_track_id]
    variables["shower_energy"][0] = most_energetic_shower
    variables["shower_theta"][0] = root_chain.shower_theta[most_energetic_shower_id]
    variables["shower_phi"][0] = root_chain.shower_phi[most_energetic_shower_id]
    variables["shower_z"][0] = root_chain.shower_dir_z[most_energetic_shower_id]
    variables["shower_distance"][0] = shower_vertex_d
    variables["reco_energy"][0] = root_chain.E
    variables["category"][0] = root_chain.category
    variables["event_weight"][0] = weight
    variables["pt"][0] = pt_plot(root_chain)

    variables["n_tracks"][0] = root_chain.n_tracks
    variables["n_showers"][0] = root_chain.n_showers
    variables["track_id"][0] = root_chain.track_id[longest_track_id]

    variables["track_shower_angle"][0] = costheta_shower_track

    variables["event"][0] = root_chain.event
    variables["run"][0] = root_chain.run
    variables["subrun"][0] = root_chain.subrun


# In[6]:

def kin_plot(root_chain, weight):
    h_n_trackss[root_chain.category].Fill(root_chain.n_tracks)
    h_n_showerss[root_chain.category].Fill(root_chain.n_showers)

    longest_track = max([root_chain.track_len[itrk] for itrk in range(root_chain.n_tracks)])
    most_z_track = max([root_chain.track_dir_z[itrk] for itrk in range(root_chain.n_tracks)])
    most_z_shower = max([root_chain.shower_dir_z[ish] for ish in range(root_chain.n_showers)])


    #if longest_track < 50 and (most_z_track < -0.8 or most_z_track > 0.2) and (most_z_shower < -0.8 or most_z_shower > 0.2):
    h_energies[root_chain.category].Fill(root_chain.E, weight)

    if root_chain.n_tracks == 1:
        v_track = TVector3(root_chain.track_dir_x[0],root_chain.track_dir_y[0],root_chain.track_dir_z[0])
        for ish in range(root_chain.n_showers):
            v_shower = TVector3(root_chain.shower_dir_x[ish],root_chain.shower_dir_y[ish],root_chain.shower_dir_z[ish])
            costheta_shower_track = v_track.Dot(v_shower)/(v_track.Mag()*v_shower.Mag())
            h_shower_track_angles[root_chain.category].Fill(costheta_shower_track,weight)


    h_vys[root_chain.category].Fill(root_chain.vy, weight)
    h_pts[root_chain.category].Fill(pt_plot(root_chain), weight)

    costheta_plot(root_chain, weight)

    for itrk in range(root_chain.n_tracks):
        h_track_thetas[root_chain.category].Fill(math.degrees(root_chain.track_theta[itrk]), weight)
        h_track_phis[root_chain.category].Fill(math.degrees(root_chain.track_phi[itrk]), weight)
        h_trklens[root_chain.category].Fill(root_chain.track_len[itrk], weight)
        h_trkzs[root_chain.category].Fill(root_chain.track_dir_z[itrk], weight)

    neutrino_vertex = [root_chain.vx,root_chain.vy,root_chain.vz]

    for ish in range(root_chain.n_showers):
        h_thetas[root_chain.category].Fill(math.degrees(root_chain.shower_theta[ish]), weight)
        h_phis[root_chain.category].Fill(math.degrees(root_chain.shower_phi[ish]), weight)
        h_shwrzs[root_chain.category].Fill(root_chain.shower_dir_z[ish], weight)
        shower_vertex = [root_chain.shower_start_x[ish],root_chain.shower_start_y[ish],root_chain.shower_start_z[ish]]
        shower_vertex_d = math.sqrt(sum([(s-n)**2 for s,n in zip(shower_vertex,neutrino_vertex)]))
        h_distances[root_chain.category].Fill(shower_vertex_d)

    h_shower_es[root_chain.category].Fill(sum([root_chain.shower_energy[ish] for ish in range(root_chain.n_showers)]), weight)


# In[7]:


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


# In[8]:


def pt_plot(root_chain):
    p_showers = []
    for ish in range(root_chain.n_showers):
        if root_chain.shower_energy[ish] > 0:
            p_vector = TVector3(root_chain.shower_dir_x[ish],root_chain.shower_dir_y[ish],root_chain.shower_dir_z[ish])
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


# In[9]:

is_signal = array("f", [ 0 ] )
reco_energy = array("f", [ 0 ] )
track_length = array("f", [ 0 ] )
track_theta = array("f", [ 0 ] )
track_phi = array("f", [ 0 ] )
track_z = array("f", [ 0 ] )
shower_theta = array("f", [ 0 ] )
shower_phi = array("f", [ 0 ] )
shower_energy = array("f", [ 0 ] )
shower_z = array("f", [ 0 ] )
event_weight = array("f", [ 0 ] )
category = array("f", [ 0 ] )
pt = array("f", [ 0 ] )
n_tracks = array("f", [ 0 ] )
n_showers = array("f", [ 0 ] )
track_shower_angle = array("f", [ 0 ] )
track_id = array("f", [ 0 ] )
event = array("f", [0])
run = array("f", [0])
subrun = array("f", [0])

shower_distance = array("f", [0])

variables = {"reco_energy":reco_energy, "track_length":track_length, "track_theta":track_theta, "track_phi":track_phi,
            "shower_theta":shower_theta,"shower_phi":shower_phi,"shower_energy":shower_energy,"shower_z":shower_z,
             "track_z":track_z,"event_weight":event_weight, "category":category,"is_signal":is_signal,"pt":pt,
             "n_tracks":n_tracks, "n_showers":n_showers, "track_shower_angle":track_shower_angle,"track_id":track_id,
             "event":event,"shower_distance":shower_distance,"run":run,"subrun":subrun}
kin_tree = TTree("kin_tree","kin_tree")

for n,b in variables.items():
    kin_tree.Branch(n,b,n+"/f")



# In[ ]:

# NU_E INTRINSIC + COSMIC SAMPLE
for i in range(chain_nue.GetEntries()):
    chain_nue.GetEntry(i)

    protons = sum(1 for i,pdg in enumerate(chain_nue.nu_daughters_pdg) if abs(pdg) == 2212)
    electrons = sum(1 for i,pdg in enumerate(chain_nue.nu_daughters_pdg) if abs(pdg) == 11)
    photons = sum(1 for i in chain_nue.nu_daughters_pdg if i == 22)
    pions = sum(1 for i in chain_nue.nu_daughters_pdg if abs(i) == 211 or abs(i) == 111)

    # MEASURE EFFICIENCY
    if electrons > 0 and photons == 0 and pions == 0 and protons > 0:
        no_photons+=1
        p = False
        p_track = False
        p_shower = False

        if protons == chain_nue.n_tracks and chain_nue.n_tracks == chain_nue.nu_matched_tracks: p_track = True

        if electrons == chain_nue.n_showers and chain_nue.n_showers == chain_nue.nu_matched_showers: p_shower = True

        if p_track and p_shower:
            p = True

        if chain_nue.true_nu_is_fiducial:
            total += 1

            if chain_nue.flash_passed:
                flash_passed += 1

            if chain_nue.shower_passed > 0:
                evt_shower += 1
                if chain_nue.track_passed > 0:
                    evt_track += 1

            if p_shower:
                perfect_showers += 1
            if p_track:
                perfect_tracks += 1

            passed+=1
            if p: perfect += 1

            e_energy_track.Fill(chain_nue.passed and p_track, chain_nue.nu_E)
            e_energy_shower.Fill(chain_nue.passed and p_shower, chain_nue.nu_E)
            e_energy.Fill(chain_nue.passed, chain_nue.nu_E)
            e_energy_perfect.Fill(chain_nue.passed and p, chain_nue.nu_E)

    if chain_nue.passed:
        h_vx_diff.Fill(chain_nue.vx-chain_nue.true_vx)
        h_vy_diff.Fill(chain_nue.vy-chain_nue.true_vy)
        h_vz_diff.Fill(chain_nue.vz-chain_nue.true_vz)
        h_vx_diff_cut.Fill(chain_nue.vx-chain_nue.true_vx_sce)
        h_vy_diff_cut.Fill(chain_nue.vy-chain_nue.true_vy_sce)
        h_vz_diff_cut.Fill(chain_nue.vz-chain_nue.true_vz_sce)

        if chain_nue.category != 1:
            h_diff.Fill((chain_nue.nu_E-chain_nue.E)/chain_nue.nu_E)
            h_true_reco_e.Fill(chain_nue.nu_E,chain_nue.E)

        proton = True
        for itrk in range(chain_nue.n_tracks):
            if chain_nue.predict_p[itrk] < 0.5:
                proton = False
                break

        if proton:
            kin_plot(chain_nue, total_pot/total_nue_pot)
            fill_kin_branches(chain_nue, total_pot/total_nue_pot, variables)
            kin_tree.Fill()


# In[ ]:

# BNB + COSMIC SAMPLE
for i in range(chain.GetEntries()):
    chain.GetEntry(i)
    numu_passed = False

    if chain.event in numu_selected_events:
        numu_passed = numu_selected_events[chain.event]



    if chain.passed and not numu_passed and abs(chain.nu_pdg) != 12:
        proton = True
        for itrk in range(chain.n_tracks):
            if chain.predict_p[itrk] < 0.5:
                proton = False
                break

        if proton:
            kin_plot(chain, total_pot/total_bnb_pot)
            fill_kin_branches(chain, total_pot/total_bnb_pot, variables)
            kin_tree.Fill()


# In[ ]:

kin_file = TFile("kin_file.root", "RECREATE")
kin_tree.Write()
kin_file.Close()

print(no_photons)
print("Entries", chain_nue.GetEntries())
print("Is fiducial", total)
print("Flash passed", total)
print("At least 1 shower", evt_shower)
print("At least 1 track", evt_track)
print("Perfect showers", perfect_showers)
print("Perfect tracks", perfect_tracks)

print (passed, perfect, total)


# In[ ]:


for i,histo in enumerate(h_energies):
    histo.SetLineColor(1)
    histo.SetLineWidth(2)
    histo.SetFillColor(colors[i])
    h_energy.Add(histo)
    l_energy.AddEntry(histo, "%s: %.0f events" % (description[i], histo.Integral(3,20)), "f")



# In[ ]:

bkg = sum([h.Integral() for h in h_energies if h != h_nu_e])
print("Signal/background",h_nu_e.Integral()/bkg)


# In[ ]:


for i,kind in enumerate(kinematic_plots):
    for j,histo in enumerate(kind):
        histo.SetLineWidth(2)
        histo.SetLineColor(colors[j])
        if j == 2:
            histo.SetFillStyle(3002)
            histo.SetFillColor(colors[j])

        if i == 0 and j == 2:
            l_plots.AddEntry(histo, "%s" % description[j], "f")
        elif i == 0 and j != 5:
            l_plots.AddEntry(histo, "%s" % description[j], "l")



# In[ ]:


pt = TPaveText(0.1,0.91,0.60,0.97, "ndc")
pt.AddText("MicroBooNE Preliminary 6.6e20 POT")
pt.SetFillColor(0)
pt.SetBorderSize(0)
pt.SetShadowColor(0)

pt2 = TPaveText(0.1,0.91,0.45,0.97, "ndc")
pt2.AddText("MicroBooNE Preliminary")
pt2.SetFillColor(0)
pt2.SetBorderSize(0)
pt2.SetShadowColor(0)


# In[ ]:


c_energy = TCanvas("c_energy")
h_energy.Draw("hist")
h_energy.GetXaxis().SetRangeUser(0.2,2)
l_energy.Draw()
c_energy.Update()
pt.Draw()
c_energy.SaveAs("plots/energy.pdf")
c_energy.Draw()


# In[ ]:


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


# In[ ]:


c_plots=[]
for i,kind in enumerate(kinematic_plots):
    c_plots.append(TCanvas("c_plot"+str(i)))


# In[ ]:

for i,kind in enumerate(kinematic_plots):
    c_plots[i].cd()
    kind[0].Draw("hist")
    for h in kind[1:5]:
        if h.Integral():
            h.Scale(kind[0].Integral()/h.Integral())
            h.Draw("hist same")
    l_plots.Draw()
    pt2.Draw()
    c_plots[i].Update()
    c_plots[i].SaveAs("plots/c%i.pdf"%i)
    c_plots[i].Draw()

input()

# In[ ]:




# In[ ]:
