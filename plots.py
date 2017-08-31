#!/usr/bin/env python3.4

from ROOT import TChain, TH1F, TH2F, gStyle, TCanvas, THStack, TLegend, gPad, TCut
from ROOT import TLine, TVector3, TPaveText, TTree, TFile
from ROOT import kRed, kGreen, kBlue, kOrange, kGray, kWhite
from array import array
from glob import glob
import math
from bdt_common import total_pot, variables, spectators, x_start, x_end, y_start, y_end, z_start, z_end
import random

gStyle.SetOptStat(0)
gStyle.SetPalette(87)
gStyle.SetNumberContours(99)

bnb_cosmic = glob("nu_files_6_42_energy/*/*.root")
nue_cosmic = glob("nue_files_6_42_energy/*/*.root")
data_bnb = glob("data_files_bnb_6_42_energy/*/*.root")
data_bnbext = glob("data_files_bnbext_6_42_energy/*/*.root")

chain = TChain("robertoana/pandoratree")
chain_pot = TChain("robertoana/pot")
chain_numu = TChain("UBXSec/tree")

chain_nue = TChain("robertoana/pandoratree")
chain_nue_pot = TChain("robertoana/pot")
chain_nue_numu = TChain("UBXSec/tree")

chain_data_bnb = TChain("robertoana/pandoratree")
chain_data_bnb_numu = TChain("UBXSec/tree")

chain_data_bnbext = TChain("robertoana/pandoratree")
chain_data_bnbext_numu = TChain("UBXSec/tree")

def is_fiducial(point):
    fidvol_cut = 50
    ok_y = point[1] > y_start+fidvol_cut and point[1] < y_end-fidvol_cut
    ok_x = point[0] > x_start+fidvol_cut and point[0] < x_end-fidvol_cut
    ok_z = point[2] > z_start+fidvol_cut and point[2] < z_end-fidvol_cut
    return ok_y and ok_x and ok_z

for f in bnb_cosmic:
    chain.Add(f)
    chain_pot.Add(f)
    chain_numu.Add(f)

for f in nue_cosmic:
    chain_nue.Add(f)
    chain_nue_numu.Add(f)
    chain_nue_pot.Add(f)

total_bnb_pot = 0
for i in range(chain_pot.GetEntries()):
    chain_pot.GetEntry(i)
    total_bnb_pot += chain_pot.pot
print("Total POT BNB {0:.2e}".format(total_bnb_pot))

total_nue_pot = 0
for i in range(chain_nue_pot.GetEntries()):
    chain_nue_pot.GetEntry(i)
    total_nue_pot += chain_nue_pot.pot
print("Total POT v_e {0:.2e}".format(total_nue_pot))

total_data_bnb_pot = 3.54e19
print("Total data POT BNB {0:.2e}".format(total_data_bnb_pot))

data_ext_scaling_factor = 1.279
print("Data EXT scaling factor {0:.2f}".format(data_ext_scaling_factor))

h_n_candidates = TH1F("h_n_candidates",";# candidates;N. Entries / 1",5,0,5)
h_true_reco_e = TH2F("h_true_reco_e",";True energy [GeV];Reco. energy [GeV]",18,0.2,2,18,0.2,2)
h_diff = TH1F("h_diff",";(True energy - Reco. energy)/True energy;N.Entries / 0.1 GeV", 30,-1,2)

colors = [kGray+2, kRed - 3, kGreen - 2, kBlue - 5, kBlue - 9, kOrange+3, kWhite, kWhite, kRed-3]

def fill_kin_branches(root_chain, weight, variables):
    longest_track = 0
    longest_track_id = 0
    most_proton_track_id = 0
    most_proton_track = 1
    most_proton_track_length = 0
    for itrk in range(root_chain.n_tracks):
        if root_chain.predict_p[itrk] < most_proton_track:
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

    track_vertex = [root_chain.track_start_x[most_proton_track_id],root_chain.track_start_y[most_proton_track_id],root_chain.track_start_z[most_proton_track_id]]
    track_end = [root_chain.track_end_x[most_proton_track_id],root_chain.track_end_y[most_proton_track_id],root_chain.track_end_z[most_proton_track_id]]
    shower_vertex = [root_chain.shower_start_x[most_energetic_shower_id],root_chain.shower_start_y[most_energetic_shower_id],root_chain.shower_start_z[most_energetic_shower_id]]
    neutrino_vertex = [root_chain.vx,root_chain.vy,root_chain.vz]
    true_neutrino_vertex = [root_chain.true_vx_sce,root_chain.true_vy_sce,root_chain.true_vz_sce]

    shower_vertex_d = math.sqrt(sum([(s-n)**2 for s,n in zip(shower_vertex,neutrino_vertex)]))
    track_vertex_d = math.sqrt(sum([(t-n)**2 for t,n in zip(track_vertex,neutrino_vertex)]))
    track_end_d = math.sqrt(sum([(t-n)**2 for t,n in zip(track_end,neutrino_vertex)]))

    track_shower_d = math.sqrt(sum([(s-t)**2 for s,t in zip(shower_vertex,track_vertex)]))
    trackend_shower_d = math.sqrt(sum([(s-t)**2 for s,t in zip(shower_vertex,track_end)]))

    direction = 1
    if trackend_shower_d < track_shower_d: direction = -1
    theta = math.acos(direction*root_chain.track_dir_z[most_proton_track_id])

    variables["is_signal"][0] = signal
    variables["track_length"][0] = root_chain.track_energy[most_proton_track_id]
    variables["track_phi"][0] = math.degrees(root_chain.track_phi[most_proton_track_id])
    variables["track_theta"][0] = math.degrees(theta)
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

    variables["shower_open_angle"][0] = math.degrees(root_chain.shower_open_angle[most_energetic_shower_id])

    dedx = root_chain.shower_dEdx[most_energetic_shower_id][2]

    variables["dedx"][0] = dedx


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

variables = dict(variables+spectators)

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

    if chain_nue.passed and chain_numu.muon_is_reco != 1:

        if chain_nue.category != 1:
            h_diff.Fill((chain_nue.nu_E-chain_nue.E)/chain_nue.nu_E)
            h_true_reco_e.Fill(chain_nue.nu_E,chain_nue.E)

        tracks_contained = min([chain_nue.track_is_fiducial[i] for i in range(chain_nue.n_tracks)])
        showers_contained = min([chain_nue.shower_is_fiducial[i] for i in range(chain_nue.n_showers)])

        track_fidvol = True
        for i in range(chain_nue.n_tracks):
            track_start = [chain_nue.track_start_x[i], chain_nue.track_start_y[i], chain_nue.track_start_z[i]]
            track_end = [chain_nue.track_end_x[i], chain_nue.track_end_y[i], chain_nue.track_end_z[i]]
            track_fidvol = track_fidvol and is_fiducial(track_start) and is_fiducial(track_end)

        shower_fidvol = True
        for i in range(chain_nue.n_showers):
            shower_start = [chain_nue.shower_start_x[i], chain_nue.shower_start_y[i], chain_nue.shower_start_z[i]]
            shower_fidvol = shower_fidvol and is_fiducial(shower_start)

        if track_fidvol and shower_fidvol:
            fill_kin_branches(chain_nue, total_pot/total_nue_pot*chain_nue.bnbweight, variables)
            kin_tree.Fill()
            h_n_candidates.Fill(chain_nue.n_candidates)



# BNB + COSMIC SAMPLE
for i in range(chain.GetEntries()):
    chain.GetEntry(i)
    chain_numu.GetEntry(i)

    if chain.passed and abs(chain.nu_pdg) != 12 and chain_numu.muon_is_reco != 1:
        track_fidvol = True
        for i in range(chain.n_tracks):
            track_start = [chain.track_start_x[i], chain.track_start_y[i], chain.track_start_z[i]]
            track_end = [chain.track_end_x[i], chain.track_end_y[i], chain.track_end_z[i]]
            track_fidvol = track_fidvol and is_fiducial(track_start) and is_fiducial(track_end)

        shower_fidvol = True
        for i in range(chain.n_showers):
            shower_start = [chain.shower_start_x[i], chain.shower_start_y[i], chain.shower_start_z[i]]
            shower_fidvol = shower_fidvol and is_fiducial(shower_start)

        if track_fidvol and shower_fidvol:
            fill_kin_branches(chain, total_pot/total_bnb_pot, variables)
            kin_tree.Fill()
            h_n_candidates.Fill(chain.n_candidates)

run_subrun_list = []
run_subrun_file = open("run_subrun_bnb.txt","w")
for i in range(chain_data_bnb.GetEntries()):
    chain_data_bnb.GetEntry(i)

    run_subrun = "%i %i" % (chain_data_bnb.run, chain_data_bnb.subrun)
    if run_subrun not in run_subrun_list:
        run_subrun_list.append(run_subrun)
        print(run_subrun, file=run_subrun_file)

    if chain_data_bnb.passed and chain_data_bnb_numu.muon_is_reco != 1:
        track_fidvol = True
        for i in range(chain_data_bnb.n_tracks):
            track_start = [chain_data_bnb.track_start_x[i], chain_data_bnb.track_start_y[i], chain_data_bnb.track_start_z[i]]
            track_end = [chain_data_bnb.track_end_x[i], chain_data_bnb.track_end_y[i], chain_data_bnb.track_end_z[i]]
            track_fidvol = track_fidvol and is_fiducial(track_start) and is_fiducial(track_end)

        shower_fidvol = True
        for i in range(chain_data_bnb.n_showers):
            shower_start = [chain_data_bnb.shower_start_x[i], chain_data_bnb.shower_start_y[i], chain_data_bnb.shower_start_z[i]]
            shower_fidvol = shower_fidvol and is_fiducial(shower_start)

        if track_fidvol and shower_fidvol:
            fill_kin_branches(chain_data_bnb,total_pot/total_data_bnb_pot,variables)
            bnb_tree.Fill()

run_subrun_list_ext = []
run_subrun_file_ext = open("run_subrun_ext.txt","w")
for i in range(chain_data_bnbext.GetEntries()):
    chain_data_bnbext.GetEntry(i)

    run_subrun = "%i %i" % (chain_data_bnbext.run, chain_data_bnbext.subrun)
    if run_subrun not in run_subrun_list_ext:
        run_subrun_list_ext.append(run_subrun)
        print(run_subrun, file=run_subrun_file_ext)

    if chain_data_bnbext.passed and chain_data_bnbext_numu.muon_is_reco != 1:
        track_fidvol = True
        for i in range(chain_data_bnbext.n_tracks):
            track_start = [chain_data_bnbext.track_start_x[i], chain_data_bnbext.track_start_y[i], chain_data_bnbext.track_start_z[i]]
            track_end = [chain_data_bnbext.track_end_x[i], chain_data_bnbext.track_end_y[i], chain_data_bnbext.track_end_z[i]]
            track_fidvol = track_fidvol and is_fiducial(track_start) and is_fiducial(track_end)

        shower_fidvol = True
        for i in range(chain_data_bnbext.n_showers):
            shower_start = [chain_data_bnbext.shower_start_x[i], chain_data_bnbext.shower_start_y[i], chain_data_bnbext.shower_start_z[i]]
            shower_fidvol = shower_fidvol and is_fiducial(shower_start)

        if track_fidvol and shower_fidvol:
            fill_kin_branches(chain_data_bnbext,data_ext_scaling_factor*total_pot/total_data_bnb_pot,variables)
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
