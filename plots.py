#!/usr/bin/env python3.4

from ROOT import TChain, TTree, TFile, TVector3
from array import array
from glob import glob
import math
from bdt_common import total_pot, variables, spectators, x_start, x_end, y_start, y_end, z_start, z_end

def is_fiducial(point):
    ok_y = point[1] > y_start+20 and point[1] < y_end-20
    ok_x = point[0] > x_start+10 and point[0] < x_end-10
    ok_z = point[2] > z_start+10 and point[2] < z_end-50
    return ok_y and ok_x and ok_z

def is_active(point):
    ok_y = point[1] > y_start and point[1] < y_end
    ok_x = point[0] > x_start and point[0] < x_end
    ok_z = point[2] > z_start and point[2] < z_end
    return ok_y and ok_x and ok_z

def choose_shower(root_chain):
    most_energetic_shower = 0
    most_energetic_shower_id = 0
    for ish in range(root_chain.n_showers):
        if root_chain.shower_energy[ish] < 3:
            if root_chain.shower_energy[ish] > most_energetic_shower:
                most_energetic_shower = root_chain.shower_energy[ish]
                most_energetic_shower_id = ish
    return most_energetic_shower_id


def fill_kin_branches(root_chain, numu_chain, weight, variables):
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

    most_energetic_shower_id = choose_shower(root_chain)

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
    variables["shower_energy"][0] = root_chain.shower_energy[most_energetic_shower_id]
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

    variables["numu_score"][0] = numu_selection(numu_chain)

    dedx = root_chain.shower_dEdx[most_energetic_shower_id][2]

    if dedx < 0: dedx = 0
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

def numu_selection(mychain):
    for i in range(mychain.nslices):
        flashmatch_cut = not (mychain.slc_flsmatch_qllx[i] - mychain.slc_flsmatch_tpcx[i] > 20)
        dist_cut = False
        if len(mychain.beamfls_z) > 0: dist_cut = mychain.slc_flsmatch_hypoz[i] - mychain.beamfls_z[0] < 100

        broken_tracks_cut = not (mychain.slc_vtxcheck_angle[i] > 2.9) and not (mychain.slc_vtxcheck_angle[i] < 0.05 and mychain.slc_vtxcheck_angle[i] !=-9999 )
        # one_track_cut = not (mychain.slc_ntrack[i] == 0)
        quality_cut = mychain.slc_passed_min_track_quality[i]
        if flashmatch_cut and broken_tracks_cut and quality_cut and dist_cut:
            return mychain.slc_flsmatch_score[i]
    return -99999


def fill_tree(chain, chain_numu, weight, tree, option=""):
    total_events = 0

    for i in range(chain.GetEntries()):
        chain.GetEntry(i)
        chain_numu.GetEntry(i)
        if chain.passed:
            neutrino_vertex = [chain.vx, chain.vy, chain.vz]

            track_fidvol = True
            for i in range(chain.n_tracks):
                track_start = [chain.track_start_x[i], chain.track_start_y[i], chain.track_start_z[i]]
                track_end = [chain.track_end_x[i], chain.track_end_y[i], chain.track_end_z[i]]
                track_fidvol = track_fidvol and is_active(track_start) and is_active(track_end)

            shower_fidvol = True
            for i in range(chain.n_showers):
                shower_start = [chain.shower_start_x[i], chain.shower_start_y[i], chain.shower_start_z[i]]
                shower_fidvol = shower_fidvol and is_active(shower_start)

            dedx = chain.shower_dEdx[choose_shower(chain)][2]

            option_check = True
            event_weight = weight
            if option == "bnb":
                option_check = abs(chain.nu_pdg) != 12
            if option == "nue":
                event_weight = weight*chain.bnbweight
                option_check = abs(chain.nu_pdg) == 12

            most_energetic_shower_id = choose_shower(chain)
            dedx = chain.shower_dEdx[most_energetic_shower_id][2]

            if is_active(neutrino_vertex) and track_fidvol and shower_fidvol and option_check:
                total_events += event_weight
                fill_kin_branches(chain,chain_numu,event_weight,variables)
                tree.Fill()


    return total_events

cosmic_mc = glob("cosmic_only/*/*.root")

# nue_cosmic = glob("nue_files_6_42_energy/*/*.root")
# bnb_cosmic = glob("nu_files_6_42_energy/*/*.root")
# data_bnb = glob("data_files_bnb_6_42_energy/*/*.root")
# data_bnbext = glob("data_files_bnbext_6_42_energy/*/*.root")
# data_ext_scaling_factor = 1.2640

data_ext_scaling_factor = 1.2693
nue_cosmic = glob("mc_nue_nofidvol/*/*.root")
bnb_cosmic = glob("mc_bnb_nofidvol/*/*.root")
data_bnb = glob("data_bnb_nofidvol/*/*.root")
data_bnbext = glob("data_ext_nofidvol/*/*.root")

chain_cosmic_mc = TChain("robertoana/pandoratree")
chain_cosmic_mc_numu = TChain("UBXSec/tree")

chain = TChain("robertoana/pandoratree")
chain_pot = TChain("robertoana/pot")
chain_numu = TChain("UBXSec/tree")

chain_nue = TChain("robertoana/pandoratree")
chain_nue_pot = TChain("robertoana/pot")
chain_nue_numu = TChain("UBXSec/tree")

chain_data_bnb = TChain("robertoana/pandoratree")
chain_data_bnb_pot = TChain("robertoana/pot")
chain_data_bnb_numu = TChain("UBXSec/tree")

chain_data_bnbext = TChain("robertoana/pandoratree")
chain_data_bnbext_pot = TChain("robertoana/pot")
chain_data_bnbext_numu = TChain("UBXSec/tree")

for f in cosmic_mc:
    chain_cosmic_mc.Add(f)
    chain_cosmic_mc_numu.Add(f)

for f in nue_cosmic:
    chain_nue.Add(f)
    chain_nue_pot.Add(f)
    chain_nue_numu.Add(f)

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

run_subrun_list = []

total_nue_pot = 0
for i in range(chain_nue_pot.GetEntries()):
    chain_nue_pot.GetEntry(i)
    total_nue_pot += chain_nue_pot.pot

total_bnb_pot = 0
for i in range(chain_pot.GetEntries()):
    chain_pot.GetEntry(i)
    total_bnb_pot += chain_pot.pot

print("Total POT BNB {0:.2e}".format(total_bnb_pot))

run_subrun_bnb = open("run_subrun_bnb.txt","w")
for i in range(chain_data_bnb_pot.GetEntries()):
    chain_data_bnb_pot.GetEntry(i)
    run_subrun = "%i %i" % (chain_data_bnb_pot.run, chain_data_bnb_pot.subrun)
    run_subrun_list.append(run_subrun)
    print(run_subrun, file=run_subrun_bnb)
run_subrun_bnb.close()

total_data_bnb_pot = 3.787e19
print("Total data POT BNB {0:.2e}".format(total_data_bnb_pot))

run_subrun_ext = open("run_subrun_ext.txt","w")
for i in range(chain_data_bnbext_pot.GetEntries()):
    chain_data_bnbext_pot.GetEntry(i)
    run_subrun = "%i %i" % (chain_data_bnbext_pot.run, chain_data_bnbext_pot.subrun)
    run_subrun_list.append(run_subrun)
    print(run_subrun, file=run_subrun_ext)
run_subrun_ext.close()

print("Data EXT scaling factor {0:.2f}".format(data_ext_scaling_factor))

variables = dict(variables+spectators)

cosmic_mc_tree = TTree("cosmic_mc_tree","cosmic_mc_tree")
mc_tree = TTree("mc_tree","mc_tree")
bnb_tree = TTree("bnb_tree","bnb_tree")
bnbext_tree = TTree("bnbext_tree","bnbext_tree")

for n,b in variables.items():
    cosmic_mc_tree.Branch(n,b,n+"/f")
    mc_tree.Branch(n,b,n+"/f")
    bnb_tree.Branch(n,b,n+"/f")
    bnbext_tree.Branch(n,b,n+"/f")

print ("*** MC cosmic sample ***")
total_cosmic_mc = fill_tree(chain_cosmic_mc, chain_cosmic_mc_numu, 1, cosmic_mc_tree)
print("MC cosmic {} events".format(total_cosmic_mc))

print ("*** MC BNB + cosmic sample ***")
total_mc = fill_tree(chain, chain_numu, total_pot/total_bnb_pot, mc_tree, "bnb")
print("MC {0:.0f} events".format(total_mc))

print ("*** MC nu_e + cosmic sample ***")
total_nu_e = fill_tree(chain_nue, chain_nue_numu, total_pot/total_nue_pot, mc_tree, "nue")
print("MC nu_e {0:.0f} events".format(total_nu_e))

print ("*** Data BNB sample ***")
total_data_bnb = fill_tree(chain_data_bnb, chain_data_bnb_numu, total_pot/total_data_bnb_pot, bnb_tree)
print("Data BNB {0:.0f} events".format(total_data_bnb))

print ("*** Data EXT sample ***")
total_data_ext = fill_tree(chain_data_bnbext, chain_data_bnbext_numu, data_ext_scaling_factor*total_pot/total_data_bnb_pot, bnbext_tree)
print("Data EXT {0:.0f} events".format(total_data_ext))
print("Data BNB-EXT {0:.0f} events".format(total_data_bnb-total_data_ext))

print("Ratio (BNB-EXT)/MC {0:.2f}".format((total_data_bnb-total_data_ext)/total_mc))

cosmic_mc_file = TFile("cosmic_mc_file.root", "RECREATE")
cosmic_mc_tree.Write()
cosmic_mc_file.Close()

mc_file = TFile("mc_file.root", "RECREATE")
mc_tree.Write()
mc_file.Close()

bnb_file = TFile("bnb_file.root", "RECREATE")
bnb_tree.Write()
bnb_file.Close()

bnbext_file = TFile("bnbext_file.root", "RECREATE")
bnbext_tree.Write()
bnbext_file.Close()
