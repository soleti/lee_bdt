#!/usr/local/bin/python3

from ROOT import TChain, TTree, TFile, TVector3
from glob import glob
import math
from bdt_common import total_pot, variables, spectators, total_data_bnb_pot
from bdt_common import x_start, x_end, y_start, y_end, z_start, z_end
from random import random


def is_fiducial(point):
    ok_y = y_start + 20 < point[1] < y_end - 20
    ok_x = x_start + 10 < point[0] < x_end - 10
    ok_z = z_start + 10 < point[2] < z_end - 50
    return ok_y and ok_x and ok_z


def is_active(point):
    ok_y = y_start < point[1] < y_end
    ok_x = x_start < point[0] < x_end
    ok_z = z_start < point[2] < z_end
    return ok_y and ok_x and ok_z


def choose_shower(root_chain):
    most_energetic_shower = 0
    shower_id = 0
    for ish in range(root_chain.n_showers):
        if root_chain.shower_energy[ish] < 3:
            if root_chain.shower_energy[ish] > most_energetic_shower:
                most_energetic_shower = root_chain.shower_energy[ish]
                shower_id = ish
    return shower_id


def choose_track(root_chain):
    track_id = 0
    least_proton_track = 1

    for itrk in range(root_chain.n_tracks):
        if root_chain.predict_p[itrk] < least_proton_track:
            least_proton_track = root_chain.predict_p[itrk]
            track_id = itrk
    return track_id


def fill_kin_branches(root_chain, numu_chain, weight, variables, option=""):
    longest_track = 0
    longest_track_id = 0
    track_id = 0
    most_proton_track = 1
    most_proton_track_length = 0
    for itrk in range(root_chain.n_tracks):
        if root_chain.predict_p[itrk] < most_proton_track:
            most_proton_track = root_chain.predict_p[itrk]
            track_id = itrk
            most_proton_track_length = root_chain.track_len[itrk]
        if root_chain.track_len[itrk] > longest_track:
            longest_track = root_chain.track_len[itrk]
            longest_track_id = itrk

    shower_id = choose_shower(root_chain)

    v_track = TVector3(
        root_chain.track_dir_x[longest_track_id],
        root_chain.track_dir_y[longest_track_id],
        root_chain.track_dir_z[longest_track_id])
    v_shower = TVector3(
        root_chain.shower_dir_x[shower_id],
        root_chain.shower_dir_y[shower_id],
        root_chain.shower_dir_z[shower_id])
    costheta_shower_track = v_track.Dot(
        v_shower) / (v_track.Mag() * v_shower.Mag())

    signal = 0
    if root_chain.category == 2:
        signal = 1

    track_vertex = [
        root_chain.track_start_x[track_id],
        root_chain.track_start_y[track_id],
        root_chain.track_start_z[track_id]]
    track_end = [
        root_chain.track_end_x[track_id],
        root_chain.track_end_y[track_id],
        root_chain.track_end_z[track_id]]
    shower_vertex = [
        root_chain.shower_start_x[shower_id],
        root_chain.shower_start_y[shower_id],
        root_chain.shower_start_z[shower_id]]
    neutrino_vertex = [root_chain.vx, root_chain.vy, root_chain.vz]
    true_neutrino_vertex = [
        root_chain.true_vx_sce,
        root_chain.true_vy_sce,
        root_chain.true_vz_sce]

    shower_vertex_d = math.sqrt(
        sum([(s - n)**2 for s, n in zip(shower_vertex, neutrino_vertex)]))
    track_vertex_d = math.sqrt(
        sum([(t - n)**2 for t, n in zip(track_vertex, neutrino_vertex)]))
    track_end_d = math.sqrt(
        sum([(t - n)**2 for t, n in zip(track_end, neutrino_vertex)]))

    track_shower_d = math.sqrt(
        sum([(s - t)**2 for s, t in zip(shower_vertex, track_vertex)]))
    trackend_shower_d = math.sqrt(
        sum([(s - t)**2 for s, t in zip(shower_vertex, track_end)]))

    direction = 1
    if trackend_shower_d < track_shower_d:
        direction = -1
    theta = math.acos(direction * root_chain.track_dir_z[track_id])

    variables["is_signal"][0] = signal
    variables["track_length"][0] = root_chain.track_len[track_id]
    variables["track_phi"][0] = math.degrees(
        root_chain.track_phi[track_id])
    variables["track_theta"][0] = math.degrees(theta)
    variables["shower_energy"][0] = root_chain.shower_energy[shower_id]

    total_shower_energy = sum(root_chain.shower_energy)

    variables["total_track_energy"][0] = root_chain.E - total_shower_energy
    variables["total_shower_energy"][0] = total_shower_energy

    variables["shower_theta"][0] = math.degrees(
        root_chain.shower_theta[shower_id])
    variables["shower_phi"][0] = math.degrees(
        root_chain.shower_phi[shower_id])
    variables["shower_distance"][0] = shower_vertex_d
    variables["track_distance"][0] = track_vertex_d

    variables["track_start_x"][0] = root_chain.track_start_x[track_id]
    variables["track_end_x"][0] = root_chain.track_end_x[track_id]

    variables["track_start_y"][0] = root_chain.track_start_y[track_id]
    variables["track_end_y"][0] = root_chain.track_end_y[track_id]

    variables["track_start_z"][0] = root_chain.track_start_z[track_id]
    variables["track_end_z"][0] = root_chain.track_end_z[track_id]

    variables["shower_start_x"][0] = root_chain.shower_start_x[shower_id]
    # variables["shower_end_x"][0] = root_chain.shower_end_x[shower_id]

    variables["shower_start_y"][0] = root_chain.shower_start_y[shower_id]
    # variables["shower_end_y"][0] = root_chain.shower_end_y[shower_id]

    variables["shower_start_z"][0] = root_chain.shower_start_z[shower_id]
    # variables["shower_end_z"][0] = root_chain.shower_end_z[shower_id]

    variables["reco_energy"][0] = root_chain.E

    if option == "cosmic_mc" or option == "ext_data":
        variables["category"][0] = 0
    elif option == "lee":
        variables["category"][0] = 8
    else:
        variables["category"][0] = root_chain.category

    # if option == "bnb_data" and 0.005 < total_shower_energy < 0.01:
    #     print("data", root_chain.run, root_chain.subrun, root_chain.event)
    #     print(root_chain.n_tracks, root_chain.n_showers, shower_vertex_d, root_chain.track_start_z[track_id])

    variables["event_weight"][0] = weight
    variables["pt"][0] = pt_plot(root_chain)

    variables["n_tracks"][0] = root_chain.n_tracks
    variables["n_showers"][0] = root_chain.n_showers
    variables["track_shower_angle"][0] = costheta_shower_track

    variables["event"][0] = root_chain.event
    variables["run"][0] = root_chain.run
    variables["subrun"][0] = root_chain.subrun
    variables["proton_score"][0] = max(0, root_chain.predict_p[track_id])
    variables["interaction_type"][0] = root_chain.interaction_type
    variables["shower_open_angle"][0] = math.degrees(
        root_chain.shower_open_angle[shower_id])

    try:
        variables["shower_pca"][0] = max(0, root_chain.shower_pca[shower_id])
        variables["track_pca"][0] = max(0, root_chain.track_pca[track_id])
    except AttributeError:
        variables["shower_pca"][0] = 0
        variables["track_pca"][0] = 0

    if numu_selection(numu_chain) < 1 and numu_selection(numu_chain) > 0:
        variables["numu_score"][0] = numu_selection(numu_chain)
    else:
        variables["numu_score"][0] = 0
    # variables["numu_score"][0] = 0
    dedx = root_chain.shower_dEdx[shower_id][2]

    variables["dedx"][0] = max(0, dedx)


def pt_plot(root_chain):
    p_showers = []
    for ish in range(root_chain.n_showers):
        if root_chain.shower_energy[ish] > 0:
            p_vector = TVector3(
                root_chain.shower_dir_x[ish],
                root_chain.shower_dir_y[ish],
                root_chain.shower_dir_z[ish])
            if root_chain.shower_energy[ish] < 10:
                p_vector.SetMag(
                    math.sqrt(
                        (root_chain.shower_energy[ish] + 0.052)**2 - 0.052**2))
            p_showers.append(p_vector)

    p_tracks = []
    for itr in range(root_chain.n_tracks):
        if root_chain.track_energy[itr] > 0:
            p_vector = TVector3(
                root_chain.track_dir_x[itr],
                root_chain.track_dir_y[itr],
                root_chain.track_dir_z[itr])
            p_vector.SetMag(
                math.sqrt(
                    (root_chain.track_energy[itr] +
                     0.938)**2 -
                    0.938**2))
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

    pt = (p_track_sum + p_shower_sum).Perp()
    return pt


def numu_selection(mychain):
    if mychain.GetEntries() == 0: return 0
    for i in range(mychain.nslices):
        flashmatch_cut = not (
            mychain.slc_flsmatch_qllx[i] -
            mychain.slc_flsmatch_tpcx[i] > 20)
        dist_cut = False
        if len(mychain.beamfls_z) > 0:
            dist_cut = mychain.slc_flsmatch_hypoz[i] - \
                mychain.beamfls_z[0] < 100

        broken_tracks_cut = not (mychain.slc_vtxcheck_angle[i] > 2.9) and not (
            mychain.slc_vtxcheck_angle[i] < 0.05 and mychain.slc_vtxcheck_angle[i] != -9999)
        # one_track_cut = not (mychain.slc_ntrack[i] == 0)
        quality_cut = mychain.slc_passed_min_track_quality[i]

        # print(mychain.slc_flsmatch_score[i])
        if flashmatch_cut and broken_tracks_cut and quality_cut and dist_cut and mychain.slc_flsmatch_score[i] > 0:
            return mychain.slc_flsmatch_score[i]

    return 0


def fill_tree(chain, chain_numu, weight, tree, option=""):
    total_events = 0

    for i in range(chain.GetEntries()):
        chain.GetEntry(i)
        chain_numu.GetEntry(i)

        if chain.passed:
            track_fidvol = True
            for i in range(chain.n_tracks):
                track_start = [
                    chain.track_start_x[i],
                    chain.track_start_y[i],
                    chain.track_start_z[i]]
                track_end = [
                    chain.track_end_x[i],
                    chain.track_end_y[i],
                    chain.track_end_z[i]]
                track_fidvol = track_fidvol and is_fiducial(track_start) and is_fiducial(track_end)

            shower_fidvol = True
            for i in range(chain.n_showers):
                shower_start = [
                    chain.shower_start_x[i],
                    chain.shower_start_y[i],
                    chain.shower_start_z[i]]
                shower_fidvol = shower_fidvol and is_fiducial(shower_start)


            option_check = True
            event_weight = weight

            if option == "bnb":
                option_check = abs(chain.nu_pdg) != 12 # and 111 not in chain.nu_daughters_pdg
            if option == "nue":
                event_weight = weight * chain.bnbweight
                option_check = abs(chain.nu_pdg) == 12

            shower_id = choose_shower(chain)
            track_id = choose_track(chain)

            track_vertex = [
                chain.track_start_x[track_id],
                chain.track_start_y[track_id],
                chain.track_start_z[track_id]]

            shower_vertex = [
                chain.shower_start_x[shower_id],
                chain.shower_start_y[shower_id],
                chain.shower_start_z[shower_id]]

            neutrino_vertex = [chain.vx, chain.vy, chain.vz]

            shower_vertex_d = math.sqrt(
                sum([(s - n)**2 for s, n in
                     zip(shower_vertex, neutrino_vertex)]))

            track_vertex_d = math.sqrt(
                sum([(t - n)**2 for t, n in
                     zip(track_vertex, neutrino_vertex)]))

            pions = 0

            if option_check and is_fiducial(neutrino_vertex) and track_fidvol and shower_fidvol:
                total_events += event_weight

                if option == "nue" and chain.category == 2:
                    for i in range(len(chain.nu_daughters_pdg)):
                        if chain.nu_daughters_pdg[i] == 111 or chain.nu_daughters_pdg[i] == 211:
                            pions += event_weight

                fill_kin_branches(chain, chain_numu, event_weight, variables, option)
                tree.Fill()

    print(total_events, pions)
    return total_events


data_ext_scaling_factor = 1.299
samples = ["pi0", "cosmic_mc", "bnb", "nue", "bnb_data", "ext_data", "lee"]

tree_files = [glob("pi0/*/*.root"),
              glob("cosmic_intime_dedx/*/*.root"),
              glob("mc_bnb_slimmed/*/*.root"),
              glob("mc_nue_pca/*/*.root"),
              glob("data_bnb_dedx/*/*.root"),
              glob("data_ext_dedx/*/*.root"),
              glob("lee/*/*.root")]

# tree_files = [glob("pi0/*/*.root"),
#               glob("wouter_files/soft_intime_83.root"),
#               glob("wouter_files/soft_nu_84.root"),
#               glob("wouter_files/soft_nue_84.root"),
#               glob("data_bnb_dedx/*/*.root"),
#               glob("wouter_files/soft_extbnb_84.root"),
#               glob("lee/*/*.root")]
#
# tree_files = [glob("pi0/*/*.root"),
#               glob("wouter_files/new_intime_83.root"),
#               glob("wouter_files/new_nu_83.root"),
#               glob("wouter_files/new_nue_84.root"),
#               glob("data_bnb_dedx/*/*.root"),
#               glob("wouter_files/new_extbnb_84.root"),
#               glob("lee/*/*.root")]
#
# tree_files = [glob("pi0/*/*.root"),
#               glob("hardcuts/tree_strong_intime_83.root"),
#               glob("hardcuts/tree_strong_nu_83.root"),
#               glob("hardcuts/tree_strong_nue_84.root"),
#               glob("data_bnb_dedx/*/*.root"),
#               glob("hardcuts/tree_strong_extbnb_84.root")]

chains = []
chains_numu = []
chains_pot = []
for i, files in enumerate(tree_files):
    chains.append(TChain("robertoana/pandoratree"))
    chains_numu.append(TChain("UBXSec/tree"))
    chains_pot.append(TChain("robertoana/pot"))

    for j, f in enumerate(files):
        chains[i].Add(f)
        # if i != 1:
        #     chains_numu[i].Add(f)
        chains_pot[i].Add(f)

pots = []
for c in chains_pot:
    total_pot_file = 0
    for i in range(c.GetEntries()):
        c.GetEntry(i)
        total_pot_file += c.pot

    pots.append(total_pot_file)

pots_dict = dict(zip(samples, pots))
chains_dict = dict(zip(samples, chains))
chains_numu_dict = dict(zip(samples, chains_numu))
chains_pot_dict = dict(zip(samples, chains_pot))
variables = dict(variables + spectators)
wouter_scaling = 1.2244
roberto_scaling = 1.350

weights = [total_pot / pots_dict["pi0"],
           data_ext_scaling_factor * total_pot / total_data_bnb_pot * roberto_scaling *
           chains_dict["ext_data"].GetEntries() / chains_dict["cosmic_mc"].GetEntries(),
           total_pot / pots_dict["bnb"],
           total_pot / pots_dict["nue"],
           1,
           data_ext_scaling_factor,
           total_pot / pots_dict["lee"]]

files = ["pi0_file.root", "cosmic_mc_file.root", "mc_file.root",
         "nue_file.root", "bnb_file.root", "bnbext_file.root", "lee_file.root"]
tree_names = ["pi0_tree", "cosmic_mc_tree", "mc_tree",
              "nue_tree", "bnb_tree", "bnbext_tree", "lee_tree"]

trees = []

for t in tree_names:
    trees.append(TTree(t, t))

for n, b in variables.items():
    for t in trees:
        t.Branch(n, b, n + "/f")

samples = ["pi0", "cosmic_mc", "bnb", "nue", "bnb_data", "ext_data", "lee"]
print(chains[0].GetEntries(), pots_dict["pi0"])

for i, s in enumerate(samples):
    print(s)
    print("Weight", weights[i])
    print("Events", fill_tree(chains[i], chains_numu[i],
                              weights[i], trees[i], s))

for f, t in zip(files, trees):
    tfile = TFile(f, "RECREATE")
    t.Write()
    tfile.Close()

run_subrun_bnb = open("run_subrun_bnb.txt", "w")

for i in range(chains_pot_dict["bnb_data"].GetEntries()):
    chains_pot_dict["bnb_data"].GetEntry(i)
    run_subrun = "%i %i" % (chains_pot_dict["bnb_data"].run,
                            chains_pot_dict["bnb_data"].subrun)
    print(run_subrun, file=run_subrun_bnb)

run_subrun_bnb.close()

run_subrun_ext = open("run_subrun_ext.txt", "w")
for i in range(chains_pot_dict["ext_data"].GetEntries()):
    chains_pot_dict["ext_data"].GetEntry(i)
    run_subrun = "%i %i" % (chains_pot_dict["ext_data"].run,
                            chains_pot_dict["ext_data"].subrun)
    print(run_subrun, file=run_subrun_ext)

run_subrun_ext.close()
