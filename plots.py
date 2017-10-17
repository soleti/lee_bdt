#!/usr/local/bin/python3

from ROOT import TChain, TTree, TFile, TVector3
from glob import glob
import math
from bdt_common import total_pot, variables, spectators
from bdt_common import x_start, x_end, y_start, y_end, z_start, z_end


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


def fill_kin_branches(root_chain, numu_chain, weight, variables):
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
    if root_chain.category == 2 and root_chain.E < 1:
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
    variables["category"][0] = root_chain.category
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

    # if numu_selection(numu_chain) < 1 and numu_selection(numu_chain) > 0:
    #     variables["numu_score"][0] = numu_selection(numu_chain)
    # else:
    #     variables["numu_score"][0] = 0

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
        # if flashmatch_cut and broken_tracks_cut and quality_cut and dist_cut and mychain.slc_flsmatch_score[i] > 0:
        if quality_cut:
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
                option_check = abs(chain.nu_pdg) != 12
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
                sum([(s - n)**2 for s, n in zip(shower_vertex, neutrino_vertex)]))

            track_vertex_d = math.sqrt(
                sum([(t - n)**2 for t, n in zip(track_vertex, neutrino_vertex)]))

            dedx = chain.shower_dEdx[shower_id][2] < 3.5 and chain.shower_dEdx[shower_id][2] > 1
            openangle = math.degrees(chain.shower_open_angle[shower_id]) < 15 and math.degrees(chain.shower_open_angle[shower_id]) > 1
            shower_angle = math.degrees(chain.shower_theta[shower_id]) < 80
            track_distance = track_vertex_d < 5
            shower_distance = shower_vertex_d < 7

            if option_check and is_fiducial(neutrino_vertex) and shower_fidvol and track_fidvol:# and numu_score <  0.5 and dedx and openangle and shower_angle and track_distance and shower_distance:  # :
                total_events += event_weight
                fill_kin_branches(chain, chain_numu, event_weight, variables)
                tree.Fill()

    return total_events


cosmic_mc = glob("cosmic_intime/*/*.root")

# nue_cosmic = glob("nue_files_6_42_energy/*/*.root")
# bnb_cosmic = glob("nu_files_6_42_energy/*/*.root")
# data_bnb = glob("data_files_bnb_6_42_energy/*/*.root")
# data_bnbext = glob("data_files_bnbext_6_42_energy/*/*.root")
# data_ext_scaling_factor = 1.2640

data_ext_scaling_factor = 1.299
nue_cosmic = glob("mc_nue_old/*/*.root")
bnb_cosmic = glob("mc_bnb_mcc83/*/*.root")
data_bnb = glob("data_bnb_mcc83/*/*.root")
data_bnbext = glob("data_ext_mcc83/*/*.root")

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
    #chain_cosmic_mc_numu.Add(f)

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

run_subrun_bnb = open("run_subrun_bnb.txt", "w")

for i in range(chain_data_bnb_pot.GetEntries()):
    chain_data_bnb_pot.GetEntry(i)
    run_subrun = "%i %i" % (chain_data_bnb_pot.run, chain_data_bnb_pot.subrun)
    run_subrun_list.append(run_subrun)
    print(run_subrun, file=run_subrun_bnb)

run_subrun_bnb.close()

total_data_bnb_pot = 4.758e19
print("Total data POT BNB {0:.2e}".format(total_data_bnb_pot))

run_subrun_ext = open("run_subrun_ext.txt", "w")
for i in range(chain_data_bnbext_pot.GetEntries()):
    chain_data_bnbext_pot.GetEntry(i)
    run_subrun = "%i %i" % (chain_data_bnbext_pot.run,
                            chain_data_bnbext_pot.subrun)
    run_subrun_list.append(run_subrun)
    print(run_subrun, file=run_subrun_ext)
run_subrun_ext.close()

print("Data EXT scaling factor {0:.2f}".format(data_ext_scaling_factor))

variables = dict(variables + spectators)

cosmic_mc_tree = TTree("cosmic_mc_tree", "cosmic_mc_tree")
mc_tree = TTree("mc_tree", "mc_tree")
bnb_tree = TTree("bnb_tree", "bnb_tree")
bnbext_tree = TTree("bnbext_tree", "bnbext_tree")

for n, b in variables.items():
    cosmic_mc_tree.Branch(n, b, n + "/f")
    mc_tree.Branch(n, b, n + "/f")
    bnb_tree.Branch(n, b, n + "/f")
    bnbext_tree.Branch(n, b, n + "/f")

print("*** MC cosmic sample ***")
total_cosmic_mc = fill_tree(chain_cosmic_mc, chain_cosmic_mc_numu, 1,
                            cosmic_mc_tree, "cosmic")
print("MC cosmic {} events".format(total_cosmic_mc))

print("*** MC BNB + cosmic sample ***")
print("Weight {0:.2f}".format(total_pot / total_bnb_pot))
total_mc = fill_tree(chain, chain_numu, total_pot / total_bnb_pot, mc_tree,
                     "bnb")
print("MC {0:.0f} events".format(total_mc))

print("*** MC nu_e + cosmic sample ***")
print("Weight {0:.2f}".format(total_pot / total_nue_pot))
total_nu_e = fill_tree(chain_nue, chain_nue_numu, total_pot / total_nue_pot,
                       mc_tree, "nue")
print("MC nu_e {0:.0f} events".format(total_nu_e))

print("*** Data BNB sample ***")
print("Weight {0:.2f}".format(total_pot / total_data_bnb_pot))
total_data_bnb = fill_tree(chain_data_bnb, chain_data_bnb_numu,
                           total_pot / total_data_bnb_pot, bnb_tree)
print("Data BNB {0:.0f} events".format(total_data_bnb))

print("*** Data EXT sample ***")
print("Weight {0:.2f}".format(data_ext_scaling_factor * total_pot / total_data_bnb_pot))
total_data_ext = fill_tree(
    chain_data_bnbext,
    chain_data_bnbext_numu,
    data_ext_scaling_factor * total_pot / total_data_bnb_pot,
    bnbext_tree)
print("Data EXT {0:.0f} events".format(total_data_ext))

bnb_ext = total_data_bnb - total_data_ext
print("Data BNB-EXT {0:.0f} events".format(bnb_ext))

print("Ratio (BNB-EXT)/MC {0:.2f}".format(bnb_ext / total_mc))

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
