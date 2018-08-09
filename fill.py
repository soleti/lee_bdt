#!/usr/local/bin/python3

from ROOT import TChain, TTree, TFile, TVector3, TMVA
from glob import glob
from array import array
import math
from bdt_common import total_pot, variables, spectators, choose_shower, printProgressBar
from bdt_common import is_active, is_fiducial, distance
from bdt_common import ELECTRON_MASS, PROTON_MASS, MAX_N_SHOWERS, MAX_N_TRACKS
from bdt_common import PROTON_THRESHOLD, ELECTRON_THRESHOLD
from proton_energy import length2energy
import time
import statistics
import sys


TMVA.Tools.Instance()
reader_dqdx = TMVA.Reader(":".join([
    "!V",
    "!Silent",
    "Color"]))
dqdx = array("f", [0])
length = array("f", [0])
reader_dqdx.AddVariable("dqdx", dqdx)
reader_dqdx.AddVariable("len", length)
reader_dqdx.BookMVA("dqdx BDT",
                    "dataset/weights/TMVAClassification_dqdx BDT.weights.xml")

reader_dedx = TMVA.Reader(":".join([
    "!V",
    "!Silent",
    "Color"]))
dedx = array("f", [0])
nhits = array("f", [0])
reader_dedx.AddVariable("dedx", dedx)
reader_dedx.AddVariable("nhits", nhits)
reader_dedx.BookMVA("dedx BDT",
                    "dataset/weights/TMVAClassification_dedx BDT.weights.xml")

reader_reclass = TMVA.Reader(":".join([
    "!V",
    "!Silent",
    "Color"]))
angle = array("f", [0])
pca = array("f", [0])
res = array("f", [0])
open_angle = array("f", [0])
n_hits = array("f", [0])
ratio = array("f", [0])

reader_reclass.AddVariable("angle", angle)
reader_reclass.AddVariable("pca", pca)
reader_reclass.AddVariable("res", res)
reader_reclass.AddVariable("open_angle", open_angle)
reader_reclass.AddVariable("ratio", ratio)
reader_reclass.AddVariable("n_hits", n_hits)

reader_reclass.BookMVA("shower_ntuple BDT",
                       "dataset/weights/TMVAClassification_shower_ntuple BDT.weights.xml")

def min_dedx(root_chain):
    min_score = 1

    for i_sh in range(root_chain.n_showers):
        shower_dedx = root_chain.shower_dEdx[i_sh][2]

        shower_hits = root_chain.shower_nhits[i_sh][2]
        score = dedx_hits(shower_dedx, shower_hits)
        if score < min_score:
            min_score = score

    return min_score


def max_dqdx(root_chain):
    max_score = 0

    for i_tr in range(root_chain.n_tracks):
        track_dqdx = root_chain.track_dQdx[i_tr][2]
        if root_chain.category == 0 or root_chain.category == 6:
            track_dqdx *= 1.2

        track_length = root_chain.track_len[i_tr]
        score = dqdx_length(track_dqdx, track_length)
        if score > max_score:
            max_score = score

    return max_score


def min_dqdx(root_chain):
    min_score = 1

    for i_tr in range(root_chain.n_tracks):
        track_dqdx = root_chain.track_dQdx[i_tr][2]
        if root_chain.category == 0 or root_chain.category == 6:
            track_dqdx *= 1.2

        track_length = root_chain.track_len[i_tr]
        score = dqdx_length(track_dqdx, track_length)
        if score < min_score:
            min_score = score

    return min_score


def dqdx_length(v_dqdx, v_length):

    if v_dqdx < 0:
        return -1

    dqdx[0] = v_dqdx
    length[0] = v_length

    BDT_response = reader_dqdx.EvaluateMVA("dqdx BDT")
    return BDT_response


def shower_score(v_angle, v_pca, v_res, v_open_angle, v_n_hits, v_ratio):
    if v_angle > 0 and v_pca > 0 and v_res > 0 and v_open_angle > 0 and v_n_hits > 0 and v_ratio > 0:
        angle[0] = v_angle
        pca[0] = v_pca
        res[0] = v_res
        open_angle[0] = v_open_angle
        n_hits[0] = v_n_hits
        ratio[0] = v_ratio
        BDT_response = reader_reclass.EvaluateMVA("shower_ntuple BDT")
        return BDT_response
    else:
        return -1

def dedx_hits(v_dedx, v_hits):
    if v_dedx < 0:
        return -1

    dedx[0] = v_dedx
    nhits[0] = v_hits
    BDT_response = reader_dedx.EvaluateMVA("dedx BDT")
    return BDT_response

def pi0_mass(root_chain):
    if root_chain.n_showers < 2:
        return -1

    shower_energies = sorted(
        [e[2] for e in root_chain.shower_energy], reverse=True)
    shower_ids = []

    for i_sh in range(root_chain.n_showers):
        if root_chain.shower_energy[i_sh][2] in shower_energies:
            shower_ids.append(i_sh)

    v_1 = TVector3(
        root_chain.shower_dir_x[shower_ids[0]],
        root_chain.shower_dir_y[shower_ids[0]],
        root_chain.shower_dir_z[shower_ids[0]])

    v_2 = TVector3(
        root_chain.shower_dir_x[shower_ids[1]],
        root_chain.shower_dir_y[shower_ids[1]],
        root_chain.shower_dir_z[shower_ids[1]])

    cos = v_1.Dot(v_2) / (v_1.Mag() * v_2.Mag())
    angle = math.acos(cos)

    e1 = root_chain.shower_energy[shower_ids[0]][2]
    e2 = root_chain.shower_energy[shower_ids[1]][2]

    for i_sh in range(root_chain.n_showers):

        if i_sh in shower_ids:
            continue

        v_x = TVector3(
            root_chain.shower_dir_x[i_sh],
            root_chain.shower_dir_y[i_sh],
            root_chain.shower_dir_z[i_sh])

        cos_1 = v_x.Dot(v_1) / (v_x.Mag() * v_1.Mag())
        cos_2 = v_x.Dot(v_2) / (v_x.Mag() * v_2.Mag())
        if math.acos(cos_1) < math.acos(cos_2):
            e1 += root_chain.shower_energy[i_sh][2]
        else:
            e2 += root_chain.shower_energy[i_sh][2]

    pi0_mass = math.sqrt(4 * e1 * e2 * (math.sin(angle / 2)**2))

    return pi0_mass


def choose_track(root_chain):
    max_score = 0
    chosen_track = 0

    for i_tr in range(root_chain.n_tracks):
        track_dqdx = root_chain.track_dQdx[i_tr][2]
        if root_chain.category == 0 or root_chain.category == 6:
            track_dqdx *= 1.2

        track_length = root_chain.track_len[i_tr]
        score = dqdx_length(track_dqdx, track_length)
        if score > max_score:
            max_score = score
            chosen_track = i_tr

    return chosen_track


def choose_plane(root_chain):
    total_hits = [0, 0, 0]
    shower_hits = [0, 0, 0]
    track_hits = [0, 0, 0]

    for i_sh in range(root_chain.n_showers):
        for i_plane in range(len(root_chain.shower_nhits[i_sh])):
            total_hits[i_plane] += root_chain.shower_nhits[i_sh][i_plane]
            shower_hits[i_plane] += root_chain.shower_nhits[i_sh][i_plane]

    for i_tr in range(root_chain.n_tracks):
        for i_plane in range(len(root_chain.track_nhits[i_tr])):
            total_hits[i_plane] += root_chain.track_nhits[i_tr][i_plane]
            track_hits[i_plane] += root_chain.track_nhits[i_tr][i_plane]

    product = [s for t, s in zip(track_hits, shower_hits)]

    return product.index(max(product))

def fill_kin_branches(root_chain, weight, variables, option=""):
    no_tracks = False
    hit_index = 2
    shower_id = choose_shower(root_chain, hit_index)
    track_id = choose_track(root_chain)

    variables["E_dep"][0] = 0

    if "nue" in option:
        for i, energy in enumerate(root_chain.nu_daughters_E):
            if root_chain.nu_daughters_pdg[i] == 2212:
                if energy - 0.938 > PROTON_THRESHOLD:
                    variables["E_dep"][0] += energy - 0.938

            if root_chain.nu_daughters_pdg[i] == 11:
                if energy - 0.51e-3 > ELECTRON_THRESHOLD:
                    variables["E_dep"][0] += energy - 0.51e-3

    track_like_shower_id = -1
    if root_chain.n_tracks == 0 and root_chain.n_showers > 1:
        shower_res_std_list = list(root_chain.shower_res_std)
        track_like_shower_id = shower_res_std_list.index(min(shower_res_std_list))
        if track_like_shower_id == shower_id:
            track_like_shower_id = abs(shower_id - 1)
        no_tracks = True

    vec_shower = TVector3(
        root_chain.shower_dir_x[shower_id],
        root_chain.shower_dir_y[shower_id],
        root_chain.shower_dir_z[shower_id])

    true_vertex = [root_chain.true_vx_sce, root_chain.true_vy_sce, root_chain.true_vz_sce]
    neutrino_vertex = [root_chain.vx, root_chain.vy, root_chain.vz]

    variables["is_signal"][0] = int(root_chain.category == 2)
    variables["true_nu_is_fidvol"][0] = int(is_fiducial(true_vertex))
    variables["n_objects"][0] = root_chain.n_tracks + root_chain.n_showers
    variables["no_tracks"][0] = int(no_tracks)
    variables["nu_E"][0] = max(-999, root_chain.nu_E)

    y_shower = sum([root_chain.shower_nhits[i_sh][2] for i_sh in range(root_chain.n_showers)])
    y_track = sum([root_chain.track_nhits[i_sh][2]for i_sh in range(root_chain.n_tracks)])

    u_shower = sum([root_chain.shower_nhits[i_sh][0] for i_sh in range(root_chain.n_showers)])
    u_track = sum([root_chain.track_nhits[i_sh][0] for i_sh in range(root_chain.n_tracks)])

    v_shower = sum([root_chain.shower_nhits[i_sh][1] for i_sh in range(root_chain.n_showers)])
    v_trackh = sum([root_chain.track_nhits[i_sh][1] for i_sh in range(root_chain.n_tracks)])

    variables["total_hits"][0] = u_shower + u_track + v_shower + v_trackh + y_shower + y_track

    variables["total_hits_u"][0] = u_shower + u_track
    variables["total_hits_v"][0] = v_shower + v_trackh
    variables["total_hits_y"][0] = y_shower + y_track

    if option == "cosmic_mc" or option == "ext_data":
        variables["category"][0] = 0
    elif option == "lee":
        variables["category"][0] = 10
    elif option == "nue_cc":
        variables["category"][0] = 8
    else:
        variables["category"][0] = root_chain.category

    variables["event_weight"][0] = weight

    if option == "nue":
        if root_chain.cosmic_fraction > 0.5 and root_chain.category == 7:
            variables["category"][0] = 2

    # variables["pt"][0] = pt_plot(root_chain, 2)
    variables["event"][0] = root_chain.event
    variables["run"][0] = root_chain.run
    variables["subrun"][0] = root_chain.subrun
    variables["interaction_type"][0] = root_chain.interaction_type

    total_shower_energy = 0
    total_shower_energy_cali = 0
    total_track_energy_length = 0

    for i_tr in range(MAX_N_TRACKS):
        variables["track_pdg"][i_tr] = -999
        variables["track_phi"][i_tr] = -999
        variables["track_theta"][i_tr] = -999
        variables["track_length"][i_tr] = -999
        variables["track_start_x"][i_tr] = -999
        variables["track_start_y"][i_tr] = -999
        variables["track_start_z"][i_tr] = -999
        variables["track_res_mean"][i_tr] = -999
        variables["track_res_std"][i_tr] = -999
        variables["track_pca"][i_tr] = -999
        variables["track_end_x"][i_tr] = -999
        variables["track_end_y"][i_tr] = -999
        variables["track_end_z"][i_tr] = -999
        variables["track_shower_angle"][i_tr] = -999
        variables["track_distance"][i_tr] = -999
        variables["track_pidchipr"][i_tr] = -999
        variables["track_likelihood"][i_tr] = -999
        variables["track_mip_likelihood"][i_tr] = -999
        variables["track_p_likelihood"][i_tr] = -999
        variables["track_dqdx"][i_tr] = -999

    for i_sh in range(MAX_N_SHOWERS):
        variables["shower_pdg"][i_sh] = -999
        variables["shower_start_x"][i_sh] = -999
        variables["shower_start_y"][i_sh] = -999
        variables["shower_start_z"][i_sh] = -999
        variables["shower_energy"][i_sh] = -999
        variables["shower_pca"][i_sh] = -999
        variables["shower_open_angle"][i_sh] = -999
        variables["shower_dedx"][i_sh] = -999
        variables["shower_dedx_u"][i_sh] = -999
        variables["shower_dedx_v"][i_sh] = -999
        variables["shower_dedx_cali"][i_sh] = -999
        variables["shower_res_mean"][i_sh] = -999
        variables["shower_res_std"][i_sh] = -999
        variables["shower_phi"][i_sh] = -999
        variables["shower_theta"][i_sh] = -999
        variables["shower_distance"][i_sh] = -999

    variables["n_showers"][0] = 0
    variables["n_tracks"][0] = 0
    variables["shower_id"][0] = shower_id
    variables["track_id"][0] = track_id

    total_track_nhits = 0
    total_shower_nhits = 0

    for i_sh in range(root_chain.n_showers):
        shower_v = [
            root_chain.shower_start_x[i_sh],
            root_chain.shower_start_y[i_sh],
            root_chain.shower_start_z[i_sh]
        ]

        v_sh = TVector3(
            root_chain.shower_dir_x[i_sh],
            root_chain.shower_dir_y[i_sh],
            root_chain.shower_dir_z[i_sh]
        )

        cos = v_sh.Dot(vec_shower) / (v_sh.Mag() * vec_shower.Mag())
        shower_angle = math.degrees(math.acos(min(cos, 1)))
        shower_pca = root_chain.shower_pca[i_sh][2]
        converted_showers = 0
        shower_res_std = root_chain.shower_res_std[i_sh]
        shower_open_angle = root_chain.shower_open_angle[i_sh]
        shower_n_hits = root_chain.shower_nhits[i_sh][2]
        length = root_chain.shower_length[i_sh]
        ratio = shower_n_hits / length
        score = shower_score(shower_angle,
                             shower_pca,
                             shower_res_std,
                             shower_open_angle,
                             shower_n_hits,
                             ratio)
        shower_end = [s+d*length for s, d in zip(shower_v, v_sh)]

        if (i_sh == track_like_shower_id or score > 0) and i_sh != shower_id:
            variables["n_tracks"][0] += 1
            variables["track_pdg"][root_chain.n_tracks + converted_showers] = root_chain.matched_showers[i_sh]
            variables["track_start_x"][root_chain.n_tracks + converted_showers] = root_chain.shower_start_x[i_sh]
            variables["track_start_y"][root_chain.n_tracks + converted_showers] = root_chain.shower_start_y[i_sh]
            variables["track_start_z"][root_chain.n_tracks + converted_showers] = root_chain.shower_start_z[i_sh]
            variables["track_phi"][root_chain.n_tracks + converted_showers] = math.degrees(root_chain.shower_phi[i_sh])
            variables["track_theta"][root_chain.n_tracks + converted_showers] = math.degrees(root_chain.shower_theta[i_sh])
            variables["track_length"][root_chain.n_tracks + converted_showers] = root_chain.shower_length[i_sh]
            variables["track_pca"][root_chain.n_tracks + converted_showers] = max(-999, root_chain.shower_pca[i_sh][2])
            variables["track_res_mean"][root_chain.n_tracks + converted_showers] = max(-999, root_chain.shower_res_mean[i_sh])
            variables["track_res_std"][root_chain.n_tracks + converted_showers] = max(-999, root_chain.shower_res_std[i_sh])
            variables["track_pidchipr"][root_chain.n_tracks + converted_showers] = -999
            variables["track_likelihood"][root_chain.n_tracks + converted_showers] = -999
            variables["track_mip_likelihood"][root_chain.n_tracks + converted_showers] = -999
            variables["track_p_likelihood"][root_chain.n_tracks + converted_showers] = -999
            variables["track_dqdx"][root_chain.n_tracks + converted_showers] = -999
            variables["track_p_likelihood"][root_chain.n_tracks + converted_showers] = -999
            variables["track_end_x"][root_chain.n_tracks + converted_showers] = shower_end[0]
            variables["track_end_y"][root_chain.n_tracks + converted_showers] = shower_end[1]
            variables["track_end_z"][root_chain.n_tracks + converted_showers] = shower_end[2]
            length_e = length2energy(root_chain.shower_length[i_sh])
            total_track_energy_length += length_e
            converted_showers += 1
            total_track_nhits += root_chain.shower_nhits[i_sh][2]
            costheta_shower_track = v_sh.Dot(vec_shower) / (v_sh.Mag() * vec_shower.Mag())
            variables["track_shower_angle"][root_chain.n_tracks + converted_showers] = costheta_shower_track
            track_vertex_d = math.sqrt(sum([(t - n)**2 for t, n in zip(shower_v, neutrino_vertex)]))
            variables["track_distance"][root_chain.n_tracks + converted_showers] = track_vertex_d
        else:
            variables["n_showers"][0] += 1
            variables["shower_pdg"][i_sh] = root_chain.matched_showers[i_sh]
            variables["shower_start_x"][i_sh] = root_chain.shower_start_x[i_sh]
            variables["shower_start_y"][i_sh] = root_chain.shower_start_y[i_sh]
            variables["shower_start_z"][i_sh] = root_chain.shower_start_z[i_sh]
            variables["shower_theta"][i_sh] = math.degrees(root_chain.shower_theta[i_sh])
            variables["shower_phi"][i_sh] = math.degrees(root_chain.shower_phi[i_sh])
            shower_energy_cali = root_chain.shower_energy[i_sh][2] * root_chain.shower_energy_cali[i_sh][2]
            variables["shower_energy"][i_sh] = max(-999, shower_energy_cali)
            dedx = root_chain.shower_dEdx[i_sh][2]
            variables["shower_dedx"][i_sh] = max(-999, dedx)
            dedx_cali = dedx * root_chain.shower_dQdx_cali[i_sh][2]
            variables["shower_dedx_cali"][i_sh] = max(-999, dedx_cali)
            dedx_u = root_chain.shower_dEdx[shower_id][0] * root_chain.shower_dQdx_cali[shower_id][0]
            variables["shower_dedx_u"][i_sh] = max(-999, dedx_u)
            dedx_v = root_chain.shower_dEdx[shower_id][1] * root_chain.shower_dQdx_cali[shower_id][1]
            variables["shower_dedx_v"][i_sh] = max(-999, dedx_v)
            variables["shower_res_mean"][i_sh] = max(-999, root_chain.shower_res_mean[i_sh])
            variables["shower_res_std"][i_sh] = max(-999, root_chain.shower_res_std[i_sh])
            variables["shower_open_angle"][i_sh] = math.degrees(root_chain.shower_open_angle[i_sh])
            variables["shower_pca"][i_sh] = max(0, root_chain.shower_pca[i_sh][2])
            total_shower_energy += root_chain.shower_energy[i_sh][hit_index]
            total_shower_energy_cali += shower_energy_cali
            total_shower_nhits += root_chain.shower_nhits[i_sh][2]
            shower_vertex_d = math.sqrt(sum([(t - n)**2 for t, n in zip(shower_v, neutrino_vertex)]))
            variables["shower_distance"][i_sh] = shower_vertex_d

    for i_tr in range(root_chain.n_tracks):
        v_track = TVector3(
            root_chain.track_dir_x[i_tr],
            root_chain.track_dir_y[i_tr],
            root_chain.track_dir_z[i_tr]
        )

        track_vertex = [
            root_chain.track_start_x[i_tr],
            root_chain.track_start_y[i_tr],
            root_chain.track_start_z[i_tr]
        ]

        total_track_nhits += root_chain.track_nhits[i_tr][2]
        variables["n_tracks"][0] += 1
        variables["track_pdg"][i_tr] = root_chain.matched_tracks[i_tr]
        variables["track_pidchipr"][i_tr] = root_chain.track_pid_chipr[i_tr]
        if root_chain.track_bragg_p[i_tr] != 0:
            variables["track_likelihood"][i_tr] = math.log(root_chain.track_bragg_mip[i_tr]/root_chain.track_bragg_p[i_tr])
        else:
            variables["track_likelihood"][i_tr] = -999
        variables["track_mip_likelihood"][i_tr] = root_chain.track_bragg_mip[i_tr]
        variables["track_p_likelihood"][i_tr] = root_chain.track_bragg_p[i_tr]
        variables["track_dqdx"][i_tr] = max(-999, root_chain.track_dQdx[i_tr][2])
        variables["track_phi"][i_tr] = math.degrees(root_chain.track_phi[i_tr])
        variables["track_theta"][i_tr] = math.degrees(root_chain.track_theta[i_tr])
        variables["track_length"][i_tr] = root_chain.track_len[i_tr]
        variables["track_start_x"][i_tr] = root_chain.track_start_x[i_tr]
        variables["track_start_y"][i_tr] = root_chain.track_start_y[i_tr]
        variables["track_start_z"][i_tr] = root_chain.track_start_z[i_tr]
        variables["track_pca"][i_tr] = max(-999, root_chain.track_pca[i_tr][2])
        variables["track_res_mean"][i_tr] = max(-999, root_chain.track_res_mean[i_tr])
        variables["track_res_std"][i_tr] = max(-999, root_chain.track_res_std[i_tr])
        variables["track_end_x"][i_tr] = root_chain.track_end_x[i_tr]
        variables["track_end_y"][i_tr] = root_chain.track_end_y[i_tr]
        variables["track_end_z"][i_tr] = root_chain.track_end_z[i_tr]
        length_e = length2energy(root_chain.track_len[i_tr])
        total_track_energy_length += length_e
        costheta_shower_track = v_track.Dot(vec_shower) / (v_track.Mag() * vec_shower.Mag())
        variables["track_shower_angle"][i_tr] = costheta_shower_track
        track_vertex_d = math.sqrt(
            sum([(t - n)**2 for t, n in zip(track_vertex, neutrino_vertex)]))
        variables["track_distance"][i_tr] = track_vertex_d

    variables["track_hits"][0] = total_track_nhits
    variables["shower_hits"][0] = total_shower_nhits
    variables["hits_ratio"][0] = total_shower_nhits/(total_track_nhits+total_shower_nhits)

    variables["total_shower_energy"][0] = total_shower_energy
    variables["total_shower_energy_cali"][0] = total_shower_energy_cali
    variables["total_track_energy_length"][0] = total_track_energy_length
    variables["reco_energy"][0] = ((total_shower_energy_cali + 9.85322e-03) / 7.82554e-01) + total_track_energy_length
    # variables["reco_energy"][0] = ((total_shower_energy_cali + 1.36881e-02) / 7.69908e-01) + total_track_energy_length
    variables["numu_score"][0] = root_chain.numu_cuts


def pt_plot(root_chain, plane):
    p_showers = []
    for ish in range(root_chain.n_showers):
        if root_chain.shower_energy[ish][plane] > 0:
            p_vector = TVector3(
                root_chain.shower_dir_x[ish],
                root_chain.shower_dir_y[ish],
                root_chain.shower_dir_z[ish])
            if root_chain.shower_energy[ish][plane] < 10:
                p_vector.SetMag(
                    math.sqrt(
                        (root_chain.shower_energy[ish][plane] + ELECTRON_MASS)**2 - ELECTRON_MASS**2))
            p_showers.append(p_vector)

    p_tracks = []
    for itr in range(root_chain.n_tracks):
        if root_chain.track_energy_hits[itr][plane] > 0:
            p_vector = TVector3(
                root_chain.track_dir_x[itr],
                root_chain.track_dir_y[itr],
                root_chain.track_dir_z[itr])
            p_vector.SetMag(
                math.sqrt(
                    (root_chain.track_energy_hits[itr][plane] +
                     PROTON_MASS)**2 -
                    PROTON_MASS**2))
            p_tracks.append(p_vector)

    p_track_sum = TVector3()
    if p_tracks:
        p_track_sum = p_tracks[0]
        for itr in p_tracks[1:]:
            p_track_sum += itr

    p_shower_sum = TVector3()
    if p_showers:
        p_shower_sum = p_showers[0]
        for ish in p_showers[1:]:
            p_shower_sum += ish

    pt = (p_track_sum + p_shower_sum).Perp()
    return pt


def fill_tree(chain, weight, tree, option=""):
    total_events = 0
    total_entries = int(chain.GetEntries() / 1)

    for ievt in range(total_entries):
        chain.GetEntry(ievt)
        printProgressBar(ievt, total_entries, prefix="Progress:", suffix="Complete", length=20)

        if not chain.passed:
            continue

        track_fidvol = True

        for i_tr in range(chain.n_tracks):
            track_start = [
                chain.track_start_x[i_tr],
                chain.track_start_y[i_tr],
                chain.track_start_z[i_tr]
            ]
            track_end = [
                chain.track_end_x[i_tr],
                chain.track_end_y[i_tr],
                chain.track_end_z[i_tr]
            ]

            track_fidvol = track_fidvol and is_fiducial(track_start) and is_fiducial(track_end)
            if not track_fidvol:
                break

        if not track_fidvol:
            continue

        shower_fidvol = True

        for i_sh in range(chain.n_showers):
            shower_start = [
                chain.shower_start_x[i_sh],
                chain.shower_start_y[i_sh],
                chain.shower_start_z[i_sh]
            ]

            shower_fidvol = shower_fidvol and is_fiducial(shower_start)
            if not shower_fidvol:
                break

        if not shower_fidvol:
            continue

        option_check = True
        event_weight = weight

        if option == "bnb":
            option_check = abs(chain.nu_pdg) != 12
        if abs(chain.nu_pdg) == 12:
            event_weight = weight * chain.bnbweight
        if "nue" in option:
            option_check = abs(chain.nu_pdg) == 12

        if not option_check:
            continue

        neutrino_vertex = [chain.vx, chain.vy, chain.vz]

        if not is_fiducial(neutrino_vertex):
            continue

        # If there are no tracks we require at least two showers
        showers_2_tracks_0 = True
        if chain.n_tracks == 0 and chain.n_showers == 1:
            showers_2_tracks_0 = False
        if not showers_2_tracks_0:
            continue

        total_events += event_weight

        if "nue" in option and chain.category == 2:
            if 2212 not in chain.nu_daughters_pdg or abs(111) in chain.nu_daughters_pdg or abs(211) in chain.nu_daughters_pdg:
                option = "nue_cc"
            else:
                option = "nue"

        fill_kin_branches(chain, event_weight, variables, option)
        tree.Fill()

    return total_events

begin_time = time.time()
# To be obtained with Zarko's POT tool
data_ext_scaling_factor = 0.1327933846  # Sample with remapped PMTs

samples = ["nue", "bnb", "bnb_data", "ext_data"]

tree_files = [glob("data_files/mc_nue_pid/*.root"),
              glob("data_files/mc_bnb_pid/*/*.root"),
              glob("data_files/data_bnb_pid/*/*.root"),
              glob("data_files/data_ext_pid/*/*.root")]

chains = []
chains_pot = []
for i, files in enumerate(tree_files):
    chains.append(TChain("robertoana/pandoratree"))
    chains_pot.append(TChain("robertoana/pot"))

    for j, f in enumerate(files):
        chains[i].Add(f)
        chains_pot[i].Add(f)

pots = []
for k, c in enumerate(chains_pot):
    total_pot_file = 0
    for i in range(c.GetEntries()):
        c.GetEntry(i)
        total_pot_file += c.pot

    pots.append(total_pot_file)

pots_dict = dict(zip(samples, pots))
chains_dict = dict(zip(samples, chains))
chains_pot_dict = dict(zip(samples, chains_pot))
variables = dict(variables + spectators)

weights = [total_pot / pots_dict["nue"] * 1.028, # Position of the detector is slightly wrong
           total_pot / pots_dict["bnb"] * 1.028,
           1,
           data_ext_scaling_factor]

print(weights)
files = ["nue_file.root", "mc_file.root",
         "bnb_file.root", "bnbext_file.root"]
tree_names = ["nue_tree", "mc_tree",
              "bnb_tree", "bnbext_tree"]

trees = []

for t in tree_names:
    trees.append(TTree(t, t))

for n, b in variables.items():
    for t in trees:
        if len(b) > 1:
            if "track" in n:
                t.Branch(n, b, "%s[%i]/%s" % (n, MAX_N_TRACKS, b.typecode))
            elif "shower" in n:
                t.Branch(n, b, "%s[%i]/%s" % (n, MAX_N_SHOWERS, b.typecode))
        else:
            t.Branch(n, b, n + "/%s" % b.typecode)

samples = ["nue", "bnb", "bnb_data", "ext_data"]

for i, s in enumerate(samples):
    start_time = time.time()
    print("******************************")
    print("Sample", s)
    print("Weight", weights[i])
    events = fill_tree(chains[i], weights[i], trees[i], s)
    print("\nEvents", events)
    print("Time to fill %.1f s"  % (time.time() - start_time))

for f, t in zip(files, trees):
    tfile = TFile("root_files/" + f, "RECREATE")
    t.Write()
    tfile.Close()

print("Total time %.2f" % (time.time() - begin_time))
