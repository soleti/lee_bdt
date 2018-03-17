#!/usr/local/bin/python3

from ROOT import TChain, TTree, TFile, TVector3, TH1F, TCanvas
from ROOT import kRed, kBlue
from ROOT import gPad, TGaxis
from glob import glob
import math
from bdt_common import total_pot, variables, spectators, total_data_bnb_pot
from bdt_common import x_start, x_end, y_start, y_end, z_start, z_end, is_fiducial
from bdt_common import ELECTRON_MASS, PROTON_MASS, SIGNAL_INTERVAL, printProgressBar
from proton_energy import length2energy
import time
import statistics


h_angle_proton = TH1F("h_angle_proton", ";Angle from main shower [#circ];N. Entries / 5#circ", 36, 0, 180)
h_angle_nonproton = TH1F("h_angle_nonproton", ";Angle from main shower [#circ];N. Entries / 5#circ", 36, 0, 180)

h_pca_proton = TH1F(
    "h_pca_proton", ";PCA;N. Entries / 5#circ", 50, 0.9, 1)
h_pca_nonproton = TH1F(
    "h_pca_nonproton", ";PCA;N. Entries / 5#circ", 50, 0.9, 1)


def is_active(point):
    ok_y = y_start < point[1] < y_end
    ok_x = x_start < point[0] < x_end
    ok_z = z_start < point[2] < z_end
    return ok_y and ok_x and ok_z


def choose_shower(root_chain, plane):
    most_energetic_shower = 0
    shower_id = 0
    for ish in range(root_chain.n_showers):
        if root_chain.shower_energy[ish][plane] < 3:
            if root_chain.shower_energy[ish][plane] > most_energetic_shower:
                most_energetic_shower = root_chain.shower_energy[ish][plane]
                shower_id = ish
    return shower_id


def choose_track(root_chain):
    if root_chain.n_tracks:
        return -1

    track_id = 0
    least_proton_track = 1

    for itrk in range(root_chain.n_tracks):
        if root_chain.predict_p[itrk] < least_proton_track:
            least_proton_track = root_chain.predict_p[itrk]
            track_id = itrk
    return track_id


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

def distance(p1, p2):
    return math.sqrt(
        sum([(t - n)**2 for t, n in zip(p1, p2)]))

def fill_kin_branches(root_chain, weight, variables, option=""):
    no_tracks = False
    longest_track = 0
    longest_track_id = 0
    hit_index = 2#choose_plane(root_chain)
    shower_id = choose_shower(root_chain, hit_index)
    track_id = choose_track(root_chain)


    track_like_shower_id = -1
    if root_chain.n_tracks == 0 and root_chain.n_showers > 1:
        max_pca = 0
        for i, pca in enumerate(root_chain.shower_pca):
            if pca > max_pca and i != shower_id:
                max_pca = pca
                track_like_shower_id = i
        no_tracks = True

    for itrk in range(root_chain.n_tracks):
        if root_chain.track_len[itrk] > longest_track:
            longest_track = root_chain.track_len[itrk]
            longest_track_id = itrk

    if no_tracks:
        v_track = TVector3(
            root_chain.shower_dir_x[track_like_shower_id],
            root_chain.shower_dir_y[track_like_shower_id],
            root_chain.shower_dir_z[track_like_shower_id])
    else:
        v_track = TVector3(
            root_chain.track_dir_x[longest_track_id],
            root_chain.track_dir_y[longest_track_id],
            root_chain.track_dir_z[longest_track_id])

    vec_shower = TVector3(
        root_chain.shower_dir_x[shower_id],
        root_chain.shower_dir_y[shower_id],
        root_chain.shower_dir_z[shower_id])
    costheta_shower_track = v_track.Dot(
        vec_shower) / (v_track.Mag() * vec_shower.Mag())

    signal = 0
    if root_chain.category == 2:
        signal = 1

    if no_tracks:
        track_vertex = [
            root_chain.shower_start_x[track_like_shower_id],
            root_chain.shower_start_y[track_like_shower_id],
            root_chain.shower_start_z[track_like_shower_id]]
    else:
        track_vertex = [
            root_chain.track_start_x[track_id],
            root_chain.track_start_y[track_id],
            root_chain.track_start_z[track_id]]

    if no_tracks:
        length = root_chain.shower_length[track_like_shower_id]
        shower_end_x = root_chain.shower_start_x[track_like_shower_id] + length * root_chain.shower_dir_x[track_like_shower_id]
        shower_end_y = root_chain.shower_start_y[track_like_shower_id] + length * root_chain.shower_dir_y[track_like_shower_id]
        shower_end_z = root_chain.shower_start_z[track_like_shower_id] + length * root_chain.shower_dir_z[track_like_shower_id]

        track_end= [
            shower_end_x,
            shower_end_y,
            shower_end_z]
    else:
        track_end = [
            root_chain.track_end_x[track_id],
            root_chain.track_end_y[track_id],
            root_chain.track_end_z[track_id]]
    shower_vertex = [
        root_chain.shower_start_x[shower_id],
        root_chain.shower_start_y[shower_id],
        root_chain.shower_start_z[shower_id]]

    neutrino_vertex = [root_chain.vx, root_chain.vy, root_chain.vz]

    shower_vertex_d = math.sqrt(
        sum([(s - n)**2 for s, n in zip(shower_vertex, neutrino_vertex)]))
    track_vertex_d = math.sqrt(
        sum([(t - n)**2 for t, n in zip(track_vertex, neutrino_vertex)]))

    variables["is_signal"][0] = signal
    # total_shower_energy = sum([root_chain.shower_energy[i_sh][hit_index] for i_sh in range(root_chain.n_showers)])
    total_shower_energy_y = sum([root_chain.shower_energy[i_sh][2] for i_sh in range(root_chain.n_showers)])

    total_shower_nhits = sum([root_chain.shower_nhits[i_sh][hit_index] for i_sh in range(root_chain.n_showers)])
    total_shower_nhits_y = sum([root_chain.shower_nhits[i_sh][2] for i_sh in range(root_chain.n_showers)])

    variables["n_objects"][0] = root_chain.n_tracks + root_chain.n_showers
    variables["no_tracks"][0] = int(no_tracks)

    if no_tracks:
        variables["track_length"][0] = root_chain.shower_length[track_like_shower_id]
        variables["track_phi"][0] = math.degrees(root_chain.shower_phi[track_like_shower_id])
        variables["track_theta"][0] = math.degrees(root_chain.shower_theta[track_like_shower_id])
        variables["track_start_x"][0] = root_chain.shower_start_x[track_id]
        variables["track_start_y"][0] = root_chain.shower_start_y[track_id]
        variables["track_start_z"][0] = root_chain.shower_start_z[track_id]
        total_track_energy = root_chain.shower_energy[track_like_shower_id][2]
        total_track_energy_y = root_chain.shower_energy[track_like_shower_id][2]
        total_track_energy_dedx = total_track_energy
        total_track_nhits = root_chain.shower_nhits[track_like_shower_id][2]
        total_shower_nhits_y -= root_chain.shower_nhits[track_like_shower_id][2]
        # total_shower_energy -= total_track_energy
        total_shower_energy_y -= total_track_energy_y
        total_shower_nhits -= total_track_nhits
        # total_track_energy_length = length2energy(root_chain.shower_length[track_like_shower_id])
        variables["track_energy"][0] = total_track_energy
        variables["track_energy_dedx"][0] = total_track_energy
        variables["proton_score"][0] = 1
        variables["track_pca"][0] = max(0, root_chain.shower_pca[track_like_shower_id])
        variables["n_tracks"][0] = 1
        variables["n_showers"][0] = root_chain.n_showers - 1
        track_dedx = root_chain.shower_dEdx[track_like_shower_id][2]
        variables["track_dedx"][0] = max(0, track_dedx)
        variables["track_energy_length"][0] = length2energy(root_chain.shower_length[track_like_shower_id])
    else:
        variables["track_length"][0] = sum(root_chain.track_len)
        variables["track_phi"][0] = math.degrees(root_chain.track_phi[track_id])
        variables["track_theta"][0] = math.degrees(root_chain.track_theta[track_id])
        variables["track_start_x"][0] = root_chain.track_start_x[track_id]
        variables["track_start_y"][0] = root_chain.track_start_y[track_id]
        variables["track_start_z"][0] = root_chain.track_start_z[track_id]
        # total_track_energy_length = sum(
        #     [length2energy(root_chain.track_len[i_tr]) for i_tr in range(root_chain.n_tracks)])
        total_track_energy = sum(
            [root_chain.track_energy_hits[i_tr][2] for i_tr in range(root_chain.n_tracks)])
        total_track_energy_dedx = sum(
            [root_chain.track_energy_dedx[i_tr] for i_tr in range(root_chain.n_tracks)])
        total_track_nhits = sum(
            [root_chain.track_nhits[i_tr][2] for i_tr in range(root_chain.n_tracks)])
        variables["track_energy"][0] = root_chain.track_energy_hits[track_id][2]
        variables["track_energy_dedx"][0] = root_chain.track_energy_dedx[track_id]
        variables["proton_score"][0] = max(0, root_chain.predict_p[track_id])
        variables["track_pca"][0] = max(0, root_chain.track_pca[track_id])
        variables["n_tracks"][0] = root_chain.n_tracks
        variables["n_showers"][0] = root_chain.n_showers
        track_dedx = root_chain.track_dEdx[track_id][2]
        variables["track_dedx"][0] = max(0, track_dedx)
        variables["track_energy_length"][0] = length2energy(root_chain.track_len[track_id])
    
    variables["shower_length"][0] = root_chain.shower_length[shower_id]
    variables["track_end_x"][0] = track_end[0]
    variables["track_end_y"][0] = track_end[1]
    variables["track_end_z"][0] = track_end[2]
    variables["track_distance"][0] = track_vertex_d

    variables["shower_start_x"][0] = root_chain.shower_start_x[shower_id]
    # variables["shower_end_x"][0] = root_chain.shower_end_x[shower_id]
    variables["shower_start_y"][0] = root_chain.shower_start_y[shower_id]
    # variables["shower_end_y"][0] = root_chain.shower_end_y[shower_id]
    variables["shower_start_z"][0] = root_chain.shower_start_z[shower_id]
    # variables["shower_end_z"][0] = root_chain.shower_end_z[shower_id]

    total_dedx_hits = []
    for ish in range(root_chain.n_showers):
        if ish != track_like_shower_id:
            total_dedx_hits += root_chain.shower_dEdx_hits[ish]

    dedx_y = 0
    if len(total_dedx_hits):
        dedx_y = statistics.median(total_dedx_hits)

    # variables["total_shower_energy"][0] = total_shower_energy
    variables["total_shower_energy_y"][0] = total_shower_energy_y
    variables["total_track_energy"][0] = total_track_energy
    # variables["total_track_energy_length"][0] = total_track_energy_length
    variables["total_track_energy_dedx"][0] = total_track_energy_dedx
    variables["track_hits"][0] = total_track_nhits
    variables["shower_hits"][0] = total_shower_nhits
    variables["shower_hits_y"][0] = total_shower_nhits_y
    variables["shower_energy"][0] = root_chain.shower_energy[shower_id][hit_index]
    variables["shower_energy_y"][0] = root_chain.shower_energy[shower_id][2]
    variables["shower_theta"][0] = math.degrees(root_chain.shower_theta[shower_id])
    variables["shower_phi"][0] = math.degrees(root_chain.shower_phi[shower_id])
    variables["shower_distance"][0] = shower_vertex_d

    y_shower = sum([root_chain.shower_nhits[i_sh][2] for i_sh in range(root_chain.n_showers)])
    y_track = sum([root_chain.track_nhits[i_sh][2]for i_sh in range(root_chain.n_tracks)])

    u_shower = sum([root_chain.shower_nhits[i_sh][0] for i_sh in range(root_chain.n_showers)])
    u_track = sum([root_chain.track_nhits[i_sh][0] for i_sh in range(root_chain.n_tracks)])

    v_shower = sum([root_chain.shower_nhits[i_sh][1] for i_sh in range(root_chain.n_showers)])
    v_track = sum([root_chain.track_nhits[i_sh][1] for i_sh in range(root_chain.n_tracks)])

    variables["total_hits"][0] = y_shower + y_track + u_shower + u_track + v_shower + v_track
    variables["total_hits_u"][0] = u_shower + u_track
    variables["total_hits_v"][0] = v_shower + v_track
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
    variables["pt"][0] = pt_plot(root_chain, 2)
    variables["track_shower_angle"][0] = costheta_shower_track
    variables["event"][0] = root_chain.event
    variables["run"][0] = root_chain.run
    variables["subrun"][0] = root_chain.subrun
    variables["interaction_type"][0] = root_chain.interaction_type
    variables["shower_open_angle"][0] = math.degrees(root_chain.shower_open_angle[shower_id])
    variables["shower_pca"][0] = max(0, root_chain.shower_pca[shower_id])
    variables["reco_energy"][0] = 0

    total_shower_energy = 0
    total_track_energy_length = 0
    max_pca = 0
    for i_sh in range(root_chain.n_showers):
        shower_vertex = [root_chain.shower_start_x[i_sh],
                         root_chain.shower_start_y[i_sh],
                         root_chain.shower_start_z[i_sh]]

        shower_d = distance(shower_vertex, neutrino_vertex)

        v_sh = TVector3(
            root_chain.shower_dir_x[i_sh],
            root_chain.shower_dir_y[i_sh],
            root_chain.shower_dir_z[i_sh])

        cos = v_sh.Dot(vec_shower) / (v_sh.Mag() * vec_shower.Mag())
        shower_angle = math.degrees(math.acos(min(cos, 1)))

        if root_chain.shower_pca[i_sh] > max_pca:
            max_pca = root_chain.shower_pca[i_sh]

        if i_sh != shower_id and (root_chain.matched_showers[i_sh] == 2212 or root_chain.matched_showers[i_sh] == 13):
            h_angle_proton.Fill(shower_angle)
            # h_pca_proton.Fill(root_chain.track_pca[i_tr])
        elif i_sh != shower_id:
            h_angle_nonproton.Fill(shower_angle)
            # h_pca_nonproton.Fill(root_chain.track_pca[i_tr])


        if i_sh != shower_id and shower_angle > 15:
            length_e = length2energy(root_chain.shower_length[i_sh])
            variables["reco_energy"][0] += length_e
            total_track_energy_length += length_e
        else:
            variables["reco_energy"][0] += root_chain.shower_energy[i_sh][hit_index]
            total_shower_energy += root_chain.shower_energy[i_sh][hit_index]

    for i_tr in range(root_chain.n_tracks):
        if root_chain.track_pca[i_tr] < max_pca:
            if root_chain.matched_tracks[i_tr] == 2212 or root_chain.matched_tracks[i_tr] == 13:
                h_pca_proton.Fill(root_chain.track_pca[i_tr])
            else:
                h_pca_nonproton.Fill(root_chain.track_pca[i_tr])


        if root_chain.track_pca[i_tr] < max_pca:
            variables["reco_energy"][0] += root_chain.track_energy_hits[i_tr][hit_index]
            total_shower_energy += root_chain.track_energy_hits[i_tr][hit_index]
        else:
            length_e = length2energy(root_chain.track_len[i_tr])
            variables["reco_energy"][0] += length_e
            total_track_energy_length += length_e

    variables["total_shower_energy"][0] = total_shower_energy
    variables["total_track_energy_length"][0] = total_track_energy_length

    variables["numu_score"][0] = root_chain.numu_cuts

    # dE/dx for the collection plane only
    dedx = root_chain.shower_dEdx[shower_id][hit_index]
    variables["dedx"][0] = max(0, dedx)
    variables["dedx_y"][0] = dedx_y


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


def fill_tree(chain, weight, tree, option=""):
    total_events = 0
    total_entries = int(chain.GetEntries() / 1)

    for i in range(total_entries):
        chain.GetEntry(i)
        printProgressBar(i, total_entries, prefix = 'Progress:', suffix = 'Complete', length = 20)

        if chain.passed:
            hit_index = choose_plane(chain)

            track_fidvol = True

            for i_tr in range(chain.n_tracks):
                track_start = [
                    chain.track_start_x[i_tr],
                    chain.track_start_y[i_tr],
                    chain.track_start_z[i_tr]]
                track_end = [
                    chain.track_end_x[i_tr],
                    chain.track_end_y[i_tr],
                    chain.track_end_z[i_tr]]

                track_fidvol = track_fidvol and is_fiducial(track_start) and is_fiducial(track_end)              

            shower_fidvol = True

            for i_sh in range(chain.n_showers):
                shower_start = [
                    chain.shower_start_x[i_sh],
                    chain.shower_start_y[i_sh],
                    chain.shower_start_z[i_sh]]

                shower_fidvol = shower_fidvol and is_fiducial(shower_start)
                

            option_check = True
            event_weight = weight

            if option == "bnb":
                option_check = abs(chain.nu_pdg) != 12
            if abs(chain.nu_pdg) == 12:
                event_weight = weight * chain.bnbweight
            if "nue" in option:
                option_check = abs(chain.nu_pdg) == 12

            neutrino_vertex = [chain.vx, chain.vy, chain.vz]
            
            # If there are no tracks we require at least two showers
            showers_2_tracks_0 = True
            if chain.n_tracks == 0 and chain.n_showers == 1:
                showers_2_tracks_0 = False

            if option_check and is_fiducial(neutrino_vertex) and track_fidvol and shower_fidvol and showers_2_tracks_0:
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
data_ext_scaling_factor = 0.142

samples = ["cosmic_mc", "nue", "bnb", "bnb_data", "ext_data", "lee"]

tree_files = [glob("cosmic_intime_mcc86/*/a*.root"),
              glob("mc_nue_ubxsec/*.root"),
              glob("mc_bnb_ubxsec/*/*.root"),
              glob("data_bnb_ubxsec/*/*.root"),
              glob("data_ext_ubxsec/*/*.root"),
              glob("lee/*.root")]

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
# pots_dict["nue"] = 1
# pots_dict["bnb"] = 1

weights = [1,#data_ext_scaling_factor * total_pot / total_data_bnb_pot * ext_mc_scaling *
           #chains_dict["ext_data"].GetEntries() / chains_dict["cosmic_mc"].GetEntries(),
           total_pot / pots_dict["nue"],
           total_pot / pots_dict["bnb"],
           1,
           data_ext_scaling_factor,
           total_pot / pots_dict["lee"]]

files = ["cosmic_mc_file.root", "nue_file.root",
         "mc_file.root", "bnb_file.root", "bnbext_file.root", "lee_file.root"]
tree_names = ["cosmic_mc_tree", "nue_tree",
              "mc_tree", "bnb_tree", "bnbext_tree", "lee_tree"]

trees = []

for t in tree_names:
    trees.append(TTree(t, t))

for n, b in variables.items():
    for t in trees:
        t.Branch(n, b, n + "/f")

samples = ["cosmic_mc", "nue", "bnb", "bnb_data", "ext_data", "lee"]


for i, s in enumerate(samples):
    start_time = time.time()
    print("******************************")
    print("Sample", s)
    print("Weight", weights[i])
    events = fill_tree(chains[i], weights[i], trees[i], s)
    print("\nEvents", events)
    print("Time to fill %.1f s"  % (time.time() - start_time))

for f, t in zip(files, trees):
    tfile = TFile(f, "RECREATE")
    t.Write()
    tfile.Close()


# Files needed for Zarko's POT counting tool
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
print("Total time %.2f" % (time.time() - begin_time))

c_angle = TCanvas("c_angle")

h_ratio = h_angle_proton.Clone()
h_den = h_angle_nonproton.Clone()
h_den.Add(h_angle_proton)
h_ratio.Divide(h_den)

h_angle_nonproton.Draw("hist")
h_angle_proton.Draw("hist same")
h_angle_proton.SetLineColor(kRed + 1)
h_angle_nonproton.SetLineColor(kBlue + 1)

c_angle.Update()

c_pca = TCanvas("c_pca")
h_pca_nonproton.Draw("hist")
h_pca_proton.Draw("hist same")
c_pca.Update()

input()
