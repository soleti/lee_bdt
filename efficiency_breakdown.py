#!/usr/bin/env python3.4

import ROOT
from glob import glob

nue_cosmic = sorted(glob("nue_efficiency/*/Pandora*.root"))[:1]
chain = ROOT.TChain("robertoana/pandoratree")
print(nue_cosmic)
for f in nue_cosmic:
    chain.Add(f)
entries = chain.GetEntries()
total = 0
splitted_event = 0
incomplete_event = 0
perfect_event = 0
wrong_event = 0
track_ok_shower_mis = 0
track_ok_shower_no = 0
shower_ok_track_mis = 0
shower_ok_track_no = 0
track_no_shower_no = 0
flash_not_passed = 0
for evt in range(entries):
    chain.GetEntry(evt)

    p = 0
    e = 0
    photons = 0
    pions = 0
    for i, energy in enumerate(chain.nu_daughters_E):
        if abs(chain.nu_daughters_pdg[i]) == 2212:
            if energy - 0.938 > 0.05:
                p += 1

        if abs(chain.nu_daughters_pdg[i]) == 11:
            if energy > 0.03:
                e += 1

        if chain.nu_daughters_pdg[i] == 22:
            # if energy > 0.035:
            photons += 1

        if chain.nu_daughters_pdg[i] == 111:
            # if energy > 0.06:
            pions += 1

    eNp = e == 1 and photons == 0 and pions == 0 and p > 0

    if eNp and chain.nu_E > 0.1:

        if chain.true_nu_is_fiducial:
            total += 1
            primary_indexes = []
            shower_passed = []
            track_passed = []
            flash_passed = []

            for i in range(chain.n_primaries):
                primary_indexes.append(chain.primary_indexes[i])
                shower_passed.append(chain.shower_passed[i])
                track_passed.append(chain.track_passed[i])
                flash_passed.append(chain.flash_passed[i] + 1)

            if chain.passed:
                candidate_id = primary_indexes.index(chain.chosen_candidate)

                chosen_showers = shower_passed[candidate_id]
                chosen_tracks = track_passed[candidate_id]

                if chosen_showers == e and chosen_tracks == p:
                    perfect_event += 1
                elif chosen_showers > 0 or chosen_tracks > 0:
                    if chosen_showers < e or chosen_tracks < p:
                        incomplete_event += 1
                    elif chosen_showers >= e or chosen_tracks >= p:
                        splitted_event += 1
                else:
                    wrong_event += 1
            else:
                find_track = False
                find_shower = False
                if 1 in flash_passed:
                    for i in range(chain.n_primaries):
                        if track_passed[i] > 0:
                            find_track = True
                            if track_passed[i] > 1:
                                track_ok_shower_mis += 1
                                break
                            else:
                                track_ok_shower_no += 1
                                break
                        if shower_passed[i] > 0:
                            find_shower = True
                            if shower_passed[i] > 1:
                                shower_ok_track_mis += 1
                                break
                            else:
                                shower_ok_track_no += 1
                                break

                    if not find_track and not find_shower:
                        track_no_shower_no += 1

                else:
                    flash_not_passed += 1


print("Passed event, perfect {:.1f} %".format(perfect_event / total * 100))
print("Passed event, incomplete {:.1f} %"
      .format(incomplete_event / total * 100))
print("Passed event, splitted {:.1f} %".format(splitted_event / total * 100))
print("Passed event, wrong {:.1f} %".format(wrong_event / total * 100))
print("Not passed event, flash not passed {:.1f}% "
      .format(flash_not_passed / total * 100))
print("Not passed event, flash passed, track ok, shower misid. {:.1f}% "
      .format(track_ok_shower_mis / total * 100))
print("Not passed event, flash passed, shower ok, track misid. {:.1f}% "
      .format(shower_ok_track_mis / total * 100))
print("Not passed event, flash passed, track ok, no shower {:.1f}% "
      .format(track_ok_shower_no / total * 100))
print("Not passed event, flash passed, shower ok, no track {:.1f}% "
      .format(shower_ok_track_no / total * 100))
print("Not passed event, flash passed, no track, no shower {:.1f}% "
      .format(track_no_shower_no / total * 100))
