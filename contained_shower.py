#!/usr/local/bin/python3

import math
import ROOT
from glob import glob
from array import array
import statistics
from bdt_common import x_start, y_start, z_start, x_end, y_end, z_end, is_fiducial

ROOT.gStyle.SetOptFit(0)

def profile_style(profile):
    profile.GetYaxis().SetRangeUser(0, 1)
    # profile.Fit("pol9")
    profile.SetFillStyle(3002)
    profile.SetFillColor(ROOT.kBlack)
    profile.SetLineColor(1)


ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetPalette(87)
ROOT.gStyle.SetNumberContours(99)


nue_cosmic = glob("mc_nue_mcc86_2/*/Pandora*.root")

chain_nue = ROOT.TChain("robertoana/pandoratree")

for f in nue_cosmic:
    chain_nue.Add(f)

entries = chain_nue.GetEntries()
print("entries", entries)
ELECTRON_THRESHOLD = 0.035
ELECTRON_MASS = 0.51e-3
E_C = 0.035
X_0 = 14
h_angle_e = ROOT.TH2F("h_angle_e","",40,0,1,40,0,1)
p_angle_e = ROOT.TProfile("p_angle_e","",20,0,1,0,1)
total = 0
outside = 0
h_e_diff = ROOT.TH1F(
    "h_e_diff", ";E_{reco}/E_{true};Fraction", 28, 0, 1.4)
h_e_diff_e = ROOT.TH2F(
    "h_e_diff_e", ";E_{reco}/E_{true};Fraction", 40, 0, 2, 20, 0, 1.4)
p_e_diff_e = ROOT.TProfile(
    "p_e_diff_e", ";E_{true} [GeV];E_{reco}/E_{true}", 20, 0, 2, 0, 1.4, "s")
p_e_diff_x = ROOT.TProfile(
    "p_e_diff_x", ";#Deltax [cm];E_{reco}/E_{true}", 20, 0, x_end / 2, 0, 1.4, "s")
p_e_diff_y = ROOT.TProfile(
    "p_e_diff_y", ";#Deltay [cm];E_{reco}/E_{true}", 20, 0, y_end, 0, 1.4, "s")
p_e_diff_z = ROOT.TProfile(
    "p_e_diff_z", ";#Deltaz [cm];E_{reco}/E_{true}", 20, 0, z_end / 2, 0, 1.4, "s")
p_e_diff_d = ROOT.TProfile(
    "p_e_diff_d", ";#Deltad [cm];E_{reco}/E_{true}", 20, 0, 600, 0, 1.4, "s")

h_x_y_z = ROOT.TH3F("h_x_y_z", ";#Deltax [cm];#Deltay [cm];#Deltaz [cm]", 5, 0, x_end/2, 5, 0, y_end, 5, 0, z_end / 2)
h_x_y_z_entries = ROOT.TH3F("h_x_y_z_entries", "", 5, 0, x_end / 2,
                    5, 0, y_end, 5, 0, z_end / 2)

e_containment = ROOT.TEfficiency("e_containment", ";e^{-} energy [GeV];Contained events", 40, 0.035, 1)

for i in range(entries):
    if i % 1000 == 0: print(i)
    chain_nue.GetEntry(i)
    electrons = 0
    electron_index = 0
    electron_energy = 0
    for i_pdg, energy in enumerate(chain_nue.nu_daughters_E):

        if chain_nue.nu_daughters_pdg[i_pdg] == 11:
            electron_energy = energy
            if energy - 0.51e-3 > ELECTRON_THRESHOLD:
                electrons += 1
                electron_index = i_pdg
                break


    if electrons > 0 and chain_nue.passed and chain_nue.category == 2:

        electron_starting_point = [chain_nue.nu_daughters_vx[electron_index],
                                   chain_nue.nu_daughters_vy[electron_index],
                                   chain_nue.nu_daughters_vz[electron_index]]

        if is_fiducial(electron_starting_point):
            neutrino_vertex = [chain_nue.vx, chain_nue.vy, chain_nue.vz]

            p = math.sqrt(electron_energy**2 - ELECTRON_MASS**2)
            # print(p, chain_nue.nu_daughters_px[electron_index])
            electron_direction = [chain_nue.nu_daughters_px[electron_index]/p,
                                    chain_nue.nu_daughters_py[electron_index]/p,
                                    chain_nue.nu_daughters_pz[electron_index]/p]

            t_max = math.log(electron_energy/E_C) - 1
            length = t_max * X_0
            electron_end_point = [s+length*d for s, d in zip(electron_starting_point, electron_direction)]
            # print(electron_end_point)
            total += 1
            # print(electron_energy, length, is_fiducial(electron_end_point))
            shower_fidvol = True
            electron_shower = False

            shower_vertex = [
                chain_nue.shower_start_x[0],
                chain_nue.shower_start_y[0],
                chain_nue.shower_start_z[0]]
            shower_dir = [
                chain_nue.shower_dir_x[0],
                chain_nue.shower_dir_y[0],
                chain_nue.shower_dir_z[0]
            ]

            shower_end = [s + d * chain_nue.shower_length[0]
                            for s, d in zip(shower_vertex, shower_dir)]

            if chain_nue.n_showers == 1 and chain_nue.matched_showers[0] == 11 and chain_nue.shower_nhits[0][2] > 5 and chain_nue.shower_energy[0][2] > 0.02 and chain_nue.category == 2:

                closest_x, closest_y, closest_z = x_end, y_end, z_end

                if abs(shower_end[0] - x_start) <= abs(shower_end[0] - x_end):
                    closest_x = x_start
                if abs(shower_end[1] - y_start) <= abs(shower_end[1] - y_end):
                    closest_x = y_start
                if abs(shower_end[2] - z_start) <= abs(shower_end[2] - z_end):
                    closest_x = z_start

                distance_from_x = abs(shower_vertex[0] - closest_x)
                distance_from_y = abs(shower_vertex[1] - closest_y)
                distance_from_z = abs(shower_vertex[2] - closest_z)



                contained_fraction = 1 - \
                    (electron_energy - chain_nue.shower_energy[0][2]) / electron_energy

                distance_walls = math.sqrt(distance_from_x**2+distance_from_y**2+distance_from_z**2)
                
                shower_vertex_d = math.sqrt(
                    sum([(s - n)**2 for s, n in zip(shower_vertex, neutrino_vertex)]))

                shower_vertex_d2 = math.sqrt(
                    sum([(s - n)**2 for s, n in zip(shower_end, neutrino_vertex)]))
                    

                h_e_diff.Fill(contained_fraction)
                h_e_diff_e.Fill(electron_energy, contained_fraction)
                p_e_diff_e.Fill(electron_energy, contained_fraction)
                p_e_diff_x.Fill(distance_from_x, contained_fraction)
                p_e_diff_y.Fill(distance_from_y, contained_fraction)
                p_e_diff_z.Fill(distance_from_z, contained_fraction)
                p_e_diff_d.Fill(distance_walls, contained_fraction)
                h_x_y_z.Fill(distance_from_x, distance_from_y, distance_from_z, contained_fraction)
                h_x_y_z_entries.Fill(distance_from_x, distance_from_y, distance_from_z)

            for i in range(chain_nue.n_showers):
                shower_start = [
                    chain_nue.shower_start_x[i],
                    chain_nue.shower_start_y[i],
                    chain_nue.shower_start_z[i]]
                shower_dir = [
                    chain_nue.shower_dir_x[i],
                    chain_nue.shower_dir_y[i],
                    chain_nue.shower_dir_z[i]
                ]

                shower_pdg = chain_nue.matched_showers[i]
                shower_end = [s+d*chain_nue.shower_length[i] for s, d in zip(shower_start, shower_dir)]

                # Check if the electron-matched shower is within the fiducial volume
                if shower_pdg == 11:
                    electron_shower = True
                    shower_fidvol = shower_fidvol and is_fiducial(shower_start)# and is_fiducial(shower_end)

            # Plot percentage of contained shower as a function of the energy)
            e_containment.Fill(not(not is_fiducial(electron_end_point) and shower_fidvol and electron_shower), electron_energy)
            if not is_fiducial(electron_end_point) and shower_fidvol and electron_shower:
                outside += 1

h_x_y_z.Divide(h_x_y_z_entries)
c_containment = ROOT.TCanvas("c_containment")
h_e_diff.Scale(1/h_e_diff.Integral())
h_e_diff.Draw("hist")
h_e_diff.SetLineColor(ROOT.kBlue + 1)
h_e_diff.SetFillColor(ROOT.kAzure + 1)
h_e_diff.GetXaxis().SetTitleOffset(1.5)
c_containment.SetLeftMargin(0.155)

c_containment.Update()
print(outside/total)



c_diff_e = ROOT.TCanvas("c_diff_e")
profile_style(p_e_diff_e)
h_e_diff_e.Draw("colz")
p_e_diff_e.Draw("e2 same")
p_e_diff_e2 = p_e_diff_e.Clone()
p_e_diff_e2.SetFillStyle(0)
p_e_diff_e2.Draw("hist same")

l = ROOT.TLegend(0.696, 0.911, 0.953, 0.985)
l.AddEntry(p_e_diff_e, "1#sigma interval", "fl")
l.Draw()
c_diff_e.SetLeftMargin(0.165)
c_diff_e.Update()

c_diffx = ROOT.TCanvas("c_diffx")
profile_style(p_e_diff_x)
p_e_diff_x.Draw("e2")
p_e_diff_x2 = p_e_diff_x.Clone()
p_e_diff_x2.SetFillStyle(0)
p_e_diff_x2.Draw("hist same")
c_diffx.SetLeftMargin(0.165)
l.Draw()
c_diffx.Update()

c_diffy = ROOT.TCanvas("c_diffy")
profile_style(p_e_diff_y)
p_e_diff_y.Draw("e2")
p_e_diff_y2 = p_e_diff_y.Clone()
p_e_diff_y2.SetFillStyle(0)
p_e_diff_y2.Draw("hist same")
c_diffy.SetLeftMargin(0.165)
l.Draw()
c_diffy.Update()

c_diffz = ROOT.TCanvas("c_diffz")
profile_style(p_e_diff_z)
p_e_diff_z.Draw("e2")
p_e_diff_z2 = p_e_diff_z.Clone()
p_e_diff_z2.SetFillStyle(0)
p_e_diff_z2.Draw("hist same")
l.Draw()
c_diffz.SetLeftMargin(0.165)
c_diffz.Update()

c_diffd = ROOT.TCanvas("c_diffd")
profile_style(p_e_diff_d)
p_e_diff_d.Draw("e2")
p_e_diff_d2 = p_e_diff_d.Clone()
p_e_diff_d2.SetFillStyle(0)
p_e_diff_d2.Draw("hist same")

c_diffd.SetLeftMargin(0.165)
c_diffd.Update()

c_3d = ROOT.TCanvas("c_3d")
h_x_y_z.SetLineColor(1)
h_x_y_z.SetLineWidth(1)
h_x_y_z.SetBinContent(h_x_y_z.GetBin(3, 3, 3), 1)
h_x_y_z.SetBinContent(h_x_y_z.GetBin(4, 4, 4), 0)
h_x_y_z.GetXaxis().SetTitleOffset(1.3)
h_x_y_z.Draw("BOX2")
c_3d.Update()
input()
