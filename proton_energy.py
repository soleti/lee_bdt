#!/usr/local/bin/python3

import ROOT
from array import array

ARGON_DENSITY = 1.4


def length2energy(length):
    """Return the energy (in GeV) of the track based on its length and
    the stopping power in the argon"""
    return g_proton_range.Eval(length * ARGON_DENSITY, 0, "S") / 1000


with open("config_files/proton_range.txt", "r") as table_file:
    lines = table_file.readlines()
    a_range = array("f", [])
    a_energy = array("f", [])

    for line in lines:
        values = line.split()
        a_range.append(float(values[1]))
        a_energy.append(float(values[0]))

    g_proton_range = ROOT.TGraph(len(lines), a_range, a_energy)


if __name__ == "__main__":
    c_proton_range = ROOT.TCanvas("c_proton_range")
    g_proton_range.Draw("AL")
    g_proton_range.GetYaxis().SetTitle("Kinetic energy [GeV]")
    g_proton_range.GetXaxis().SetTitle("Range [g/cm^{2}]")
    g_proton_range.GetXaxis().SetRangeUser(0, 100)
    g_proton_range.GetYaxis().SetRangeUser(0, 1000)

    c_proton_range.Update()
    input()
