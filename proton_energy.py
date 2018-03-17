#!/usr/local/bin/python3

import ROOT
from array import array

ARGON_DENSITY = 1.379


with open("proton_range.txt", "r") as table_file:
    lines = table_file.readlines()
    a_range = array("f", [])
    a_energy = array("f", [])

    for line in lines:
        values = line.split()
        a_range.append(float(values[1]))
        a_energy.append(float(values[0]))

    g_proton_range = ROOT.TGraph(len(lines), a_range, a_energy)


def length2energy(length):
    """Return the energy (in GeV) of the track based on its length and
    the stopping power in the argon"""
    return g_proton_range.Eval(length * ARGON_DENSITY, 0, "S") / 1000
    
