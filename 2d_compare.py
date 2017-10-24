import ROOT
ROOT.gStyle.SetOptStat(0)
mc = ROOT.TFile("2d_mc.root")
h_xy_track_start_mc = mc.Get("h_xy_track_start")

data = ROOT.TFile("2d_data.root")
h_xy_track_start_data = data.Get("h_xy_track_start")

c_xy = ROOT.TCanvas("c_xy")
h_xy_track_start_mc.Scale(h_xy_track_start_data.Integral() /
                          h_xy_track_start_mc.Integral())
h_xy_track_start_data.Divide(h_xy_track_start_mc)
h_xy_track_start_data.Draw("colz texte")
h_xy_track_start_data.GetZaxis().SetRangeUser(0, 2)
c_xy.Update()
input()
