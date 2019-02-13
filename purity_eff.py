from array import array
import ROOT

ROOT.gStyle.SetOptStat(0)

tot = array("f", [3259, 1877, 1123, 683, 205, 180, 75, 17.2, 16, 15.2])
nue = array("f", [15.2, 14.6, 13.7, 12.5, 6.4*3.3/2.9, 5.9*3.3/2.9, 4.2*3.3/2.9, 3.1*3.3/2.9, 3.0*3.3/2.9, 2.9*3.3/2.9])
x_bins = array("f", [0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5])
purity = array("f", [])
eff = array("f", [])
descriptions = ["Selection", "N. of hits", "Shower energy", "Hits ratio", "Shower dE/dx",
                "Track d.", "Shower d.", "Proton #chi^{2}", "Track-shower angle", "Track length"]


for i, n in enumerate(nue):
    t = tot[i]
    purity.append(n/t)
    print(n/t, descriptions[i])
    eff.append(n/33.8)


h_axis = ROOT.TH1F("h","",10,0,10)
h_axis.GetYaxis().SetTitle("Fraction")
for i in range(1, h_axis.GetNbinsX()+1):
    h_axis.GetXaxis().SetBinLabel(i, descriptions[i-1])
g_purity = ROOT.TGraph(len(tot), x_bins, purity)
g_eff = ROOT.TGraph(len(tot), x_bins, eff)

c_effpurity = ROOT.TCanvas("c_effpurity")
h_axis.Draw()
g_purity.Draw("LP SAME")
g_purity.GetYaxis().SetTitle("Fraction")
h_axis.GetYaxis().SetTitleOffset(0.8)
g_eff.Draw("LP SAME")

g_purity.SetLineWidth(3)
g_eff.SetLineWidth(3)
g_purity.SetLineColor(ROOT.kBlue+1)
g_eff.SetLineColor(ROOT.kRed+1)
h_axis.GetYaxis().SetRangeUser(0.002, 1)

leg = ROOT.TLegend(0.6, 0.75, 0.88, 0.86)
leg.AddEntry(g_eff, "Efficiency", "l")
leg.AddEntry(g_purity, "Purity", "l")
leg.Draw()
c_effpurity.SetLogy()
c_effpurity.Update()
c_effpurity.SaveAs("plots/purity_eff.pdf")
input()
