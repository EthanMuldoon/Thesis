#Standard setup stuff
import ROOT as R
RF = R.RooFit
import sys
import os
sys.path.append(os.environ["NEUTRONS"]+'/RooFit')

from pyUtils import *
outputDir = "../figures/"
#fPDF = pdf(outputDir+"kinetics_wall_electron_neutron_comparison.pdf")

from Kinetics import Kinetics
k = Kinetics("k")
k.SetTitle("Muon Kinetics Signals")
k.frame.GetXaxis().SetTitle("time [ns]")
k.frame.GetYaxis().SetTitle("signal")

can = R.TCanvas("can")
k.DrawComponent(k.electron, [RF.LineColor(R.kBlue), RF.Name("el_d")])
k.DrawComponent(k.capture, [RF.LineColor(R.kRed), RF.Name("neu_d")])
k.L_d = 0.01
k.L_q = 0.01
k.DrawComponent(k.electron, [RF.LineColor(R.kBlue), RF.LineStyle(2), RF.Name("el_w")])
k.DrawComponent(k.capture, [RF.LineColor(R.kRed), RF.LineStyle(2), RF.Name("neu_w")])

can.SetLogy()
k.frame.GetXaxis().SetRangeUser(0,4000)
k.frame.GetYaxis().SetRangeUser(1e-7,1e-2)

leg = R.TLegend(0.5,0.7,0.9,0.9)
leg.AddEntry(k.frame.findObject("el_d"), "Decay Electrons (Deuterium)", "l")
leg.AddEntry(k.frame.findObject("neu_d"), "Capture Neutrons (Deuterium)", "l")
leg.AddEntry(k.frame.findObject("el_w"), "Decay Electrons (Wall Stop)", "l")
leg.AddEntry(k.frame.findObject("neu_w"), "Capture Neutrons (Wall Stop)", "l")
leg.Draw()

can.Print(outputDir+"kinetics_wall_electron_neutron_comparison.pdf")

