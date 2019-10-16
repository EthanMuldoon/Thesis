#Standard setup stuff
import ROOT as R
RF = R.RooFit
import sys
import os
sys.path.append(os.environ["NEUTRONS"]+'/RooFit')

from pyUtils import *
outputDir = "../figures/"
can = R.TCanvas("can")
can.SetLogy()
tRange = [-10000,10000]

from Neutrons import Neutrons
m = Neutrons("m")

m.neu.ch0.clearFrame()
m.neu.ch0.capture.plotOn(m.neu.ch0.frame, R.RooFit.LineColor(R.kBlack), R.RooFit.Precision(1e-5))
m.neu.ch0.frame.SetTitle("Capture Neutron Signal")
m.neu.ch0.frame.GetXaxis().SetTitle("time [ns]")
m.neu.ch0.frame.GetYaxis().SetTitle("")
m.neu.ch0.frame.GetXaxis().SetRangeUser(tRange[0],tRange[1])
m.neu.ch0.frame.GetYaxis().SetRangeUser(5e-5,4e-3)
m.neu.ch0.frame.Draw()
can.Print(outputDir+"shape_capture.pdf")

m.neu.ch0.clearFrame()
m.neu.ch0.fusion.plotOn(m.neu.ch0.frame, R.RooFit.LineColor(R.kBlack), R.RooFit.Precision(1e-5))
m.neu.ch0.frame.SetTitle("Fusion Neutron Signal")
m.neu.ch0.frame.GetXaxis().SetTitle("time [ns]")
m.neu.ch0.frame.GetYaxis().SetTitle("")
m.neu.ch0.frame.GetXaxis().SetRangeUser(tRange[0],tRange[1])
m.neu.ch0.frame.GetYaxis().SetRangeUser(1e-5,3e-2)
m.neu.ch0.frame.Draw()
can.Print(outputDir+"shape_fusion.pdf")


hists1a = R.TFile("/scratch/ethan/mta/pass1809/dsn1a_hists_el.root")
norm1a = R.TFile("/scratch/ethan/mta/pass1809/dsn1a.root")
hNorm = norm1a.Get("hNeuBiasStops")

print("Analyzing Fiducial Stops")
m.SetNormalization(hNorm, 1, 4, 15, 55, 2, 7)
h3d = hists1a.Get("hNeuTimeVsStopY_HighE")
h3d.GetZaxis().SetRange(h3d.GetZaxis().FindBin(15),h3d.GetZaxis().FindBin(56)-1)
hNeu = h3d.Project3D("yx")
h3d = hists1a.Get("hGammaTimeVsStopY")
h3d.GetZaxis().SetRange(h3d.GetZaxis().FindBin(15),h3d.GetZaxis().FindBin(56)-1)
hGamma = h3d.Project3D("yx")

hRF = hists1a.Get("hGammaRF")
m.fitRF(hRF)

m.SetTitle("Fiducial Volume Fit")
m.Free(["eff_background","eff_signal"])
m.fitTo(hNeu, hGamma)

m.Free(["amp_kick"])
m.doFit()

m.SetRangeX(-40,100)
m.gamma.Free(["time_start", "res_expo"])
m.doFit()

m.setBinned(False)
can.cd()
can.SetLogy(False)
m.neu.ch0.clearFrame()
m.neu.ch0.background.plotOn(m.neu.ch0.frame, R.RooFit.LineColor(R.kBlack), R.RooFit.Precision(1e-5))
m.neu.ch0.frame.SetTitle("Background Neutron Signal")
m.neu.ch0.frame.GetXaxis().SetTitle("time [ns]")
m.neu.ch0.frame.GetYaxis().SetTitle("")
m.neu.ch0.frame.GetXaxis().SetRangeUser(tRange[0],tRange[1])
m.neu.ch0.frame.GetYaxis().SetRangeUser(0.001,0.0021)
m.neu.ch0.frame.Draw()
can.Print(outputDir+"shape_neuBG.pdf")

rfAmp = m.gamma.ch0.rf_amplitude.getVal()
m.gamma.ch0.rf_amplitude = 0
m.neu.ch0.clearFrame()
m.gamma.ch0.background.plotOn(m.neu.ch0.frame, R.RooFit.LineColor(R.kBlack), R.RooFit.Precision(1e-5))
m.neu.ch0.frame.SetTitle("Background Gamma Signal")
m.neu.ch0.frame.GetXaxis().SetTitle("time [ns]")
m.neu.ch0.frame.GetYaxis().SetTitle("")
m.neu.ch0.frame.GetXaxis().SetRangeUser(tRange[0],tRange[1])
m.neu.ch0.frame.GetYaxis().SetRangeUser(0.0014,0.00155)
m.neu.ch0.frame.Draw()
can.Print(outputDir+"shape_gammaBG.pdf")

