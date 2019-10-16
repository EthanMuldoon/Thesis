import ROOT as R
can = R.TCanvas()

sim1k = R.TFile("long1k.root")
hKickRealKick = sim1k.Get("hKickRealKick")
hKickRealKick.Rebin(10)
hKickAutoCorr = R.TH1D("hKickAutoCorr", "Kicker Autocorrelation (Simulation)", hKickRealKick.GetNbinsX(), hKickRealKick.GetXaxis().GetXmin()*0.001, hKickRealKick.GetXaxis().GetXmax()*0.001)
for b in range(1,hKickRealKick.GetNbinsX()+1):
  hKickAutoCorr.SetBinContent(b, hKickRealKick.GetBinContent(b))

hKickAutoCorr.SetLineColor(R.kBlack)
hKickAutoCorr.SetLineWidth(3)
hKickAutoCorr.GetXaxis().SetRangeUser(0,150)
hKickAutoCorr.GetXaxis().SetTitle("time [#mus]")
hKickAutoCorr.GetXaxis().SetNdivisions(6,False)
hKickAutoCorr.SetStats(False)

hEntranceStatisticsEM = sim1k.Get("hEntranceStatisticsEM")
counts = hEntranceStatisticsEM.GetBinContent(2) * 0.1 # bin width
hKickAutoCorr.Scale(1.0/counts)
hKickAutoCorr.GetYaxis().SetTitle("rate [#mus^{-1}]")

hKickAutoCorr.Draw("hist")
can.Print("../figures/KickAutoCorrSim.pdf")
