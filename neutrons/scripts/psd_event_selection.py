import ROOT
can = ROOT.TCanvas()
can.SetLogz()

channel = 0
xMin = 8
xMax = 12.2
yMin = -0.2
yMax = 0.5

def getPSD(fIn, hName, channel, xMin,xMax, yMin,yMax):
  h3d = getattr(fIn, hName)
  h3d.GetZaxis().SetRange(channel+1)
  h2d = h3d.Project3D("yx")
  h2d.GetXaxis().SetRangeUser(xMin,xMax)
  h2d.GetYaxis().SetRangeUser(yMin,yMax)
  h2d.GetYaxis().SetTitle("Tail Integral / Total Integral")
  return h2d

fIn = ROOT.TFile("/scratch/ethan/mta/pass1806/psds/psds_dsn1a.root")
hGamma = getPSD(fIn,"hNeuPSDRawLog_PromptMu", channel, xMin,xMax, yMin,yMax)
hNeu = getPSD(fIn,"hNeuPSDRawLog_WallStop", channel, xMin,xMax, yMin,yMax)
hFusion = getPSD(fIn,"hNeuPSDRawLog_Fusion", channel, xMin,xMax, yMin,yMax)
hNeu.Add(hFusion)

hGamma.SetTitle("Gamma-Selected PSD (Atomic Capture)")
hGamma.Draw("colz")
can.Print("../figures/PSD_gammma_atomic.pdf")

hNeu.SetTitle("Neutron-Selected PSD (Fusion & Wall Capture)")
hNeu.Draw("colz")
can.Print("../figures/PSD_neu_fusion_wall.pdf")

