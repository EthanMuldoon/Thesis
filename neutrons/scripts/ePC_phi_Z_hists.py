import ROOT
can = ROOT.TCanvas()

fAllHits = ROOT.TFile("/scratch/ethan/hist/stampede/sum.root")
hAllHits = getattr(fAllHits,"hepc1AnodeInnerPhiZ")
hAllHits.RebinX(4)
b1 = hAllHits.ProjectionY().FindFirstBinAbove(0)
b2 = hAllHits.ProjectionY().FindLastBinAbove(0)
hAllHits.GetYaxis().SetRange(b1,b2)
hAllHits.GetZaxis().SetRangeUser(0,1.6e6)
hAllHits.SetTitle("R2015 ePC #phi vs Z")
hAllHits.Draw("colz")
can.Print("../figures/ePC_phi_Z_All.pdf")

fNeuSel = ROOT.TFile("/scratch/ethan/mta/run9/pass4/mta_silver.root")
hNeuSel = getattr(fNeuSel, "hNeuEPC1PhiZ")
hNeuSel.RebinX(4)
hNeuSel.GetYaxis().SetRange(b1,b2)
hNeuSel.GetZaxis().SetRangeUser(0,400)
hNeuSel.SetTitle("R2016 ePC #phi vs Z with Coincident Neutron")
hNeuSel.Draw("colz")
can.Print("../figures/ePC_phi_Z_Neu.pdf")