import ROOT as R
import sys
import os
ran = R.TRandom3()
sys.path.append(os.environ["NEUTRONS"]+'/RooFit')
pdfDir = "/home/ethan/Nextcloud/neutrons/"
#fPDF = pdf(pdfDir+"LT_calib.pdf")

binWidth  = 40
plotRange = [0,30000]
fitRange  = [1000,24000]

nEvents   = 1e10
doRand    = False

rBG       = 3.5e-5
lMu       = 0.00045
lWall     = 1.1e-2 #silver

nBins = int(round((plotRange[1]-plotRange[0])/binWidth))
binWidth = (plotRange[1]-plotRange[0])/nBins
hLT = R.TH1D("hLT","Lifetime Fit", nBins,plotRange[0],plotRange[1])
can = R.TCanvas()

def doFit(alpha = 0.0):
  hLT.SetTitle("Lifetime Fit (#alpha=%f)" % alpha)

  vD = [1,1]
  vW = [1,1]
  nD = 0
  nW = 0
  nBG = nEvents*rBG
  for b in range(1,nBins+1):
    t = binWidth * b
    if alpha < 1:
      vD[0] = vD[1]
      vD[1] = R.TMath.Exp(-lMu*t)
      nD = nEvents*(1-alpha)*(vD[0]-vD[1])
    if alpha > 0:
      vW[0] = vW[1]
      vW[1] = R.TMath.Exp(-(lMu+lWall)*t)
      nW = nEvents*(alpha)*(vW[0]-vW[1])
    nTot = nD + nW + nBG
    if doRand:
      hLT.SetBinContent(b, ran.PoissonD(nTot))
    else:
      hLT.SetBinContent(b, nTot)
  
  can.SetLogy()
  hLT.Draw()
  fExpo1 = R.TF1("expo1", "expo(0)+pol0(2)", fitRange[0], fitRange[1])
  fExpo1.SetParameters(R.TMath.Log(nEvents*(1-alpha)/binWidth), -lMu, nBG)
  res = hLT.Fit(fExpo1, "RSQ", "")
  can.Update()
  #fPDF.Print(can)
  shift = (-fExpo1.GetParameter(1)-lMu)*1e9
  shiftErr = fExpo1.GetParError(1)*1e9
  if doRand:
    print( "LT Shift = %f +/- %f" % ( shift, shiftErr ) )
  else:
    print( "LT Shift = %f" % shift )
  return [[fExpo1.GetParameter(0),fExpo1.GetParError(0)],
          [shift,shiftErr],
          [fExpo1.GetParameter(2),fExpo1.GetParError(2)],
          [res.Chi2()/res.Ndf()]]


def doAlphaScan(aMax = 0.5, aSteps = 50):
  print("Doing Scan")
  hists = []
  aWidth = aMax / aSteps
  hists.append(R.TH1D("hN","N;#alpha",aSteps,-0.5*aWidth,aMax+0.5*aWidth))
  hists.append(R.TH1D("hR","R;#alpha",aSteps,-0.5*aWidth,aMax+0.5*aWidth))
  hists.append(R.TH1D("hBG","BG;#alpha",aSteps,-0.5*aWidth,aMax+0.5*aWidth))
  hists.append(R.TH1D("hChi2","Chi2;#alpha",aSteps,-0.5*aWidth,aMax+0.5*aWidth))
  
  for a in range(aSteps+1):
    alpha = aMax*a/aSteps
    res = doFit(alpha)
    for i in range(3):
      hists[i].SetBinContent(a+1, res[i][0])
      if doRand:
        hists[i].SetBinError(a+1, res[i][1])
    hists[3].SetBinContent(a+1, res[3][0])
  
  can.SetLogy(False)
  hists[1].Draw("c")
  can.Update()
  return hists


def findAlpha(target = 1.0, err = 0.001):
  alphas = [0,0.01]
  shifts = [doFit(alphas[0])[1][0], doFit(alphas[1])[1][0]]
  shift0 = shifts[0]
  while abs(shifts[1]-target) > err:
    slope = (shifts[1]-shifts[0])/(alphas[1]-alphas[0])
    alphas[0] = alphas[1]
    shifts[0] = shifts[1]
    alphas[1] = alphas[0] + (target-shifts[0])/slope
    print("alpha = %e" % alphas[1])
    if alphas[1] < 0:
      alphas[1] = 0
      shifts[1] = shift0
      print( "LT Shift = %f" % shift0 )
    elif alphas[1] >= 1:
      alphas[1] = 0.5+0.5*alphas[0]
      shifts[1] = doFit(alphas[1])[1][0]
    else:
      shifts[1] = doFit(alphas[1])[1][0]
    #raw_input()
  return alphas[1]


def doStartScan(target = 1.0, err = 0.001):
  hAlpha = R.TH1D("hAlpha", "Silver Wall Stop Fraction Required for a %.1f s^{-1} Rate Shift;Start Time [ns]" % target, 22,140,1020)
  for b in range(1,22+1):
    fitRange[0] = hAlpha.GetBinCenter(b)
    print( "Start Time = %d ns" % fitRange[0])
    hAlpha.SetBinContent(b, findAlpha(target, err))
  hAlpha.Draw("c")
  can.Update()
  return hAlpha

hAlpha = doStartScan()
hAlpha.GetYaxis().SetRangeUser(1e-5,1e-1)
hAlpha.SetLineWidth(3)
hAlpha.Draw("c")
can.Update()
can.Print("../figures/lifetime_calib.pdf")
raw_input()

