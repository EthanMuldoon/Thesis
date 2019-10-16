#include "TH1.h"

#define nBuffer 1000

TRandom3 ran;

double tKickStart = 800;
double tKickStop = 25000;
double beamRate = 0.0000415;
double decayRate = 0.00045;
//double tWindow = 35000;
double tWindow = 1000000;
double tWindowExt = tWindow + tKickStop;
double tPileup = 25000;
double tBin = 10;
double tClock;
double tBeam = 1.0/beamRate;
double fExtinction = 100;

double tRoll = 1e9;

vector<double> vElDecay;

TH1I* hMuReal;
TH1I* hMuRealKick;
TH1I* hMuRealKickPP;
TH1I* hMuClock;
TH1I* hMuClockKick;
TH1I* hMuClockKickPP;

TH1I* hKickReal;
TH1I* hKickRealKick;
TH1I* hKickRealKickPP;
TH1I* hKickClock;
TH1I* hKickClockKick;
TH1I* hKickClockKickPP;

TH1D* hEntranceStatisticsEM;

struct muon{
  double time;
  bool isClock;
  bool isKick;
  bool isPileup;
  bool isDone;
};

muon muons[nBuffer] = {0};
int idx = 0;

void printBuffer() {
  cout << "idx = " << idx << endl;
  for(int i=0; i<nBuffer; i++){
    cout << " " << i << ": " << muons[i].time << " " << muons[i].isClock << " " << muons[i].isKick << endl;
  }
}

void fillMuHists(muon mu1, muon mu2) {
  if(mu2.isClock) return;
  double dt = mu2.time - mu1.time;
  if(mu1.isClock) {
    hMuClock->Fill(dt);
    if(mu1.isKick) {
      hMuClockKick->Fill(dt);
      if(!mu1.isPileup) {
        hMuClockKickPP->Fill(dt);
      }
    }
  } else {
    hMuReal->Fill(dt);
    if(mu1.isKick) {
      hMuRealKick->Fill(dt);
      if(!mu1.isPileup) {
        hMuRealKickPP->Fill(dt);
      }
    }
  }
}

void fillKickHists(muon mu1, double dt) {
  if(mu1.isClock) {
    hEntranceStatisticsEM->Fill(3);
    hKickClock->Fill(dt);
    if(mu1.isKick) {
      hEntranceStatisticsEM->Fill(4);
      hKickClockKick->Fill(dt);
      if(!mu1.isPileup) {
        hEntranceStatisticsEM->Fill(5);
        hKickClockKickPP->Fill(dt);
      }
    }
  } else {
    hEntranceStatisticsEM->Fill(0);
    hKickReal->Fill(dt);
    if(mu1.isKick) {
      hEntranceStatisticsEM->Fill(1);
      hKickRealKick->Fill(dt);
      if(!mu1.isPileup) {
        hEntranceStatisticsEM->Fill(2);
        hKickRealKickPP->Fill(dt);
      }
    }
  }
}

void fillHists(int j) {
  if(muons[j].time < tWindowExt) return;
  bool foundKick = false;
  for(int k = j + nBuffer - 1; k>j; k--) {
    double dt = muons[j].time - muons[k%nBuffer].time;
    if(dt > tWindowExt) break;
    if(dt <= 0) {
      cout << "Warning: Buffer Overflow 1! (j=" << j << ",k=" << k << ")" << endl;
      cout << " dt = " << dt << endl;
      printBuffer();
      break;
    }
    fillMuHists(muons[j], muons[k%nBuffer]);
    if(muons[k%nBuffer].isKick) {
      fillKickHists(muons[j], -dt);
      /*if(!foundKick) {
        fillAdjKickHists(muons[j], -dt);
        foundKick = true;
      }*/
    }
  }
  for(int k = j + 1; k<j+nBuffer; k++) {
    double dt = muons[j].time - muons[k%nBuffer].time;
    if(-dt > tWindow) break;
    if(dt >= 0) {
      cout << "Warning: Buffer Overflow 2! (j=" << j << ",k=" << k << ")" << endl;
      cout << " dt = " << dt << endl;
      printBuffer();
      break;
    }
    fillMuHists(muons[j], muons[k%nBuffer]);
    if(muons[k%nBuffer].isKick) {
      fillKickHists(muons[j], -dt);
    }
  }
}

void push(muon mu) {
  /*if(mu.time > tRoll1) {
    mu.time -= tRoll2;
    for(int i=0; i<nBuffer; i++) {
      muons[i].time -= tRoll2;
    }
  }*/
  muons[idx] = mu;
  int j;
  for(j=idx + nBuffer - 1; j>=idx; j--){
    if(mu.time - muons[j%nBuffer].time < tPileup) {
      //mu.isPileup = true;
      muons[idx].isPileup = true;
      muons[j%nBuffer].isPileup = true;
    }
    if(mu.time - muons[j%nBuffer].time > tWindow) {
      if(muons[j%nBuffer].isDone) {
        break;
      } else {
        muons[j%nBuffer].isDone = true;
        fillHists(j%nBuffer);
      }
    }
    if(muons[j%nBuffer].time == 0) break;
  }
  if(j < idx) {
    cout << "Warning: Buffer Overflow 3!" << endl;
    printBuffer();
  }
  idx = (idx+1)%nBuffer;
}



double nextClock(double time) {
  return ( 1 + (int)(time / tClock) ) * tClock;
}

double nextMu(double time, bool kick=false) {
  if(kick) {
    return time + fExtinction*ran.Exp(tBeam);
  }
  return time + ran.Exp(tBeam);
}

void event(double& tEnt, bool& clock, double& tNextMu, double& tNextClock) {
  while(tNextMu < tEnt + tKickStart) {
    if(tNextClock < tNextMu) {
      //fill(tNextClock, true, false);
      push({tNextClock, true, false, false, false});
      tNextClock += tClock;
    }
    //fill(tNextMu, false, false);
    push({tNextMu, false, false, false, false});
    tNextMu = nextMu(tNextMu);
  }
  
  tNextMu = (tNextMu - (tEnt + tKickStart))*fExtinction + (tEnt + tKickStart);
  
  while(tNextMu < tEnt + tKickStop) {
    if(tNextClock < tNextMu) {
      //fill(tNextClock, true, false);
      push({tNextClock, true, false, false, false});
      tNextClock += tClock;
    }
    //fill(tNextMu, false, false);
    push({tNextMu, false, false, false, false});
    tNextMu = nextMu(tNextMu, true);
  }
  
  tNextMu = (tNextMu - (tEnt + tKickStop))/fExtinction + (tEnt + tKickStop);
  
  //tNextMu = nextMu(tEnt + tKickStop);
  
  if(tNextClock < tEnt + tKickStop) {
    //fill(tNextClock, true, false);
    push({tNextClock, true, false, false, false});
    tNextClock += tClock;
  }
  
  if(tNextMu < tNextClock) {
    clock = false;
    tEnt = tNextMu;
    tNextMu = nextMu(tEnt);
  } else {
    clock = true;
    tEnt = tNextClock;
    tNextClock += tClock;
  }
  //fill(tEnt, clock, true);
  push({tEnt, clock, true, false, false});
}

void backgroundMC(double clockRate, long count, TString outFile="background.root") {
  TFile* fOut = new TFile(outFile.Data(), "RECREATE");
  hMuReal = new TH1I("hMuReal", "Real Muon Entrance Autocorrelation", 2*tWindow/tBin,-tWindow,tWindow);
  hMuRealKick = new TH1I("hMuRealKick", "Real Muon Entrance Autocorrelation with Kick", 2*tWindow/tBin,-tWindow,tWindow);
  hMuRealKickPP = new TH1I("hMuRealKickPP", "Real Muon Entrance Autocorrelation with Kick and PP", 2*tWindow/tBin,-tWindow,tWindow);
  hMuClock = new TH1I("hMuClock", "Clock Muon Entrance Autocorrelation", 2*tWindow/tBin,-tWindow,tWindow);
  hMuClockKick = new TH1I("hMuClockKick", "Clock Muon Entrance Autocorrelation with Kick", 2*tWindow/tBin,-tWindow,tWindow);
  hMuClockKickPP = new TH1I("hMuClockKickPP", "Clock Muon Entrance Autocorrelation with Kick and PP", 2*tWindow/tBin,-tWindow,tWindow);

  hKickReal = new TH1I("hKickReal", "Kicks vs Real Muon", 2*tWindow/tBin,-tWindow,tWindow);
  hKickRealKick = new TH1I("hKickRealKick", "Kicks vs Real Muon with Kick", 2*tWindow/tBin,-tWindow,tWindow);
  hKickRealKickPP = new TH1I("hKickRealKickPP", "Kicks vs Real Muon with Kick and PP", 2*tWindow/tBin,-tWindow,tWindow);
  hKickClock = new TH1I("hKickClock", "Kicks vs Clock Muon", 2*tWindow/tBin,-tWindow,tWindow);
  hKickClockKick = new TH1I("hKickClockKick", "Kicks vs Clock Muon with Kick", 2*tWindow/tBin,-tWindow,tWindow);
  hKickClockKickPP = new TH1I("hKickClockKickPP", "Kicks vs Clock Muon with Kick and PP", 2*tWindow/tBin,-tWindow,tWindow);
  
  hEntranceStatisticsEM = new TH1D("hEntranceStatisticsEM", "Entrance Statistics", 6,0,6);
  hEntranceStatisticsEM->GetXaxis()->SetBinLabel(1, "Real");
  hEntranceStatisticsEM->GetXaxis()->SetBinLabel(2, "Real Kick");
  hEntranceStatisticsEM->GetXaxis()->SetBinLabel(3, "Real Kick PP");
  hEntranceStatisticsEM->GetXaxis()->SetBinLabel(4, "Clk");
  hEntranceStatisticsEM->GetXaxis()->SetBinLabel(5, "Clk Kick");
  hEntranceStatisticsEM->GetXaxis()->SetBinLabel(6, "Clk Kick PP");
  
  tClock = 1000000/clockRate;
  double tEnt = 0;
  //fill(tEnt, true, true);
  push({tEnt, true, true, false, false});
  double tNextMu = tEnt + tKickStop + ran.Exp(tBeam);
  double tNextClock = tEnt + tClock;
  bool clock;
  
  for(long i=1; i<count; i++) {
    if(i % 1000000 == 0) {
      cout << "Processed " << i/1000000 << " M events" << endl;
    }
    event(tEnt, clock, tNextMu, tNextClock);
    if(tEnt > tRoll) {
      double dt = 10*tWindowExt - tRoll;
      tEnt += dt;
      tNextMu += dt;
      tNextClock += dt;
      for(int j=0; j<nBuffer; j++) {
        muons[j].time  += dt;
      }
    }
  }
  
  fOut->Write();
  fOut->Close();
  
}
