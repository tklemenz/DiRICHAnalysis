#include <TH2.h>
#include <TH1.h>
#include <TTree.h>
#include <TFile.h>
#include <TNtupleD.h>
#include <TCanvas.h>
#include <Rtypes.h>

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <unistd.h>

#include "CTSEvent.h"
#include "Utility.h"
#include "Constants.h"

///< usage: ./convertToCTSEvents -i inputfile -o outputfileName -n numberOfSignalsToBeProcessed
///< n = -1 by default which means the whole file is processed

extern char* optarg;

static const Int_t timeWindow = 5; // time window for coincidence
static const Int_t fiberWindow = 1; // fiber window for coincidence

void convertToCTSEvents(const char *inputFile, const char *outputFile, ULong_t procNr)
{
  TFile* f = TFile::Open(inputFile);

  if (f->IsOpen()==kFALSE){
    printf("\n\n%s%sOpen input file failed%s\n\n",text::BOLD,text::RED,text::RESET);
    exit(1);
  }

  TNtupleD *signals = (TNtupleD*)f->Get("Signals");

  Double_t eventNr(-1), chID(0), TDC(-1), layer(-1), x(-1), y(-1), signalNr(-1), timeStamp(-1), ToT(-1), padiwaConfig(-1), refTime(-1), prevEventNr(-1), calibTime(0);
  Int_t prevSigNr(0), prevCh(-1), firstCounter(0), secondCounter(0), signalCounter(0); 

  ULong_t nSignals = procNr;

  signals->SetBranchAddress("EventNr",      &eventNr);
  signals->SetBranchAddress("timeStamp",    &timeStamp);    // in seconds
  signals->SetBranchAddress("ToT",          &ToT);          // in seconds
  signals->SetBranchAddress("chID",         &chID);
  signals->SetBranchAddress("TDC",          &TDC);          // 0 = DiRICH 1200, 1 = DiRICH 1201 etc
  signals->SetBranchAddress("layer",        &layer);        // 1-8
  signals->SetBranchAddress("x",            &x);            // odd layers have  x != 0
  signals->SetBranchAddress("y",            &y);            // even layers have y != 0
  signals->SetBranchAddress("signalNr",     &signalNr);     // Nth Signal per channel and event
  signals->SetBranchAddress("padiwaConfig", &padiwaConfig);
  signals->SetBranchAddress("refTime",      &refTime);

  TH1D* hNSignalsEvent    = new TH1D("hNSignalsEvent","n Signals in Event;n signals; counts", 200,0,200);
  TH1D* hToTAll           = new TH1D("hToTAll","Overall ToT distribution;ToT;counts",200,0,200);
  TH1D* hToTOneSigInEvent = new TH1D("hToTOneSigInEvent","ToT distribution for single signal events;ToT;counts",200,0,200);

  TFile *fout = new TFile(Form("%s",outputFile),"recreate");
  TTree *tree = new TTree("dummy","RadMap data in fancy objects -> CTSEvents");

  Signal signal               = Signal();
  std::vector<Module> modules = std::vector<Module>{ Module(), Module()};
  CTSEvent *event = new CTSEvent();

  // for coincidence
  std::vector<Signal> signalsInEvent{};
  std::vector<std::vector<TH2D*>> coincidenceVec{};
  std::vector<std::vector<TH2D*>> coincidenceVecFiber{};
  for(Int_t layer=0; layer<8; layer++) {
    coincidenceVec.emplace_back(std::vector<TH2D*>{});
    coincidenceVecFiber.emplace_back(std::vector<TH2D*>{});
    for(Int_t coi=0; coi<8; coi++) {
      coincidenceVec.back().emplace_back(new TH2D(Form("hCoiToTL%i%i",layer+1,coi+1),Form("Signal coincidence in L%i_%i, within %d ns;ToT L%i;ToT L %i",layer+1, coi+1, timeWindow, layer+1, coi+1),200,0,50,200,0,50));
      coincidenceVecFiber.back().emplace_back(new TH2D(Form("hCoiFiberL%i%i",layer+1,coi+1),Form("Signal coincidence in L%i_%i, within %d ns;fiber L%i;fiber L %i",layer+1, coi+1, timeWindow, layer+1, coi+1),33,0,33,33,0,33));
    }
  }

  TH2D* hCoiToTL13fibRange = new TH2D("hCoiToTL13fibRange","Signal coincidence in L1_3 if signals in fiber range +-1;ToT L1;ToT L3",200,0,50,200,0,50);
  TH2D* hCoiToTL35fibRange = new TH2D("hCoiToTL35fibRange","Signal coincidence in L3_5 if signals in fiber range +-1;ToT L3;ToT L5",200,0,50,200,0,50);
  TH2D* hCoiToTL57fibRange = new TH2D("hCoiToTL57fibRange","Signal coincidence in L5_7 if signals in fiber range +-1;ToT L5;ToT L7",200,0,50,200,0,50);
  TH2D* hCoiToTL24fibRange = new TH2D("hCoiToTL24fibRange","Signal coincidence in L2_4 if signals in fiber range +-1;ToT L2;ToT L4",200,0,50,200,0,50);
  TH2D* hCoiToTL46fibRange = new TH2D("hCoiToTL46fibRange","Signal coincidence in L4_6 if signals in fiber range +-1;ToT L4;ToT L6",200,0,50,200,0,50);
  TH2D* hCoiToTL68fibRange = new TH2D("hCoiToTL68fibRange","Signal coincidence in L6_8 if signals in fiber range +-1;ToT L6;ToT L8",200,0,50,200,0,50);

  // for multiplicity (how many layers were hit)
  std::vector<bool> layerHit{ false, false, false, false, false, false, false, false };
  Int_t hitLayerCounter = 0;
  TH1D* nHitLayers = new TH1D("hNHitLayers","n hit layers in Event; layers with hits; counts",9,0,9);

  tree->Branch("Events","CTSEvent",&event,32000,1);

  if ((nSignals == -1) || (nSignals > signals->GetEntries())) { nSignals = signals->GetEntries(); }

  printf("signals to process: %lu\t %.1f%% of the file\n", nSignals, Float_t(100*nSignals)/Float_t(signals->GetEntries()));

  for (ULong_t entry = 0; entry < nSignals; entry++) {
    if ((((entry+1)%10000) == 0) || (entry == (nSignals-1))) {
      printf("\rprocessing signal %lu...", entry+1);
      fflush(stdout);
      std::cout<<std::setw(5)<<std::setiosflags(std::ios::fixed)<<std::setprecision(1)<<" "<<(100.*(entry+1))/nSignals<<" % done\r"<<std::flush;
    }

    signals->GetEntry(entry);

    calibTime = (timeStamp-refTime)*1e9 - constants::dirichTimeCorr.at(TDC).at(chID);
    
    if ((ULong_t(eventNr) != prevEventNr) && (prevEventNr !=1)) {
      for (auto& module : modules) {
        module.removeEmpty();
      }
      hNSignalsEvent->Fill(signalCounter);
      event->setNSignals(signalCounter);
      event->setModules(modules);
      tree->Fill();
      signalCounter = 0;
      for (auto& module : modules) {
        module.reset();
      }

      if (signalsInEvent.size() == 1) {
        hToTOneSigInEvent->Fill(signalsInEvent.at(0).getToT());
      }

      // check coincidences
      for (auto& signal : signalsInEvent) {
        Int_t sigLayer = signal.getLayer();
        Float_t sigFiberNr = mapping::getFiberNr(signal.getConfiguration(), signal.getChannelID(), signal.getTDCID());
        for (auto& other : signalsInEvent) {
          Float_t otherFiberNr = mapping::getFiberNr(other.getConfiguration(), other.getChannelID(), other.getTDCID());
          if (std::abs(signal.getTimeStamp()-other.getTimeStamp() <= timeWindow)) {
            //printf("%s%s %sAccessing coincidence histos at Layer %d and Layer %d...%s\n", text::LBLU, process.c_str(), text::YEL, cluster.getLayer(), other.getLayer(), text::RESET);
            coincidenceVec.at(signal.getLayer()-1).at(other.getLayer()-1)->Fill(signal.getToT(), other.getToT());
            coincidenceVecFiber.at(signal.getLayer()-1).at(other.getLayer()-1)->Fill(sigFiberNr, otherFiberNr);
            if (std::abs(sigFiberNr - otherFiberNr)>fiberWindow) { continue; }

            switch(sigLayer) {
              case 1:
                if (other.getLayer() == 3) { hCoiToTL13fibRange->Fill(signal.getToT(), other.getToT()); }
                break;
              case 2:
                if (other.getLayer() == 4) { hCoiToTL24fibRange->Fill(signal.getToT(), other.getToT()); }
                break;
              case 3:
                if (other.getLayer() == 5) { hCoiToTL35fibRange->Fill(signal.getToT(), other.getToT()); }
                break;
              case 4:
                if (other.getLayer() == 6) { hCoiToTL46fibRange->Fill(signal.getToT(), other.getToT()); }
                break;
              case 5:
                if (other.getLayer() == 7) { hCoiToTL57fibRange->Fill(signal.getToT(), other.getToT()); }
                break;
              case 6:
                if (other.getLayer() == 8) { hCoiToTL68fibRange->Fill(signal.getToT(), other.getToT()); }
                break;
              default:
                break;
            }
          }
        }
        // check n hit layers
        layerHit.at(signal.getLayer()-1) = true;
      }
      for (Int_t i=0; i<layerHit.size(); i++) {
        if (layerHit.at(i) == true) { hitLayerCounter++; }
      }
      nHitLayers->Fill(hitLayerCounter);

      hitLayerCounter = 0;
      for (Int_t i=0; i<layerHit.size(); i++) { layerHit.at(i) = false; }
      signalsInEvent.clear();
    }

    if (TDC == 5 || TDC == 8 || TDC == 10) {
      ToT = ToT*1e9 -10;
    }
    else {
      ToT = ToT*1e9;
    }

    signal = Signal(ToT,calibTime,signalNr,chID,layer,TDC,padiwaConfig);

    modules.at(mapping::getModule(padiwaConfig, TDC)).addSignal(signal);

    signalsInEvent.emplace_back(signal);

    event->setEventNr(eventNr);
    event->setPadiwaConfig(UShort_t(padiwaConfig));

    prevEventNr = eventNr;

    if (signalNr == 1) { signalCounter++; }
  }

  tree->Write("data");
  fout->WriteObject(hNSignalsEvent, hNSignalsEvent->GetName());
  fout->WriteObject(nHitLayers, nHitLayers->GetName());
  fout->WriteObject(hToTOneSigInEvent, hToTOneSigInEvent->GetName());

  for (auto& vec : coincidenceVec) {
    for (auto& hist : vec) {
      if(hist->GetEntries() != 0) { fout->WriteObject(hist, hist->GetName()); }
    }
  }
  for (auto& vec : coincidenceVecFiber) {
    for (auto& hist : vec) {
      if(hist->GetEntries() != 0) { fout->WriteObject(hist, hist->GetName()); }
    }
  }
  fout->WriteObject(hCoiToTL13fibRange, hCoiToTL13fibRange->GetName());
  fout->WriteObject(hCoiToTL24fibRange, hCoiToTL24fibRange->GetName());
  fout->WriteObject(hCoiToTL35fibRange, hCoiToTL35fibRange->GetName());
  fout->WriteObject(hCoiToTL46fibRange, hCoiToTL46fibRange->GetName());
  fout->WriteObject(hCoiToTL57fibRange, hCoiToTL57fibRange->GetName());
  fout->WriteObject(hCoiToTL68fibRange, hCoiToTL68fibRange->GetName());
  fout->Close();

  delete event;
  event=nullptr;
}

int main(int argc, char** argv)
{
  char    inputFile[512]="";
  char    outputFile[512]="convertToCTSEvents_output.root";
  ULong_t procNr=-1;

  int argsforloop;
  while ((argsforloop = getopt(argc, argv, "hi:o:n:")) != -1) {
    switch (argsforloop) {
      case '?':
        ///TODO: write usage function
        exit(EXIT_FAILURE);
      case 'i':
        strncpy(inputFile, optarg, 512);
        break;
      case 'o':
        strncpy(outputFile, optarg, 512);
        break;
      case 'n':
        procNr = std::atoi(optarg);
        break;
      default:
        printf("\n\n%s%sdefault case%s\n\n",text::BOLD,text::RED,text::RESET);
        exit(EXIT_FAILURE);
    }
  }

  printf("\n\n%sRunning convertToCTSEvents%s\n\n",text::BOLD,text::RESET);
  
  convertToCTSEvents(inputFile,outputFile,procNr);

  printf("\n\n%s%sDONE!%s\n\n",text::BOLD,text::GRN,text::RESET);
}
