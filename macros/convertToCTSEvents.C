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

  TH1D* hNSignalsEvent  = new TH1D("hNSignalsEvent","n Signals in Event",200,0,200);
  TH1D* hToTAll         = new TH1D("hToTAll","Overall ToT distribution",200,0,200);

  TFile *fout = new TFile(Form("%s",outputFile),"recreate");
  TTree *tree = new TTree("dummy","RadMap data in fancy objects -> CTSEvents");

  Signal signal               = Signal();
  std::vector<Module> modules = std::vector<Module>{ Module(), Module()};
  CTSEvent *event = new CTSEvent();

  // for coincidence
  std::vector<Signal> signalsInEvent{};
  std::vector<std::vector<TH2D*>> coincidenceVec{};
  for(Int_t layer=0; layer<8; layer++) {
    coincidenceVec.emplace_back(std::vector<TH2D*>{});
    for(Int_t coi=0; coi<8; coi++) {
      coincidenceVec.back().emplace_back(new TH2D(Form("hCoiToTL%i%i",layer+1,coi+1),Form("Signal coincidence in L%i%i;ToT L%i;ToT L %i",layer+1, coi+1, layer+1, coi+1),500,0,50,500,0,50));
    }
  }

  // for multiplicity (how many layers were hit)
  std::vector<bool> layerHit{ false, false, false, false, false, false, false, false };
  Int_t hitLayerCounter = 0;
  TH1D* nHitLayers = new TH1D("hNHitLayers","n hit layers in Event",8,0,8);

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

      // check coincidences
      for (auto& signal : signalsInEvent) {
        for (auto& other : signalsInEvent) {
          if (std::abs(signal.getTimeStamp()-other.getTimeStamp() <= timeWindow)) {
            //printf("%s%s %sAccessing coincidence histos at Layer %d and Layer %d...%s\n", text::LBLU, process.c_str(), text::YEL, cluster.getLayer(), other.getLayer(), text::RESET);
            coincidenceVec.at(signal.getLayer()-1).at(other.getLayer()-1)->Fill(signal.getToT(), other.getToT());
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
  for (auto& vec : coincidenceVec) {
    for (auto& hist : vec) {
      if(hist->GetEntries() != 0) { fout->WriteObject(hist, hist->GetName()); }
    }
  }
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
