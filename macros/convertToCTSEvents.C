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

///< usage: ./convertToCTSEvents -i inputfile -o outputfileName -n numberOfSignalsToBeProcessed
///< n = -1 by default which means the whole file is processed

extern char* optarg;

void convertToCTSEvents(const char *inputFile, const char *outputFile, ULong_t procNr)
{
  TFile* f = TFile::Open(inputFile);

  if (f->IsOpen()==kFALSE){
    printf("\n\n%s%sOpen input file failed%s\n\n",text::BOLD,text::RED,text::RESET);
    exit(1);
  }

  TNtupleD *signals = (TNtupleD*)f->Get("Signals");

  Double_t eventNr(-1), chID(0), TDC(-1), layer(-1), x(-1), y(-1), signalNr(-1), timeStamp(-1), ToT(-1), padiwaConfig(-1), refTime(-1), prevEventNr(-1);
  Int_t prevSigNr(0), prevCh(-1), firstCounter(0), secondCounter(0), signalCounter(0);

  bool module1(false), module2(false);

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
  std::vector<Module> modules = std::vector<Module>{};
  CTSEvent *event = new CTSEvent();

  tree->Branch("Events","CTSEvent",&event,32000,1);

  if ((nSignals == -1) || (nSignals > signals->GetEntries())) { nSignals = signals->GetEntries(); }

  printf("signals to process: %lu\t %.1f%% of the file\n", nSignals, Float_t(100*nSignals)/Float_t(signals->GetEntries()));

  for (ULong_t entry = 0; entry < nSignals; entry++) {
    if ((((entry+1)%100000) == 0) || (entry == (nSignals-1))) {
      printf("\rprocessing signal %lu...", entry+1);
      fflush(stdout);
      std::cout<<std::setw(5)<<std::setiosflags(std::ios::fixed)<<std::setprecision(1)<<" "<<(100.*(entry+1))/nSignals<<" % done\r"<<std::flush;
    }

    signals->GetEntry(entry);

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
    }

    if (TDC == 5 || TDC == 8 || TDC == 10) {
      ToT = ToT*1e9 -10;
    }
    else {
      ToT = ToT*1e9;
    }

    signal = Signal(ToT,(timeStamp-refTime)*1e9,signalNr,chID,layer,TDC,padiwaConfig);

    //printf("module Nr: %d\n", mapping::getModule(padiwaConfig, TDC)-1);
    if (mapping::getModule(padiwaConfig, TDC) == 0 && !module1) {
      modules.emplace_back(Module());
      module1 = true;
    }
    if (mapping::getModule(padiwaConfig, TDC) == 1 && !module2) {
      modules.emplace_back(Module());
      module2 = true;
    }
    modules.at(mapping::getModule(padiwaConfig, TDC)).addSignal(signal);

    event->setEventNr(eventNr);
    event->setPadiwaConfig(UShort_t(padiwaConfig));

    prevEventNr = eventNr;

    if (signalNr == 1) { signalCounter++; }
  }

  tree->Write("data");
  fout->WriteObject(hNSignalsEvent, hNSignalsEvent->GetName());
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
