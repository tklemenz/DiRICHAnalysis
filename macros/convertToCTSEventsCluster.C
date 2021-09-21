#include <TTree.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <Rtypes.h>
#include <TCanvas.h>

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <unistd.h>

#include "CTSEvent.h"
#include "CTSEventClusters.h"
#include "Utility.h"

///< usage: ./convertToCTSEventsCluster -i inputfile -o outputfile -n numberOfEventsToBeProcessed
///< n = -1 by default which means the whole file is processed

extern char* optarg;

void convertToCTSEventsCluster(const char *inputFile, const char *outputFile, ULong_t procNr)
{
  Int_t fullTrackCounter = 0;
  Int_t fullTracksTmp = 0;

  TFile* f = TFile::Open(inputFile);

  if (f->IsOpen()==kFALSE){
    printf("\n\n%s%sOpen input file failed%s\n\n",text::BOLD,text::RED,text::RESET);
    exit(1);
  }

  TTree *data = (TTree*)f->Get("data");

  ULong_t nEvents = procNr;
  if ((nEvents == -1) || (nEvents > data->GetEntries())) { nEvents = data->GetEntries(); }

  Float_t eventNr      = -1;
  Int_t   padiwaConfig = -1;
  std::vector<Module> modules{};
  Clusterer clusterer = Clusterer();

  std::vector<Cluster> clusters; 
  std::vector<Fiber> fibers;

  Int_t fiberMult(0), layer(-1), x(-1), y(-1), clusterCounter(0);

  TFile *fout = new TFile(Form("%s",outputFile),"recreate");
  TTree *treeout = new TTree("dummy","RadMap data in CTSEvents -> CTSEventClusters");

  CTSEvent *event;
  CTSEvent eventBuffer;
  event = new CTSEvent();
  eventBuffer = CTSEvent();
  CTSEventClusters *ctsEventCluster;
  ctsEventCluster = new CTSEventClusters();

  data->SetBranchAddress("Events", &event);
  treeout->Branch("CTSEventsCluster","ctsEventCluster",ctsEventCluster,32000,1);

  printf("events to process: %lu\t %.1f%% of the file\n", nEvents, Float_t(100*nEvents)/Float_t(data->GetEntries()));


  for (ULong_t entry = 0; entry < nEvents; entry++) {
    if ((((entry+1)%1000) == 0) || (entry == (nEvents-1))) {
      printf("\rprocessing event %lu...", entry+1);
      fflush(stdout);
      std::cout<<std::setw(5)<<std::setiosflags(std::ios::fixed)<<std::setprecision(1)<<" "<<(100.*(entry+1))/nEvents<<" % done\r"<<std::flush;
    }
    clusterCounter = 0;
    data->GetEntry(entry);
    eventNr       = event->getEventNr();
    padiwaConfig  = event->getPadiwaConfig();
    modules       = event->getModules();
    eventBuffer.setModules(modules);

    fullTracksTmp = clusterer.findClusters(eventBuffer, ParticleType::Cosmic);
    clusters = clusterer.getClusters();

    for(auto& cluster : clusters){
      if(cluster.getNSignals() > clusterCounter) { clusterCounter = cluster.getNSignals();}
      //ctsEventCluster->addCluster(cluster);
    }
    ctsEventCluster->addCluster(clusters);
    treeout->Fill();
    for (auto& module : modules) {
      module.reset();
    }
    clusterer.reset();
    fullTrackCounter+=fullTracksTmp;
    fullTracksTmp = 0;
    
  } /// loop over file

  treeout->Write("data");
  //fout->WriteObject(clusterer.MhTimeDiff,"MhTimeDiff");
  //fout->WriteObject(clusterer.MhSpaceDiff,"MhSpaceDiff");
  fout->WriteObject(clusterer.mHnSignal,clusterer.mHnSignal->GetName());

  fout->Close();

  delete ctsEventCluster;
  ctsEventCluster=nullptr;

  //printf("\n\n High multiplicity cluster counter:%i\n\n", clusterCounter);
  //printf("\n\n Full track counter:%i\n\n", fullTrackCounter);

}

int main(int argc, char** argv)
{
  char    inputFile[512]="";
  char    outputFile[512]="convertToCTSEventsCluster_output.root";
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

  printf("\n\n%sRunning convertToCTSEventsCluster%s\n\n",text::BOLD,text::RESET);
  
  convertToCTSEventsCluster(inputFile,outputFile,procNr);

  printf("\n\n%s%sDONE!%s\n\n",text::BOLD,text::GRN,text::RESET);
}