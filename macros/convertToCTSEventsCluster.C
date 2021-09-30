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

#include "CTSEventClusters.h"
#include "Clusterer.h"
#include "Utility.h"

///< usage: ./convertToCTSEventsCluster -i inputfile -o outputfile -n numberOfEventsToBeProcessed
///< n = -1 by default which means the whole file is processed

extern char* optarg;

void convertToCTSEventsCluster(const char *inputFile, const char *outputFile, ULong_t procNr, const bool& debugClusterer)
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

  std::vector<TH2D*> clusterQtotDists{};
  std::vector<TH2D*> clusterQmaxDists{};
  for(Int_t i=0; i<16; i++) {
    clusterQtotDists.emplace_back(new TH2D(Form("hToTtotalL%i",i+1),Form("ToT total distribution of clusters vs fiber in L%i;fiber;ToT",i+1),33,0,33,500,0,50));
    clusterQmaxDists.emplace_back(new TH2D(Form("hToTmaxL%i",i+1),Form("max ToT distribution of clusters vs fiber in L%i;fiber;ToT",i+1),33,0,33,500,0,50));
  }

  TH1D* hNSignalsCluster  = new TH1D("hNSignalsCluster","n Signals in Cluster",20,0,20);
  TH1D* hNClusters = new TH1D("hNClusters","n Clusters in Event",20,0,20);

  Int_t fiberMult(0), layer(-1), x(-1), y(-1), clusterCounter(0);

  TFile *fout = new TFile(Form("%s",outputFile),"recreate");
  TTree *treeout = new TTree("dummy","RadMap data in CTSEvents -> CTSEventClusters");

  CTSEvent *event;
  //CTSEvent eventBuffer;
  event = new CTSEvent();
  //eventBuffer = CTSEvent();
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
    //eventBuffer.setModules(modules);

    clusterer.findClusters(*event, ParticleType::Cosmic, debugClusterer);
    clusters = clusterer.getClusters();

    for(auto& cluster : clusters){
      clusterQtotDists.at(cluster.getLayer()-1)->Fill(Int_t(std::round(cluster.getMeanFiber())), cluster.getQTot());
      clusterQmaxDists.at(cluster.getLayer()-1)->Fill(Int_t(std::round(cluster.getMeanFiber())), cluster.getQMax());
      hNSignalsCluster->Fill(cluster.getNSignals());
    }
    hNClusters->Fill(clusters.size());

    ctsEventCluster->setClusters(clusters);
    ctsEventCluster->setEventNr(eventNr);
    ctsEventCluster->setPadiwaConfig(padiwaConfig);
    treeout->Fill();
    for (auto& module : modules) {
      module.reset();
    }
    clusterer.reset();
  } /// loop over file

  treeout->Write("data");
  fout->WriteObject(clusterer.hNSignalsEvent, clusterer.hNSignalsEvent->GetName());
  fout->WriteObject(clusterer.hNGoodSignalsEvent, clusterer.hNGoodSignalsEvent->GetName());
  fout->WriteObject(hNClusters, hNClusters->GetName());
  fout->WriteObject(hNSignalsCluster,hNSignalsCluster->GetName());

  Int_t histCounter = 0;
  for(auto& hist : clusterQtotDists)  { if(hist->GetEntries() != 0) { fout->WriteObject(hist, hist->GetName()); histCounter++; } }
  for(auto& hist : clusterQmaxDists)  { if(hist->GetEntries() != 0) { fout->WriteObject(hist, hist->GetName()); } }

  TCanvas *c1 = new TCanvas("ctotToTDists","ctotToTDists");
  c1->DivideSquare(histCounter);
  Int_t padIter = 1;
  for(auto& hist : clusterQtotDists) {
    if(hist->GetEntries() == 0) { continue; }
    c1->cd(padIter);
    gPad->SetLogz();
    hist->Draw("COLZ");
    padIter++;
  }

  TCanvas*c2 = new TCanvas("cmaxToTDists","cmaxToTDists");
  c2->DivideSquare(histCounter);
  padIter = 1;
  for(auto& hist : clusterQmaxDists) {
  if(hist->GetEntries() == 0) { continue; }
    c2->cd(padIter);
    gPad->SetLogz();
    hist->Draw("COLZ");
    padIter++; 
  }

  fout->WriteObject(c1, c1->GetName());
  fout->WriteObject(c2, c2->GetName());

  fout->Close();

  delete ctsEventCluster;
  ctsEventCluster=nullptr;
}

int main(int argc, char** argv)
{
  char    inputFile[512]="";
  char    outputFile[512]="convertToCTSEventsCluster_output.root";
  ULong_t procNr=-1;
  bool    debugClusterer = false;

  int argsforloop;
  while ((argsforloop = getopt(argc, argv, "hi:o:n:d:")) != -1) {
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
      case 'd':
        if (std::atoi(optarg) > 0) { debugClusterer = true; }
        else if (std::atoi(optarg) == 0) { debugClusterer = false; }
        else { printf("Valid debug options: '1', '0' (default)"); }
        break;
      default:
        printf("\n\n%s%sdefault case%s\n\n",text::BOLD,text::RED,text::RESET);
        exit(EXIT_FAILURE);
    }
  }

  printf("\n\n%sRunning convertToCTSEventsCluster%s\n\n",text::BOLD,text::RESET);
  
  convertToCTSEventsCluster(inputFile,outputFile,procNr,debugClusterer);

  printf("\n\n%s%sDONE!%s\n\n",text::BOLD,text::GRN,text::RESET);
}