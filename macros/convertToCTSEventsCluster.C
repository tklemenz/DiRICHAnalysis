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

static const std::string process = "[MACRO]";
static const Int_t timeWindow = 5; // time window for correlation plots

void convertToCTSEventsCluster(const char *inputFile, const char *outputFile, ULong_t procNr, const bool& debugClusterer)
{
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

  std::vector<std::vector<TH2D*>> coincidenceVecQMax{};
  std::vector<std::vector<TH2D*>> coincidenceVecQTot{};
  for(Int_t layer=0; layer<8; layer++) {
    coincidenceVecQMax.emplace_back(std::vector<TH2D*>{});
    coincidenceVecQTot.emplace_back(std::vector<TH2D*>{});
    for(Int_t coi=0; coi<8; coi++) {
      coincidenceVecQMax.back().emplace_back(new TH2D(Form("hCoiMaxToTL%i%i",layer+1,coi+1),Form("Cluster coincidence in L%i%i;max ToT L%i;max ToT L %i",layer+1, coi+1, layer+1, coi+1),500,0,50,500,0,50));
      coincidenceVecQTot.back().emplace_back(new TH2D(Form("hCoiTotToTL%i%i",layer+1,coi+1),Form("Cluster coincidence in L%i%i;tot ToT L%i;tot ToT L %i",layer+1, coi+1, layer+1, coi+1),500,0,50,500,0,50));
    }
  }

  TH1D* hNSignalsCluster  = new TH1D("hNSignalsCluster","n Signals in Cluster;n signals;counts",20,0,20);
  TH1D* hNClusters = new TH1D("hNClusters","n Clusters in Event;n Clusters;counts",20,0,20);
  TH1D* hQmax = new TH1D("hQmax","max Cluster ToT; max ToT [ns];counts",200,0,100);
  TH1D* hQtot = new TH1D("hQtot","total Cluster ToT; total ToT [ns]",200,0,100);

  Int_t fiberMult(0), layer(-1), x(-1), y(-1);

  TFile *fout = new TFile(Form("%s",outputFile),"recreate");
  TTree *treeout = new TTree("dummy","RadMap data in CTSEvents -> CTSEventClusters");

  CTSEvent *event = new CTSEvent();
  CTSEventClusters *ctsEventCluster = new CTSEventClusters();

  data->SetBranchAddress("Events", &event);
  treeout->Branch("CTSEventsCluster","ctsEventCluster",ctsEventCluster,32000,1);

  printf("events to process: %lu\t %.1f%% of the file\n", nEvents, Float_t(100*nEvents)/Float_t(data->GetEntries()));

  for (ULong_t entry = 0; entry < nEvents; entry++) {
    if ((((entry+1)%1000) == 0) || (entry == (nEvents-1))) {
      printf("\rprocessing event %lu...", entry+1);
      fflush(stdout);
      std::cout<<std::setw(5)<<std::setiosflags(std::ios::fixed)<<std::setprecision(1)<<" "<<(100.*(entry+1))/nEvents<<" % done\r"<<std::flush;
    }
    data->GetEntry(entry);
    eventNr       = event->getEventNr();
    padiwaConfig  = event->getPadiwaConfig();
    modules       = event->getModules();

    clusterer.findClusters(*event, ParticleType::Cosmic, debugClusterer);
    clusters = clusterer.getClusters();

    // check coincidences
    for (auto& cluster : clusters) {
      for (auto& other : clusters) {
        if (std::abs(cluster.getMeanTimeStamp()-other.getMeanTimeStamp() <= timeWindow)) {
          //printf("%s%s %sAccessing coincidence histos at Layer %d and Layer %d...%s\n", text::LBLU, process.c_str(), text::YEL, cluster.getLayer(), other.getLayer(), text::RESET);
          coincidenceVecQMax.at(cluster.getLayer()-1).at(other.getLayer()-1)->Fill(cluster.getQMax(), other.getQMax());
          coincidenceVecQTot.at(cluster.getLayer()-1).at(other.getLayer()-1)->Fill(cluster.getQTot(), other.getQTot());
        }
      }
    }

    for(auto& cluster : clusters){
      clusterQtotDists.at(cluster.getLayer()-1)->Fill(cluster.getMeanFiber(), cluster.getQTot());
      clusterQmaxDists.at(cluster.getLayer()-1)->Fill(Int_t(std::round(cluster.getMeanFiber())), cluster.getQMax());
      hNSignalsCluster->Fill(cluster.getNSignals());
      hQmax->Fill(cluster.getQMax());
      hQtot->Fill(cluster.getQTot());
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
  fout->WriteObject(hQmax,hQmax->GetName());
  fout->WriteObject(hQtot,hQtot->GetName());

  Int_t histCounter = 0;
  for(auto& hist : clusterQtotDists)  { if(hist->GetEntries() != 0) { fout->WriteObject(hist, hist->GetName()); histCounter++; } }
  for(auto& hist : clusterQmaxDists)  { if(hist->GetEntries() != 0) { fout->WriteObject(hist, hist->GetName()); } }

  for (auto& vec : coincidenceVecQMax) {
    for (auto& hist : vec) {
      if(hist->GetEntries() != 0) {
        fout->WriteObject(hist, hist->GetName());
      }
    }
  }

  for (auto& vec : coincidenceVecQTot) {
    for (auto& hist : vec) {
      if(hist->GetEntries() != 0) {
        fout->WriteObject(hist, hist->GetName());
      }
    }
  }

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