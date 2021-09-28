#include <TTree.h>
#include <TFile.h>
//#include <TH1.h>
//#include <TH2.h>
#include <Rtypes.h>
//#include <TCanvas.h>

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <unistd.h>

#include "CTSEventClusters.h"
#include "CTSEventTracks.h"
#include "Tracker.h"
#include "Utility.h"

///< usage: ./convertToCTSEventsTrack -i inputfile -o outputfile -n numberOfEventsToBeProcessed
///< n = -1 by default which means the whole file is processed

extern char* optarg;

void convertToCTSEventsTrack(const char *inputFile, const char *outputFile, ULong_t procNr, const bool& debugTracker)
{
  TFile* f = TFile::Open(inputFile);

  if (f->IsOpen()==kFALSE){
    printf("\n\n%s%sOpen input file failed%s\n\n",text::BOLD,text::RED,text::RESET);
    exit(1);
  }

  TTree *data = (TTree*)f->Get("data");

  ULong_t nEvents = procNr;
  if ((nEvents == -1) || (nEvents > data->GetEntries())) { nEvents = data->GetEntries(); }

  Tracker tracker = Tracker();

  Float_t eventNr = -1;
  Int_t padiwaConfig = -1;
  std::vector<Track> tracks; 

  TFile *fout = new TFile(Form("%s",outputFile),"recreate");
  TTree *treeout = new TTree("dummy","RadMap data in CTSEventClusters -> CTSEventTracks");

  CTSEventClusters *event;
  CTSEventClusters eventBuffer;
  event = new CTSEventClusters();
  eventBuffer = CTSEventClusters();
  CTSEventTracks *ctsEventTracks;
  ctsEventTracks = new CTSEventTracks();

  data->SetBranchAddress("CTSEventsCluster", &event);
  treeout->Branch("CTSEventsTracks","ctsEventTracks",ctsEventTracks,32000,1);

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
    printf("%s[MACRO] %s(Event %g)%s nClusters: %lu\n", text::LBLU, text::LYEL, eventNr, text::RESET, event->getClusters().size());
    eventBuffer.setClusters(event->getClusters());

    tracker.run(eventBuffer, ParticleType::Cosmic, debugTracker);

    tracks = tracker.getTracks();

    for(auto& track : tracks){
      printf("%s[MACRO]%s Number of clusters in track: %lu\n", text::LBLU, text::RESET, track.getClusters().size());
    }

    ctsEventTracks->setTracks(tracks);
    ctsEventTracks->setEventNr(eventNr);
    ctsEventTracks->setPadiwaConfig(padiwaConfig);

    treeout->Fill();

    tracker.reset();
    
  } /// loop over file

  treeout->Write("data");

  fout->Close();

  delete ctsEventTracks;
  ctsEventTracks=nullptr;

}

int main(int argc, char** argv)
{
  char    inputFile[512]="";
  char    outputFile[512]="convertToCTSEventsTrack_output.root";
  ULong_t procNr=-1;
  bool    debugTracker = false;

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
        if (std::atoi(optarg) > 0) { debugTracker = true; }
        else { printf("Valid debug options: '1', '0' (default)"); }
        break;
      default:
        printf("\n\n%s%sdefault case%s\n\n",text::BOLD,text::RED,text::RESET);
        exit(EXIT_FAILURE);
    }
  }

  printf("\n\n%sRunning convertToCTSEventsTrack%s\n\n",text::BOLD,text::RESET);
  
  convertToCTSEventsTrack(inputFile,outputFile,procNr,debugTracker);

  printf("\n\n%s%sDONE!%s\n\n",text::BOLD,text::GRN,text::RESET);
}