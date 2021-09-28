#include <TTree.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
#include <Rtypes.h>
#include <TCanvas.h>
#include <TChain.h>

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <unistd.h>

#include "CTSEvent.h"
#include "Utility.h"

//#include <fmt/format.h>
//#include <boost/log/trivial.hpp>

///< usage: ./fitCosmics -i inputfile -o outputfile -n numberOfEventsToBeProcessed
///< n = -1 by default which means the whole file is processed

extern char* optarg;


  //void histoSettings (TH1D hist) {
	  //hist->SetFillColor(kAzure-9);
	  //hist->GetXaxis()->SetLabelSize(0.03);
	  //hist->GetYaxis()->SetLabelSize(0.03);
	  //hist->GetXaxis()->SetTitleSize(0.04);
	  //hist->GetYaxis()->SetTitleSize(0.04);	  
  //}
void fitCosmicss(const TString inputFiles, const char *outputFile, ULong_t procNr)
{

  TChain chain("data", "data");
  fileHandling::makeChain(chain, inputFiles);
  TObjArray* files = chain.GetListOfFiles();
  printf("%sFiles to be processed:%s\n", text::BOLD, text::RESET);
  for (int ifile=0; ifile<files->GetEntriesFast(); ++ifile){ printf("%s\n", files->At(ifile)->GetTitle()); }

  /* Define variables
  ==========================================================
  ==========================================================*/

  ULong_t nEvents    = -1;
  Int_t   layer      = -1;
  Int_t   x          = -1;
  Int_t   y          = -1;
  std::string outputName = std::string();

  /* Define histograms and other useful containers
  ==========================================================
  ==========================================================*/
  CTSEvent *event = nullptr;
  Float_t eventNr      = -1;
  Int_t   padiwaConfig = -1;
  std::vector<Module> modules{};
  std::vector<Fiber> fibers;
  
  std::vector<TH1D*> totFiberLy1Vec{}; 
  std::vector<TH1D*> totFiberLy2Vec{}; 
  std::vector<TH1D*> totFiberLy3Vec{}; 
  std::vector<TH1D*> totFiberLy4Vec{}; 
  std::vector<TH1D*> totFiberLy5Vec{}; 
  std::vector<TH1D*> totFiberLy6Vec{}; 
  std::vector<TH1D*> totFiberLy7Vec{}; 
  std::vector<TH1D*> totFiberLy8Vec{}; 
  
  
  
  //std::vector<std::vector<TH1D>> totFiberVec{};
 
 
  for(Int_t i=0; i<32; i++) {
	totFiberLy1Vec.emplace_back(new TH1D(Form("hToTLy1F%i", i+1),Form("ToT distribustions of first Signals in Layer 1 Fiber %i;ToT;Counts", i+1),400,0,40)); 
	totFiberLy2Vec.emplace_back(new TH1D(Form("hToTLy2F%i", i+1),Form("ToT distribustions of first Signals in Layer 2 Fiber %i;ToT;Counts", i+1),400,0,40)); 
    totFiberLy3Vec.emplace_back(new TH1D(Form("hToTLy3F%i", i+1),Form("ToT distribustions of first Signals in Layer 3 Fiber %i;ToT;Counts", i+1),400,0,40)); 
	totFiberLy4Vec.emplace_back(new TH1D(Form("hToTLy4F%i", i+1),Form("ToT distribustions of first Signals in Layer 4 Fiber %i;ToT;Counts", i+1),400,0,40)); 
	totFiberLy5Vec.emplace_back(new TH1D(Form("hToTLy5F%i", i+1),Form("ToT distribustions of first Signals in Layer 5 Fiber %i;ToT;Counts", i+1),400,0,40)); 
	totFiberLy6Vec.emplace_back(new TH1D(Form("hToTLy6F%i", i+1),Form("ToT distribustions of first Signals in Layer 6 Fiber %i;ToT;Counts", i+1),400,0,40)); 
	totFiberLy7Vec.emplace_back(new TH1D(Form("hToTLy7F%i", i+1),Form("ToT distribustions of first Signals in Layer 7 Fiber %i;ToT;Counts", i+1),400,0,40)); 
	totFiberLy8Vec.emplace_back(new TH1D(Form("hToTLy8F%i", i+1),Form("ToT distribustions of first Signals in Layer 8 Fiber %i;ToT;Counts", i+1),400,0,40)); 
	}
	//totFiberVec.insert(std::end(totFiberVec), std::begin(totFiberLy1Vec), std::end(totFiberLy1Vec));  
   
  //for(Int_t i = 0; i<8; i++) {
	  //totFiberVec.insert(std::end(totFiberVec), std::begin(Form(totFiberLy"%i"Vec,i+1), std::end(Form(totFiberLy"%i",i+1));  
  //}
  
  /*========================================================
  ==========================================================*/

  if (fileHandling::splitString(inputFiles.Data(), ",").size() == 1) {
    outputName = fileHandling::splitString(fileHandling::splitString(outputFile).back().data(), ".").front().data();
    outputName.append("_");
    outputName.append(fileHandling::splitString(fileHandling::splitString(inputFiles.Data()).back().data(), ".").front().data());
    outputName.append(".root");
  }
  else { 
    outputName = fileHandling::splitString(fileHandling::splitString(outputFile).back().data(), ".").front().data();
    outputName.append(".root");
  }

  TFile *fout = new TFile(Form("%s", outputName.c_str()), "recreate");

  for (Int_t ifile=0; ifile<files->GetEntriesFast(); ++ifile){
    TFile *EventFile = new TFile(Form("%s", files->At(ifile)->GetTitle()));
    TTree *data = (TTree*)EventFile->Get("data");
    data->SetBranchAddress("Events", &event);
    nEvents = procNr;
    if ((nEvents == -1) || (nEvents > data->GetEntries())) { nEvents = data->GetEntries(); }
    printf("events to process: %lu\t %.1f%% of the file\n", nEvents, Float_t(100*nEvents)/Float_t(data->GetEntries()));

    for (ULong_t entry = 0; entry < nEvents; entry++) {
      if ((((entry+1)%100000) == 0) || (entry == (nEvents-1))) {
        printf("\rprocessing event %lu...", entry+1);
        fflush(stdout);
        std::cout<<std::setw(5)<<std::setiosflags(std::ios::fixed)<<std::setprecision(1)<<" "<<(100.*(entry+1))/nEvents<<" % done\r"<<std::flush;
      }

      data->GetEntry(entry);

      for(Int_t module=0; module < event->getModules().size(); module++) {
        modules.emplace_back(event->getModules().at(module));
      };

      eventNr       = event->getEventNr();
      padiwaConfig  = event->getPadiwaConfig();
      modules       = event->getModules();

      for (auto& module : modules) {
        for (auto& fiber : module.getFibers()) {
          if(fiber.getSignals().size() > 0) {
            layer = fiber.getLayer();
          }
          for(auto& signal : fiber.getSignals()) {
            if(signal.getSignalNr() == 1) {
              Float_t fiberNr = mapping::getFiberNr(signal.getConfiguration(),signal.getChannelID(),signal.getTDCID());
              if(layer == 1) {//std::cout<<layer<<std::endl;
				  totFiberLy1Vec.at(fiberNr-1)->Fill(signal.getToT());
			  }
			  else if (layer == 2) {
				  totFiberLy2Vec.at(fiberNr-1)->Fill(signal.getToT());
			  }
			  else if (layer == 3) {
				  totFiberLy2Vec.at(fiberNr-1)->Fill(signal.getToT());
			  }
			  else if (layer == 4) {
				  totFiberLy2Vec.at(fiberNr-1)->Fill(signal.getToT());
			  }
			  else if (layer == 5) {
				  totFiberLy2Vec.at(fiberNr-1)->Fill(signal.getToT());
			  }
			  else if (layer == 6) {
				  totFiberLy2Vec.at(fiberNr-1)->Fill(signal.getToT());
			  }
			  else if (layer == 7) {
				  totFiberLy2Vec.at(fiberNr-1)->Fill(signal.getToT());
			  }
			  else if (layer == 8) {
				  totFiberLy2Vec.at(fiberNr-1)->Fill(signal.getToT());
			  }
			  else {printf("Wrong Layer Nr \n");}
			  
            }
          } // loop over signals in fiber
        } // loop over fibers in module
      } // loop over modules
      modules.clear();
    } // loop over file
  } // loop over all files
  
  
  TF1 *fit  = new TF1("fit", "expo(pol2)", 0, 10); 
  fit->SetParameters(1,1,0,10000);
  TF1 *fit2 = new TF1("fit2", "gaus",  10, 30);
  
  //TF1 *total = new TF1("total", "[0]*10^([1]*x² + [2]*x +[3] + gaus", 0,30);
  TCanvas *c1 = new TCanvas("combined_fit","combined_fit");
  
  c1->Divide(3);
  
  c1->cd(1);
  gPad->SetLogy();
  totFiberLy1Vec.at(1)->Fit("fit" , "R");
  totFiberLy1Vec.at(1)->SetFillColor(kAzure-9);
  totFiberLy1Vec.at(1)->GetXaxis()->SetLabelSize(0.03);
  totFiberLy1Vec.at(1)->GetYaxis()->SetLabelSize(0.03);
  totFiberLy1Vec.at(1)->GetXaxis()->SetTitleSize(0.04);
  totFiberLy1Vec.at(1)->GetYaxis()->SetTitleSize(0.04);	  
  totFiberLy1Vec.at(1)->Draw("");
  
  totFiberLy1Vec.at(1)->GetListOfFunctions()->Remove(totFiberLy1Vec.at(1)->GetFunction("expo(pol2)"));
  
  c1->cd(2);
  gPad->SetLogy();
  totFiberLy1Vec.at(1)->Fit("fit2" , "R+");
  totFiberLy1Vec.at(1)->SetFillColor(kAzure-9);
  totFiberLy1Vec.at(1)->GetXaxis()->SetLabelSize(0.03);
  totFiberLy1Vec.at(1)->GetYaxis()->SetLabelSize(0.03);
  totFiberLy1Vec.at(1)->GetXaxis()->SetTitleSize(0.04);
  totFiberLy1Vec.at(1)->GetYaxis()->SetTitleSize(0.04);	 
  totFiberLy1Vec.at(1)->Draw("");
  
  totFiberLy1Vec.at(1)->GetListOfFunctions()->Remove(totFiberLy1Vec.at(1)->GetFunction("gaus"));
  
  c1->cd(3);
  gPad->SetLogy();
  totFiberLy1Vec.at(1)->Fit("fit" , "R");
  totFiberLy1Vec.at(1)->Fit("fit2" , "R+");
  totFiberLy1Vec.at(1)->SetFillColor(kAzure-9);
  totFiberLy1Vec.at(1)->GetXaxis()->SetLabelSize(0.03);
  totFiberLy1Vec.at(1)->GetYaxis()->SetLabelSize(0.03);
  totFiberLy1Vec.at(1)->GetXaxis()->SetTitleSize(0.04);
  totFiberLy1Vec.at(1)->GetYaxis()->SetTitleSize(0.04);	 
  totFiberLy1Vec.at(1)->Draw("");
	
 fout->WriteObject(c1, c1->GetName());
  
  
  
  for(auto& hist : totFiberLy1Vec)  { 
	  if(hist->GetEntries() != 0) { 
	  hist->SetFillColor(kAzure-9);
	  hist->GetXaxis()->SetLabelSize(0.03);
	  hist->GetYaxis()->SetLabelSize(0.03);
	  hist->GetXaxis()->SetTitleSize(0.04);
	  hist->GetYaxis()->SetTitleSize(0.04);	  
	  fout->WriteObject(hist, hist->GetName());
	  } 
  }
 
  fout->Close();
}

int main(int argc, char** argv)
{
  char    inputFile[512]="";
  char    outputFile[512]="fitCosmicss_output.root";
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

  printf("\n\n%sRunning fitCosmicss%s\n\n",text::BOLD,text::RESET);
  
  fitCosmicss(inputFile,outputFile,procNr);

  printf("\n\n%s%sDONE!%s\n\n",text::BOLD,text::GRN,text::RESET);
}
