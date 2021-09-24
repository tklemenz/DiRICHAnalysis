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
#include <fstream>
#include <iomanip>
#include <unistd.h>

#include "CTSEvent.h"
#include "Utility.h"

//#include <fmt/format.h>
//#include <boost/log/trivial.hpp>

///< usage: ./timeCalibration -i inputfile -o outputfile -n numberOfEventsToBeProcessed
///< n = -1 by default which means the whole file is processed

extern char* optarg;

void timeCalibration(const TString inputFiles, const char *outputFile, ULong_t procNr)
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
  Int_t   chID       = -1;
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
  std::vector<Fiber>  fibersMod0; 
  std::vector<Fiber>  fibersMod1; 
  std::vector<Fiber>  fibers;
  std::vector<Fiber>  signals; //basically just to rename fibers to signals, both are vector containing vector of fibers
  std::vector<Signal> SignalsInCurrentChannel;
  std::vector<Signal> signalsRefChannel;
  
  std::vector<TH1D*> timeChannelVec1200{};
  std::vector<TH1D*> timeChannelVec1201{};
  std::vector<TH1D*> timeChannelVec1202{};
  std::vector<TH1D*> timeChannelVec1203{};
  std::vector<TH1D*> timeChannelVec1204{};
  std::vector<TH1D*> timeChannelVec1205{};
  std::vector<TH1D*> timeChannelVec1206{};
  std::vector<TH1D*> timeChannelVec1207{};
  std::vector<TH1D*> timeChannelVec1208{};
  std::vector<TH1D*> timeChannelVec1209{};
  std::vector<TH1D*> timeChannelVec120A{};
  std::vector<TH1D*> timeChannelVec120B{};
  
  
  for(Int_t i=0; i<32; i++) {
	timeChannelVec1200.emplace_back(new TH1D(Form("hTime1200Ch%i",i+1),Form("TimeStamp distribution in TDC 1200 Channel %i", i+1), 1000,-20000,100));
	timeChannelVec1201.emplace_back(new TH1D(Form("hTime1201Ch%i",i+1),Form("TimeStamp distribution in TDC 1201 Channel %i", i+1), 1000,-20000,100));
	timeChannelVec1202.emplace_back(new TH1D(Form("hTime1202Ch%i",i+1),Form("TimeStamp distribution in TDC 1202 Channel %i", i+1), 2000,-100,100));
	timeChannelVec1203.emplace_back(new TH1D(Form("hTime1203Ch%i",i+1),Form("TimeStamp distribution in TDC 1203 Channel %i", i+1), 2000,-100,100));
	timeChannelVec1204.emplace_back(new TH1D(Form("hTime1204Ch%i",i+1),Form("TimeStamp distribution in TDC 1204 Channel %i", i+1), 2000,-100,100));
	timeChannelVec1205.emplace_back(new TH1D(Form("hTime1205Ch%i",i+1),Form("TimeStamp distribution in TDC 1205 Channel %i", i+1), 2000,-100,100));
	timeChannelVec1206.emplace_back(new TH1D(Form("hTime1206Ch%i",i+1),Form("TimeStamp distribution in TDC 1206 Channel %i", i+1), 2000,-100,100));
	timeChannelVec1207.emplace_back(new TH1D(Form("hTime1207Ch%i",i+1),Form("TimeStamp distribution in TDC 1207 Channel %i", i+1), 2000,-100,100));
	timeChannelVec1208.emplace_back(new TH1D(Form("hTime1208Ch%i",i+1),Form("TimeStamp distribution in TDC 1208 Channel %i", i+1), 2000,-100,100));
	timeChannelVec1209.emplace_back(new TH1D(Form("hTime1209Ch%i",i+1),Form("TimeStamp distribution in TDC 1209 Channel %i", i+1), 2000,-100,100));
	timeChannelVec120A.emplace_back(new TH1D(Form("hTime120ACh%i",i+1),Form("TimeStamp distribution in TDC 120A Channel %i", i+1), 2000,-100,100));
	timeChannelVec120B.emplace_back(new TH1D(Form("hTime120BCh%i",i+1),Form("TimeStamp distribution in TDC 120B Channel %i", i+1), 2000,-100,100));
	}
  /*========================================================
  ==========================================================*/
  //finds and sets name for output file
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
  //creates/overwrites output file
  TFile *fout = new TFile(Form("%s", outputName.c_str()), "recreate");
  
  //gets acess to file from convertToCTSEvent or TChain
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
 
       
      eventNr      = event->getEventNr();
      padiwaConfig = event->getPadiwaConfig();
      //modules      = event->getModules();
      Double_t refTimeStamp;
      fibersMod0 = modules.at(0).getFibers();
      signalsRefChannel = fibersMod0.at(0).getSignals();
    
      //printf("Modules size:"); std::cout<<modules.size()<<std::endl;
      ////fibersMod0 = modules.at(0).getFibers();
      //Int_t i = 0;
      //for (auto& module : modules){ 
		  
		  //if (i=0) { fibersMod0 = module.getFibers(); }
		  //else break;	  
		  //i++;
	  //}
	    //printf("Signalss size:"); std::cout<<fibersMod0.size()<<std::endl;
	  //i = 0;
	  //for(auto& fiber : fibersMod0) {
		 
		  //if (i=0) { signalsRefChannel = fiber.getSignals();}
		  //else break;
		  //i++;
	  //}
     //printf("Signal size:"); std::cout<<signalsRefChannel.size()<<std::endl;
      
       //printf("SIgnal size:"); std::cout<<signalsRefChannel.size()<<std::endl;
       ////printf("ElementNr:"); std::cout<<signalsRefChannel.at(i).getTimeStamp()<<std::endl;
       Int_t i = 0;
      for(auto& signal : signalsRefChannel) {
		   //printf("SIgnal size:"); std::cout<<signalsRefChannel.size()<<std::endl;
		   //printf("ElementNr:"); std::cout<<signalsRefChannel.at(i).getTimeStamp()<<std::endl;
		   
		  if (signalsRefChannel.at(i).getSignalNr() == 1) {
		    refTimeStamp = signalsRefChannel.at(i).getTimeStamp();
		  }
		i++;
	  }
  
      
      
      //fibersMod0   = modules.at(0).getFibers();
      //fibersMod1   = modules.at(1).getFibers();
  
      //fibers.resize(fibersMod0.size() + fibersMod1.size());
      //std::copy(fibersMod0.begin(), fibersMod0.end(), fibers.begin());
      //std::copy(fibersMod1.begin(), fibersMod1.end(), fibers.begin() + fibersMod0.size());
     
      Signal signal;
      
      for(auto& module : modules) {//printf("modules size"); std::cout<<modules.size()<<std::endl;
	    fibers       = module.getFibers();
	    signals      = fibers;
        for(Int_t i = 0; i < signals.size(); i++) { 
	      SignalsInCurrentChannel=signals.at(i).getSignals();
		  for(Int_t j =0; j < SignalsInCurrentChannel.size(); j++) {
		    if (SignalsInCurrentChannel.at(j).getSignalNr() == 1){
			  signal=SignalsInCurrentChannel.at(j);
			}
	      }
				
			if (signal.getTDCID() == 0) {
				timeChannelVec1200.at(signal.getChannelID()-1)->Fill(signal.getTimeStamp());
				}
				
			else if (signal.getTDCID() == 1)  { 
				timeChannelVec1201.at(signal.getChannelID()-1)->Fill(signal.getTimeStamp());
				}
				
			else if (signal.getTDCID() == 2)  { 
				timeChannelVec1202.at(signal.getChannelID()-1)->Fill(signal.getTimeStamp());
				}
				
			else if (signal.getTDCID() == 3)  { 
				timeChannelVec1203.at(signal.getChannelID()-1)->Fill(signal.getTimeStamp());
			    }
			    
			else if (signal.getTDCID() == 4)  {
				timeChannelVec1204.at(signal.getChannelID()-1)->Fill(signal.getTimeStamp());
				}
				
			else if (signal.getTDCID() == 5)  { 
				timeChannelVec1205.at(signal.getChannelID()-1)->Fill(signal.getTimeStamp());
				}
				
			else if (signal.getTDCID() == 6)  { 
				timeChannelVec1206.at(signal.getChannelID()-1)->Fill(signal.getTimeStamp());
				}
				
			else if (signal.getTDCID() == 7)  { 
				timeChannelVec1207.at(signal.getChannelID()-1)->Fill(signal.getTimeStamp());
				}
				
			else if (signal.getTDCID() == 8)  { 
				timeChannelVec1208.at(signal.getChannelID()-1)->Fill(signal.getTimeStamp());
				}
				
			else if (signal.getTDCID() == 9)  { 
				timeChannelVec1209.at(signal.getChannelID()-1)->Fill(signal.getTimeStamp());
				}
				
			else if (signal.getTDCID() == 10) { 
				timeChannelVec120A.at(signal.getChannelID()-1)->Fill(signal.getTimeStamp());
				}
				
			else if (signal.getTDCID() == 11) {  
				timeChannelVec120A.at(signal.getChannelID()-1)->Fill(signal.getTimeStamp());
				}
           
                 
        }//loop over all channels and signals
      } //loop over all modules
    } // loop over file
  } // loop over all files

   //TH1D *histRef = new TH1D("histRef", "histRef",2000000,-1000000000,100);
   //TH1D *histTest = new TH1D("histTest", "histTest",200,-100,100);
   //histRef = timeChannelVec1200.at(0);
   //histTest = timeChannelVec1201.at(5);
   //histTest->Add(histRef,-1);
   //fout->WriteObject(histTest, histTest->GetName());    

  //makes gaussian fit for the time stamps
  //Double_t mean            = -1;
  //Double_t sigma           = -1;
  //Double_t meanRefChannel =  -1;
 
  
  //TF1 *fit = new TF1("fit", "gaus");
  
  //timeChannelVec1200.at(0)->Fit("fit", "R");
  //meanRefChannel = fit->GetParameter(1);
  
  //std::fstream calibrationValues;
  //calibrationValues.open("calibrationValues.txt",std::ios::out);
  
  //for (Int_t i=0; i<32; i++) {
  //for(auto& hist : timeChannelVec1200) { hist->Fit("fit", "R"); mean = fit->GetParameter(1); calibrationValues << Form("TDC 1200 Channel %i",i+1) << std::setw(5)<< mean-meanRefChannel << "\n";}
  //for(auto& hist : timeChannelVec1201) { hist->Fit("fit", "R"); mean = fit->GetParameter(1); calibrationValues << Form("TDC 1201 Channel %i",i+1) << std::setw(5)<< mean-meanRefChannel << "\n";}
  //for(auto& hist : timeChannelVec1202) { hist->Fit("fit", "R"); mean = fit->GetParameter(1); calibrationValues << Form("TDC 1202 Channel %i",i+1) << std::setw(5)<< mean-meanRefChannel << "\n";}
  //for(auto& hist : timeChannelVec1203) { hist->Fit("fit", "R"); mean = fit->GetParameter(1); calibrationValues << Form("TDC 1203 Channel %i",i+1) << std::setw(5)<< mean-meanRefChannel << "\n";}
  //for(auto& hist : timeChannelVec1204) { hist->Fit("fit", "R"); mean = fit->GetParameter(1); calibrationValues << Form("TDC 1204 Channel %i",i+1) << std::setw(5)<< mean-meanRefChannel << "\n";}
  //for(auto& hist : timeChannelVec1205) { hist->Fit("fit", "R"); mean = fit->GetParameter(1); calibrationValues << Form("TDC 1205 Channel %i",i+1) << std::setw(5)<< mean-meanRefChannel << "\n";}
  //for(auto& hist : timeChannelVec1206) { hist->Fit("fit", "R"); mean = fit->GetParameter(1); calibrationValues << Form("TDC 1206 Channel %i",i+1) << std::setw(5)<< mean-meanRefChannel << "\n";}
  //for(auto& hist : timeChannelVec1207) { hist->Fit("fit", "R"); mean = fit->GetParameter(1); calibrationValues << Form("TDC 1207 Channel %i",i+1) << std::setw(5)<< mean-meanRefChannel << "\n";}
  //for(auto& hist : timeChannelVec1208) { hist->Fit("fit", "R"); mean = fit->GetParameter(1); calibrationValues << Form("TDC 1208 Channel %i",i+1) << std::setw(5)<< mean-meanRefChannel << "\n";}
  //for(auto& hist : timeChannelVec1209) { hist->Fit("fit", "R"); mean = fit->GetParameter(1); calibrationValues << Form("TDC 1209 Channel %i",i+1) << std::setw(5)<< mean-meanRefChannel << "\n";}
  //for(auto& hist : timeChannelVec120A) { hist->Fit("fit", "R"); mean = fit->GetParameter(1); calibrationValues << Form("TDC 120A Channel %i",i+1) << std::setw(5)<< mean-meanRefChannel << "\n";}
  //for(auto& hist : timeChannelVec120B) { hist->Fit("fit", "R"); mean = fit->GetParameter(1); calibrationValues << Form("TDC 120B Channel %i",i+1) << std::setw(5)<< mean-meanRefChannel << "\n";}
  //}
   
   //calibrationValues.close();
   
  Int_t histCounter = 0;
  for(auto& hist : timeChannelVec1200) { if(hist->GetEntries() != 0) { fout->WriteObject(hist, hist->GetName()); histCounter++; } }
  for(auto& hist : timeChannelVec1201) { if(hist->GetEntries() != 0) { fout->WriteObject(hist, hist->GetName()); histCounter++; } }
  for(auto& hist : timeChannelVec1202) { if(hist->GetEntries() != 0) { fout->WriteObject(hist, hist->GetName()); histCounter++; } }
  for(auto& hist : timeChannelVec1203) { if(hist->GetEntries() != 0) { fout->WriteObject(hist, hist->GetName()); histCounter++; } }
  for(auto& hist : timeChannelVec1204) { if(hist->GetEntries() != 0) { fout->WriteObject(hist, hist->GetName()); histCounter++; } }
  for(auto& hist : timeChannelVec1205) { if(hist->GetEntries() != 0) { fout->WriteObject(hist, hist->GetName()); histCounter++; } }
  for(auto& hist : timeChannelVec1206) { if(hist->GetEntries() != 0) { fout->WriteObject(hist, hist->GetName()); histCounter++; } }
  for(auto& hist : timeChannelVec1207) { if(hist->GetEntries() != 0) { fout->WriteObject(hist, hist->GetName()); histCounter++; } }
  for(auto& hist : timeChannelVec1208) { if(hist->GetEntries() != 0) { fout->WriteObject(hist, hist->GetName()); histCounter++; } }
  for(auto& hist : timeChannelVec1209) { if(hist->GetEntries() != 0) { fout->WriteObject(hist, hist->GetName()); histCounter++; } }
  for(auto& hist : timeChannelVec120A) { if(hist->GetEntries() != 0) { fout->WriteObject(hist, hist->GetName()); histCounter++; } }
  for(auto& hist : timeChannelVec120B) { if(hist->GetEntries() != 0) { fout->WriteObject(hist, hist->GetName()); histCounter++; } } 
 
  
  //makes canvas with "nicer" plots, not yet implemented for time calibration
  //TCanvas *c1 = new TCanvas("cToTDists","cToTDists");
  //c1->DivideSquare(histCounter);

  //Int_t padIter = 1;
  //for(auto& hist : timeChannelVec1209) {
    //if(hist->GetEntries() == 0) { continue; }
    //c1->cd(padIter);
    //gPad->SetLogz();
    //hist->Draw("COLZ");
    //padIter++;
  //}

  //fout->WriteObject(c1, c1->GetName());

  fout->Close();
}

int main(int argc, char** argv)
{
  char    inputFile[512]="";
  char    outputFile[512]="timeCalibration_output.root";
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

  printf("\n\n%sRunning timeCalibration%s\n\n",text::BOLD,text::RESET);
  
  timeCalibration(inputFile,outputFile,procNr);

  printf("\n\n%s%sDONE!%s\n\n",text::BOLD,text::GRN,text::RESET);
}
