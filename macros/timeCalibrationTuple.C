#include <TTree.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
#include <TNtupleD.h>
#include <Rtypes.h>
#include <TCanvas.h>


#include <cstdlib>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <unistd.h>



#include "Utility.h"

///< usage: ./timeCalibrationTuple -i inputfile -o outputfileName -n numberOfSignalsToBeProcessed
///< n = -1 by default which means the whole file is processed

extern char* optarg;

void timeCalibrationTuple(const char *inputFile, const char *outputFile, ULong_t procNr)
{  
   
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
  
  
  //std::vector<TH2D*> tdcTimeCalibVec{};
  std::vector<TH2D*> tdcTimeUnCalibVec{};
  for (Int_t i=0; i<10;i++) {
	  tdcTimeUnCalibVec.emplace_back(new TH2D(Form("hTimeUnCalib120%i",i),Form("TDC 120%i Time uncalibrated; chID; TimeDiff",i),33,0,33,2000,-10,10));
	}
	
	  tdcTimeUnCalibVec.emplace_back(new TH2D(Form("hTimeUnCalib120A"),Form("TDC 120A Time uncalibrated; chID; TimeDiff"),33,0,33,2000,-10,10));
	  tdcTimeUnCalibVec.emplace_back(new TH2D(Form("hTimeUnCalib120B"),Form("TDC 120B Time uncalibrated; chID; TimeDiff"),33,0,33,2000,-10,10));
  
  
  //define histograms
  for(Int_t i=0; i<32; i++) {
	timeChannelVec1200.emplace_back(new TH1D(Form("hTime1200Ch%i",i+1),Form("TimeStamp distribution in TDC 1200 Channel %i; TimeDiff; Counts", i+1), 2000,-10,10));
	timeChannelVec1201.emplace_back(new TH1D(Form("hTime1201Ch%i",i+1),Form("TimeStamp distribution in TDC 1201 Channel %i; TimeDiff; Counts", i+1), 2000,-10,10));
	timeChannelVec1202.emplace_back(new TH1D(Form("hTime1202Ch%i",i+1),Form("TimeStamp distribution in TDC 1202 Channel %i; TimeDiff; Counts", i+1), 2000,-10,10));
	timeChannelVec1203.emplace_back(new TH1D(Form("hTime1203Ch%i",i+1),Form("TimeStamp distribution in TDC 1203 Channel %i; TimeDiff; Counts", i+1), 2000,-10,10));
	timeChannelVec1204.emplace_back(new TH1D(Form("hTime1204Ch%i",i+1),Form("TimeStamp distribution in TDC 1204 Channel %i; TimeDiff; Counts", i+1), 2000,-10,10));
	timeChannelVec1205.emplace_back(new TH1D(Form("hTime1205Ch%i",i+1),Form("TimeStamp distribution in TDC 1205 Channel %i; TimeDiff; Counts", i+1), 2000,-10,10));
	timeChannelVec1206.emplace_back(new TH1D(Form("hTime1206Ch%i",i+1),Form("TimeStamp distribution in TDC 1206 Channel %i; TimeDiff; Counts", i+1), 2000,-10,10));
	timeChannelVec1207.emplace_back(new TH1D(Form("hTime1207Ch%i",i+1),Form("TimeStamp distribution in TDC 1207 Channel %i; TimeDiff; Counts", i+1), 2000,-10,10));
	timeChannelVec1208.emplace_back(new TH1D(Form("hTime1208Ch%i",i+1),Form("TimeStamp distribution in TDC 1208 Channel %i; TimeDiff; Counts", i+1), 2000,-10,10));
	timeChannelVec1209.emplace_back(new TH1D(Form("hTime1209Ch%i",i+1),Form("TimeStamp distribution in TDC 1209 Channel %i; TimeDiff; Counts", i+1), 2000,-10,10));
	timeChannelVec120A.emplace_back(new TH1D(Form("hTime120ACh%i",i+1),Form("TimeStamp distribution in TDC 120A Channel %i; TimeDiff; Counts", i+1), 2000,-10,10));
	timeChannelVec120B.emplace_back(new TH1D(Form("hTime120BCh%i",i+1),Form("TimeStamp distribution in TDC 120B Channel %i; TimeDiff; Counts", i+1), 2000,-10,10));
	}
	
	
  TFile* f = TFile::Open(inputFile);

  if (f->IsOpen()==kFALSE){
    printf("\n\n%s%sOpen input file failed%s\n\n",text::BOLD,text::RED,text::RESET);
    exit(1);
  }

  TNtupleD *signals = (TNtupleD*)f->Get("Signals");

  Double_t eventNr(-1), chID(0), TDC(-1), layer(-1), x(-1), y(-1), signalNr(-1), timeStamp(-1), ToT(-1), padiwaConfig(-1), refTime(-1);
  Int_t prevSigNr(0), prevCh(-1), prevEventNr(-1), firstCounter(0), secondCounter(0);


  ULong_t nSignals = procNr;

  signals->SetBranchAddress("EventNr",      &eventNr);
  signals->SetBranchAddress("timeStamp",    &timeStamp);    // in seconds
  signals->SetBranchAddress("ToT",          &ToT);          // in seconds
  signals->SetBranchAddress("chID",         &chID);
  signals->SetBranchAddress("TDC",          &TDC);          // 0 = TDC1500, 1= 1510 etc
  signals->SetBranchAddress("layer",        &layer);        // 1-8
  signals->SetBranchAddress("x",            &x);            // odd layers have  x != 0
  signals->SetBranchAddress("y",            &y);            // even layers have y != 0
  signals->SetBranchAddress("signalNr",     &signalNr);     // Nth Signal per channel and event
  signals->SetBranchAddress("padiwaConfig", &padiwaConfig);
  signals->SetBranchAddress("refTime",      &refTime);

  TFile *fout = new TFile(Form("%s",outputFile),"recreate");
  

  if ((nSignals == -1) || (nSignals > signals->GetEntries())) { nSignals = signals->GetEntries(); }

  printf("signals to process: %lu\t %.1f%% of the file\n", nSignals, Float_t(100*nSignals)/Float_t(signals->GetEntries()));
  
  std::vector<Double_t> RefTimeStampVec;
  std::vector<Double_t> RefEventNrVec;
  Int_t refchID = 2, refTDC=0;
  for (ULong_t entry = 0; entry < nSignals; entry++) {
	  signals->GetEntry(entry);
	  if(TDC == 0 && chID == 2 && signalNr ==1) {
		  RefTimeStampVec.emplace_back(timeStamp-refTime);
		  RefEventNrVec.emplace_back(eventNr);
	  }
  }
   
  
  for (ULong_t entry = 0; entry < nSignals; entry++) { //loop over all Signals
    if ((((entry+1)%100000) == 0) || (entry == (nSignals-1))) {
      printf("\rprocessing signal %lu...", entry+1);
      fflush(stdout);
      std::cout<<std::setw(5)<<std::setiosflags(std::ios::fixed)<<std::setprecision(1)<<" "<<(100.*(entry+1))/nSignals<<" % done\r"<<std::flush;
    }

    signals->GetEntry(entry);
    
    Double_t timeRefChannel = 0;
	    
	if (signalNr == 1) {
	  for (Int_t i = 0; i<RefEventNrVec.size();i++) {
	    if (RefEventNrVec.at(i) == eventNr) {
		  timeRefChannel = RefTimeStampVec.at(i);
		  break;
		}
     }
				
			if (TDC == 0 ) {
				timeChannelVec1200.at(chID-1)->Fill((timeStamp-refTime-timeRefChannel)*1e9); 
				}
				
			else if (TDC == 1)  { 
				timeChannelVec1201.at(chID-1)->Fill((timeStamp-refTime-timeRefChannel)*1e9);
				}
				
			else if (TDC == 2)  { 
				timeChannelVec1202.at(chID-1)->Fill((timeStamp-refTime-timeRefChannel)*1e9);
				}
				
			else if (TDC == 3)  { 
				timeChannelVec1203.at(chID-1)->Fill((timeStamp-refTime-timeRefChannel)*1e9);
			    }
			    
			else if (TDC == 4)  {
				timeChannelVec1204.at(chID-1)->Fill((timeStamp-refTime-timeRefChannel)*1e9);
				}
				
			else if (TDC == 5)  { 
				timeChannelVec1205.at(chID-1)->Fill((timeStamp-refTime-timeRefChannel)*1e9);
				}
				
			else if (TDC == 6)  { 
				timeChannelVec1206.at(chID-1)->Fill((timeStamp-refTime-timeRefChannel)*1e9);
				}
				
			else if (TDC == 7)  { 
				timeChannelVec1207.at(chID-1)->Fill((timeStamp-refTime-timeRefChannel)*1e9);
				}
				
			else if (TDC == 8)  { 
				timeChannelVec1208.at(chID-1)->Fill((timeStamp-refTime-timeRefChannel)*1e9);
				}
				
			else if (TDC == 9)  { 
				timeChannelVec1209.at(chID-1)->Fill((timeStamp-refTime-timeRefChannel)*1e9);
				}
				
			else if (TDC == 10) { 
				timeChannelVec120A.at(chID-1)->Fill((timeStamp-refTime-timeRefChannel)*1e9);
				}
				
			else if (TDC == 11) {  
				timeChannelVec120B.at(chID-1)->Fill((timeStamp-refTime-timeRefChannel)*1e9);
				}
				
				
		    tdcTimeUnCalibVec.at(TDC)->Fill(chID, (timeStamp-refTime-timeRefChannel)*1e9);
			}
			
		 //std::cout<<tdcTimeUnCalibVec.size()<<std::endl;
      
	    //tdcTimeUnCalibVec.at(TDC)->Fill(chID, (timeStamp-refTime-timeRefChannel)*1e9);
		
			
		}
		//tdcTimeUnCalibVec.at(TDC)->Fill(chID, (timeStamp-refTime-timeRefChannel)*1e9);
		//copy(tdcTimeUnCalibVec.begin(), tdcTimeUnCalibVec.end(), back_inserter(tdcTimeCalibVec));
	
  //makes gaussian fit for the time diff
  Double_t mean            = -1;
  Double_t sigma           = -1;
  
  TF1 *fit = new TF1("fit", "gaus", -10,10);
    
  std::fstream calibrationValues;
  calibrationValues.open("calibrationValues.txt",std::ios::out);
  
  
  for(Int_t i=0;i<timeChannelVec1200.size();i++) { 
	  if (i != refchID -1 ){
		timeChannelVec1200.at(i)->Fit("fit", "R"); 
		mean = fit->GetParameter(1); 
	  }
	  else {mean = 0;}
	  if (i == 0) { calibrationValues << "\n TDC 1200 \n";}
	  calibrationValues << mean << ", ";
  }
  
  for(Int_t i=0;i<timeChannelVec1201.size();i++) { timeChannelVec1201.at(i)->Fit("fit", "R"); mean = fit->GetParameter(1); if (i == 0) { calibrationValues << "\n TDC 1201 \n";} calibrationValues << mean << ", ";}
  for(Int_t i=0;i<timeChannelVec1202.size();i++) { timeChannelVec1202.at(i)->Fit("fit", "R"); mean = fit->GetParameter(1); if (i == 0) { calibrationValues << "\n TDC 1202 \n";} calibrationValues << mean << ", ";}
  for(Int_t i=0;i<timeChannelVec1203.size();i++) { timeChannelVec1203.at(i)->Fit("fit", "R"); mean = fit->GetParameter(1); if (i == 0) { calibrationValues << "\n TDC 1203 \n";} calibrationValues << mean << ", ";}
  for(Int_t i=0;i<timeChannelVec1204.size();i++) { timeChannelVec1204.at(i)->Fit("fit", "R"); mean = fit->GetParameter(1); if (i == 0) { calibrationValues << "\n TDC 1204 \n";} calibrationValues << mean << ", ";}
  for(Int_t i=0;i<timeChannelVec1205.size();i++) { timeChannelVec1205.at(i)->Fit("fit", "R"); mean = fit->GetParameter(1); if (i == 0) { calibrationValues << "\n TDC 1205 \n";} calibrationValues << mean << ", ";}
  for(Int_t i=0;i<timeChannelVec1206.size();i++) { timeChannelVec1206.at(i)->Fit("fit", "R"); mean = fit->GetParameter(1); if (i == 0) { calibrationValues << "\n TDC 1206 \n";} calibrationValues << mean << ", ";}
  for(Int_t i=0;i<timeChannelVec1207.size();i++) { timeChannelVec1207.at(i)->Fit("fit", "R"); mean = fit->GetParameter(1); if (i == 0) { calibrationValues << "\n TDC 1207 \n";} calibrationValues << mean << ", ";}
  for(Int_t i=0;i<timeChannelVec1208.size();i++) { timeChannelVec1208.at(i)->Fit("fit", "R"); mean = fit->GetParameter(1); if (i == 0) { calibrationValues << "\n TDC 1208 \n";} calibrationValues << mean << ", ";}
  for(Int_t i=0;i<timeChannelVec1209.size();i++) { timeChannelVec1209.at(i)->Fit("fit", "R"); mean = fit->GetParameter(1); if (i == 0) { calibrationValues << "\n TDC 1209 \n";} calibrationValues << mean << ", ";}
  for(Int_t i=0;i<timeChannelVec120A.size();i++) { timeChannelVec120A.at(i)->Fit("fit", "R"); mean = fit->GetParameter(1); if (i == 0) { calibrationValues << "\n TDC 120A \n";} calibrationValues << mean << ", ";}
  for(Int_t i=0;i<timeChannelVec120B.size();i++) { timeChannelVec120B.at(i)->Fit("fit", "R"); mean = fit->GetParameter(1); if (i == 0) { calibrationValues << "\n TDC 120B \n";} calibrationValues << mean << ", ";}
  
   
   calibrationValues.close();

   
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
  
  //for(auto& hist : tdcTimeCalibVec)   { fout->WriteObject(hist, hist->GetName());}
  for(auto& hist : tdcTimeUnCalibVec) { 
	  fout->WriteObject(hist, hist->GetName());
	  hist->GetXaxis()->SetTitleSize(0.04);
	  hist->GetYaxis()->SetTitleSize(0.04);
	  gPad->SetLogy();	  
	  }
  
  
  fout->Close();

}

int main(int argc, char** argv)
{
  char    inputFile[512]="";
  char    outputFile[512]="timeCalibrationTuple_output.root";
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

  printf("\n\n%sRunning timeCalibrationTuple%s\n\n",text::BOLD,text::RESET);
  
  timeCalibrationTuple(inputFile,outputFile,procNr);

  printf("\n\n%s%sDONE!%s\n\n",text::BOLD,text::GRN,text::RESET);
}
