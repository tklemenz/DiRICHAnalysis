#include <TH2.h>
#include <TH1.h>
#include <TTree.h>
#include <TFile.h>
#include <TNtupleD.h>
#include <TCanvas.h>
#include <Rtypes.h>

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <unistd.h>

#include "Utility.h"

///< usage: ./plotTimeCalib -i inputfile -o outputfileName -n numberOfSignalsToBeProcessed
///< n = -1 by default which means the whole file is processed

//rudimentary, only works with "calibrationValues.txt" in same folder
extern char* optarg;

void plotTimeCalib(const char *inputFile, const char *outputFile, ULong_t procNr)
{
  TFile* f = TFile::Open(inputFile);

  if (f->IsOpen()==kFALSE){
    printf("\n\n%s%sOpen input file failed%s\n\n",text::BOLD,text::RED,text::RESET);
    exit(1);
  }
  
  std::vector<TH2D*> tdcTimeCalibVec{};
  std::vector<TH2D*> tdcTimeUnCalibVec{};
  
      tdcTimeUnCalibVec.emplace_back(new TH2D("hTimeUnCalib1200","TDC 1200 Time uncalibrated; chID; TimeDiff",33,0,33,2000,-10,10));
	  tdcTimeUnCalibVec.emplace_back(new TH2D("hTimeUnCalib1201","TDC 1201 Time uncalibrated; chID; TimeDiff",33,0,33,2000,-10,10));
	  tdcTimeUnCalibVec.emplace_back(new TH2D("hTimeUnCalib1202","TDC 1202 Time uncalibrated; chID; TimeDiff",33,0,33,2000,-10,10));
	  tdcTimeUnCalibVec.emplace_back(new TH2D("hTimeUnCalib1203","TDC 1203 Time uncalibrated; chID; TimeDiff",33,0,33,2000,-10,10));	
	  tdcTimeUnCalibVec.emplace_back(new TH2D("hTimeUnCalib1204","TDC 1204 Time uncalibrated; chID; TimeDiff",33,0,33,2000,-10,10));
	  tdcTimeUnCalibVec.emplace_back(new TH2D("hTimeUnCalib1205","TDC 1205 Time uncalibrated; chID; TimeDiff",33,0,33,2000,-10,10));
	  tdcTimeUnCalibVec.emplace_back(new TH2D("hTimeUnCalib1206","TDC 1206 Time uncalibrated; chID; TimeDiff",33,0,33,2000,-10,10));
	  tdcTimeUnCalibVec.emplace_back(new TH2D("hTimeUnCalib1207","TDC 1207 Time uncalibrated; chID; TimeDiff",33,0,33,2000,-10,10));
	  tdcTimeUnCalibVec.emplace_back(new TH2D("hTimeUnCalib1208","TDC 1208 Time uncalibrated; chID; TimeDiff",33,0,33,2000,-10,10));
	  tdcTimeUnCalibVec.emplace_back(new TH2D("hTimeUnCalib1209","TDC 1209 Time uncalibrated; chID; TimeDiff",33,0,33,2000,-10,10));
	  tdcTimeUnCalibVec.emplace_back(new TH2D("hTimeUnCalib120A","TDC 120A Time uncalibrated; chID; TimeDiff",33,0,33,2000,-10,10));
	  tdcTimeUnCalibVec.emplace_back(new TH2D("hTimeUnCalib120B","TDC 120B Time uncalibrated; chID; TimeDiff",33,0,33,2000,-10,10));
  

  //for (Int_t i=0; i<10;i++) {
	  //tdcTimeCalibVec.emplace_back(new TH2D(Form("hTimeCalib120%i",i),Form("TDC 120%i Time calibrated; chID; TimeDiff",i),33,0,33,2000,-10,10));//gives warning and messes up the Histo title
	//}
	
      	

      tdcTimeCalibVec.emplace_back(new TH2D("hTimeCalib1200","TDC 1200 Time calibrated; chID; TimeDiff",33,0,33,2000,-10,10));
	  tdcTimeCalibVec.emplace_back(new TH2D("hTimeCalib1201","TDC 1201 Time calibrated; chID; TimeDiff",33,0,33,2000,-10,10));
	  tdcTimeCalibVec.emplace_back(new TH2D("hTimeCalib1202","TDC 1202 Time calibrated; chID; TimeDiff",33,0,33,2000,-10,10));
	  tdcTimeCalibVec.emplace_back(new TH2D("hTimeCalib1203","TDC 1203 Time calibrated; chID; TimeDiff",33,0,33,2000,-10,10));	
	  tdcTimeCalibVec.emplace_back(new TH2D("hTimeCalib1204","TDC 1204 Time calibrated; chID; TimeDiff",33,0,33,2000,-10,10));
	  tdcTimeCalibVec.emplace_back(new TH2D("hTimeCalib1205","TDC 1205 Time calibrated; chID; TimeDiff",33,0,33,2000,-10,10));
	  tdcTimeCalibVec.emplace_back(new TH2D("hTimeCalib1206","TDC 1206 Time calibrated; chID; TimeDiff",33,0,33,2000,-10,10));
	  tdcTimeCalibVec.emplace_back(new TH2D("hTimeCalib1207","TDC 1207 Time calibrated; chID; TimeDiff",33,0,33,2000,-10,10));
	  tdcTimeCalibVec.emplace_back(new TH2D("hTimeCalib1208","TDC 1208 Time calibrated; chID; TimeDiff",33,0,33,2000,-10,10));
	  tdcTimeCalibVec.emplace_back(new TH2D("hTimeCalib1209","TDC 1209 Time calibrated; chID; TimeDiff",33,0,33,2000,-10,10));
	  tdcTimeCalibVec.emplace_back(new TH2D("hTimeCalib120A","TDC 120A Time calibrated; chID; TimeDiff",33,0,33,2000,-10,10));
	  tdcTimeCalibVec.emplace_back(new TH2D("hTimeCalib120B","TDC 120B Time calibrated; chID; TimeDiff",33,0,33,2000,-10,10));
	  
  

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
   
   //read input file
   std::vector<Double_t> CalibValVec;
   std::ifstream file ("calibrationValues.txt");
   Double_t inputVal;
   if(file.is_open()){ std::cout<<"yuhu"<<std::endl;}
   else { std::cout<<"I didnt work"<<std::endl;}
   while(file >> inputVal){
	   CalibValVec.emplace_back(inputVal);
   }
   
   
  
   //for(Int_t i=0;i<272;i++) { CalibValVec.emplace_back(0);}//whack fix for out of range bcs file not long enough
    //for(auto& val : CalibValVec){std::cout<<val<<std::endl;} 
   
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


  for (ULong_t entry = 0; entry < nSignals; entry++) {
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
     
    tdcTimeUnCalibVec.at(TDC)->Fill(chID, (timeStamp-refTime-timeRefChannel)*1e9); 
    tdcTimeCalibVec.at(TDC)->Fill(chID, (timeStamp-refTime-timeRefChannel)*1e9-CalibValVec.at(TDC*32+chID-1)); 
  
  }
  
  }
  for(auto& hist : tdcTimeUnCalibVec)  { 
	 hist->GetXaxis()->SetTitleSize(0.04);
	 hist->GetYaxis()->SetTitleSize(0.04);
	 //gPad->SetLogz();	  
	 fout->WriteObject(hist, hist->GetName());
       
	  }  
  
  
  for(auto& hist : tdcTimeCalibVec)  { 
	 hist->GetXaxis()->SetTitleSize(0.04);
	 hist->GetYaxis()->SetTitleSize(0.04);
	 //gPad->SetLogz();	  
	 fout->WriteObject(hist, hist->GetName());
       
	  }  
  fout->Close();

}

int main(int argc, char** argv)
{
  char    inputFile[512]="";
  char    outputFile[512]="plotTimeCalib_output.root";
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

  printf("\n\n%sRunning plotTimeCalib%s\n\n",text::BOLD,text::RESET);
  
  plotTimeCalib(inputFile,outputFile,procNr);

  printf("\n\n%s%sDONE!%s\n\n",text::BOLD,text::GRN,text::RESET);
}
