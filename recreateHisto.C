// Recreate Histograms
// Recreates the MIDAS histogram file from a fragment tree, as the histogram files are used for
// calibration purposes.  There is space to insert a function to check the validity of the wave
// --------------
// Usage::
// root -l
// .x recreateHisto.C++("filename.root")
// ---------------
// 30 Sep 14 :: evitts  :: Created file
// 06 Oct 14 :: moukaddam  :: updated file to produce other histograms from the fragment tree (ChargeCal, time...)
// ---------------
// Variable nomenclature (though sometimes subjective)::
// gName = Created in global
// fName = Created in main function
// iName = Created for purposes of iteration (e.g. in a for loop)
// sName = Created in a smaller scope e.g. if or while
// ------------------------------------------------------------------------------------------------

#include <iostream>
#include <vector>

#include "TFile.h"
#include "TFolder.h"
#include "TChain.h"
#include "TBranch.h"
#include "TH1F.h"
#include "TString.h"

using namespace std;

// ------------------------------------------------------------------------------------------------

Int_t gChannelNumber;
Int_t gCharge;
Int_t gChargeCal;
Int_t gCfd;
Int_t gLed;
Int_t gTime;
vector<int> gWaveBuffer;

Double_t gIntegration = 125.0;

TBranch *b_TTigFragment_ChannelNumber;
TBranch *b_TTigFragment_Charge;
TBranch *b_TTigFragment_wavebuffer;
TBranch *b_TTigFragment_ChargeCal;
TBranch *b_TTigFragment_Cfd;
TBranch *b_TTigFragment_Led;
TBranch *b_TTigFragment_Time;

Int_t fProcessWaveform(TString fDetectorType);
Int_t fProcessHistogram();

  // Create histograms
Int_t gNbrHistos = 1211;
TH1F *gHistogram[1211];
TH1F *gHistogramCal[1211];
TH1F *gHistogramCfd[1211];
TH1F *gHistogramLed[1211];
TH1F *gHistogramTime[1211];


// ------------------------------------------------------------------------------------------------

void recreateHisto(TString gInputFile) {

  // Open file
  TChain *fChain = new TChain("FragmentTree");
  fChain->Add(gInputFile);
  fChain->SetMakeClass(1);
  
  // Create output root file
  TFile *fRootOutput = new TFile("hisXXXXX.root", "recreate");
	TFolder *fRootHistoDir = new TFolder("histos", "histos");
	TFolder *fRootHistoCalDir = new TFolder("histosCal", "histosCal");
	TFolder *fRootHistoCfdDir = new TFolder("histosCfd", "histosCfd");
	TFolder *fRootHistoLedDir = new TFolder("histosLed", "histosLed");
	TFolder *fRootHistoTimeDir = new TFolder("histosTime", "histosTime");
	fRootOutput->Add(fRootHistoDir);
	fRootOutput->Add(fRootHistoCalDir);   
	fRootOutput->Add(fRootHistoCfdDir);  
	fRootOutput->Add(fRootHistoLedDir);  
	fRootOutput->Add(fRootHistoTimeDir);  
    
  for (Int_t iHistoNbr=0; iHistoNbr<gNbrHistos; iHistoNbr++) {
    TString sHistoTitle = Form("Chrg%04i", iHistoNbr);
    gHistogram[iHistoNbr] = new TH1F(sHistoTitle, sHistoTitle, 65536, 0, 65536);
    sHistoTitle = Form("ChrgCal%04i", iHistoNbr);
    gHistogramCal[iHistoNbr] = new TH1F(sHistoTitle, sHistoTitle, 65536, 0, 65536); // limits not accurate
    sHistoTitle = Form("Cfd%04i", iHistoNbr);
    gHistogramCfd[iHistoNbr] = new TH1F(sHistoTitle, sHistoTitle, 65536, 0, 65536); // limits not accurate
    sHistoTitle = Form("Led%04i", iHistoNbr);
    gHistogramLed[iHistoNbr] = new TH1F(sHistoTitle, sHistoTitle, 65536, 0, 65536); // limits not accurate
    sHistoTitle = Form("Time%04i", iHistoNbr);
    gHistogramTime[iHistoNbr] = new TH1F(sHistoTitle, sHistoTitle, 65536, 0, 65536); // limits not accurate
  }
  
  // Assign branch addresses to globals
  fChain->SetBranchAddress("ChannelNumber", &gChannelNumber, &b_TTigFragment_ChannelNumber);
  fChain->SetBranchAddress("Charge", &gCharge, &b_TTigFragment_Charge);
  fChain->SetBranchAddress("ChargeCal", &gChargeCal, &b_TTigFragment_ChargeCal);
  fChain->SetBranchAddress("Cfd", &gCfd, &b_TTigFragment_Cfd);
  fChain->SetBranchAddress("Led", &gLed, &b_TTigFragment_Led);
  fChain->SetBranchAddress("TimeToTrig", &gTime, &b_TTigFragment_Time);
  fChain->SetBranchAddress("wavebuffer", &gWaveBuffer, &b_TTigFragment_wavebuffer);
  
  Long64_t fNbrEntries = fChain->GetEntries();
   cout << " Nb of Entries = " << fNbrEntries << endl ; 
   
  // For each entry ..
  for (Long64_t iEntry=0; iEntry<fNbrEntries; iEntry++) {
  	if ( iEntry % (5000) == 0 || iEntry==fNbrEntries) printf(" Filling... %2.0f \% \r",100.*iEntry/fNbrEntries);
    fChain->GetEntry(iEntry);
    Bool_t sGoodWaveform = "False";
    if (gChannelNumber >= 1000 && gChannelNumber < 1060)
      sGoodWaveform = fProcessWaveform("S3");
    else if (gChannelNumber >= 1060 && gChannelNumber <= 1180)
      sGoodWaveform = fProcessWaveform("SPICE");
    if (sGoodWaveform) fProcessHistogram();  
  }
  
  // Write histograms to output root file
  for (Int_t iHistoNbr=0; iHistoNbr<gNbrHistos; iHistoNbr++){
	if ( iHistoNbr % (5) == 0 || iHistoNbr==gNbrHistos) printf("Adding histos of Charge to Folder... %2.0f \% \r",100.*iHistoNbr/gNbrHistos);
    fRootHistoDir->Add(gHistogram[iHistoNbr]);
    }
    printf("\nWriting histos... \n");
  fRootHistoDir->Write();
    
  for (Int_t iHistoNbr=0; iHistoNbr<gNbrHistos; iHistoNbr++){
	if ( iHistoNbr % (5) == 0 || iHistoNbr==gNbrHistos) printf("Adding histos of ChargeCal to Folder... %2.0f \% \r",100.*iHistoNbr/gNbrHistos);
    fRootHistoCalDir->Add(gHistogramCal[iHistoNbr]);
    }
    printf("\nWriting histos... \n");
  fRootHistoCalDir->Write();
   
  for (Int_t iHistoNbr=0; iHistoNbr<gNbrHistos; iHistoNbr++){
	if ( iHistoNbr % (5) == 0 || iHistoNbr==gNbrHistos) printf("Adding histos of Cfd to Folder... %2.0f \% \r",100.*iHistoNbr/gNbrHistos);
    fRootHistoCfdDir->Add(gHistogramCfd[iHistoNbr]);    
    } 
    printf("\nWriting histos... \n");
  fRootHistoCfdDir->Write();

  for (Int_t iHistoNbr=0; iHistoNbr<gNbrHistos; iHistoNbr++){
	if ( iHistoNbr % (5) == 0 || iHistoNbr==gNbrHistos) printf("Adding histos of Led to Folder... %2.0f \% \r",100.*iHistoNbr/gNbrHistos); 
    fRootHistoLedDir->Add(gHistogramLed[iHistoNbr]); 
    } 
    printf("\nWriting histos... \n");
  	fRootHistoLedDir->Write();

  for (Int_t iHistoNbr=0; iHistoNbr<gNbrHistos; iHistoNbr++){
	if ( iHistoNbr % (5) == 0 || iHistoNbr==gNbrHistos) printf("Adding histos of timeToTrig to Folder... %2.0f \% \r",100.*iHistoNbr/gNbrHistos);
    fRootHistoTimeDir->Add(gHistogramTime[iHistoNbr]); 
    } 
    printf("\nWriting histos... \n");
  	fRootHistoTimeDir->Write();

//Close       
  fRootOutput->Close();

}

Int_t fProcessWaveform(TString fDetectorType) {
  return 1;
}

Int_t fProcessHistogram() {
  gCharge /= gIntegration;
  gHistogram[gChannelNumber]->Fill(gCharge);
  gHistogramCal[gChannelNumber]->Fill(gChargeCal);
  gHistogramCfd[gChannelNumber]->Fill(gCfd);
  gHistogramLed[gChannelNumber]->Fill(gLed);
  gHistogramTime[gChannelNumber]->Fill(gTime);
  return 0;
}
