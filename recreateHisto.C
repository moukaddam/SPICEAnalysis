// Recreate Histograms
// Recreates the MIDAS histogram file from a fragment tree, as the histogram files are used for
// calibration purposes.  There is space to insert a function to check the validity of the wave
// --------------
// Usage::
// root -l
// .x recreateHisto.C++("filename.root")
// ---------------
// 30 Sep 14 :: evitts  :: Created file
//
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
vector<int> gWaveBuffer;

Double_t gIntegration = 125.0;

TBranch *b_TTigFragment_ChannelNumber;
TBranch *b_TTigFragment_Charge;
TBranch *b_TTigFragment_wavebuffer;

Int_t fProcessWaveform(TString fDetectorType);
Int_t fProcessHistogram(TH1F **fHistogram);

// ------------------------------------------------------------------------------------------------

void recreateHisto(TString gInputFile) {

  // Open file
  TChain *fChain = new TChain("FragmentTree");
  fChain->Add(gInputFile);
  fChain->SetMakeClass(1);
  
  // Create output root file
  TFile *fRootOutput = new TFile("hisXXXXX.root", "recreate");
  TFolder *fRootHistoDir = new TFolder("histos", "histos");
  fRootOutput->Add(fRootHistoDir);
  
  // Create histograms
  Int_t fNbrHistos = 1211;
  TH1F *fHistogram[fNbrHistos];
  for (Int_t iHistoNbr=0; iHistoNbr<fNbrHistos; iHistoNbr++) {
    TString sHistoTitle = Form("Chrg%04i", iHistoNbr);
    fHistogram[iHistoNbr] = new TH1F(sHistoTitle, sHistoTitle, 65536, 0, 65536);
  }
  
  // Assign branch addresses to globals
  fChain->SetBranchAddress("ChannelNumber", &gChannelNumber, &b_TTigFragment_ChannelNumber);
  fChain->SetBranchAddress("Charge", &gCharge, &b_TTigFragment_Charge);
  fChain->SetBranchAddress("wavebuffer", &gWaveBuffer, &b_TTigFragment_wavebuffer);
  
  Long64_t fNbrEntries = fChain->GetEntries();
   
  // For each entry ..
  for (Long64_t iEntry=0; iEntry<fNbrEntries; iEntry++) {
    fChain->GetEntry(iEntry);
    Bool_t sGoodWaveform = "False";
    if (gChannelNumber >= 1000 && gChannelNumber < 1060)
      sGoodWaveform = fProcessWaveform("S3");
    else if (gChannelNumber >= 1060 && gChannelNumber <= 1180)
      sGoodWaveform = fProcessWaveform("SPICE");
    if (sGoodWaveform) fProcessHistogram(fHistogram);  
  }
  
  // Write histograms to output root file
  for (Int_t iHistoNbr=0; iHistoNbr<fNbrHistos; iHistoNbr++) 
    fRootHistoDir->Add(fHistogram[iHistoNbr]);
  
  fRootHistoDir->Write();
  fRootOutput->Close();

}

Int_t fProcessWaveform(TString fDetectorType) {
  return 1;
}

Int_t fProcessHistogram(TH1F **fHistogram) {
  gCharge /= gIntegration;
  fHistogram[gChannelNumber]->Fill(gCharge);
  return 0;
}
