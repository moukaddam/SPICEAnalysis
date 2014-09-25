// S3 Calibration
// Calibrates each S3 channel in the fragment tree according to either the default triple-alpha
// or a user-specified source (in format energy error\n).  Outputs channel and energy to a file 
// which another program (today named vectorFit.C) will fit a function.
// --------------
// Usage::
// root -l
// .x S3Cal.C++
// .x S3Cal.C++("filename")
// 
// First option runs the default triple-alpha calibration
// Second option takes a user input source file to calibrate with
// ---------------
// 25 Sep 14 :: evitts  :: Created file
//
// ---------------
// TODO: User input filename
// ---------------
// Variable nomenclature (though sometimes subjective)::
// gName = Created in global
// fName = Created in main function
// iName = Created for purposes of iteration (e.g. in a for loop)
// sName = Created in a small scope e.g. for loop
// ------------------------------------------------------------------------------------------------

#include <iostream>
#include <fstream>
#include <vector>

#include "TError.h"
#include "TFile.h"
#include "TChain.h"
#include "TSpectrum.h"
#include "TString.h"
#include "TH1F.h"
#include "TF1.h"

using namespace std;

// ------------------------------------------------------------------------------------------------

// Globals, output files & root file structure
Bool_t gAlphaOption = "True";

TFile *gRootOutput = new TFile("S3PeakSearch.root", "recreate");
TDirectory *gRootHistoDir = gRootOutput->mkdir("1. Histograms");
TDirectory *gRootSearchDir = gRootOutput->mkdir("2. Peak Search");
TDirectory *gRootGaussDir = gRootOutput->mkdir("3. Peak Fit");

vector<Double_t> gEnergies;
vector<Double_t> gEnergyErrors;

// Initialize functions
Int_t fPeakSearch(Int_t iChannel, TH1F* fHistogram);

TString fGetArgument(TString Var, TString histName) {
	TString Argument = Var+">>"+histName;
	return Argument; 
}

// ------------------------------------------------------------------------------------------------

void S3Cal(TString fCalOption="") {

  // This is used to hide warnings produced by root
  // (warnings present because it wants to find more peaks than we desire)
  gErrorIgnoreLevel = 3000;
  
  // Create output calibration data file
  ofstream fOutputCal("S3.fitme");
  fOutputCal.close();

  // Set custom energies or triple alpha energies
  if (fCalOption!="") {
    gAlphaOption = "False";
    ifstream sInputCal(fCalOption);
    Double_t sE, sEd;
    while (true) {
      sInputCal >> sE >> sEd;
      if (sInputCal.eof()) break;
      gEnergies.push_back( sE );
      gEnergyErrors.push_back( sEd );  
    }
    sInputCal.close();
    cout << "User: custom energies read from file" << endl;
  } else {
    cout << "Default: triple-alpha calibration" << endl;
    gEnergies.push_back( 5244.50 ); // 239Pu
    gEnergyErrors.push_back( 0.23 );
    gEnergies.push_back( 5637.81 ); // 241Am
    gEnergyErrors.push_back( 0.12 );
    gEnergies.push_back( 5901.61 ); // 244Cm
    gEnergyErrors.push_back( 0.05 );
  }

  // Open fragment tree
  TChain *fChain = new TChain("FragmentTree");
  fChain->Add("/data2/evitts/SpiceTestSep2014/FragmentTrees/fragment30126_000.root"); // TODO user option
  
  // S3 Details
  Int_t fNbrSectors = 32;
  Int_t fNbrRings = 24;
  Int_t fStartChannel = 1000;
  Int_t fFinalChannel = fStartChannel + fNbrSectors + fNbrRings;
  
  // For each channel ..
  for (Int_t iChannel=fStartChannel; iChannel<=fFinalChannel; iChannel++) {
    // .. draw spectrum and output to histogram
    TString sCut = Form("ChannelNumber==%i",iChannel);
    TH1F *fHistogram = new TH1F("fHistogram", sCut, 200, 40000, 80000);
    fChain->Draw( fGetArgument("Charge",fHistogram->GetName() ) ,sCut,"");
    // .. add histo to relevant directory in root file
    TString fHistogramName = Form("Chrg%i",iChannel);
    fHistogram->SetName(fHistogramName);
    gRootHistoDir->WriteTObject(fHistogram);
    // .. peak search function
    fPeakSearch(iChannel, fHistogram);
    // .. delete pointers
    fHistogram->Delete();
  } // END for (Int_t sChannel ...)

} // END S3Cal()

Int_t fPeakSearch(Int_t iChannel, TH1F* fHistogram) {

  // Number of peaks to find, depends on what calibration source is
  Int_t sPeaks = 3;
  if (!gAlphaOption) sPeaks = gEnergies.size();
  
  // Find peaks 
  TSpectrum *sSpectrum = new TSpectrum( sPeaks, 1.0); // (maxpeaks, resolution)
  Int_t sPeaksFound = sSpectrum->Search(fHistogram, 1.0, "", 0.05); // (TH1* hist, sigma, option, threshold)
  Float_t *sPeaksPosition = sSpectrum->GetPositionX();
  
  // If number of found peaks is not expected, discontinue
  if (sPeaksFound != sPeaks) {
    cout << fHistogram->GetName() << " has only " << sPeaksFound << " peaks, skipped." << endl;
    return 0;
  }

  // Re-order peaks in increasing numerical order
  for (Int_t i=0; i<sPeaksFound; i++) {
    for (Int_t j=0; j<(sPeaksFound-1); j++) {
      if (sPeaksPosition[j]>sPeaksPosition[j+1]) {
        Float_t sStore = sPeaksPosition[j];
        sPeaksPosition[j] = sPeaksPosition[j+1];
        sPeaksPosition[j+1] = sStore;
      }
    }
  }

  // If triple-alpha calibration, check peaks are suitable distance apart
  Double_t sFitWidth = 150.;
  vector<Double_t> sPeakConstant;
  vector<Double_t> sPeakCentroid;
  vector<Double_t> sPeakSigma;
  if (gAlphaOption) {
    Double_t sLeftGap = sPeaksPosition[1] - sPeaksPosition[0];
    Double_t sRightGap = sPeaksPosition[2] - sPeaksPosition[1];
    if (abs(sLeftGap-sRightGap) > 1000) {
      cout << fHistogram->GetName() << " has failed to do a reasonable peak search, skipped." << endl;
      return 0;
    }
    if (sLeftGap < sRightGap) sFitWidth = sLeftGap/2.5; else sFitWidth = sRightGap/2.5;
  }
  
  gRootSearchDir->WriteTObject(fHistogram);
  
  // Then, for each peak ..
  for (Int_t iPeak=0; iPeak<sPeaksFound; iPeak++) {
    // .. fit a function
    TF1 *sGaussFit = new TF1("sGaussFit", "gaus", sPeaksPosition[iPeak]-sFitWidth, sPeaksPosition[iPeak]+sFitWidth);
    fHistogram->Fit(sGaussFit, "RQM+");
    // .. and retrieve the parameters of the fit
    sPeakConstant.push_back( sGaussFit->GetParameter(0) );
    sPeakCentroid.push_back( sGaussFit->GetParameter(1) );
    sPeakSigma.push_back( sGaussFit->GetParameter(2) );
  
  } // END for(Int_t iPeak)

  gRootGaussDir->WriteTObject(fHistogram);

  // Write to file for fitting  
  ofstream sOutputCal("S3.fitme", ios::app);
  sOutputCal << "Segment " << iChannel << endl;
  for (Int_t iPeak=0; iPeak<sPeaksFound; iPeak++) {
    sOutputCal << sPeakCentroid[iPeak] << "\t" << sPeakSigma[iPeak] << "\t";
    sOutputCal << gEnergies[iPeak] << "\t" << gEnergyErrors[iPeak] << endl;
  }
  sOutputCal << "============" << endl;

  return 0;
  
} // END fPeakSearch()
