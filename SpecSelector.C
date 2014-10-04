//c++
#include <iostream>
#include <fstream>
#include <istream>
#include <sstream>
#include <vector>
#include <iomanip>      // std::setprecision
using namespace std ; 

// Root
#include "TMath.h"
#include "TFile.h"
#include "TString.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TColor.h"
#include "TGraphErrors.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TPolyLine.h"
#include "TVector2.h"
#include "TDirectory.h"
#include "TH2Poly.h"
#include "TKey.h"
#include "TStyle.h"
#include "TBenchmark.h"
#include "TRandom3.h"


// functions 
void Stop();
 
#define STOP Stop();


//MAIN FUNCTION 
void SpecSelector() {
   
	TCanvas *c1 = new TCanvas("SiLi_Segments","SiLi_Segments",10,10,700,700);	
   c1->ToggleEventStatus();
   c1->SetGridx();
   c1->SetGridy();

   TFile *f;
   f = TFile::Open("SiLiSegmentsGraph.root","READ");

   TH2Poly *p = new TH2Poly("SiLi","SiLi (something in title...)",-20,+20,-20,+20);
   //p->GetXaxis()->SetNdivisions(520);
   p->GetXaxis()->SetTitle("X");
   p->GetYaxis()->SetTitle("Y");

   p->SetContour(100);

   TMultiGraph *mg;
   TKey *key;
	TObject* obj;
   TIter nextkey(f->GetListOfKeys());
   while (key = (TKey*)nextkey()) {
      obj = key->ReadObj();
      if (obj->InheritsFrom("TGraph")) {
         mg = (TMultiGraph*)obj;
			if (mg->GetName()=="all") continue ; 
         p->AddBin(mg);
      }
   }


   gBenchmark->Start("Partitioning");
   p->ChangePartition(100, 100);
   gBenchmark->Show("Partitioning");

   TRandom3 r;
   gBenchmark->Start("Filling");
   for (int i=0; i<5000; i++) {
      float x = r.Uniform(-15,+15);
      float y  = r.Uniform(-15,+15);
      p->Fill(x,y);
   }
   gBenchmark->Show("Filling");

   gStyle->SetOptStat(0000);
   gStyle->SetPalette(1);
   p->Draw("COL");
			
}




void Stop() {
	char c; 
	c = getchar(); 
	if (c=='q') exit(-1) ;
	}

