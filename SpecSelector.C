//c++
#include <iostream>
#include <fstream>
#include <istream>
#include <sstream>
#include <vector>
#include <iomanip>      // std::setprecision
using namespace std ; 

// Root
#include "TROOT.h"
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
#include "TVirtualX.h"



// functions 
void Stop();
void SegmentClicked();
#define STOP Stop();

//Global variables 
	TCanvas *gCanvas ;	

//MAIN FUNCTION 
void SpecSelector() {
   
	gCanvas = new TCanvas("SiLi_Segments","SiLi_Segments",10,10,700,700);	
   gCanvas->ToggleEventStatus();
   //gCanvas->SetGridx();
   //gCanvas->SetGridy();

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

gCanvas->AddExec("ex","SegmentClicked()");
			
}


 

void SegmentClicked() {
	//this action function is called whenever you move the mouse
	//it prints the informations of the segments 
	/*

	I still have a problem with the markers...
	*/
	TObject *select = gPad->GetSelected();
	if(!select) return;
	if (!select->InheritsFrom(TH2Poly::Class())) {gPad->SetUniqueID(0); return;}
	TH2 *h2poly = (TH2Poly*)select;
	gPad->GetCanvas()->FeedbackMode(kTRUE);

	//in case many pads 
	//TPad *pad = (TPad *) gCanvas->GetSelectedPad();
	//if (!pad) return;

	//Get the event
	int event = gPad->GetEvent();

	//Get the abscissa  
	int px = gPad->GetEventX();
	int py = gPad->GetEventY();

	float x = gPad->AbsPixeltoX(px);
	float y = gPad->AbsPixeltoY(py);

	if (event == 51) 
		{

      // Get the sector and the ring 
		//GetRing(x,y) ;  
		//GetSector(x,y) ;

		//erase old position and draw a line at current position
		int pyold = gPad->GetUniqueID();
		float uxmin = gPad->GetUxmin();
		float uxmax = gPad->GetUxmax();
		//cout << "uxmin, uxmax  " << uxmin << " " << uxmax << endl;

		int pxmin = gPad->XtoAbsPixel(uxmin);
		int pxmax = gPad->XtoAbsPixel(uxmax);
		//cout << "pxmin, pxmax  " << pxmin << " " << pxmax << endl;

		if(pyold) gVirtualX->DrawLine(pxmin,pyold,pxmax,pyold);
		gVirtualX->DrawLine(pxmin,y,pxmax,y);
		gPad->SetUniqueID(py);
		//Float_t upy = gPad->AbsPixeltoY(py);
		//Float_t y = gPad->PadtoY(upy);

//TLine *line = new TLine(0,0,0,c2->GetUymax());
//line->Draw();

//create or set the new canvas c2
		TVirtualPad *padsav = gPad;
		TCanvas *c2 = (TCanvas*)gROOT->GetListOfCanvases()->FindObject("c2");
		if(c2) delete c2->GetPrimitive("Projection"); // Delete the (primitive) subpad of the object named "Projection" 
		else   c2 = new TCanvas("c2","Projection Canvas",710,10,400,400);
		c2->SetGrid();
		c2->cd();

//Draw Histogram corresponding to mouse positionon the new Canvas	
		
		//  find bin 
		//  cout<< "bin #" <<h2poly->FindBin(x,y) <<endl ;
		//  Int_t biny = h->GetYaxis()->FindBin(y);

		//call the right histogram for now its dummy
		TH1F *hp = new TH1F("Projection","hist title",100,-20,20);
		TRandom ran;
		ran.SetSeed(x*y) ;
		for (int i = 0; i<1000; i++) {
		hp->Fill(ran.Gaus(0,3));
		}

		// Modify the histogram
		hp->SetFillColor(38);
		hp->Draw();

		//update canvas	
		//gCanvas->Update();					
		c2->Update();
		padsav->cd();

		//return ;
	}

/*
		if (event == 12) {
		cout<< "Middle click" <<endl ;

		//TH2Poly *h2poly = (TH2Poly*)select;
		cout<< "x and y :" <<x << "\t"<< y <<endl ;
		cout<< "bin #" <<h2poly->FindBin(x,y) <<endl ;
      int binxy = h2poly->FindBin(x,y) ;
		cout<< "bin content" <<h2poly->GetBinContent( binxy ) <<endl ;
		}
*/

		return ;

}



void Stop() {
	char c; 
	c = getchar(); 
	if (c=='q') exit(-1) ;
	}


/*
		if (event == 1) {
		cout<< "Right click" <<endl ;
		}
		if (event == 11) {
		cout<< " Right Release" <<endl ;
		}

		if (event == 61) {
		cout<< "Right double click" <<endl ;
		}
*/
