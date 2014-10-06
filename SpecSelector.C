//c++
#include <iostream>
#include <fstream>
#include <istream>
#include <sstream>
#include <vector>
#include <iomanip>      // std::setprecision
using namespace std ; 

// Root
#include "TBenchmark.h"
#include "TCanvas.h"
#include "TColor.h"
#include "TDirectory.h"
#include "TEllipse.h"
#include "TGraphErrors.h"
#include "TGraph.h"
#include "TH2Poly.h"
#include "TH1F.h"
#include "TFile.h"
#include "TFolder.h"
#include "TKey.h"
#include "TLine.h"
#include "TMath.h"
#include "TMultiGraph.h"
#include "TPolyLine.h"
#include "TRandom3.h"
#include "TROOT.h"
#include "TString.h"
#include "TStyle.h"
#include "TVirtualX.h"
#include "TVector2.h"


// functions 
void Stop();
void SegmentClicked();
void CreateSiLiMap();
void OpenDataFile();
int  OpenDataThisFile(TString filename);
void GetHistogram(int indicator);
#define STOP Stop();

//Global variables 
TCanvas *gCanvas ;
TCanvas *gCanvasHist ;	
TFolder *gFolderHistos ; 


//MAIN FUNCTION 
void SpecSelector() {

//Global settings
   gStyle->SetOptStat(0000);
   gStyle->SetPalette(1);   

//Open File for the histigram selection 
	OpenDataFile();
// Set the Mapping parameters, basically it will be an ascii file (in case the mapping is right during the experiment, this can be done through the mnemonics.	
	//SetSiLiMap();
	
//Create the main map used by the user	
	CreateSiLiMap();

}


void CreateSiLiMap(){


	if(!gCanvas) {
		//gCanvas->Delete();
		TString title = "Bi207" ;
		gCanvas = new TCanvas("SiLi_Segments","SiLi_Segments",10,10,700,700);
		gCanvas->ToggleEventStatus();
		//gCanvas->SetGridx();
		//gCanvas->SetGridy();
		gCanvas->AddExec("ex","SegmentClicked()");
		//gCanvas->Divide(1,3) ; 
		}
	
			
//Create The h2 poly graph 
   TH2Poly *H2Poly = new TH2Poly("SiLi","SiLi (something in title...)",-20,+20,-20,+20);
   //H2Poly->GetXaxis()->SetNdivisions(520);
   H2Poly->GetXaxis()->SetTitle("X");
   H2Poly->GetYaxis()->SetTitle("Y");
   H2Poly->SetContour(100);


//Add the bin map 
   TFile *f = TFile::Open("SiLiSegmentsGraph.root","READ");
   TMultiGraph *mg;
   TKey *key;
	TObject* obj;
   TIter nextkey(f->GetListOfKeys());
   while (key = (TKey*)nextkey()) {
      obj = key->ReadObj();
      if (obj->InheritsFrom("TGraph")) {
			mg = (TMultiGraph*)obj;
			if (mg->GetName()=="all") continue ; 
			H2Poly->AddBin(mg);
      	}
   }


   //gBenchmark->Start("Partitioning");
   H2Poly->ChangePartition(100, 100);
   //gBenchmark->Show("Partitioning");

//Fill the map (I nthe future fill it withthe number f content from a specefic range of energy, useful if it's calibrated')
   TRandom3 r;
   //gBenchmark->Start("Filling");
   for (int i=0; i<10000; i++) {
      float x = r.Uniform(-15,+15);
      float y  = r.Uniform(-15,+15);
      H2Poly->Fill(x,y);
   }
  //gBenchmark->Show("Filling");

   H2Poly->Draw("COL");	
	
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

		//if(pyold) gVirtualX->DrawLine(pxmin,pyold,pxmax,pyold);
		//gVirtualX->DrawLine(pxmin,y,pxmax,y);
		//gPad->SetUniqueID(py);
		//Float_t upy = gPad->AbsPixeltoY(py);
		//Float_t y = gPad->PadtoY(upy);

//detecotr parameters
		float InnerRadius = 5 ; 
		float OuterRadius = 15 ; 
		int   NumberOfRings = 10 ; 		
		//TLine *line = new TLine(0,0,x,y);
		//line->Draw();
		
		float theta = TMath::ATan(y/x) * TMath::RadToDeg(); // in Radians 
		if (x < 0 ) theta = theta + 180 ;  
		float angle1 = theta-(15);
		float angle2 = angle1+(30);

		TEllipse  *lsector = new TEllipse(0,0,OuterRadius,0,angle1,angle2);
		lsector->SetFillStyle(0);
		lsector->SetLineWidth(2);
   		lsector->Draw();
   		 
   		 float RingInnerRadius = TMath::Sqrt(x*x+y*y) - 0.5*(OuterRadius-InnerRadius)/NumberOfRings; 
   		 float RingOuterRadius = RingInnerRadius + (OuterRadius-InnerRadius)/NumberOfRings;  
		 TEllipse  *lring1 = new TEllipse(0.0,0.0,RingInnerRadius,RingInnerRadius);
		 TEllipse  *lring2 = new TEllipse(0.0,0.0,RingOuterRadius,RingOuterRadius);
		 lring1->SetFillStyle(0);
		 lring1->SetLineColor(1);
		 lring2->SetLineWidth(2);
   		 lring2->SetFillStyle(lring1->GetFillStyle());
   		 lring2->SetLineColor(lring1->GetLineColor());
   		 lring1->SetLineWidth(lring2->GetLineWidth());
   		 lring1->Draw();
   		 lring2->Draw();
  

		gCanvas->Update();
		//line->Delete();
		lring1->Delete();
		lring2->Delete();
		lsector->Delete();
		//return ;
	}


		if (event == 12) {
		cout<< "Middle click" <<endl ;

		//TH2Poly *h2poly = (TH2Poly*)select;
		cout<< "x and y :" <<x << "\t"<< y <<endl ;
		cout<< "bin #" <<h2poly->FindBin(x,y) <<endl ;
        int binxy = h2poly->FindBin(x,y) ;
		cout<< "bin content" <<h2poly->GetBinContent( binxy ) <<endl ;
		
		//Draw Histogram corresponding to mouse positionon the new Canvas	
		
		GetHistogram( binxy ) ; 

		//update canvas	
		//gCanvas->Update();					
		gCanvasHist->Update();
		//padsav->cd();
		}


		return ;

}


//_______________________________________
 void OpenDataFile() {
 
 // add command asking for the name/path...
 TString fname = "his30093"; 
 
 OpenDataThisFile(fname) ; 
 }
 
//_______________________________________ 
 int OpenDataThisFile(TString fname ) {
 
 	int success = 0 ;
 	
 	//TString fname = "his30093" ; 
       // open the ROOT file to process

   TString path ="/data1/moukaddam/SpiceTestSep2014/Calibration/Files/";
   TString extension =  ".root" ;
   TString inFileName = fname+extension;
   TFile *inFile = new TFile(path + inFileName);
   
   if ( ! inFile->IsOpen() )  { //try present directory
   path ="./";
   inFile->SetName(path + inFileName);
   }
   else {
   success = 1 ;
   cout << "file is opened "<< endl; 
    inFile->ls();   	
   } 

	   if (gDebug) { 
	   cout << "Opening the root file and grabing histograms " << inFile->GetName() << endl ;  
	   inFile->ls();
	   getchar();
	   }
	   
	 gFolderHistos = (TFolder*)(inFile->FindObjectAny("histos"));
	 return success ;
   }  
   
  
  //___________________________
void GetHistogram(int channel  ) {

//make sure file is opened 
	if (!gFolderHistos) {
		cout << " You need to open a histogram folder first " <<endl ; 
		return ; 
	}
	

	//create or set the new canvas gCanvasHist
		//TVirtualPad *padsav = gPad;
		gCanvasHist = (TCanvas*)gROOT->GetListOfCanvases()->FindObject("gCanvasHist");
		if(gCanvasHist) delete gCanvasHist->GetPrimitive("Projection"); // Delete the (primitive) subpad of the object named "Projection" 
		else   gCanvasHist = new TCanvas("gCanvasHist","Projection Canvas",710,10,700,700);
		//gCanvasHist->SetGrid();
		gCanvasHist->SetFrameFillColor(kBlack);
		gCanvasHist->SetCrosshair(2);
		gCanvasHist->cd();
		
       
	TString hname = Form("Chrg%d", channel+1060); 
	TH1F* hist = (TH1F*)(gFolderHistos->FindObjectAny(hname));
	cout << "=========================================";
	cout << " Retreiving Channel : " << channel ;
	cout << " ========================================="<< endl;

	//gCanvas->cd(1);
	hist->GetXaxis()->SetRangeUser(0,3500);
	hist->GetXaxis()->SetTitle("Energy (keV)");
	hist->GetXaxis()->SetAxisColor(kRed);
	hist->GetXaxis()->SetTickLength(-0.01);
	
	hist->GetYaxis()->SetTitle("Counts ");
	hist->GetYaxis()->SetAxisColor(kRed);
	hist->GetYaxis()->SetTickLength(-0.01);
	
	hist->SetLineColor(kYellow);
	
	//Draw on the other Canvas
	gCanvasHist->cd() ;
	hist->Draw(); 
	
	gCanvasHist->ToggleEventStatus();
	gCanvasHist->Update(); 
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
