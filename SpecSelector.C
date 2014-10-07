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
void PQSelector(int select );
#define STOP Stop();

//Global variables 
TCanvas *gCanvas ;
TCanvas *gCanvasHist ;
TFile *gInFile;	
TFolder *gFolderHistos ; 
TString gPhysQuantity="Chrg";
TH1F* gHist;


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
		
		//Draw Histogram corresponding to mouse position the new Canvas			
		GetHistogram( binxy ) ; 

		//update canvas				
		gCanvasHist->Update();
		//padsav->cd();
		}


		return ;

}


//_______________________________________
 void OpenDataFile() {
 
 // add command asking for the name/path...
 TString fname = "hisXXXXX"; 
 
 OpenDataThisFile(fname) ; 
 }
 
//_______________________________________ 
 int OpenDataThisFile(TString fname ) {
 
 	int success = 0 ;
 	
 	//TString fname = "his30093" ; 
       // open the ROOT file to process

   TString path ="./";
   TString extension =  ".root" ;
   TString inFileName = fname+extension;
   gInFile = new TFile(path + inFileName);
   
   if ( ! gInFile->IsOpen() )  { //try present directory
   path ="./";
   gInFile->SetName(path + inFileName);
   }
   else {
   success = 1 ;
   cout << "file is opened "<< endl; 
    gInFile->ls();   	
   } 

	   if (gDebug) { 
	   cout << "Opening the root file and grabing histograms " << gInFile->GetName() << endl ;  
	   gInFile->ls();
	   getchar();
	   }
	   
	 gFolderHistos = (TFolder*)(gInFile->FindObjectAny("histos"));
 	 
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
		
       
	TString hname = Form(gPhysQuantity+"%d", channel+1060); 
	gHist = (TH1F*)(gFolderHistos->FindObjectAny(hname));
	cout << "=========================================";
	cout << " Retreiving Channel : " << channel ;
	cout << " ========================================="<< endl;

	//gCanvas->cd(1);
	gHist->GetXaxis()->SetRangeUser(0,3500);
	gHist->GetXaxis()->SetTitle("Energy (keV)");
	gHist->GetXaxis()->SetAxisColor(kRed);
	gHist->GetXaxis()->SetTickLength(-0.01);
	
	gHist->GetYaxis()->SetTitle("Counts ");
	gHist->GetYaxis()->SetAxisColor(kRed);
	gHist->GetYaxis()->SetTickLength(-0.01);
	
	gHist->SetLineColor(kYellow);
	
	//Draw on the other Canvas
	gCanvasHist->cd() ;
	gHist->Draw(); 
	
	gCanvasHist->ToggleEventStatus();
	gCanvasHist->Update(); 
}
 

void PQSelector(int select ) {

if ( select == 0 ) {
	gPhysQuantity="Chrg"; 
	gFolderHistos = (TFolder*)(gInFile->FindObjectAny("histos"));
	return ;
	}
	
if ( select == 1 ) { 
	gPhysQuantity="ChrgCal"; 
	gFolderHistos = (TFolder*)(gInFile->FindObjectAny("histosCal"));
	return ;
	}

if ( select == 2 ) { 
	gPhysQuantity="Cfd"; 
	gFolderHistos = (TFolder*)(gInFile->FindObjectAny("histosLed"));
	return ;
	}
	
if ( select == 3 ) { 
	gPhysQuantity="Led"; 
	gFolderHistos = (TFolder*)(gInFile->FindObjectAny("histosCfd"));
	return ;
	}
	
if ( select == 4 ) { 
	gPhysQuantity="Time"; 
	gFolderHistos = (TFolder*)(gInFile->FindObjectAny("histosTime"));
	return ;
	}
	
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
