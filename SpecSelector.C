//c++
#include <iostream>
#include <fstream>
#include <istream>
#include <sstream>
#include <vector>
#include <iomanip>      // setprecision
#include <map>
#include <stdlib.h>     // atof, atoi
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
void CreateSiLiDynamicMap();
void OpenDataFile();
int  OpenDataThisFile(TString filename);
void GetHistogram(int indicator);
void PQSelector(int select );
void ReadSiLiMap(TString CsvFileName);
#define STOP Stop();

//Global variables 
TCanvas *gCanvas ;
TCanvas *gCanvasHist ;
TFile *gInFile;	
TFolder *gFolderHistos ; 
TString gPhysQuantity="Chrg";
TH1F* gHist_current;
TH2Poly* gH2Poly ;


struct map_element_st {

	TString Mnemonic; // Mnemonic of the segment, refer to TigWiki  https://www.triumf.info/wiki/tigwiki/index.php/Detector_Nomenclature 
	int Sector; // 0-11 Number of sector as provided by Semikon
	int Ring; // 0-9 Number of ring as provided by Semikon
	TString SegIDSemikon; // Segment ID as provided by semikon map : SxRy (sector x, Ring y)  
	int SemikonPCBconnector; // Number 1-4 describing the position of the channel in one of the 4 groups on the PCB 
	TString FETBoardSlot; // Number 1-8 describing the position of the channel in one of the 8 slots on the vertical Copper plate 
	TString PreampAbsolutePosition; // Two letters describing the position of the premaps (Left/North or Right/South) and the assigned colour (Blue/Green/Red/White) : LR	
	int PreampPin; // Number 01-15 describing premap board pin number 	
	TString CloverCable; // Name of the Clover cable connected to the preamp e.g. Cassie	
	int CloverNumber;	// 01-16 Clover number from TIGRESS, corresponding to the physical position in space
	TString TIG10Cable ; // Tig10 cable (from pig-tail) colour<B,G,R,W>+number<01-16> e.g. R11 	
	int CollectorNumber	; // Number of Collector (holding several tig10 cards)
	int TIG10PortNumberOnCollector ; // Number 01-12 of port on the corresponding collector  
	TString TIG10PortNumberOnCollectorHEXA ; //	Hexa decimal Number 01-09,a,b,c of port on the corresponding collector  
	int TIG10PlugNumber; // Number 0-9 corresponding to the plug on the TIG10 card
	TString FSCP; // Word concatinating the FrontENd+Collector+Port+Plug number 
	int ChannelNumberFromOdb ; // a number of channel in the hit port pattern.
	float RadiusMidArea; // radius of the segment from the center 	
	float Theta; // Theta angle as seen from the target in Degrees	
	float Phi; // Phi angle as seen from the target	in Degrees
	int ReversedRingx12plusSector; // Number 000-119 (Reversed ring * 12) + sector, this formula is used to assign the mnemonic as described in the tigwiki 	
	int ReversedRing ; // Reversed ring	9-ring, used in constructing the ReversedRingx12plusSector parameter
	int Ringx12plusSector; // Number 000-119 (ring * 12) + sector	

};


//MAIN FUNCTION 
void SpecSelector() {

//Global settings
   gStyle->SetOptStat(0000);
   gStyle->SetPalette(1);   

//Open File for the histigram selection 
	OpenDataFile();
	
// Set the Mapping parameters, basically it will be an ascii file (in case the mapping is right during the experiment, this can be done through the mnemonics.	
	ReadSiLiMap("spicepreampmap_20141006.csv");
	
//Create the main map used by the user	
	CreateSiLiDynamicMap();

}

void ReadSiLiMap(TString CsvFileName ){
// This file is exported from a spread sheet of the SiLi Map.

	TString path ="/data1/moukaddam/SpiceTestSep2014/DetectorMapping/";

	 map<TString,map_element_st> gMap;
	 map<TString,map_element_st>::iterator it;

	int success = 1 ;
	int nb_columns = 23 ;
	string dummy_string="";
	string dummy_word[nb_columns]; // Number of columns 
	string dummy_value[nb_columns]; 

	ifstream csv_file ; 
	csv_file.open((path + CsvFileName).Data()) ; 

	if ( !csv_file.is_open() )  { //try present directory
	cout << "File doesn't exist in the directory  : " << path << endl ; 
	cout << "Trying the present working directory : ./" << endl ; 
	path ="./";
	csv_file.open((path + CsvFileName).Data()) ; 
	}
	if(csv_file.is_open()){
	cout << "File opened..." << path << endl ; 
	}
	else {
	cout << "File is not found.. EXIT!"<< endl;
	success = 0 ;
	exit(-1);
	} 

	//Read the titles, they will be used as indicators to fill the map 
	getline (csv_file,dummy_string);
	cout << dummy_string <<" "<< endl;
	istringstream ss(dummy_string);

	int anchor=-1;
	TString key_title = "SegIDSemikon" ; 
	
	for (int i=0 ; i < nb_columns ; i++){
		ss>>dummy_word[i];
		cout << " "<< dummy_word[i] <<endl;
		if( dummy_word[i] == key_title ){
			cout << " The map will be organised according to parameter "<< key_title <<endl;
			anchor = i ; 
			}
	}
	
	// in case the title is not found 
	if (anchor == -1 ) {
	key_title = dummy_word[0]; 
	cout << " Default : The map will be organised according to parameter "<< key_title <<endl;
	anchor == 0 ; 
	}

	
	while (true) {

			getline (csv_file,dummy_string);
			cout << " Reading : "<<dummy_string << endl;
	 		istringstream ss(dummy_string);
	 		
			if (csv_file.eof()) break;
	
			for (int i=0 ; i < nb_columns ; i++){
				ss>>dummy_value[i];
				cout<< " "<<dummy_value[i]<<endl ;
				if ( dummy_word[i] == "Mnemonic")  						gMap[dummy_value[anchor]].Mnemonic		        = dummy_value[i];
				if ( dummy_word[i] == "Sector")  						gMap[dummy_value[anchor]].Sector 		        = atoi( (dummy_value[i]).c_str() );
				if ( dummy_word[i] == "Ring")  							gMap[dummy_value[anchor]].Ring 			        = atoi( (dummy_value[i]).c_str() );
				if ( dummy_word[i] == "SegIDSemikon")  						gMap[dummy_value[anchor]].SegIDSemikon	 	        = dummy_value[i];
				if ( dummy_word[i] == "SemikonPCBconnector")  					gMap[dummy_value[anchor]].SemikonPCBconnector	        = atoi( (dummy_value[i]).c_str() );
				if ( dummy_word[i] == "FETBoardSlot")  						gMap[dummy_value[anchor]].FETBoardSlot	     	        = dummy_value[i];
				if ( dummy_word[i] == "PreampAbsolutePosition")  				gMap[dummy_value[anchor]].PreampAbsolutePosition         = dummy_value[i];
				if ( dummy_word[i] == "PreampPin")  						gMap[dummy_value[anchor]].PreampPin		        = atoi( (dummy_value[i]).c_str() );
				if ( dummy_word[i] == "CloverCable")  						gMap[dummy_value[anchor]].CloverCable		        = dummy_value[i];
				if ( dummy_word[i] == "CloverNumber")  						gMap[dummy_value[anchor]].CloverNumber		        = atoi( (dummy_value[i]).c_str() );
				if ( dummy_word[i] == "TIG10Cable")  						gMap[dummy_value[anchor]].TIG10Cable		        = dummy_value[i];
				if ( dummy_word[i] == "CollectorNumber")  					gMap[dummy_value[anchor]].CollectorNumber	        = atoi( (dummy_value[i]).c_str() );
				if ( dummy_word[i] == "TIG10PortNumberOnCollector")  				gMap[dummy_value[anchor]].TIG10PortNumberOnCollector     = atoi( (dummy_value[i]).c_str() );
				if ( dummy_word[i] == "TIG10PortNumberOnCollectorHEXA")  			gMap[dummy_value[anchor]].TIG10PortNumberOnCollectorHEXA = dummy_value[i];
				if ( dummy_word[i] == "TIG10PlugNumber")  					gMap[dummy_value[anchor]].TIG10PlugNumber	        = atoi( (dummy_value[i]).c_str() );
				if ( dummy_word[i] == "FSCP") 				 			gMap[dummy_value[anchor]].FSCP			        = dummy_value[i];
				if ( dummy_word[i] == "ChannelNumberFromOdb")  					gMap[dummy_value[anchor]].ChannelNumberFromOdb	        = atoi( (dummy_value[i]).c_str() );
				if ( dummy_word[i] == "RadiusMidArea")  					gMap[dummy_value[anchor]].RadiusMidArea		        = atof( (dummy_value[i]).c_str() );
				if ( dummy_word[i] == "Theta")  						gMap[dummy_value[anchor]].Theta			        = atof( (dummy_value[i]).c_str() );
				if ( dummy_word[i] == "Phi")  							gMap[dummy_value[anchor]].Phi			        = atof( (dummy_value[i]).c_str() );
				if ( dummy_word[i] == "ReversedRingx12plusSector")  				gMap[dummy_value[anchor]].ReversedRingx12plusSector      = atoi( (dummy_value[i]).c_str() );																      
				if ( dummy_word[i] == "ReversedRing")  						gMap[dummy_value[anchor]].ReversedRing		        = atoi( (dummy_value[i]).c_str() );
				if ( dummy_word[i] == "Ringx12plusSector")  					gMap[dummy_value[anchor]].Ringx12plusSector	        = atoi( (dummy_value[i]).c_str() );     
			}
		
	}

	
	if (1) { 
	  cout << " Dump the content of the map " << endl; 
	  getchar();
	  
	  for (it=gMap.begin(); it!=gMap.end(); it++) {
	  
 
	      cout << it->first << " => " << it->second.Mnemonic << '\n';
	      cout << it->first << " => " << it->second.Sector << '\n';
	      cout << it->first << " => " << it->second.Ring << '\n';
	      cout << it->first << " => " << it->second.SegIDSemikon << '\n';	      
	      cout << it->first << " => " << it->second.SemikonPCBconnector << '\n';
	      cout << it->first << " => " << it->second.FETBoardSlot << '\n';	       
	      cout << it->first << " => " << it->second.PreampAbsolutePosition << '\n';
	      cout << it->first << " => " << it->second.PreampPin << '\n';
	      cout << it->first << " => " << it->second.CloverCable << '\n';
	      cout << it->first << " => " << it->second.CloverNumber << '\n';
	      cout << it->first << " => " << it->second.TIG10Cable << '\n';
	      cout << it->first << " => " << it->second.CollectorNumber << '\n';
	      cout << it->first << " => " << it->second.TIG10PortNumberOnCollector << '\n';
	      cout << it->first << " => " << it->second.TIG10PortNumberOnCollectorHEXA << '\n';
	      cout << it->first << " => " << it->second.TIG10PlugNumber << '\n';
	      cout << it->first << " => " << it->second.FSCP << '\n';
	      cout << it->first << " => " << it->second.ChannelNumberFromOdb << '\n';
	      cout << it->first << " => " << it->second.RadiusMidArea << '\n';
	      cout << it->first << " => " << it->second.Theta << '\n';
	      cout << it->first << " => " << it->second.Phi << '\n';
	      cout << it->first << " => " << it->second.ReversedRingx12plusSector << '\n';
	      cout << it->first << " => " << it->second.ReversedRing << '\n';
	      cout << it->first << " => " << it->second.Ringx12plusSector << '\n';
          STOP
	  }

	
	}

}


void CreateSiLiDynamicMap(){


	if(!gCanvas) {
		TString title = "Bi207" ;
		gCanvas = new TCanvas("SiLi_Segments","SiLi_Segments",10,10,700,700);
		gCanvas->ToggleEventStatus();
		gCanvas->AddExec("ex","SegmentClicked()");
		}
	
//Create The h2 poly graph 
	gH2Poly = new TH2Poly("SiLi","SiLi (something in title...)",-20,+20,-20,+20);
   //gH2Poly->GetXaxis()->SetNdivisions(520);
   gH2Poly->GetXaxis()->SetTitle("X");
   gH2Poly->GetYaxis()->SetTitle("Y");
   gH2Poly->SetContour(100);


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
			gH2Poly->AddBin(mg);
	  	}
	}


   gH2Poly->ChangePartition(100, 100);

//Fill the map (I nthe future fill it withthe number f content from a specefic range of energy, useful if it's calibrated')
   TRandom3 r;
   //gBenchmark->Start("Filling");
   for (int i=0; i<10000; i++) {
      float x = r.Uniform(-15,+15);
      float y  = r.Uniform(-15,+15);
      gH2Poly->Fill(x,y);
   }
  //gBenchmark->Show("Filling");

   gH2Poly->Draw("LF COL");	
	
}

void SegmentClicked() {
	
	TObject *select = gPad->GetSelected();
	if(!select) return;
	if (!select->InheritsFrom(TH2Poly::Class())) {gPad->SetUniqueID(0); return;}
	
	gH2Poly = (TH2Poly*)select;
	gPad->GetCanvas()->FeedbackMode(kTRUE);

	//Get the event
	int event = gPad->GetEvent();

	//Get the abscissa  
	int px = gPad->GetEventX();
	int py = gPad->GetEventY();

	float x = gPad->AbsPixeltoX(px);
	float y = gPad->AbsPixeltoY(py);


	if (event == 51) 
		{
//Detector parameters
		float InnerRadius = 5 ; 
		float OuterRadius = 15 ; 
		int   NumberOfRings = 10 ; 		

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
		lring1->Delete();
		lring2->Delete();
		lsector->Delete();
		//return ;
	}

	//this action function is called whenever you middle click on the mouse
		if (event == 12) {
		 // Get the sector and the ring 
		//GetRing(x,y) ;  
		//GetSector(x,y) ;
		//it prints the informations of the segments 
		cout<< "x and y :" <<x << "\t"<< y <<endl ;
		cout<< "bin #" <<gH2Poly->FindBin(x,y) <<endl ;
        int binxy = gH2Poly->FindBin(x,y) ;
		cout<< "bin content" <<gH2Poly->GetBinContent( binxy ) <<endl ;
		
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
		/*gCanvasHist = (TCanvas*)gROOT->GetListOfCanvases()->FindObject("gCanvasHist");
		if(gCanvasHist) delete gCanvasHist->GetPrimitive("Projection"); // Delete the (primitive) subpad of the object named "Projection" 
		else   gCanvasHist = new TCanvas("gCanvasHist","Projection Canvas",710,10,700,700);
		//gCanvasHist->SetGrid();
		gCanvasHist->SetFrameFillColor(kBlack);
		gCanvasHist->SetCrosshair(2);
		gCanvasHist->cd();
		gCanvasHist->ToggleEventStatus();
		*/
		
		if(!gCanvasHist) {
		//gCanvasHist->Delete();
		TString title = "Bi207" ;
		gCanvasHist = new TCanvas("gCanvasHist","Selected",710,10,700,700);
		gCanvasHist->SetCrosshair(2);
		gCanvasHist->SetFrameFillColor(kBlack);
		gCanvasHist->SetFillColor(kBlack);
		gCanvasHist->ToggleEventStatus();
		//gCanvasHist->AddExec("ex","RangeClicked()");
		//gCanvas->Divide(1,3) ; 
	}
	
       
	TString hname = Form(gPhysQuantity+"%d", channel+1060); 
	gHist_current = (TH1F*)(gFolderHistos->FindObjectAny(hname));
	cout << "=========================================";
	cout << " Retreiving Channel : " << channel ;
	cout << " ========================================="<< endl;

	//gCanvas->cd(1);
	
	gHist_current->GetXaxis()->SetRangeUser(0,3500);
	gHist_current->GetXaxis()->SetTitle("Energy (keV)");
	gHist_current->GetXaxis()->SetTickLength(-0.01);
	gHist_current->GetXaxis()->SetTitleColor(kRed);
	gHist_current->GetXaxis()->SetAxisColor(kRed);
	gHist_current->GetXaxis()->SetLabelColor(kRed);
		
	gHist_current->GetYaxis()->SetTitle("Counts ");
	gHist_current->GetYaxis()->SetTickLength(-0.01);
	gHist_current->GetYaxis()->SetAxisColor(kRed);
	gHist_current->GetYaxis()->SetTitleColor(kRed);
	gHist_current->GetYaxis()->SetLabelColor(kRed);
	
	gHist_current->SetLineColor(kYellow);
	
	//Draw on the other Canvas
	gCanvasHist->cd() ;
	gHist_current->Draw(); 
	
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
	gFolderHistos = (TFolder*)(gInFile->FindObjectAny("histosCfd"));
	return ;
	}
	
if ( select == 3 ) { 
	gPhysQuantity="Led"; 
	gFolderHistos = (TFolder*)(gInFile->FindObjectAny("histosLed"));
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
