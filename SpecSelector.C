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


//Global variables 

//Detector parameters
	int rings_number = 10 ; //10
	float rad_start = 8 ;
	float rad_end = 47 ; 
	float rad_pitch = (rad_end-rad_start)/rings_number ;

	int sectors_number = 12 ; //12      
	float phi_start = 0*TMath::DegToRad(); ;
	float phi_end = 360*TMath::DegToRad(); ; 
	float phi_pitch = ((phi_end-phi_start)/sectors_number); // 30 degrees
	
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
	int SemikonPCBpin; // Number 1-30 describing the position of the pin on the PCB, useful to investigate cross talk 
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
map<TString,map_element_st> gMap;

// functions 
void Stop();
void SegmentClicked();
void CreateSiLiDynamicMap();
void OpenDataFile();
int  OpenDataThisFile(TString filename);
void GetHistogram(int indicator);
void PQSelector(int select );
void ReadSiLiMap(TString CsvFileName);
void PrintInfo(const char* SegmentID) ; 
TGraph* RotateGraph(TGraph* g, float angle ) ;
void MakeSiLiGraphs(); 
void SetColor(Int_t red, Int_t green, Int_t blue); 

//Definitions 
#define DEBUG 0
#define STOP Stop();
// NODES should always be of form n+2+1 where:
// n is the number of jumps from point(node) to point(node)m this number is always even, 
// n+2 is the number of nodes
// n+2+1 the additional point correspond to the final point superimposed on the starting point(node) 
#define NODES 8+3   



//MAIN FUNCTION 
void SpecSelector() {
//Global settings
   gStyle->SetOptStat(0000);
   gStyle->SetPalette(1);   

// Set the Mapping parameters, basically it will be an ascii file (in case the mapping is right during the experiment, this can be done through the mnemonics.	
	ReadSiLiMap("spicepreampmap_20141009.csv");
	
//Make the graphs of the SiLi 
	//MakeSiLiGraphs(); 

//Open File for the histigram selection 
	OpenDataFile();

//Create the main map used by the user	
	CreateSiLiDynamicMap();
}



void ReadSiLiMap(TString CsvFileName ){
// This file is exported from a spread sheet of the SiLi Map.

	TString path ="/data1/moukaddam/SpiceTestSep2014/DetectorMapping/";
	 map<TString,map_element_st>::iterator it;

	int success = 1 ;
	int nb_columns = 30 ;
	string buffer_string="";
	string buffer_word[nb_columns]; // Number of columns 
	string buffer_value[nb_columns]; 

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
	getline (csv_file,buffer_string);
	cout << buffer_string <<" "<< endl;
	istringstream ss(buffer_string);

	int anchor=-1;
	TString key_title = "SegIDSemikon" ; 
	
	for (int i=0 ; i < nb_columns ; i++){
		ss>>buffer_word[i];
		cout << " "<< buffer_word[i] ;
		if( buffer_word[i] == key_title ){
			cout << " <<<<<< The map will be organised according to this parameter ";
			anchor = i ; 
			}
			cout << endl ; 
		}
	
	// in case the title is not found 
	if (anchor == -1 ) {
	key_title = buffer_word[0]; 
	cout << " Default : The map will be organised according to parameter "<< key_title <<endl;
	anchor == 0 ; 
	}
	
	while (true) {

            getline (csv_file,buffer_string);
			//cout << "\nReading : "<<buffer_string << endl;
	 		istringstream ss(buffer_string);
	 		
	 		//if end of file break 
			if (csv_file.eof()) break;
           
           //read all the values 	
			for (int i=0 ; i < nb_columns ; i++){
				ss>>buffer_value[i];
				}
			//Fill all the values in the right slot 	
			for (int i=0 ; i < nb_columns ; i++){
				buffer_value[i];
				//cout<< " * "<< buffer_value[anchor] << "/"<< buffer_value[i]<<" " ;
				if ( buffer_word[i] == "Mnemonic")  						gMap[buffer_value[anchor]].Mnemonic		        			= buffer_value[i];
				if ( buffer_word[i] == "Sector")  						gMap[buffer_value[anchor]].Sector 		       		 		= atoi( (buffer_value[i]).c_str() );
				if ( buffer_word[i] == "Ring")  							gMap[buffer_value[anchor]].Ring 			        			= atoi( (buffer_value[i]).c_str() );
				if ( buffer_word[i] == "SegIDSemikon")  					gMap[buffer_value[anchor]].SegIDSemikon	 	       		 	= buffer_value[i];
				if ( buffer_word[i] == "SemikonPCBconnector")  			gMap[buffer_value[anchor]].SemikonPCBconnector	        	= atoi( (buffer_value[i]).c_str() );
				if ( buffer_word[i] == "SemikonPCBpin")  			gMap[buffer_value[anchor]].SemikonPCBpin	        	= atoi( (buffer_value[i]).c_str() );
				if ( buffer_word[i] == "FETBoardSlot")  					gMap[buffer_value[anchor]].FETBoardSlot	     	        	= buffer_value[i];
				if ( buffer_word[i] == "PreampAbsolutePosition")  		gMap[buffer_value[anchor]].PreampAbsolutePosition         	= buffer_value[i];
				if ( buffer_word[i] == "PreampPin")  					gMap[buffer_value[anchor]].PreampPin		        			= atoi( (buffer_value[i]).c_str() );
				if ( buffer_word[i] == "CloverCable")  					gMap[buffer_value[anchor]].CloverCable		        		= buffer_value[i];
				if ( buffer_word[i] == "CloverNumber")  					gMap[buffer_value[anchor]].CloverNumber		        		= atoi( (buffer_value[i]).c_str() );
				if ( buffer_word[i] == "TIG10Cable")  					gMap[buffer_value[anchor]].TIG10Cable		        		= buffer_value[i];
				if ( buffer_word[i] == "CollectorNumber")  				gMap[buffer_value[anchor]].CollectorNumber	        		= atoi( (buffer_value[i]).c_str() );
				if ( buffer_word[i] == "TIG10PortNumberOnCollector")  	gMap[buffer_value[anchor]].TIG10PortNumberOnCollector     	= atoi( (buffer_value[i]).c_str() );
				if ( buffer_word[i] == "TIG10PortNumberOnCollectorHEXA") gMap[buffer_value[anchor]].TIG10PortNumberOnCollectorHEXA 	= buffer_value[i];
				if ( buffer_word[i] == "TIG10PlugNumber")  				gMap[buffer_value[anchor]].TIG10PlugNumber	        		= atoi( (buffer_value[i]).c_str() );
				if ( buffer_word[i] == "FSCP") 				 			gMap[buffer_value[anchor]].FSCP			        			= buffer_value[i];
				if ( buffer_word[i] == "ChannelNumberFromOdb")  			gMap[buffer_value[anchor]].ChannelNumberFromOdb	        	= atoi( (buffer_value[i]).c_str() );
				if ( buffer_word[i] == "RadiusMidArea")  				gMap[buffer_value[anchor]].RadiusMidArea		        		= atof( (buffer_value[i]).c_str() );
				if ( buffer_word[i] == "Theta")  						gMap[buffer_value[anchor]].Theta			        			= atof( (buffer_value[i]).c_str() );
				if ( buffer_word[i] == "Phi")  							gMap[buffer_value[anchor]].Phi			        			= atof( (buffer_value[i]).c_str() );
				if ( buffer_word[i] == "ReversedRingx12plusSector")  	gMap[buffer_value[anchor]].ReversedRingx12plusSector      	= atoi( (buffer_value[i]).c_str() );																      
				if ( buffer_word[i] == "ReversedRing")  					gMap[buffer_value[anchor]].ReversedRing		        		= atoi( (buffer_value[i]).c_str() );
				if ( buffer_word[i] == "Ringx12plusSector")  			gMap[buffer_value[anchor]].Ringx12plusSector	        		= atoi( (buffer_value[i]).c_str() );     
			//STOP
			}
		
	}

	
	if (DEBUG) { 
	  cout << " Dump the content of the map " << endl; 
	  getchar();
	  
	  for (it=gMap.begin(); it!=gMap.end(); it++) {
	      cout << it->first << " => " << it->second.Mnemonic << '\n';
	      cout << it->first << " => " << it->second.Sector << '\n';
	      cout << it->first << " => " << it->second.Ring << '\n';
	      cout << it->first << " => " << it->second.SegIDSemikon << '\n';	      
	      cout << it->first << " => " << it->second.SemikonPCBconnector << '\n';
	      cout << it->first << " => " << it->second.SemikonPCBpin << '\n';
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
	  }

	
	}

}


void CreateSiLiDynamicMap(){


	if(!gCanvas) {
		TString title = "Bi207" ;
		gCanvas = new TCanvas("SiLiDynamicMap","SiLiDynamicMap",10,10,700,700);
		gCanvas->ToggleEventStatus();
		gCanvas->SetLogz();
		gCanvas->AddExec("ex","SegmentClicked()");
		}
	
//Create The h2 poly graph 
	gH2Poly = new TH2Poly("SiLi","SPICE Si(Li) << Looking Up-Stream >> ",-50,+50,-50,+50);
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
      float x = r.Uniform(-50,+50);
      float y  = r.Uniform(-50,+50);
      gH2Poly->Fill(x,y);
   }
  //gBenchmark->Show("Filling");

   gH2Poly->Draw("LF COLZ");	
	
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
		float theta = TMath::ATan(y/x) * TMath::RadToDeg(); // from Radians to degrees
		if (x > 0 && y > 0 ) theta = theta ; // x<0 ; y<0 	
		if (x < 0 && y > 0 ) theta = theta+180 ; // x<0 ; y<0 
		if (x > 0 && y < 0 ) theta = theta-30 ; // x>0 ; y<0 
		if (x < 0 && y < 0 ) theta = theta+180 ; // x<0 ; y<0 

		int angle1 = (int)theta - (int)theta % 30; 
		int angle2 = angle1 + 30;

    //Get the color of premap 
	int binxy = gH2Poly->FindBin(x,y) ;
	TString PreampPosition = gMap[gH2Poly->GetBinTitle(binxy)].PreampAbsolutePosition;
    int color = kWhite ;  
	if (PreampPosition == "LR" || PreampPosition == "RR"  ) color = kRed ;
		else if (PreampPosition == "LW" || PreampPosition == "RW"  ) color = kBlack ; 
			else if (PreampPosition == "LG" || PreampPosition == "RG"  ) color = kGreen+2 ;
				else if (PreampPosition == "LB" || PreampPosition == "RB"  ) color = kBlue ;
   
		TEllipse  *lsector = new TEllipse(0,0,rad_end,0,angle1,angle2);
		lsector->SetFillStyle(0);
		lsector->SetLineWidth(4);
		lsector->SetLineColor(color);
		lsector->Draw();

        int radius = TMath::Sqrt(x*x+y*y); 
		float RingInnerRadius = (TMath::Floor((radius-rad_start)/rad_pitch) * rad_pitch) + rad_start;
		float RingOuterRadius = RingInnerRadius + rad_pitch;  
		TEllipse  *lring1 = new TEllipse(0.0,0.0,RingInnerRadius,RingInnerRadius);
		TEllipse  *lring2 = new TEllipse(0.0,0.0,RingOuterRadius,RingOuterRadius);
		lring1->SetFillStyle(0);
		lring1->SetLineColor(color);
		lring2->SetLineWidth(4);
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
		cout<< "x and y :" <<x << "\t"<< y <<endl ;
		cout<< "bin #" <<gH2Poly->FindBin(x,y) <<endl ;
        int binxy = gH2Poly->FindBin(x,y) ;
		//cout<< "bin content" <<gH2Poly->GetBinContent( binxy ) <<endl ;
		
		//Draw Histogram corresponding to mouse position the new Canvas
		int channel = gMap[gH2Poly->GetBinTitle(binxy)].ChannelNumberFromOdb ; 
		GetHistogram( channel ) ; 
		
		//Print the info of the selected segment on the screen
		PrintInfo( gH2Poly->GetBinTitle(binxy) ) ; 


		//update canvas				
		gCanvasHist->Update();
		//padsav->cd();
		}

		return ;

}


void  PrintInfo(const char* segment){
 
	cout << "===================================================================================="<< endl; 	     
	cout << "Semikon Segment ID :  " << gMap[segment].SegIDSemikon   << endl; 
	cout << "------------------------------------------------------------------ F E T - B O A R D"<< endl; 	     
	cout << "SemikonPCBconnector\t\t" << 		  gMap[segment].SemikonPCBconnector << endl;
    cout << "SemikonPCBpin\t\t\t" << 		  	  gMap[segment].SemikonPCBpin << endl;	  	  
	cout << "FETBoardSlot\t\t\t" <<   			  gMap[segment].FETBoardSlot << endl; 
	cout << "---------------------------------------------------------------------- P R E A M P S"<< endl; 	     
	cout << "PreampAbsolutePosition\t\t" << 	  gMap[segment].PreampAbsolutePosition	  << endl; 
	cout << "PreampPin\t\t\t" << 				  gMap[segment].PreampPin  	  << endl; 
	cout << "--------------------------------------------------------- T I G R E S S - C A B L E S"<< endl; 	   
	cout << "CloverCable\t\t\t" << 				  gMap[segment].CloverCable	<< endl;   
	cout << "CloverNumber\t\t\t" << 			  gMap[segment].CloverNumber << endl;
	cout << "TIG10Cable\t\t\t" << 				  gMap[segment].TIG10Cable 	  << endl; 
	cout << "----------------------------------------------------- S P A C I A L - P O S I T I O N"<< endl; 
	cout << "Sector\t\t\t\t" << 				  gMap[segment].Sector	<< endl;   
	cout << "Ring\t\t\t\t" << 					  gMap[segment].Ring << endl; 
	cout << "Theta\t\t\t\t" <<  				  gMap[segment].Theta	  << endl; 
	cout << "Phi\t\t\t\t" << 					  gMap[segment].Phi << endl;
	cout << "RadiusMidArea\t\t\t" <<  			  gMap[segment].RadiusMidArea	  << endl; 
	cout << "------------------------------------------------------------------- T I G 1 0 - D A Q"<< endl; 
	cout << "CollectorNumber\t\t\t" << 			  gMap[segment].CollectorNumber	  << endl; 
	//cout << "TIG10PortNumberOnCollector\t"  <<  gMap[segment].TIG10PortNumberOnCollector    << endl; 
	cout << "TIG10PortNumberOnCollectorHEXA\t"<<  gMap[segment].TIG10PortNumberOnCollectorHEXA << endl; 	  
	cout << "TIG10PlugNumber\t\t\t" << 			  gMap[segment].TIG10PlugNumber	  << endl;
	cout << "ChannelNumberFromOdb\t\t" << 		  gMap[segment].ChannelNumberFromOdb << endl;
	cout << "FSCP\t\t\t\t" << 					  gMap[segment].FSCP	<< endl;
	cout << "Mnemonic\t\t\t" << 				  gMap[segment].Mnemonic	 << endl;
	cout << "====================================================================================="<< endl; 	  
	//cout << "ReversedRingx12plusSector\t\t\t" << 	 gMap[segment].ReversedRingx12plusSector	<< endl; 	  
	//cout << "ReversedRing\t\t\t" << 				 gMap[segment].ReversedRing << endl; 
	//cout << "Ringx12plusSector\t\t\t" << 			 gMap[segment].Ringx12plusSector  	<< endl;   
 
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
void GetHistogram(int channel) {

//make sure file is opened 
	if (!gFolderHistos) {
		cout << " You need to open a histogram folder first " <<endl ; 
		return ; 
	}
	

	//create or set the new canvas gCanvasHist		
		if(!gCanvasHist) {
		TString title = "Bi207" ;
		gCanvasHist = new TCanvas("gCanvasHist","Selected",710,10,700,700);
		gCanvasHist->SetCrosshair(2);
		gCanvasHist->SetFrameFillColor(kBlack);
		gCanvasHist->SetFillColor(kBlack);
		gCanvasHist->ToggleEventStatus();
		//gCanvasHist->AddExec("ex","RangeClicked()");
		//gCanvas->Divide(1,3) ; 
	}
	
       
	TString hname = Form(gPhysQuantity+"%d", channel); 
	gHist_current = (TH1F*)(gFolderHistos->FindObjectAny(hname));
	cout << "==========================================================Retreiving Channel : " << channel << endl ;     
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



void MakeSiLiGraphs() {
   
	TFile* file = new TFile("SiLiSegmentsGraph.root","RECREATE");
	TCanvas *c1 = new TCanvas("SiLiSegments","SiLiSegments",10,10,700,700);	
	TMultiGraph *mg = new TMultiGraph("all","all");
	TGraph *gr[10];
				
	Double_t x[NODES]; 
	Double_t y[NODES];

	float phi_mini_pitch = phi_pitch/((NODES-3)/2); 
	float rad_min = 0 ;
	float rad_max = 0 ;
	float phi_min = 0 ;
	float phi_max = 0 ;
	float rad_current = 0 ;
	float phi_current = 0 ;
	  	       	   	
	for (Int_t i = 0 ; i < rings_number ; i++) { // loop on rings 
		// rad_min and rad max are the borders of the current segment
		rad_min = rad_start + i*rad_pitch ; 
		rad_max = rad_start + (i+1)*rad_pitch ;

		// phi min and phi max are the borders of the current segment	
		phi_min = phi_start ;
		phi_max = phi_min + phi_pitch ;
		
		// phi and rad are now defined for one segment
		phi_current = phi_min ; // Initiate the current phi as the minimum phi on the current segment
		rad_current = rad_min ; // Initiate the current radius as the inner radius on the current segment
		int sign = +1 ; // Initiate the current radius as the minimum radius on the current segment
		
		for (Int_t k = 0; k < NODES-1 ; k++) { // loop on the nodes, first for inner border with sign = +1
		 
			x[k] =rad_current*TMath::Cos( phi_current ) ; 
			y[k] =rad_current*TMath::Sin( phi_current ) ;
		
			if (k==0) {//Assign the last point
	 	 		x[NODES-1] = x[k]; 
	 	 		y[NODES-1] = y[k];
				}
				
			 if ( phi_current == phi_max && rad_current == rad_min) { // Jump up to radius max but keep the same phi ()
				rad_current = rad_max ; 
				sign = 0;
				}
			else if ( rad_current == rad_max && phi_current > phi_min) {
				sign = -1;
				}
		
			phi_current = phi_current + sign*phi_mini_pitch;  // loop within segment on outer border, counter clock wise
 	 	
 	 		}// finish the loop with the segment 
			
		 gr[i]= new TGraph(NODES,x,y);
					
		}
		
	file->cd() ; 
	
	for (Int_t i = 0; i < rings_number ; i++) { // loop on rings		
		for (Int_t j = 0; j < sectors_number ; j++) { // loop on sectors
			TString name = Form("S%dR%d",j,i) ;
			gr[i]->SetNameTitle(name,name);
	 		TGraph* g = RotateGraph(gr[i],j*phi_pitch) ; 
			g->SetNameTitle(name,name); 
			g->Write();  
			mg->Add(g,"lp"); 
		}						
	}

	c1->cd();
	mg->Draw("a");    
    mg->Write(); 
	file->Write(); 
}

// Rotate the graphs used to define the bins in the TH2poly graph
TGraph* RotateGraph(TGraph* g, float angle) { // radian

	Double_t* x = g->GetX(); 
	Double_t* y = g->GetY();
	Int_t n = g->GetN();
		 	 	   
	//new Values 
	Double_t* xx = new Double_t[n]; 
	Double_t* yy = new Double_t[n];

	for (int i=0 ; i< n ; i++){
		TVector2 v1(x[i],y[i]);
		TVector2 v2 = v1.Rotate(angle); //radian 
		xx[i]=v2.X(); 
		yy[i]=v2.Y();	
	}
 
	TGraph* gg = new TGraph(n,xx,yy) ; 
	gg->SetLineColor(kWhite); 
	gg->SetLineWidth(0);
	gg->SetMarkerColor(kWhite);
	gg->SetMarkerStyle(1);

	return gg ; 
}



void Stop() {
	char c; 
	c = getchar(); 
	if (c=='q') exit(-1) ;
	}


void SetColor(Int_t red, Int_t green, Int_t blue ){
   
   //Default
   if (red==0 && green==0 && blue==0 )  {
	   gStyle->SetPalette(1);
	   gCanvas->Modified();
	   gCanvas->Update();
	   return; 
	   }
   
    const Int_t NumberOfStops = 5;
    const Int_t NCont = 255;

    Double_t StopsPositions[NumberOfStops] = { 0.00, 0.25, 0.50, 0.75, 1.00 };
    Double_t Red[NumberOfStops]   = { (red*0.05)/255, (red*0.25)/255, (red*0.50)/255, (red*0.75)/255, (red*1.00)/255};
    Double_t Green[NumberOfStops] = { (green*0.05)/255, (green*0.25)/255, (green*0.50)/255, (green*0.75)/255, (green*1.00)/255 };
    Double_t Blue[NumberOfStops]  = { (blue*0.05)/255, (blue*0.25)/255, (blue*0.50)/255, (blue*0.75)/255, (green*1.00)/255 };
    
    Int_t FI = TColor::CreateGradientColorTable(NumberOfStops, StopsPositions, Red, Green, Blue, NCont);
    gStyle->SetNumberContours(NCont);
   
   gCanvas->Modified();
   gCanvas->Update();
   return; 
   
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
