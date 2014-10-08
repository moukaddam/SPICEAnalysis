/*
Description 
...

Last Update :   

How to execute :
...
*/


//c++
#include <iostream>
#include <fstream>
#include <istream>
#include <vector>
using namespace std ; 

// Root
#include "TROOT.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TMath.h"
#include "TFile.h"
#include "TString.h"
#include "TTree.h"
#include "TBranch.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TColor.h"
#include "TSpectrum.h"
#include "TGraphErrors.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TF1.h"
#include "TDirectory.h"
#include "TFolder.h"
#include "TControlBar.h"
#include "TInterpreter.h"


//___________________ Global Variables
TFolder *gFolderHistos ; 
Int_t    gDebug; 
Int_t    gChannel_offset = 1060 ; // depends on odb	
Int_t    gChannel_current = gChannel_offset-1 ; // depends on odb	
TH1F 	*gHist_current ; // histogram (original) 
Int_t 	 gBining ; // Number of times the histogram was binned, pow(2,gBinning) is the number stored 
TCanvas *gCanvas ; 
Float_t  gMin = -1  ; 
Float_t  gMax = -1  ;
ofstream  gPeaksFile;
bool      gMinimumIsSet = false ; // the selection of range go in pairs (min max)  
bool      gMaximumIsSet = false ; // the selection of range go in pairs (min max)  

//___________________ FUNCTIONS
void PrintFunctionList();
void OpenDataFile();  
int  OpenDataThisFile(TString filename);
void AssignEnergy(float energy); 
void GetNextHistogram(int a);
void RebinHistogram();
void WriteRebin();
void GetControlBar();
void RangeClicked();  
void OpenPeaksFile();
void ResetParameters();
void Close(); 



//___________________ Main
void PeakSearch_207Bi_interactive() {

	cout << endl ; 
 	PrintFunctionList();

// Get control bar, one can make different control bar for different Nucleus
 	GetControlBar(/*gNucleus*/);

// Open data file 	
	OpenDataFile();
	
// Open peaks output file 
    OpenPeaksFile(); 
	
}


//___________________
void ResetParameters() {

	gBining=0 ; // Number of times the histogram was binned, pow(2,gBinning) is the number stored 
	gMin = -1  ; 
	gMax = -1  ;

	gMinimumIsSet = false ; // the selection of range go in pairs (min max)  
	gMaximumIsSet = false ; // the selection of range go in pairs (min max) 	
}

//______________________________________
void PrintFunctionList() {
	
cout << " List of functions :"<<endl
     << "GetControlBar()"<<endl
     << "OpenDataFile()"<<endl
     << endl; 
}


//_______________________________________
 void OpenDataFile() {
 
 // add command asking for the name/path...
 //TString fname = "his30093"; 
  TString fname = "hisXXXXX"; //default
 OpenDataThisFile(fname) ; 
 }
 
//_______________________________________ 
 int OpenDataThisFile(TString fname ) {
 
 	int success = 0 ;

       // open the ROOT file to process
   TString path ="./data1/moukaddam/SpiceTestSep2014/Calibration/Files/";
   TString inFileName = fname+".root";
   TFile *inFile = new TFile(path + inFileName);
   
   if ( !inFile->IsOpen() )  { //try present directory
   cout << "File doesn't exist in the directory  : " << path << endl ; 
   cout << "Trying the present working directory : ./" << endl ; 
   path ="./";
   inFile = new TFile(path + inFileName);
   }
   
   if ( inFile->IsOpen() ) {
       success = 1 ;
   	   cout << "Opening the root file and grabing the histograms from " << inFile->GetName() << endl ; 
   	   inFile->ls(); 
   }
   else {
   cout << "File is not found.. EXIT!"<< endl;
   exit(-1);
   } 
	   
	 gFolderHistos = (TFolder*)(inFile->FindObjectAny("histos"));
	 
	 return success ;
   }  
   
//__________________________________________________________________
void GetControlBar() {

   gROOT->Reset();

   //Add the tutorials directory to the macro path
   //This is necessary in case this macro is executed from another user directory
   TString dir = gSystem->UnixPathName(gInterpreter->GetCurrentMacroName());
   dir.ReplaceAll("PeakSearch_207Bi_interactive.C","");
   dir.ReplaceAll("/./","");
   const char *current = gROOT->GetMacroPath();
   gROOT->SetMacroPath(Form("%s:%s",current,dir.Data()));

   TControlBar *bar = new TControlBar("vertical", "Controls",20,10);
   //bar->SetNumberOfColumns(3);
   //bar->SetNumberOfRows(7);
   
   bar->AddButton("Root Browser",   	"new TBrowser;",         "Start the ROOT Browser");
   bar->AddButton("Open folder",   	"OpenDataFile();",         "Open data folder, provide name on screen");
   bar->AddButton("OpenPeaksFile",     	"OpenPeaksFile()",   "Open the output file of the peaks details");
   bar->AddButton("Next histogram", 		"GetNextHistogram(+1)"," Load the next histogram");
   bar->AddButton("Prev histogram", 		"GetNextHistogram(-1)"," Load the previous histogram");
   bar->AddButton("Rebin histogram ", 		"RebinHistogram()"," Rebin the current histogram");
   bar->AddButton("Store binning ", 		"WriteRebin()"," Write the adequate bining of the histogram ");
   bar->AddButton("Close",   "Close()",          "Close the files and exit");
  
 TControlBar *bar2 = new TControlBar("vertical", "Bi207 (keV)",100,100);
   //bar2->SetNumberOfColumns(3);
   //bar->SetNumberOfRows(7);
    
   bar2->AddButton("XRay  72.8",   "AssignEnergy(72.805)",          "AssignEnergy this energy : 21.4%");
   bar2->AddButton("XRay  75.0",   "AssignEnergy(74.969)",          "AssignEnergy this energy : 35.7%");
   bar2->AddButton("XRay  84.4",   "AssignEnergy(84.450)",          "AssignEnergy this energy : 4.31%");
   bar2->AddButton("XRay  84.9",   "AssignEnergy(84.938)",          "AssignEnergy this energy : 8.27%");
   bar2->AddButton("XRay  87.3",   "AssignEnergy(87.300)",          "AssignEnergy this energy : 3.02%");

   //bar->AddSeparator() ;

//TControlBar *bar3 = new TControlBar("vertical", "Bi207 ICE (keV)",50,50);
   //bar2->SetNumberOfColumns(3);
   //bar->SetNumberOfRows(7);
    
   bar2->AddButton("K   481.7",   "AssignEnergy(481.6935)",          "AssignEnergy this energy : 1.515%");
   bar2->AddButton("L   553.8",   "AssignEnergy(553.8372)",          "AssignEnergy this energy : 0.438%");
   bar2->AddButton("M   565.8",   "AssignEnergy(565.8473)",          "AssignEnergy this energy : 0.147%");
   bar2->AddButton("K   975.6",   "AssignEnergy(975.651)",          "AssignEnergy this energy : 7.03%");
   bar2->AddButton("L  1047.9",   "AssignEnergy(1047.795)",          "AssignEnergy this energy : 1.84%");
   bar2->AddButton("M  1059.8",   "AssignEnergy(1059.805)",          "AssignEnergy this energy : 0.54%");
   
   bar->SetButtonWidth(200);
   bar->Show();
   
   bar2->SetButtonWidth(200);
   bar2->SetTextColor("blue");
   //bar2->SetFont("-adobe-helvetica-bold-r-*-*-16-*-*-*-*-*-iso8859-1");
   bar2->Show();
   
   /*
   bar3->SetButtonWidth(100);
   bar3->SetTextColor("red");
   bar3->SetFont("-adobe-helvetica-bold-r-*-*-16-*-*-*-*-*-iso8859-1");
   bar3->Show();
   */
   
   gROOT->SaveContext();
}


//___________________________
void GetNextHistogram(int a ) {

	if(!gCanvas) {
		//gCanvas->Delete();
		TString title = "Bi207" ;
		gCanvas = new TCanvas(title,title, 1200,700);
		gCanvas->SetCrosshair(2);
		gCanvas->SetFrameFillColor(kBlack);
		gCanvas->SetFillColor(kBlack);
		gCanvas->ToggleEventStatus();
		gCanvas->AddExec("ex","RangeClicked()");
		//gCanvas->Divide(1,3) ; 
	}
	
	if (!gFolderHistos) {
		cout << " You need to open a histogram folder first " <<endl ; 
		return ; 
	}
	
	if (a>0)
    	gChannel_current++;
    else 
        gChannel_current--;
        
	TString hname = Form("Chrg%d", gChannel_current); 
	gHist_current = (TH1F*)(gFolderHistos->FindObjectAny(hname));
	cout << "=========================================";
	cout << " Retreiving Channel : " << gChannel_current ;
	cout << " ========================================="<< endl;
	
	//Write the channel number in file
	gPeaksFile << "=========================================\n";
	gPeaksFile << "# " <<gChannel_current << "\n" ;  

	// Reset Parameters
	ResetParameters() ; 
	
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
	gHist_current->Draw(); 
	
	gCanvas->Update(); 
}


//___________________________
void RebinHistogram() {

	if (gHist_current)	{
	gBining++ ; 
	gHist_current->Rebin() ; 	
	gHist_current->GetXaxis()->SetRangeUser(0,3500);
	gCanvas->Update();
	//gCanvas->AddExec("ex","RangeClicked()"); 
	cout << " Histogram rebinned - " << gBining << " - time(s)\n" ;  
	}
	else cout << " No histogram is loaded\n Press [Open Folder] followed by [Next] \n" ; 	
}

//___________________________  29 Sep 2014, MHD : Must find a better solution 
void WriteRebin() {
	gPeaksFile << "Rebin " <<pow(2,gBining) << "\n" ;
	cout << "Store binning " <<pow(2,gBining) << "\n" ;
}

//___________________________ 
void RangeClicked() {

	TPad *pad = (TPad *) gCanvas->GetSelectedPad();
	if (!pad) return;

	int eventtype = pad->GetEvent(); 
	
	if ( pad->GetEvent() == 12 && !gMinimumIsSet && !gMaximumIsSet) {   //   kButton1Down   //kButton1Double=61    //middle button 1 click = 12

		//Get the abscissa  
		gMin = pad->GetEventX();
		gMin = pad->AbsPixeltoX(gMin);
		printf(" [ %.3f , \n",gMin);
		gMinimumIsSet = true ;
		gMaximumIsSet = false ; 
		return ; // after this return the condition of the block below will be satisfied  
		} 
		
	if ( pad->GetEvent() == 12 && gMinimumIsSet && !gMaximumIsSet) {

		//Get the abscissa  
		gMax = pad->GetEventX();
		gMax = pad->AbsPixeltoX(gMax);
		
		//print the values
		if(gMin>=gMax) {	cout << " WARNING : Min >= Max , "<< " choose a value > "<< gMin << " \n" << endl ; return ; }
		else printf(" [ %.3f , %.3f ]\n",gMin,gMax);
		
		//reset 
		gMinimumIsSet=false ;
		gMaximumIsSet=false ;
		return ; 
		}
	
	return ;

}


//______________________________________
void OpenPeaksFile() { //Function to be called once

	if (!gPeaksFile.is_open()) {
	cout << "provide the name of the peaks output-file" <<endl ;
	string filename = "peaks_file.dat"; 
	//cin>>filename; 
	cout << "Writing the peaks details in file : " << filename <<endl ;
	gPeaksFile.open ("peaks_file.dat", std::ofstream::out);
	gPeaksFile << "Channel#\t\tmin\t\tmax\t\tenergy\t\t" << endl ;
	}
	else
	cout << "File already opened" <<endl ;

}

//______________________________________
void AssignEnergy(float energy){

	if (gPeaksFile.is_open()) {
		      cout << "        \t\t" << gMin << "\t\t" << gMax << "\t\t" << energy << "\n" ; // print on screen 
	    gPeaksFile << "        \t\t" << gMin << "\t\t" << gMax << "\t\t" << energy << "\n" ;
	    //reset 
	    gMin = -1 ; 
	    gMax = -1 ; 
		}
	else {
		cout << " The peaks file is not open " << endl ; 
		} 	
	
}

//______________________________________
void Close() {
gPeaksFile.close(); 
exit(-1); 
}

