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

// functions 
void Stop();
TGraph* RotateGraph(TGraph* g, float angle ) ;

//USER definitions
// it should always be of form n+2+1 where:
// n is the number of jumps from point(node) to point(node)m this number is always even, 
// n+2 is the number of nodes
// n+2+1 the additional point correspond to the final point superimposed on the starting point(node) 
#define NODES 8+3   
#define STOP Stop();




//MAIN FUNCTION 
void SpecSelector() {
   
	TFile* file = new TFile("SiLiSegmentsGraph.root","RECREATE");
	TCanvas *c1 = new TCanvas("SiLi_Segments","SiLi_Segments",10,10,700,700);	
	TMultiGraph *mg = new TMultiGraph("all","all");
	TGraph *gr[10];
				
	int rings_number = 10 ; //10
	float rad_start = 3 ;
	float rad_end = 15 ; 
	float rad_pitch = (rad_end-rad_start)/rings_number ;

	int sectors_number = 12 ; //12      
	float phi_start = 0*TMath::DegToRad(); ;
	float phi_end = 360*TMath::DegToRad(); ; 
	float phi_pitch = ((phi_end-phi_start)/sectors_number); // 30 degrees

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
  //c1->AddExec("ex","TriangleClicked()");



void Stop() {
	char c; 
	c = getchar(); 
	if (c=='q') exit(-1) ;
	}

TGraph* RotateGraph(TGraph* g, float angle ) { // radian

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
	gg->SetLineColor(2);
	gg->SetLineWidth(2);
	gg->SetMarkerColor(4);
	gg->SetMarkerStyle(1);

	return gg ; 
}



