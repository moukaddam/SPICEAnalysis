

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
void Stop() {
	char c; 
	c = getchar(); 
	if (c=='q') exit(-1) ;
	}


//USER definitions
#define STOP Stop();
#define NODES 9   // it should always be even n+2+1 where n is a multiplet of 6 
#define THIS if (1)//if(j==2&&i==0)






TGraph* RotateGraph(TGraph* g, float angle ) { // radian

cout << "angle " << angle*TMath::RadToDeg() << endl  ; 
Double_t* x = g->GetX(); 
Double_t* y = g->GetY();
Int_t n = g->GetN();
	 	 	   
//new Values 
Double_t* xx = new Double_t[n]; 
Double_t* yy = new Double_t[n];

for (int i=0 ; i< n ; i++){
	TVector2 v1(x[i],y[i]);
	v1.Dump(); 
	TVector2 v2 = v1.Rotate(angle); //radian 
	v2.Dump(); 
	xx[i]=v2.X(); 
	yy[i]=v2.Y();	
	 cout << "xx " << xx[i] << "  x  " << x[i] << endl  ; 
	 cout << "yy " << yy[i] << "  y  " << y[i] << endl  ;
}

	 	 	  // for (int i = 0 ; i< n ; i++){
	 	 	  // cout << "xx " << xx[i] << "  x  " << x[i] << endl  ; 
			  // cout << "yx " << yy[i] << "  y  " << y[i] << endl  ;
	 	 	  //  }

TGraph* gg = new TGraph(n,xx,yy) ; 
//gg= new TPolyLine(NODES,x,y);    
//gg->SetUniqueID(i*12+j);
//gg->SetFillColor(kRed+i);
//gg->SetFillColor(38);
gg->SetLineColor(2);
gg->SetLineWidth(2);
gg->SetMarkerColor(4);
gg->SetMarkerStyle(21);


//STOP

return gg ; 
}





void SpecSelector() {
   
	TCanvas *c1 = new TCanvas("SiLi_Segments","SiLi_Segments",10,10,1200,1200);
	//c1->Divide(12,10);
	//c1->Draw();
	
	
	 TMultiGraph *mg = new TMultiGraph("all","all");
     
	TGraph *pl[10];
	//TPolyLine *pl[120];
				
	int rings_number = 10 ; //10
	float rad_start = 5 ;
	float rad_end = 15 ; 
	float rad_pitch = (rad_end-rad_start)/rings_number ;

	int sectors_number = 12 ; //12      
	float phi_start = 0*TMath::DegToRad(); ;
	float phi_end = 360*TMath::DegToRad(); ; 
	float phi_pitch = ((phi_end-phi_start)/sectors_number); // 30 degrees
   cout << phi_pitch << endl ; 
   //STOP
	Double_t x[NODES]; // the additional point correspond to the final point superimposed on the starting point 
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
		
		cout << "rad min " << rad_min << endl  ; 
		cout << "rad max " << rad_max << endl  ;
		  
		for (Int_t j = 0; j < 1 /*3*/ ; j++) { // loop on sectors
		
		//if (j==2 && i>0) continue ; 
		
		// phi min and phi max are the borders of the current segment	
			THIS phi_min = phi_start + j*phi_pitch ;
			THIS phi_max = phi_start + (j+1)*phi_pitch ;
		
			THIS cout << "phi min " << phi_min*TMath::RadToDeg() << endl  ; 
			THIS cout << "phi max " << phi_max*TMath::RadToDeg() << endl  ;
			
			// phi and rad are now defined for one segment
			phi_current = phi_min ; // Initiate the current phi as the minimum phi on the current segment
			rad_current = rad_min ; // Initiate the current radius as the inner radius on the current segment
			int sign = +1 ; // Initiate the current radius as the minimum radius on the current segment
			
			for (Int_t k = 0; k < NODES-1 ; k++) { // loop on the nodes, first for inner border with sign = +1
			 	
			 	THIS cout << " k =  " << k  <<endl ;  
				x[k] =rad_current*TMath::Cos( phi_current ) ; 
				y[k] =rad_current*TMath::Sin( phi_current ) ;
				
				THIS cout << "rad current " << rad_current << endl  ; 
				THIS cout << "phi current " << phi_current*TMath::RadToDeg() << endl  ;
			    THIS cout << "x " << x[k] << endl  ; 
			    THIS cout << "y " << y[k] << endl  ;

				if (k==0) {	//Assign the last point
		 	 		THIS cout <<" Assign the last point for the node " << NODES <<endl ; 
		 	 		x[NODES-1] = x[k]; 
		 	 		y[NODES-1] = y[k];
					}
					
				 if ( phi_current == phi_max && rad_current == rad_min) { // Jump up to radius max but keep the same phi ()
					THIS cout << " rad_current now is max, phi will be the same "<<endl ;  
					rad_current = rad_max ; 
					sign = 0;
					}
				else if ( rad_current == rad_max && phi_current > phi_min) {
						THIS cout << " rad_current now is still max, phi will decrease"<<endl ;   
						sign = -1;
						}
			
				phi_current = phi_current + sign*phi_mini_pitch;  // loop within segment on outer border, counter clock wise
	 	 	
	 	 		}// finish the loop with the segment 

	 	 	   // for (int i = 0 ; i< NODES ; i++){
	 	 	   //cout << "x  " << x[i] << endl  ; 
			   //cout << "y  " << y[i] << endl  ; 
	 	 	   // }
	 	 	    
	 	 	    // Draw
				//c1->cd(i*12+j);

				//pl[i*12+j]= new TPolyLine(NODES,x,y);    
				//pl[i*12+j]->SetUniqueID(i*12+j);
				//pl[i*12+j]->SetFillColor(kRed+i);
				//pl[i*12+j]->SetFillColor(38);
				//pl[i*12+j]->SetLineColor(2);
				//pl[i*12+j]->SetLineWidth(2);
				//pl[i*12+j]->SetMarkerColor(4);
				//pl[i*12+j]->SetMarkerStyle(21);
				//pl[i*12+j]->Draw("ALP");
				//pl[i*12+j]->Draw();
				//pl[i*12+j]->Draw("Af");

			//STOP	
			}
			
		 pl[i]= new TGraph(NODES,x,y);
							
		}
		
for (Int_t i = 0; i < rings_number ; i++) { // loop on rings		
	for (Int_t j = 0; j < sectors_number ; j++) { // loop on sectors
	//RotateGraph(pl[j],j*phi_pitch) ; 
   cout << j*phi_pitch << endl ; 
   //STOP
		mg->Add(RotateGraph(pl[i],j*phi_pitch),"lp"); 
	}						
}

	 //c1->cd(i*12+j);
     mg->Draw("a");
     
}
  //c1->AddExec("ex","TriangleClicked()");
