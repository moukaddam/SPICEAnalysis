/*

This code is used to search for peaks in from a spectrum of a Bi207 electron source
The background is subtracted before the fit is made on the different peaks

In principle the mnemonic could be found by hist->GetTitle(), But the mapping was wrong 
so we should use the GetMnemonic function.

The code spits out two files, one root file for inspection of fit and one ascii file
to calibrate with, for now they are called :

inputData_his30093
segments_his30093.root

NB : Dont forget to change 
#define ODB_CHANNELS 1
to 
#define ODB_CHANNELS 120


Process Particle Energy 	Intensity
Auger K e  	56.7 keV 		2.9%     Low in energy
CE K e 		481.6935 keV 	1.515%
CE L e 		553.8372 keV 	0.438%
CE M e 		565.8473 keV 	0.147%
CE K e 		975.651 keV 	7.03%
CE L e 		1047.795 keV 	1.84%
CE M e 		1059.805 keV 	0.54%

XR       	10.6 keV 		33.2%     Low in energy 
XR       	72.805 keV 		21.4%
XR       	74.969 keV 		35.7%
XR       	84.45 keV 		4.31%
XR       	84.938 keV 		8.27%
XR       	87.3 keV 		3.02%
Gam 		569.698 keV 	97.76%   high in energyGam			1063.656 keV 	74.6%    high in energyGam			1770.228 keV 	6.87%    high in energy


*/

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

//User definitions
#define ODB_CHANNELS 120
#define DEBUG 0

int odb_offset = 1060 ; // depends on odb	

struct peak_buffer_st {

    TString mnemonic  ;
	TString address  ;
	
    int ring  ;
	int sector  ;
	int binning  ;
	
	std::vector<float> energy;
	std::vector<float> minimum;
	std::vector<float> maximum;

	// to be fitted between [minimum, maximum]	
	std::vector<float> mean; 
	std::vector<float> sigma; 
	std::vector<float> err;
	std::vector<float> ndf;
	std::vector<float> chi2; 
	std::vector<float> integral;
	std::vector<float> rms; 
};

// Global variables
TH1F *hist[ODB_CHANNELS] ; // histogram (original) 
TH1F *hist_bg[ODB_CHANNELS] ;  // histogram (with background subtracted)
TF1  *fGaus[10] ;  // could store up to 10 peaks
TCanvas * gCanvas; 
TFolder *gHistos ;
peak_buffer_st * peaker = new peak_buffer_st[ODB_CHANNELS];

//User functions 
TString GetMnemonic(int channel ) ;
TString GetAddress(int channel ) ;
int 	GetRingSector(int select, TString mnemonic ) ;
int     ReadInput(TString peaks_fname) ;
int     GrabHistos(TString path, TString fname) ;
void 	SubtractBackground(int odb_channel);   

// Main function
void PeakSearch_207Bi_Simple(const char* input_path = "./../../Calibration/Files/",const char* input_file_name = "his30093",
							 const char* detail_peaks_fname = "peaks_file.dat") {

	gStyle->SetPalette(1);

   // open output root file 	
   TFile *outFile ;
   TString prefix =  "background_sub_" ;
   TString output_file_name = prefix + input_file_name + ".root";
   outFile = new TFile(output_file_name, "RECREATE");
   if (outFile->IsOpen() )  {
   cout << " Output root file is created successfully" << endl ; 
   }
   else { 
   cout << " Output root file "<< outFile->GetName() <<" could not be created, will continue..." <<endl ;
   }

	//***********************
	//users message
	//***********************   
    cout<<"==============================" <<endl;
   	cout<<"Looking in path                    : "<<input_path<<endl;
	cout<<"for root file                     : "<<input_file_name<<endl<<endl;
	
	cout<<"Output in the working directory     "<<endl;
	cout<<"Background inspection root file   : "<<output_file_name<<endl;
	cout<<"Peaks details file                : "<<detail_peaks_fname<<endl;
	cout<<"==============================" <<endl;

	//***********************
	//read the ranges
	//***********************
	ReadInput(detail_peaks_fname);

	//***********************
	//Grab the gHistos
	//***********************
	GrabHistos(input_path, input_file_name);

	//**********************************************
	// Fit with guassians, and store in a root file
	//**********************************************
	for (int odb_channel =  0; odb_channel < ODB_CHANNELS; odb_channel++) {


		SubtractBackground(odb_channel); 

		for (int i = 0 ; i < peaker[odb_channel].energy.size() ; i++) {
		
			// set the limits
			int lower = peaker[odb_channel].minimum.at(i);
			int upper = peaker[odb_channel].maximum.at(i);
			hist_bg[odb_channel]->GetXaxis()->SetRangeUser(lower,upper);

 			// fit with guassian 
			float Amplitude = hist_bg[odb_channel]->GetMaximumBin() ; 
			float Center = hist_bg[odb_channel]->GetBinCenter(Amplitude);
			float Rms = hist_bg[odb_channel]->GetRMS();
		
		    TString fGausName = Form("fGaus%d",i);
			fGaus[i] = new TF1(fGausName,"gaus",lower,upper);
			fGaus[i]->SetLineColor(kRed);
			fGaus[i]->SetLineWidth(3);
	  		fGaus[i]->SetParameters(Amplitude,Center,Rms);
			hist_bg[odb_channel]->Fit(fGaus[i],"RMQ");
			gCanvas->cd(2)->Update() ;
			
			float mean		= fGaus[i]->GetParameter(1);
			float mean_err	= fGaus[i]->GetParError(1);
			float sigma 	= fGaus[i]->GetParameter(2);
			float sigma_err = fGaus[i]->GetParError(2);

			float ndf   	= fGaus[i]->GetNDF();
			float chi2  	= fGaus[i]->GetChisquare();

			//float integral= hist_bg[odb_channel]->Integral();
			//float rms		= hist_bg[odb_channel]->GetRMS();
			float integral  = fGaus[i]->Integral(mean-3*sigma,mean+3*sigma) ; // this should scale with above
		    float rms		= fGaus[i]->Variance(mean-3*sigma,mean+3*sigma) ; // this should scale with above
		    
			peaker[odb_channel].mean.push_back(mean) ;
			peaker[odb_channel].sigma.push_back(sigma);
			peaker[odb_channel].err.push_back(mean_err);
			peaker[odb_channel].ndf.push_back(ndf);
			peaker[odb_channel].chi2.push_back(chi2);
			peaker[odb_channel].integral.push_back(integral);
			peaker[odb_channel].rms.push_back(rms);
		 			
			}// loop on peaks 
			
			// Draw the Gaussians 
			hist_bg[odb_channel]->GetXaxis()->SetRangeUser(0,4000);
			hist_bg[odb_channel]->SetStats(false) ;
			for (int i = 0 ; i < peaker[odb_channel].energy.size() ; i++) {
			fGaus[i]->DrawCopy("LSAME");
			}
			gCanvas->Update();
			
			outFile->cd() ;
			gCanvas->Write();	
	
		}//LOOP on segment
		
		
	//**********************************
	// Write the peak parameters to file
	//**********************************	
    ofstream peaks_file;
	TString peaks_file_name = Form("peaks_detail");
	peaks_file_name = peaks_file_name + "_" + input_file_name + ".dat"; 
	peaks_file.open("./"+peaks_file_name);
	//peaks_file<< "channel" << "\t" << "Address"<< "\t" << "ring"  <<"\t"<< "sector"  << "\t" << "Integral"  <<"\t"<< "mean"  <<"\t"<< "gain/ch" << "\t"<<"sigma" << "\t" <<"err"<<"\t"<<"chi2"<<endl ;

	peaks_file	<< "\t"
				<< fixed << showpoint << setw(8) << setprecision(2) 
				<< "energy" <<"\t"
				<< fixed << showpoint << setw(8) << setprecision(2) 
				<< "energy_err" <<"\t" 
				<< fixed << showpoint << setw(8) << setprecision(2) 
				<< "mean" <<"\t" 
				<< fixed << showpoint << setw(8) << setprecision(2) 
				<< "mean_err" <<"\t"
				<< fixed << showpoint << setw(8) << setprecision(2) 				
				<< "sigma"  <<"\t"  
				<< fixed << showpoint << setw(8) << setprecision(2) 				
				<< "integral" <<"\t" 
				<< fixed << showpoint << setw(8) << setprecision(2) 				
				<< "rms"  <<"\t" 
				<< fixed << showpoint << setw(8) << setprecision(2) 				
				<< "red.chi2"<<"\n" ;
 						
	for (int odb_channel =  0; odb_channel < ODB_CHANNELS; odb_channel++) {	
	  		
		peaks_file << "Chan# " << (odb_channel + odb_offset)<<endl;  
	  					
		for (int i = 0 ; i < peaker[odb_channel].energy.size() ; i++) {
			peaks_file	<< "\t"
						<< fixed << showpoint << setw(8) << setprecision(2) 
						<< peaker[odb_channel].energy.at(i) <<"\t"
						<< fixed << showpoint << setw(8) << setprecision(2) 
						<< 0.1 <<"\t"  // this value should be fixed in a proper way
						<< fixed << showpoint << setw(8) << setprecision(2) 
						<< peaker[odb_channel].mean.at(i) <<"\t" 
						<< fixed << showpoint << setw(8) << setprecision(2) 
						<< peaker[odb_channel].err.at(i)  <<"\t"
						<< fixed << showpoint << setw(8) << setprecision(2) 
						<< peaker[odb_channel].sigma.at(i)  <<"\t"  
						<< fixed << showpoint << setw(8) << setprecision(2) 
						<< peaker[odb_channel].integral.at(i) <<"\t" 
						<< fixed << showpoint << setw(8) << setprecision(2) 
						<< peaker[odb_channel].rms.at(i)  <<"\t" 
						<< fixed << showpoint << setw(8) << setprecision(2) 
 						<< peaker[odb_channel].chi2.at(i) / peaker[odb_channel].ndf.at(i)<<"\n" ;
			}// loop on peaks 
			peaks_file << "===============" << endl;
		
	}// loop on segments
			
			
		//*************************
		// Finish and close files
		//**************************	
		peaks_file.close();
		outFile->Write();
		outFile->Close(); 

}// END





 int GrabHistos(TString path, TString input_file_name) {
 
 	int success = 0 ; 
	TFile *inFile ; 
	TString extension =  ".root" ;
	TString inFileName = input_file_name+extension;
	inFile = new TFile(path + inFileName);
   
	if (DEBUG) { 
		cout << "Opening the root file and grabing histograms " << inFile->GetName() << endl ;  
		inFile->ls();
		getchar();
	}
   
	if ( inFile->IsOpen() )  {
		cout << "Histogram files was opened successfully " << endl ;
		gHistos = (TFolder*)(inFile->FindObjectAny("histos"));  
		success = 1 ; 
	}
	else success = 0 ; 

	return success ;
}  
   
   
   
   
int GetRingSector(int select, TString mnemonic ) {

	TString sub( mnemonic(7,9) ) ;
	int sector = -1 ;
	int ring = -1 ;

	int ringsector = atoi(sub.Data());
	ring = ringsector/12 ; // rings goes from (0 - 9)
	sector = (ringsector-(ring*12)); // rings goes from (0 - 11)
	ring = 9 - ring ; // the rings where inverted to make the mnemonics
	
	if (select == 0) select = ring;
	else select = sector;

	return select ; 
	
	}




int ReadInput(TString peaks_fname) {  
   
	int success = 0 ;
	int channel = -1 ; 
	int odb_channel = -1 ;  
	string dummy_string="";
	string dummy_word=""; 
	float minimum = -1 ; 
	float maximum = -1 ; 
	float energy = -1 ; 
   
	ifstream peaks_file ; 
	peaks_file.open (peaks_fname.Data()) ; 

	if (peaks_file.is_open()) {
		success = 1 ; 
	  	cout << " peaks input file is opened..." << endl ;  
		getline (peaks_file,dummy_string);
 
		while (true) {
  		getline (peaks_file,dummy_string);
  		cout << dummy_string << endl;
  		
	  	string indicator = dummy_string.substr(0, 1);// get the first letter 
	  	istringstream ss(dummy_string);
		
	  	if (peaks_file.eof()) break;
		
		if ( indicator == "=")
			continue;
		else if ( indicator == "#") {
		ss>>dummy_word>>channel ; 
		odb_channel =  channel - odb_offset ;
		} 
		else 
			if ( indicator == "%") {
				if (DEBUG) cout << " Comment : " << dummy_string << endl ; 
				continue;
			}
			else
				if ( indicator == "R") {
					int binning = 0  ; 
					ss>>dummy_word>>binning ;
					if (DEBUG) cout << " Rebin Detected : " << binning << endl  ; 
					peaker[odb_channel].binning = binning ; 
				}
				else 
				//if ( indicator == "e") // this tag could be used to tag/consider x-ray or klm-electron 
				{
						if (DEBUG) cout << " Potential point Detected  " ;  
						// ss >> dummy_word ; in case there's a tag  
						ss >> minimum >> maximum  >> energy ;
						if (minimum!=-1 && maximum!=-1) {
						if (DEBUG) cout << " / Accepted :  " << minimum << " " << maximum <<endl  ;  
						peaker[odb_channel].energy.push_back(energy); 		
						peaker[odb_channel].minimum.push_back(minimum);
						peaker[odb_channel].maximum.push_back(maximum);
						}
				}
				
		} // while loop

	}// if file is opened
	else {
	  cout << "file not found ! will exit " << endl ; 
	  success = 0 ; 
	}
	
	if (DEBUG) { 
	  cout << " Dump the content of peaker " << endl; getchar();
		for ( odb_channel = 0 ; odb_channel < ODB_CHANNELS ; odb_channel++)	{
			for (int j = 0 ; j < peaker[odb_channel].energy.size() ; j++) {
			  cout << odb_channel + odb_offset << "\t" ;
			  cout << peaker[odb_channel].energy.at(j) << "\t" ;
			  cout << peaker[odb_channel].minimum.at(j)<< "\t" ;
			  cout << peaker[odb_channel].maximum.at(j)<< "\n" ;
			}
		  cout << "+++" << "\n" ;
		}
	}
  
  return success ; 
} 


void SubtractBackground(int odb_channel) {

	 	if (DEBUG) { cout << " Subtract the Background of channel " <<odb_channel<<endl; }
       
        TString canvas_title = Form("%d",odb_offset+odb_channel);
	 	TString canvas_name = "Bi207";
		if(!gCanvas) {
		gCanvas = new TCanvas(canvas_name,canvas_title, 1000,700);
        gCanvas->Divide(1,2) ;  
        }
        gCanvas->SetTitle(canvas_title) ;  
        
        
		// get the histogram, copy it and subtract the background 
		TString hname = Form("Chrg%d", odb_offset+odb_channel);
		hist[odb_channel] = (TH1F*)(gHistos->FindObjectAny(hname)); 
		if ( peaker[odb_channel].binning>0 ) 
			hist[odb_channel]->Rebin(peaker[odb_channel].binning); 
		hist[odb_channel]->GetXaxis()->SetRangeUser(0,4000);
        hist[odb_channel]->GetXaxis()->SetLabelSize(0.04);
		hist[odb_channel]->GetXaxis()->SetNdivisions(5, 5, 0); 
		hist[odb_channel]->GetYaxis()->SetLabelSize(0.04);
		hist_bg[odb_channel] = (TH1F*)hist[odb_channel]->Clone(Form("hist_bg_%d",odb_channel+odb_offset));

		// Prepare the Canvases
        gCanvas->cd(1) ;
        hist[odb_channel]->DrawCopy("") ;  
         		
		// search for the right pad position 
		TString mnemonic = GetMnemonic(odb_offset + odb_channel);
		TString address = GetAddress(odb_offset + odb_channel);
		int sector = GetRingSector(1,mnemonic); 
		int ring =  GetRingSector(0,mnemonic);
		
		//Store in the peak container
		peaker[odb_channel].mnemonic 	= mnemonic;
		peaker[odb_channel].address 	= address;
		peaker[odb_channel].ring 		= ring;
		peaker[odb_channel].sector 		= sector;

		int S = 20 ; // smoothing parameter to optimize background subtraction
		TH1F *hbackground = (TH1F*) hist[odb_channel]->ShowBackground(S, "SAME") ; // nb of iteration = 10,   
		//hbackground->SetFillStyle(3001);
		hbackground->SetFillColor(kYellow);
		hbackground->SetLineColor(kRed);

		gCanvas->cd(2) ;
		hist_bg[odb_channel]->Add(hist[odb_channel],hbackground,+1,-1); // h = +1*h - 1*background ;
		hist_bg[odb_channel]->Draw("");
	
}

TString GetAddress(int channel ){

	TString address="";

	switch (channel){
	case 1060:       address="0x0800100"; break;
	case 1061:       address="0x0800101"; break;
	case 1062:       address="0x0800102"; break;
	case 1063:       address="0x0800103"; break;
	case 1064:       address="0x0800104"; break;
	case 1065:       address="0x0800105"; break;
	case 1066:       address="0x0800106"; break;
	case 1067:       address="0x0800107"; break;
	case 1068:       address="0x0800108"; break;
	case 1069:       address="0x0800109"; break;
	case 1070:       address="0x0800200"; break;
	case 1071:       address="0x0800201"; break;
	case 1072:       address="0x0800202"; break;
	case 1073:       address="0x0800203"; break;
	case 1074:       address="0x0800204"; break;
	case 1075:       address="0x0800205"; break;
	case 1076:       address="0x0800206"; break;
	case 1077:       address="0x0800207"; break;
	case 1078:       address="0x0800208"; break;
	case 1079:       address="0x0800209"; break;
	case 1080:       address="0x0800300"; break;
	case 1081:       address="0x0800301"; break;
	case 1082:       address="0x0800302"; break;
	case 1083:       address="0x0800303"; break;
	case 1084:       address="0x0800304"; break;
	case 1085:       address="0x0800305"; break;
	case 1086:       address="0x0800306"; break;
	case 1087:       address="0x0800307"; break;
	case 1088:       address="0x0800308"; break;
	case 1089:       address="0x0800309"; break;
	case 1090:       address="0x0800400"; break;
	case 1091:       address="0x0800401"; break;
	case 1092:       address="0x0800402"; break;
	case 1093:       address="0x0800403"; break;
	case 1094:       address="0x0800404"; break;
	case 1095:       address="0x0800405"; break;
	case 1096:       address="0x0800406"; break;
	case 1097:       address="0x0800407"; break;
	case 1098:       address="0x0800408"; break;
	case 1099:       address="0x0800409"; break;
	case 1100:       address="0x0800500"; break;
	case 1101:       address="0x0800501"; break;
	case 1102:       address="0x0800502"; break;
	case 1103:       address="0x0800503"; break;
	case 1104:       address="0x0800504"; break;
	case 1105:       address="0x0800505"; break;
	case 1106:       address="0x0800506"; break;
	case 1107:       address="0x0800507"; break;
	case 1108:       address="0x0800508"; break;
	case 1109:       address="0x0800509"; break;
	case 1110:       address="0x0800600"; break;
	case 1111:       address="0x0800601"; break;
	case 1112:       address="0x0800602"; break;
	case 1113:       address="0x0800603"; break;
	case 1114:       address="0x0800604"; break;
	case 1115:       address="0x0800605"; break;
	case 1116:       address="0x0800606"; break;
	case 1117:       address="0x0800607"; break;
	case 1118:       address="0x0800608"; break;
	case 1119:       address="0x0800609"; break;
	case 1120:       address="0x0800700"; break;
	case 1121:       address="0x0800701"; break;
	case 1122:       address="0x0800702"; break;
	case 1123:       address="0x0800703"; break;
	case 1124:       address="0x0800704"; break;
	case 1125:       address="0x0800705"; break;
	case 1126:       address="0x0800706"; break;
	case 1127:       address="0x0800707"; break;
	case 1128:       address="0x0800708"; break;
	case 1129:       address="0x0800709"; break;
	case 1130:       address="0x0800800"; break;
	case 1131:       address="0x0800801"; break;
	case 1132:       address="0x0800802"; break;
	case 1133:       address="0x0800803"; break;
	case 1134:       address="0x0800804"; break;
	case 1135:       address="0x0800805"; break;
	case 1136:       address="0x0800806"; break;
	case 1137:       address="0x0800807"; break;
	case 1138:       address="0x0800808"; break;
	case 1139:       address="0x0800809"; break;
	case 1140:       address="0x0800900"; break;
	case 1141:       address="0x0800901"; break;
	case 1142:       address="0x0800902"; break;
	case 1143:       address="0x0800903"; break;
	case 1144:       address="0x0800904"; break;
	case 1145:       address="0x0800905"; break;
	case 1146:       address="0x0800906"; break;
	case 1147:       address="0x0800907"; break;
	case 1148:       address="0x0800908"; break;
	case 1149:       address="0x0800909"; break;
	case 1150:       address="0x0800a00"; break;
	case 1151:       address="0x0800a01"; break;
	case 1152:       address="0x0800a02"; break;
	case 1153:       address="0x0800a03"; break;
	case 1154:       address="0x0800a04"; break;
	case 1155:       address="0x0800a05"; break;
	case 1156:       address="0x0800a06"; break;
	case 1157:       address="0x0800a07"; break;
	case 1158:       address="0x0800a08"; break;
	case 1159:       address="0x0800a09"; break;
	case 1160:       address="0x0800b00"; break;
	case 1161:       address="0x0800b01"; break;
	case 1162:       address="0x0800b02"; break;
	case 1163:       address="0x0800b03"; break;
	case 1164:       address="0x0800b04"; break;
	case 1165:       address="0x0800b05"; break;
	case 1166:       address="0x0800b06"; break;
	case 1167:       address="0x0800b07"; break;
	case 1168:       address="0x0800b08"; break;
	case 1169:       address="0x0800b09"; break;
	case 1170:       address="0x0800c00"; break;
	case 1171:       address="0x0800c01"; break;
	case 1172:       address="0x0800c02"; break;
	case 1173:       address="0x0800c03"; break;
	case 1174:       address="0x0800c04"; break;
	case 1175:       address="0x0800c05"; break;
	case 1176:       address="0x0800c06"; break;
	case 1177:       address="0x0800c07"; break;
	case 1178:       address="0x0800c08"; break;
	case 1179:       address="0x0800c09"; break;

    }            
return address;
}


            

TString GetMnemonic(int channel ) {
TString mnemonic="";
switch (channel) {
                                                     // used for Sep 2014          // used for Dec 2013
case 1060:       mnemonic="SPI00XN019"; break;       //  SPI00XN019		    SPI00XN100
case 1061:       mnemonic="SPI00XN008"; break;       //  SPI00XN008		    SPI00XN113
case 1062:       mnemonic="SPI00XN032"; break;       //  SPI00XN032		    SPI00XN089
case 1063:       mnemonic="SPI00XN056"; break;       //  SPI00XN056		    SPI00XN065
case 1064:       mnemonic="SPI00XN080"; break;       //  SPI00XN080		    SPI00XN041
case 1065:       mnemonic="SPI00XN104"; break;       //  SPI00XN104		    SPI00XN017
case 1066:       mnemonic="SPI00XN115"; break;       //  SPI00XN115		    SPI00XN004
case 1067:       mnemonic="SPI00XN091"; break;       //  SPI00XN091		    SPI00XN028
case 1068:       mnemonic="SPI00XN067"; break;       //  SPI00XN067		    SPI00XN052
case 1069:       mnemonic="SPI00XN043"; break;       //  SPI00XN043		    SPI00XN076
case 1070:       mnemonic="SPI00XN006"; break;       //  SPI00XN006		    SPI00XN111
case 1071:       mnemonic="SPI00XN030"; break;       //  SPI00XN030		    SPI00XN087
case 1072:       mnemonic="SPI00XN054"; break;       //  SPI00XN054		    SPI00XN063
case 1073:       mnemonic="SPI00XN078"; break;       //  SPI00XN078		    SPI00XN039
case 1074:       mnemonic="SPI00XN102"; break;       //  SPI00XN102		    SPI00XN015
case 1075:       mnemonic="SPI00XN018"; break;       //  SPI00XN018		    SPI00XN099
case 1076:       mnemonic="SPI00XN042"; break;       //  SPI00XN042		    SPI00XN075
case 1077:       mnemonic="SPI00XN066"; break;       //  SPI00XN066		    SPI00XN051
case 1078:       mnemonic="SPI00XN090"; break;       //  SPI00XN090		    SPI00XN027
case 1079:       mnemonic="SPI00XN114"; break;       //  SPI00XN114		    SPI00XN003
case 1080:       mnemonic="SPI00XN007"; break;       //  SPI00XN007		    SPI00XN112
case 1081:       mnemonic="SPI00XN020"; break;       //  SPI00XN020		    SPI00XN101
case 1082:       mnemonic="SPI00XN044"; break;       //  SPI00XN044		    SPI00XN077
case 1083:       mnemonic="SPI00XN068"; break;       //  SPI00XN068		    SPI00XN053
case 1084:       mnemonic="SPI00XN092"; break;       //  SPI00XN092		    SPI00XN029
case 1085:       mnemonic="SPI00XN116"; break;       //  SPI00XN116		    SPI00XN005
case 1086:       mnemonic="SPI00XN103"; break;       //  SPI00XN103		    SPI00XN016
case 1087:       mnemonic="SPI00XN079"; break;       //  SPI00XN079		    SPI00XN040
case 1088:       mnemonic="SPI00XN055"; break;       //  SPI00XN055		    SPI00XN064
case 1089:       mnemonic="SPI00XN031"; break;       //  SPI00XN031		    SPI00XN088
case 1090:       mnemonic="SPI00XN112"; break;       //  SPI00XN112		    SPI00XN007
case 1091:       mnemonic="SPI00XN101"; break;       //  SPI00XN101		    SPI00XN020
case 1092:       mnemonic="SPI00XN077"; break;       //  SPI00XN077		    SPI00XN044
case 1093:       mnemonic="SPI00XN053"; break;       //  SPI00XN053		    SPI00XN068
case 1094:       mnemonic="SPI00XN029"; break;       //  SPI00XN029		    SPI00XN092
case 1095:       mnemonic="SPI00XN005"; break;       //  SPI00XN005		    SPI00XN116
case 1096:       mnemonic="SPI00XN016"; break;       //  SPI00XN016		    SPI00XN103
case 1097:       mnemonic="SPI00XN040"; break;       //  SPI00XN040		    SPI00XN079
case 1098:       mnemonic="SPI00XN064"; break;       //  SPI00XN064		    SPI00XN055
case 1099:       mnemonic="SPI00XN088"; break;       //  SPI00XN088		    SPI00XN031
case 1100:       mnemonic="SPI00XN099"; break;       //  SPI00XN099		    SPI00XN018
case 1101:       mnemonic="SPI00XN075"; break;       //  SPI00XN075		    SPI00XN042
case 1102:       mnemonic="SPI00XN051"; break;       //  SPI00XN051		    SPI00XN066
case 1103:       mnemonic="SPI00XN027"; break;       //  SPI00XN027		    SPI00XN090
case 1104:       mnemonic="SPI00XN003"; break;       //  SPI00XN003		    SPI00XN114
case 1105:       mnemonic="SPI00XN111"; break;       //  SPI00XN111		    SPI00XN006
case 1106:       mnemonic="SPI00XN087"; break;       //  SPI00XN087		    SPI00XN030
case 1107:       mnemonic="SPI00XN063"; break;       //  SPI00XN063		    SPI00XN054
case 1108:       mnemonic="SPI00XN039"; break;       //  SPI00XN039		    SPI00XN078
case 1109:       mnemonic="SPI00XN015"; break;       //  SPI00XN015		    SPI00XN102
case 1110:       mnemonic="SPI00XN100"; break;       //  SPI00XN100		    SPI00XN019
case 1111:       mnemonic="SPI00XN113"; break;       //  SPI00XN113		    SPI00XN008
case 1112:       mnemonic="SPI00XN089"; break;       //  SPI00XN089		    SPI00XN032
case 1113:       mnemonic="SPI00XN065"; break;       //  SPI00XN065		    SPI00XN056
case 1114:       mnemonic="SPI00XN041"; break;       //  SPI00XN041		    SPI00XN080
case 1115:       mnemonic="SPI00XN017"; break;       //  SPI00XN017		    SPI00XN104
case 1116:       mnemonic="SPI00XN004"; break;       //  SPI00XN004		    SPI00XN115
case 1117:       mnemonic="SPI00XN028"; break;       //  SPI00XN028		    SPI00XN091
case 1118:       mnemonic="SPI00XN052"; break;       //  SPI00XN052		    SPI00XN067
case 1119:       mnemonic="SPI00XN076"; break;       //  SPI00XN076		    SPI00XN043
case 1120:       mnemonic="SPI00XN118"; break;       //  SPI00XN118		    SPI00XN001
case 1121:       mnemonic="SPI00XN107"; break;       //  SPI00XN107		    SPI00XN014
case 1122:       mnemonic="SPI00XN083"; break;       //  SPI00XN083		    SPI00XN038
case 1123:       mnemonic="SPI00XN059"; break;       //  SPI00XN059		    SPI00XN062
case 1124:       mnemonic="SPI00XN035"; break;       //  SPI00XN035		    SPI00XN086
case 1125:       mnemonic="SPI00XN011"; break;       //  SPI00XN011		    SPI00XN110
case 1126:       mnemonic="SPI00XN022"; break;       //  SPI00XN022		    SPI00XN097
case 1127:       mnemonic="SPI00XN046"; break;       //  SPI00XN046		    SPI00XN073
case 1128:       mnemonic="SPI00XN070"; break;       //  SPI00XN070		    SPI00XN049
case 1129:       mnemonic="SPI00XN094"; break;       //  SPI00XN094		    SPI00XN025
case 1130:       mnemonic="SPI00XN105"; break;       //  SPI00XN105		    SPI00XN012
case 1131:       mnemonic="SPI00XN081"; break;       //  SPI00XN081		    SPI00XN036
case 1132:       mnemonic="SPI00XN057"; break;       //  SPI00XN057		    SPI00XN060
case 1133:       mnemonic="SPI00XN033"; break;       //  SPI00XN033		    SPI00XN084
case 1134:       mnemonic="SPI00XN009"; break;       //  SPI00XN009		    SPI00XN108
case 1135:       mnemonic="SPI00XN117"; break;       //  SPI00XN117		    SPI00XN000
case 1136:       mnemonic="SPI00XN093"; break;       //  SPI00XN093		    SPI00XN024
case 1137:       mnemonic="SPI00XN069"; break;       //  SPI00XN069		    SPI00XN048
case 1138:       mnemonic="SPI00XN045"; break;       //  SPI00XN045		    SPI00XN072
case 1139:       mnemonic="SPI00XN021"; break;       //  SPI00XN021		    SPI00XN096
case 1140:       mnemonic="SPI00XN106"; break;       //  SPI00XN106		    SPI00XN013
case 1141:       mnemonic="SPI00XN119"; break;       //  SPI00XN119		    SPI00XN002
case 1142:       mnemonic="SPI00XN095"; break;       //  SPI00XN095		    SPI00XN026
case 1143:       mnemonic="SPI00XN071"; break;       //  SPI00XN071		    SPI00XN050
case 1144:       mnemonic="SPI00XN047"; break;       //  SPI00XN047		    SPI00XN074
case 1145:       mnemonic="SPI00XN023"; break;       //  SPI00XN023		    SPI00XN098
case 1146:       mnemonic="SPI00XN010"; break;       //  SPI00XN010		    SPI00XN109
case 1147:       mnemonic="SPI00XN034"; break;       //  SPI00XN034		    SPI00XN085
case 1148:       mnemonic="SPI00XN058"; break;       //  SPI00XN058		    SPI00XN061
case 1149:       mnemonic="SPI00XN082"; break;       //  SPI00XN082		    SPI00XN037
case 1150:       mnemonic="SPI00XN013"; break;       //  SPI00XN013		    SPI00XN106
case 1151:       mnemonic="SPI00XN002"; break;       //  SPI00XN002		    SPI00XN119
case 1152:       mnemonic="SPI00XN026"; break;       //  SPI00XN026		    SPI00XN095
case 1153:       mnemonic="SPI00XN050"; break;       //  SPI00XN050		    SPI00XN071
case 1154:       mnemonic="SPI00XN074"; break;       //  SPI00XN074		    SPI00XN047
case 1155:       mnemonic="SPI00XN098"; break;       //  SPI00XN098		    SPI00XN023
case 1156:       mnemonic="SPI00XN109"; break;       //  SPI00XN109		    SPI00XN010
case 1157:       mnemonic="SPI00XN085"; break;       //  SPI00XN085		    SPI00XN034
case 1158:       mnemonic="SPI00XN061"; break;       //  SPI00XN061		    SPI00XN058
case 1159:       mnemonic="SPI00XN037"; break;       //  SPI00XN037		    SPI00XN082
case 1160:       mnemonic="SPI00XN000"; break;       //  SPI00XN000		    SPI00XN117
case 1161:       mnemonic="SPI00XN024"; break;       //  SPI00XN024		    SPI00XN093
case 1162:       mnemonic="SPI00XN048"; break;       //  SPI00XN048		    SPI00XN069
case 1163:       mnemonic="SPI00XN072"; break;       //  SPI00XN072		    SPI00XN045
case 1164:       mnemonic="SPI00XN096"; break;       //  SPI00XN096		    SPI00XN021
case 1165:       mnemonic="SPI00XN012"; break;       //  SPI00XN012		    SPI00XN105
case 1166:       mnemonic="SPI00XN036"; break;       //  SPI00XN036		    SPI00XN081
case 1167:       mnemonic="SPI00XN060"; break;       //  SPI00XN060		    SPI00XN057
case 1168:       mnemonic="SPI00XN084"; break;       //  SPI00XN084		    SPI00XN033
case 1169:       mnemonic="SPI00XN108"; break;       //  SPI00XN108		    SPI00XN009
case 1170:       mnemonic="SPI00XN001"; break;       //  SPI00XN001		    SPI00XN118
case 1171:       mnemonic="SPI00XN014"; break;       //  SPI00XN014		    SPI00XN107
case 1172:       mnemonic="SPI00XN038"; break;       //  SPI00XN038		    SPI00XN083
case 1173:       mnemonic="SPI00XN062"; break;       //  SPI00XN062		    SPI00XN059
case 1174:       mnemonic="SPI00XN086"; break;       //  SPI00XN086		    SPI00XN035
case 1175:       mnemonic="SPI00XN110"; break;       //  SPI00XN110		    SPI00XN011
case 1176:       mnemonic="SPI00XN097"; break;       //  SPI00XN097		    SPI00XN022
case 1177:       mnemonic="SPI00XN073"; break;       //  SPI00XN073		    SPI00XN046
case 1178:       mnemonic="SPI00XN049"; break;       //  SPI00XN049		    SPI00XN070
case 1179:       mnemonic="SPI00XN025"; break;       //  SPI00XN025		    SPI00XN094

	}	     
return mnemonic;
}




    

   
