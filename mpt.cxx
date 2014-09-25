// g++ mpt.cxx /opt/GRSISpoon/libraries/libFormat.so  /opt/GRSISpoon/libraries/libCalManager.so /opt/GRSISpoon/libraries/libRootIOManager.so /opt/GRSISpoon/libraries/libTigress.so /opt/GRSISpoon/libraries/libS3.so  /opt/GRSISpoon/libraries/libSiLi.so /opt/GRSISpoon/libraries/libCSM.so -I/opt/GRSISpoon/include --std=c++0x -o mpt  -O2 `root-config --cflags --libs` -lTreePlayer -lgsl -lgslcblas 


#include <iostream>
//#include <fstream>
#include <unordered_set>
#include <vector>

using namespace std;

#include <stdio.h>
#include <cstdlib>

#include <TFile.h>
#include <TStopwatch.h>
#include <TTree.h>
#include <TTreeIndex.h>
#include <TChain.h>


#include "CalibrationManager.h"
#include "TTigFragment.h"
#include "TSiLi.h"
#include "TS3.h"
#include "TCSM.h"
#include "TTigress.h"
#include "TChannel.h"
//#include "TRFitter.h"

TTigFragment *f = new TTigFragment();


//                          10000000000
//const size_t MEM_SIZE = (size_t)1024 *(size_t)1024*(size_t)1024*size_t(20) ; 
const long int MEM_SIZE = 20000000000;
//ofstream logfile("out_log.txt");
						  
TFile 		*outfile 	= 0;
TTree			*outtree	= 0;
TTigress	*tigress	= 0;
TSiLi *sili = 0;
TS3 *s3 = 0;

TCSM *csm = 0;
//TRf				*rf			=	0;
CalibrationManager *cal = CalibrationManager::instance();

TStopwatch w;

float integration = 125.0;

struct Mnemonic	{
  int arrayposition;
  int	segment;
  std::string system;
  std::string subsystem;
  std::string arraysubposition;
  std::string collectedcharge;
  std::string outputsensor;
};


void ParseMnemonic(std::string *name,Mnemonic *mnemonic)	{
  std::string buf;
  mnemonic->system.assign(*name,0,2);
  mnemonic->subsystem.assign(*name,2,1);
  buf.clear();	buf.assign(*name,3,2);
  mnemonic->arrayposition = atoi(buf.c_str());
  mnemonic->arraysubposition.assign(*name,5,1);
  mnemonic->collectedcharge.assign(*name,6,1);
  buf.clear(); buf.assign(*name,7,2);
  mnemonic->segment = atoi(buf.c_str());
  mnemonic->outputsensor.assign(*name,9,1);
}

//int totalevents = 0;
//int goodevents  = 0;
//int badevents   =	0;



Mnemonic mnemonic;
unsigned short color = 5;


////////////////////////////////////////////////////////////////////////////////
void ProcessEvent(std::vector<TTigFragment> *ev, TList *outlist){
		
  for (int i=0; i < ev->size(); i++) {        	
		
	// Int_t slave = ((ev->at(i).ChannelRaw & 0x00F00000) >> 20);
	// Int_t port  = ((ev->at(i).ChannelRaw & 0x00000F00) >> 8);
	// Int_t chan  =  (ev->at(i).ChannelRaw & 0x000000FF);
    
    std::string name	=	ev->at(i).ChannelName;
    if(ev->at(i).ChannelName.length() < 10) {
	  printf("name too short!! found = %s\n",name.c_str());
	  continue;
 	}
    ParseMnemonic(&name,&mnemonic);

    float energy	= ev->at(i).ChargeCal;
    int midasid 	= ev->at(i).MidasId;
    int charge 		= ev->at(i).Charge;
    int timecfd 	= ev->at(i).Cfd;
    int timeled 	= ev->at(i).Led;
    int time 		= ev->at(i).TimeToTrig;
/*
	printf(DBLUE "Name: %s" RESET_COLOR "\n",name.c_str());

	printf(DRED "\tSystem:      %s" RESET_COLOR "\n",mnemonic.system.c_str());
	printf(DRED "\tSubsystem:   %s" RESET_COLOR "\n",mnemonic.subsystem.c_str());
	printf(DRED "\tArrayPos:    %i" RESET_COLOR "\n",mnemonic.arrayposition);
	printf(DRED "\tArraySubPos: %s" RESET_COLOR "\n",mnemonic.arraysubposition.c_str());
	printf(DRED "\tSegment:     %i" RESET_COLOR "\n",mnemonic.segment);
	printf(DRED "\tCollectChg:  %s" RESET_COLOR "\n",mnemonic.collectedcharge.c_str());
	printf(DRED "\tOutputSen:   %s" RESET_COLOR "\n",mnemonic.outputsensor.c_str());
	ev->at(i).Print();
	printf(DYELLOW "=====================================================" RESET_COLOR  "\n");
 */ 
    if(mnemonic.system.compare(0,2 ,"CS")==0) { //csm data
		unsigned short pos = 0xffff;
		if(mnemonic.arraysubposition.compare(0,1,"D")==0) {
			pos = 0;
		} else if(mnemonic.arraysubposition.compare(0,1,"E")==0) {
			pos = 1;
		} else if(mnemonic.arraysubposition.compare(0,1,"F")==0) {
			pos = 2;
		}           
		if(mnemonic.collectedcharge.compare(0,1,"N")==0 ) {  // Horizontal Strips (aka "front" strips);
			 csm->SetHorizontal(&ev->at(i),(unsigned short)mnemonic.arrayposition,pos,(unsigned short)mnemonic.segment);
			 /*
				printf(DYELLOW "\tdet = %i" RESET_COLOR "\n", (unsigned short)mnemonic.arrayposition);
				printf(DYELLOW "\tpos = %i" RESET_COLOR "\n", pos);
				printf(DYELLOW "\tstr = %i" RESET_COLOR "\n", (unsigned short)mnemonic.segment);
			 */             
		} else if(mnemonic.collectedcharge.compare(0,1,"P")==0 ) {  // Vertical Strips (aka "back" strips);
			 csm->SetVertical(&ev->at(i),(unsigned short)mnemonic.arrayposition,pos,(unsigned short)mnemonic.segment);
			 /*
				printf(DGREEN "\tdet = %i" RESET_COLOR "\n", (unsigned short)mnemonic.arrayposition);
				printf(DGREEN "\tpos = %i" RESET_COLOR "\n", pos);
				printf(DGREEN "\tstr = %i" RESET_COLOR "\n", (unsigned short)mnemonic.segment);
			 */
		}
	 } 

	else if(mnemonic.system.compare(0,2,"SP")==0)	{ //spice data
		if(mnemonic.subsystem.compare(0,1,"I")==0) { 
			//printf("I AM HERE.\n");
			int tmpseg = mnemonic.segment * 10 + atoi(mnemonic.outputsensor.c_str());
			sili->SetSiLi(tmpseg,energy,timecfd,timeled,time,charge);
		} else if (mnemonic.subsystem.compare(0,1,"E")==0) {  


			if(mnemonic.collectedcharge.compare(0,1,"N")==0) {
				s3->SetSector(mnemonic.arrayposition,mnemonic.segment,energy,timecfd,timeled,time,charge);
         } else if(mnemonic.collectedcharge.compare(0,1,"P")==0) {
				s3->SetRing(mnemonic.arrayposition,mnemonic.segment,energy,timecfd,timeled,time,charge);
         }
		}
	}








	else if(mnemonic.system.compare(0,2,"TI")==0)	{ //tigress data
		color = 5;
		char asp = mnemonic.arraysubposition.c_str()[0];
		switch(asp){
			case 'B':
				color = 0;
				break;
			case 'G':
				color = 1;
				break;					
			case 'R':
				color = 2;
				break;					
			case 'W':
				color = 3;
				break;
			default:
				break;
      };
	  if(mnemonic.subsystem.compare(0,1,"G")==0)	{
	  	if(mnemonic.segment==0 ) {//||mnemonic.segment==9)	{
			if(ev->at(i).ChannelNumber % 10 == 0) {//mnemonic.segment==0)	{ //////////////////////////////////  punting.
 				tigress->SetCore(&ev->at(i),(unsigned short)mnemonic.arrayposition,color);
			}                         //////////////////////////////////
		}
		else	{
			tigress->SetSegment(&ev->at(i),(unsigned short)mnemonic.arrayposition,(unsigned short)color,(unsigned short)mnemonic.segment);
		}
      }                
      else if(mnemonic.subsystem == 'S'){
			tigress->SetBGO(&ev->at(i),(unsigned short)mnemonic.arrayposition,(unsigned short)color,(unsigned short)mnemonic.segment);	
      }
    } 
    //else if(mnemonic.system == "RF")	{
    //  rf->SetRf(&ev->at(i));
  	//} 
  }

	//sharc->BuildHits("CLEAN");
	//csm->BuildHits();
	s3->BuildHits();
	sili->BuildHits();
	tigress->BuildHits();
//	for(int x = 0; x < tigress->GetMultiplicity();x++)	{
//		TTigressHit *thit = tigress->GetTigressHit(x);	
//		char buffer[64];
//		switch(thit->GetCrystalNumber())	{
//			case 0:
//				sprintf(buffer,"TIG%02iB",thit->GetDetectorNumber());
//				break;
//			case 1:
//				sprintf(buffer,"TIG%02iG",thit->GetDetectorNumber());
//				break;
//			case 2:
//				sprintf(buffer,"TIG%02iR",thit->GetDetectorNumber());
//				break;
//			case 3:
//				sprintf(buffer,"TIG%02iW",thit->GetDetectorNumber());
//				break;
//			default:
//				printf("unknown crystal %i fopund.\n",thit->GetCrystalNumber());
//				break;
//		}
//		TChannel *chan = cal->GetChannel(buffer);
//		if(chan)	{
			//printf("Here!!\t%s\n",buffer);
			//chan->Print();
//			thit->GetCore()->SetEnergy(chan->CalibrateENG(thit->GetCore()->GetCharge()));
//		}
//		TTigress::DopplerCorrect(thit);
//	}

	tigress->BuildAddBack();
	for(int x = 0; x <   tigress->GetAddBackMultiplicity();x++)	{
		TTigress::DopplerCorrect(tigress->GetAddBackHit(x));
	}

//	printf(DGREEN "================================================" RESET_COLOR  "\n");
//	printf(DGREEN "================================================" RESET_COLOR  "\n");
//	printf(DGREEN "================================================" RESET_COLOR  "\n");
	
	//printf("sili hits:  %i\n",sili->GetMultiplicity());
  
	outtree->Fill();
	
//	tf->Clear();
// 	rf->Clear();	
  s3->Clear();
  sili->Clear();		  
  tigress->Clear();
  //csm->Clear();
 return;
}


////////////////////////////////////////////////////////////////////////////////
void BuildEvents(TChain *chain, TList *outlist){  

  TTigFragment *pFrag = 0;

  TStopwatch w;
  w.Start();
  int ntrees = chain->GetNtrees();
  int nChainEntries = chain->GetEntries();
  int treeNumber, lastTreeNumber = -1;

  chain->SetMaxVirtualSize(MEM_SIZE);

  for(int i=0;i<nChainEntries;i++) {
    chain->LoadTree(i);
    treeNumber = chain->GetTreeNumber();

    if (treeNumber != lastTreeNumber) {
      printf("Changing to tree number %d from %d at chain entry number %d.\n",
	     treeNumber, lastTreeNumber, i);		  
	     lastTreeNumber = treeNumber;
    } else {
      continue;
    }

   TTree *tree = chain->GetTree();
//	TTree *tree = (TTree*)tree_disk->Clone(0);
	printf("tree->GetEntries() = %i\n",tree->GetEntries());
	int nentries = tree->GetEntries();


	
	if(!tree->GetTreeIndex())	{
		printf("\nTree Index not found, Building index...");
		fflush(stdout);
		tree->BuildIndex("TriggerId","FragmentId");	
		printf("  done.\n");
		fflush(stdout);
	}

	tree->SetCacheSize(MEM_SIZE);
	tree->SetCacheLearnEntries(5);	

    TBranch *branch = tree->GetBranch("TTigFragment");
    branch->SetAddress(&pFrag);
    tree->SetMaxVirtualSize(MEM_SIZE);
	 int xx = tree->LoadBaskets(MEM_SIZE);
//    tree->LoadBaskets(20000000000); 
//    int xx = branch->LoadBaskets();//MEM_SIZE); 
  	 printf("loaded %i baskets\t MEM_SIZE = %zd\n",xx,MEM_SIZE);

    std::vector<TTigFragment> evFrags;
		

	int MinTriggerId = (int) tree->GetMinimum("TriggerId");
	if(MinTriggerId<0)
		MinTriggerId =0;

    int MaxTriggerId = (int) tree->GetMaximum("TriggerId");

    for(int j=MinTriggerId;j<=MaxTriggerId; j++)	{	

      Int_t fragno = 1;
      evFrags.clear();
      
      while (tree->GetEntryWithIndex(j,fragno++) != -1) {
				//pFrag->Print();	
				evFrags.push_back(*pFrag);
      }
	  //do something with the evFrag vector (contains a built event)... ProcessEvent(evFrags); 
	  if(!evFrags.empty())	{
		ProcessEvent(&evFrags,outlist); 
	  }
//	  else	{
//		printf(DYELLOW "\n\tTrigger %i is empty!!!" RESET_COLOR  "\n",j);
//	  }
     if ( (j % 10000) == 0 || j == nChainEntries - 1) {
	  	printf("processing event %i/%i from tree %i/%i\t\ttime since start: %.02f seconds     \n",j+1,MaxTriggerId,treeNumber+1,ntrees,w.RealTime());
		printf("Reading %lld bytes in %d transactions\n",tree->GetCurrentFile()->GetBytesRead(),  tree->GetCurrentFile()->GetReadCalls());
			fflush(stdout);
	    w.Continue();
     }
    }

    branch->DropBaskets("all");
	//tree->Delete();
	if(outtree)	{
		outtree->AutoSave();
		//outtree->Print();
	}
    printf("\n");
    i += (nentries - 10);
  }


}





////////////////////////////////////////////////////////////////////////////////
int main(int argc, char **argv){
  csm = new TCSM();
  tigress	= new TTigress();
  tigress->SetBeta(0.107);
  tigress->Print();
//  rf 		= new TRf();

//  cal->SetCalFileName("/data4/tigress/S1389Work_part2/calibrations/Calibrations3Alpha/test/Jun12_Calibration.cal");

  //cal->ReadCalibrationFile();
  TChain *chain = new TChain("FragmentTree");	

  if(argc <2)	{
	printf("Try  ./mpt fragmentxxxxxx.root instead.\n");
	return 1;
  }

  std::string name = "analysis";
  std::string iname = argv[1];

  // extracts the fragment number and file type -> eg fragment25090.root becomes analysis25090.root
  name.append(iname.substr(iname.find_last_of("fragment",iname.find_last_of('.'))+1,iname.find_last_of('.')));	

  outfile = new TFile(name.c_str(),"recreate");

  outtree = new TTree("AnalysisTree","AnalysisTree");

  //outtree->Branch("TSharc",&sharc);
  outtree->Branch("TTigress",&tigress);
  outtree->Branch("TSiLi",&sili);
  outtree->Branch("TS3",&s3);
 

 // outtree->Branch("TCSM",&csm);
  //outtree->Branch("TRf",&rf);
  //outtree->Branch("TTriFoil",&tf);

  for(int i=1;i<argc;i++)	{  
    printf("Adding tree %s\n", argv[i]);
    chain->Add(argv[i]);
  }

  TList *outlist = new TList();	

  BuildEvents(chain,outlist);
  
  outfile->cd();	
  outlist->Sort();
  outlist->Write();
  outfile->Write();
  outfile->Close();
  
  return EXIT_SUCCESS;
}
