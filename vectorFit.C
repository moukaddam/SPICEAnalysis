#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>

#include "TFile.h"
#include "TAxis.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TString.h"

using namespace std;

vector<Double_t> fChannel;
vector<Double_t> fChannelError;
vector<Double_t> fEnergy;
vector<Double_t> fEnergyError;

vector<TString> fSegment;

vector<Double_t> fReducedChi;

// Declaring functions
Int_t fAddToVectors(Double_t inChannel, Double_t inChannelError, Double_t inEnergy, Double_t inEnergyError);
Int_t fClearVectors();
Int_t fPlot(TString outFileName, TString sFitFunc);
Int_t fReturnStatistics();
TString fGetMnemonic(TString Segment);
TString fGetAddress(TString Segment);

TFile outHisto("vectorHisto.root", "recreate");

// Main
void vectorFit(TString inFileName="", TString outFileName="", TString sFitFunc=""){

  // Variable Declaration
  Double_t inChannel, inChannelError, inEnergy, inEnergyError;
  string inString;

  // Create output file
  ofstream outFile(outFileName);
  outFile.close();

  // Read input file
  ifstream inFile(inFileName);  
  while (true) {
    getline(inFile, inString);
    // Search for segment start and get the segment number
    if (inString.substr(0, 7) == "Segment")
      fSegment.push_back( inString.substr(8, 4) );
    else if (inString.substr(0, 1) == "=") {
      // .. if the calibration info has ended, plot the data and clear vectors
      fPlot(outFileName, sFitFunc);  
      fClearVectors();
    } else {
      // .. read in each line of the calibration info and push to vector
      istringstream inStringStream(inString);
      inStringStream >> inChannel >> inChannelError >> inEnergy >> inEnergyError;
      if (inFile.eof()) break;
      fAddToVectors(inChannel, inChannelError, inEnergy, inEnergyError);
    }
  }
  inFile.close();  

  fReturnStatistics();

}

Int_t fAddToVectors(Double_t inChannel, Double_t inChannelError, Double_t inEnergy, Double_t inEnergyError) {

  fChannel.push_back(inChannel);
  fChannelError.push_back(inChannelError);
  fEnergy.push_back(inEnergy);
  fEnergyError.push_back(inEnergyError);

  return 0;
}

Int_t fClearVectors() {

  fChannel.clear();
  fChannelError.clear();
  fEnergy.clear();
  fEnergyError.clear();

  return 0;
}

Int_t fPlot(TString outFileName, TString sFitFunc) {

  TString sSegmentName = fSegment[fSegment.size()-1];  
  TString sMnemonic = fGetMnemonic(sSegmentName);
  TString sAddress = fGetAddress(sSegmentName);

  // Assign vectors to TGraphErrors and fit, according to user specified function
  TGraphErrors *sGraphErrors = new TGraphErrors(fChannel.size(), &fChannel[0], &fEnergy[0], &fChannelError[0], &fEnergyError[0]);
  TString sGraphErrorsName = "Chrg"+sSegmentName;
  TString sGraphErrorsTitle = sGraphErrorsName+" :: "+sMnemonic+" :: "+sAddress;
  sGraphErrors->SetName(sGraphErrorsName);
  sGraphErrors->SetTitle(sGraphErrorsTitle);
  sGraphErrors->GetXaxis()->SetTitle("Channel");
  sGraphErrors->GetYaxis()->SetTitle("Energy [keV]");
  TF1 *sFit = new TF1("sFit", sFitFunc);
  sGraphErrors->Fit(sFit, "QMN");
  outHisto.WriteTObject(sGraphErrors);

  // Retrieve parameters of fit
  vector<Double_t> sParameters;
  for (Int_t i=0; i<sFit->GetNpar(); i++) 
    sParameters.push_back(sFit->GetParameter(i));
  Double_t sReducedChi = sFit->GetChisquare()/sFit->GetNDF();
  fReducedChi.push_back(sReducedChi);

  // Output parameters to file
  ofstream outFile(outFileName, ios::app);
  outFile << sMnemonic << " {" << endl << "\t//" << sSegmentName << endl;
  outFile << "\tADDRESS:" << sAddress << endl << "\tENG_COEFF: ";
  for (Int_t i=0; i<sFit->GetNpar(); i++)
    outFile << sParameters[i] << " ";
  outFile << endl << "}" << endl << endl << endl;
  outFile.close();
  
  return 0;
}

Int_t fReturnStatistics() {

  // Get the average reduced chi-squared
  Double_t sTotalChi = 0.;
  Double_t sSegments = fReducedChi.size();
  for (Int_t i=0; i<(sSegments-1); i++)
    sTotalChi += fReducedChi[i];

  cout << endl << "=================" << endl << "= Fit Complete! = " << endl << "=================" << endl << endl;
  cout << "Segments calibrated: " << sSegments << endl;
  cout << "Average reduced-chi-squared: " << sTotalChi/fReducedChi.size() << endl;
  cout << endl;

  return 0;

}

TString fGetMnemonic(TString Segment) {

  TString sMnemonic = "";  
  int sSegment = atoi(Segment);

  switch (sSegment) {
    case 1060: sMnemonic = "SPI00XN019"; break;
		case 1061: sMnemonic = "SPI00XN008"; break;
		case 1062: sMnemonic = "SPI00XN032"; break;
		case 1063: sMnemonic = "SPI00XN056"; break;
		case 1064: sMnemonic = "SPI00XN080"; break;
		case 1065: sMnemonic = "SPI00XN104"; break;
		case 1066: sMnemonic = "SPI00XN115"; break;
		case 1067: sMnemonic = "SPI00XN091"; break;
		case 1068: sMnemonic = "SPI00XN067"; break;
		case 1069: sMnemonic = "SPI00XN043"; break;
		case 1070: sMnemonic = "SPI00XN006"; break;
		case 1071: sMnemonic = "SPI00XN030"; break;
		case 1072: sMnemonic = "SPI00XN054"; break;
		case 1073: sMnemonic = "SPI00XN078"; break;
		case 1074: sMnemonic = "SPI00XN102"; break;
		case 1075: sMnemonic = "SPI00XN018"; break;
		case 1076: sMnemonic = "SPI00XN042"; break;
		case 1077: sMnemonic = "SPI00XN066"; break;
		case 1078: sMnemonic = "SPI00XN090"; break;
		case 1079: sMnemonic = "SPI00XN114"; break;
		case 1080: sMnemonic = "SPI00XN007"; break;
		case 1081: sMnemonic = "SPI00XN020"; break;
		case 1082: sMnemonic = "SPI00XN044"; break;
		case 1083: sMnemonic = "SPI00XN068"; break;
		case 1084: sMnemonic = "SPI00XN092"; break;
		case 1085: sMnemonic = "SPI00XN116"; break;
		case 1086: sMnemonic = "SPI00XN103"; break;
		case 1087: sMnemonic = "SPI00XN079"; break;
		case 1088: sMnemonic = "SPI00XN055"; break;
		case 1089: sMnemonic = "SPI00XN031"; break;
		case 1090: sMnemonic = "SPI00XN112"; break;
		case 1091: sMnemonic = "SPI00XN101"; break;
		case 1092: sMnemonic = "SPI00XN077"; break;
		case 1093: sMnemonic = "SPI00XN053"; break;
		case 1094: sMnemonic = "SPI00XN029"; break;
		case 1095: sMnemonic = "SPI00XN005"; break;
		case 1096: sMnemonic = "SPI00XN016"; break;
		case 1097: sMnemonic = "SPI00XN040"; break;
		case 1098: sMnemonic = "SPI00XN064"; break;
		case 1099: sMnemonic = "SPI00XN088"; break;
		case 1100: sMnemonic = "SPI00XN099"; break;
		case 1101: sMnemonic = "SPI00XN075"; break;
		case 1102: sMnemonic = "SPI00XN051"; break;
		case 1103: sMnemonic = "SPI00XN027"; break;
		case 1104: sMnemonic = "SPI00XN003"; break;
		case 1105: sMnemonic = "SPI00XN111"; break;
		case 1106: sMnemonic = "SPI00XN087"; break;
		case 1107: sMnemonic = "SPI00XN063"; break;
		case 1108: sMnemonic = "SPI00XN039"; break;
		case 1109: sMnemonic = "SPI00XN015"; break;
		case 1110: sMnemonic = "SPI00XN100"; break;
		case 1111: sMnemonic = "SPI00XN113"; break;
		case 1112: sMnemonic = "SPI00XN089"; break;
		case 1113: sMnemonic = "SPI00XN065"; break;
		case 1114: sMnemonic = "SPI00XN041"; break;
		case 1115: sMnemonic = "SPI00XN017"; break;
		case 1116: sMnemonic = "SPI00XN004"; break;
		case 1117: sMnemonic = "SPI00XN028"; break;
		case 1118: sMnemonic = "SPI00XN052"; break;
		case 1119: sMnemonic = "SPI00XN076"; break;
		case 1120: sMnemonic = "SPI00XN118"; break;
		case 1121: sMnemonic = "SPI00XN107"; break;
		case 1122: sMnemonic = "SPI00XN083"; break;
		case 1123: sMnemonic = "SPI00XN059"; break;
		case 1124: sMnemonic = "SPI00XN035"; break;
		case 1125: sMnemonic = "SPI00XN011"; break;
		case 1126: sMnemonic = "SPI00XN022"; break;
		case 1127: sMnemonic = "SPI00XN046"; break;
		case 1128: sMnemonic = "SPI00XN070"; break;
		case 1129: sMnemonic = "SPI00XN094"; break;
		case 1130: sMnemonic = "SPI00XN105"; break;
		case 1131: sMnemonic = "SPI00XN081"; break;
		case 1132: sMnemonic = "SPI00XN057"; break;
		case 1133: sMnemonic = "SPI00XN033"; break;
		case 1134: sMnemonic = "SPI00XN009"; break;
		case 1135: sMnemonic = "SPI00XN117"; break;
		case 1136: sMnemonic = "SPI00XN093"; break;
		case 1137: sMnemonic = "SPI00XN069"; break;
		case 1138: sMnemonic = "SPI00XN045"; break;
		case 1139: sMnemonic = "SPI00XN021"; break;
		case 1140: sMnemonic = "SPI00XN106"; break;
		case 1141: sMnemonic = "SPI00XN119"; break;
		case 1142: sMnemonic = "SPI00XN095"; break;
		case 1143: sMnemonic = "SPI00XN071"; break;
		case 1144: sMnemonic = "SPI00XN047"; break;
		case 1145: sMnemonic = "SPI00XN023"; break;
		case 1146: sMnemonic = "SPI00XN010"; break;
		case 1147: sMnemonic = "SPI00XN034"; break;
		case 1148: sMnemonic = "SPI00XN058"; break;
		case 1149: sMnemonic = "SPI00XN082"; break;
		case 1150: sMnemonic = "SPI00XN013"; break;
		case 1151: sMnemonic = "SPI00XN002"; break;
		case 1152: sMnemonic = "SPI00XN026"; break;
		case 1153: sMnemonic = "SPI00XN050"; break;
		case 1154: sMnemonic = "SPI00XN074"; break;
		case 1155: sMnemonic = "SPI00XN098"; break;
		case 1156: sMnemonic = "SPI00XN109"; break;
		case 1157: sMnemonic = "SPI00XN085"; break;
		case 1158: sMnemonic = "SPI00XN061"; break;
		case 1159: sMnemonic = "SPI00XN037"; break;
		case 1160: sMnemonic = "SPI00XN000"; break;
		case 1161: sMnemonic = "SPI00XN024"; break;
		case 1162: sMnemonic = "SPI00XN048"; break;
		case 1163: sMnemonic = "SPI00XN072"; break;
		case 1164: sMnemonic = "SPI00XN096"; break;
		case 1165: sMnemonic = "SPI00XN012"; break;
		case 1166: sMnemonic = "SPI00XN036"; break;
		case 1167: sMnemonic = "SPI00XN060"; break;
		case 1168: sMnemonic = "SPI00XN084"; break;
		case 1169: sMnemonic = "SPI00XN108"; break;
		case 1170: sMnemonic = "SPI00XN001"; break;
		case 1171: sMnemonic = "SPI00XN014"; break;
		case 1172: sMnemonic = "SPI00XN038"; break;
		case 1173: sMnemonic = "SPI00XN062"; break;
		case 1174: sMnemonic = "SPI00XN086"; break;
		case 1175: sMnemonic = "SPI00XN110"; break;
		case 1176: sMnemonic = "SPI00XN097"; break;
		case 1177: sMnemonic = "SPI00XN073"; break;
		case 1178: sMnemonic = "SPI00XN049"; break;
		case 1179: sMnemonic = "SPI00XN025"; break;
    default: sMnemonic = "";
  }

  return sMnemonic;

}

TString fGetAddress(TString Segment) {

  TString sAddress = "";  
  int sSegment = atoi(Segment);

  switch (sSegment) {
		case 1060: sAddress = "0x0800100"; break;
		case 1061: sAddress = "0x0800101"; break;
		case 1062: sAddress = "0x0800102"; break;
		case 1063: sAddress = "0x0800103"; break;
		case 1064: sAddress = "0x0800104"; break;
		case 1065: sAddress = "0x0800105"; break;
		case 1066: sAddress = "0x0800106"; break;
		case 1067: sAddress = "0x0800107"; break;
		case 1068: sAddress = "0x0800108"; break;
		case 1069: sAddress = "0x0800109"; break;
		case 1070: sAddress = "0x0800200"; break;
		case 1071: sAddress = "0x0800201"; break;
		case 1072: sAddress = "0x0800202"; break;
		case 1073: sAddress = "0x0800203"; break;
		case 1074: sAddress = "0x0800204"; break;
		case 1075: sAddress = "0x0800205"; break;
		case 1076: sAddress = "0x0800206"; break;
		case 1077: sAddress = "0x0800207"; break;
		case 1078: sAddress = "0x0800208"; break;
		case 1079: sAddress = "0x0800209"; break;
		case 1080: sAddress = "0x0800300"; break;
		case 1081: sAddress = "0x0800301"; break;
		case 1082: sAddress = "0x0800302"; break;
		case 1083: sAddress = "0x0800303"; break;
		case 1084: sAddress = "0x0800304"; break;
		case 1085: sAddress = "0x0800305"; break;
		case 1086: sAddress = "0x0800306"; break;
		case 1087: sAddress = "0x0800307"; break;
		case 1088: sAddress = "0x0800308"; break;
		case 1089: sAddress = "0x0800309"; break;
		case 1090: sAddress = "0x0800400"; break;
		case 1091: sAddress = "0x0800401"; break;
		case 1092: sAddress = "0x0800402"; break;
		case 1093: sAddress = "0x0800403"; break;
		case 1094: sAddress = "0x0800404"; break;
		case 1095: sAddress = "0x0800405"; break;
		case 1096: sAddress = "0x0800406"; break;
		case 1097: sAddress = "0x0800407"; break;
		case 1098: sAddress = "0x0800408"; break;
		case 1099: sAddress = "0x0800409"; break;
		case 1100: sAddress = "0x0800500"; break;
		case 1101: sAddress = "0x0800501"; break;
		case 1102: sAddress = "0x0800502"; break;
		case 1103: sAddress = "0x0800503"; break;
		case 1104: sAddress = "0x0800504"; break;
		case 1105: sAddress = "0x0800505"; break;
		case 1106: sAddress = "0x0800506"; break;
		case 1107: sAddress = "0x0800507"; break;
		case 1108: sAddress = "0x0800508"; break;
		case 1109: sAddress = "0x0800509"; break;
		case 1110: sAddress = "0x0800600"; break;
		case 1111: sAddress = "0x0800601"; break;
		case 1112: sAddress = "0x0800602"; break;
		case 1113: sAddress = "0x0800603"; break;
		case 1114: sAddress = "0x0800604"; break;
		case 1115: sAddress = "0x0800605"; break;
		case 1116: sAddress = "0x0800606"; break;
		case 1117: sAddress = "0x0800607"; break;
		case 1118: sAddress = "0x0800608"; break;
		case 1119: sAddress = "0x0800609"; break;
		case 1120: sAddress = "0x0800700"; break;
		case 1121: sAddress = "0x0800701"; break;
		case 1122: sAddress = "0x0800702"; break;
		case 1123: sAddress = "0x0800703"; break;
		case 1124: sAddress = "0x0800704"; break;
		case 1125: sAddress = "0x0800705"; break;
		case 1126: sAddress = "0x0800706"; break;
		case 1127: sAddress = "0x0800707"; break;
		case 1128: sAddress = "0x0800708"; break;
		case 1129: sAddress = "0x0800709"; break;
		case 1130: sAddress = "0x0800800"; break;
		case 1131: sAddress = "0x0800801"; break;
		case 1132: sAddress = "0x0800802"; break;
		case 1133: sAddress = "0x0800803"; break;
		case 1134: sAddress = "0x0800804"; break;
		case 1135: sAddress = "0x0800805"; break;
		case 1136: sAddress = "0x0800806"; break;
		case 1137: sAddress = "0x0800807"; break;
		case 1138: sAddress = "0x0800808"; break;
		case 1139: sAddress = "0x0800809"; break;
		case 1140: sAddress = "0x0800900"; break;
		case 1141: sAddress = "0x0800901"; break;
		case 1142: sAddress = "0x0800902"; break;
		case 1143: sAddress = "0x0800903"; break;
		case 1144: sAddress = "0x0800904"; break;
		case 1145: sAddress = "0x0800905"; break;
		case 1146: sAddress = "0x0800906"; break;
		case 1147: sAddress = "0x0800907"; break;
		case 1148: sAddress = "0x0800908"; break;
		case 1149: sAddress = "0x0800909"; break;
		case 1150: sAddress = "0x0801000"; break;
		case 1151: sAddress = "0x0801001"; break;
		case 1152: sAddress = "0x0801002"; break;
		case 1153: sAddress = "0x0801003"; break;
		case 1154: sAddress = "0x0801004"; break;
		case 1155: sAddress = "0x0801005"; break;
		case 1156: sAddress = "0x0801006"; break;
		case 1157: sAddress = "0x0801007"; break;
		case 1158: sAddress = "0x0801008"; break;
		case 1159: sAddress = "0x0801009"; break;
		case 1160: sAddress = "0x0801100"; break;
		case 1161: sAddress = "0x0801101"; break;
		case 1162: sAddress = "0x0801102"; break;
		case 1163: sAddress = "0x0801103"; break;
		case 1164: sAddress = "0x0801104"; break;
		case 1165: sAddress = "0x0801105"; break;
		case 1166: sAddress = "0x0801106"; break;
		case 1167: sAddress = "0x0801107"; break;
		case 1168: sAddress = "0x0801108"; break;
		case 1169: sAddress = "0x0801109"; break;
		case 1170: sAddress = "0x0801200"; break;
		case 1171: sAddress = "0x0801201"; break;
		case 1172: sAddress = "0x0801202"; break;
		case 1173: sAddress = "0x0801203"; break;
		case 1174: sAddress = "0x0801204"; break;
		case 1175: sAddress = "0x0801205"; break;
		case 1176: sAddress = "0x0801206"; break;
		case 1177: sAddress = "0x0801207"; break;
		case 1178: sAddress = "0x0801208"; break;
		case 1179: sAddress = "0x0801209"; break;

    default: sAddress = "";
  }

  return sAddress;

}