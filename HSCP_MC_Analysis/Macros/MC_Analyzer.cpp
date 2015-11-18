#include <cctype>
#include <cstddef> //size_t
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <string>
#include <cmath>
#include <map>
#include <vector>
#include <array>
#include <iterator>
#include "format.h"
#include <TROOT.h> //for gROOT
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TH1.h>
#include <TH1F.h>
#include <TH2.h>
#include <TH2F.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TMathBase.h>
#include <TMath.h>
#include <TNetFile.h>
#include <TAuthenticate.h>
#include <TObject.h>
#include <TApplication.h>
#include <TNamed.h>
#include <TStyle.h>
#include <TAttLine.h>
#include <TAttFill.h>
#include <TAttMarker.h>
#include <TAttText.h>
#include <TGraph.h>
#include <TMultiGraph.h>
#include <TColor.h>
#include <Rtypes.h>
#include <TVectorD.h>
#include <TAxis.h>
#include <TLegend.h>
#include <THStack.h>
#include <TLine.h>
#include <TF1.h>
#include <TPad.h>

using namespace std;
  
class Key{
public:
  Key(string dist,string type,int charge,int mass){
    this->dist = dist;
    this->charge = charge;
    this->mass = mass;
    this->type = type;
    this->name = fmt::format("{}_{}_{}_{}",type,dist,charge,mass);
  }

  string dist;
  string type;
  int charge;
  int mass;
  string name;
  
  bool operator<(const Key &k) const{
    return (this->name<k.name);
  }
  
  bool operator==(const Key &rhs){
    return this->dist == rhs.dist && this->charge == rhs.charge && this->mass == rhs.mass && this->type == rhs.type;
  }

  void operator=(const Key &aKey){
    this->dist = aKey.dist;
    this->charge = aKey.charge;
    this->mass = aKey.mass;
    this->name = aKey.name;
  }
  
};

struct distObjects{
  TCanvas *canvas;
  TH1F *distribution;
  vector<float> data;
  distObjects() : canvas(NULL), distribution(NULL) {}
  distObjects(const distObjects &arg) : canvas(arg.canvas), distribution(arg.distribution) {}
  //distObjects(TCanvas *&aCanv, TH1F *&aDist) : canvas((TCanvas*)aCanv->Clone()), distribution((TH1F*)aDist->Clone()) {}
  distObjects(TCanvas *&aCanv, TH1F *&aDist) : canvas(NULL), distribution((TH1F*)aDist->Clone()) {}
};

typedef map<Key,distObjects> distMap;
typedef pair<string,string> axes;
typedef pair<double,double> limits;

struct distProp{
  limits axisLimits;
  double nBins;
  axes axisTitles;
};

struct Particle{
  float Pt,P,Theta,Eta,Phi,E,Ih,corrMass,uncorrMass,Beta,TOF;
  int Event,PDG;
  Particle() : Pt(0),P(0),Theta(0),Eta(0),Phi(0),E(0),Ih(0),corrMass(0),uncorrMass(0),Beta(0),TOF(0),Event(-1),PDG(-999) {}
  Particle(const float Pt, const float P, const float Theta, const float Eta, const float Phi, const float E, const float Ih, const float corrMass, const float uncorrMass, const float Beta, const float TOF, const int PDG, const int Event) : Pt(Pt),P(P),Theta(Theta),Eta(Eta),Phi(Phi),E(E),Ih(Ih),corrMass(corrMass),uncorrMass(uncorrMass),Beta(Beta),TOF(TOF),PDG(PDG),Event(Event) {}
  Particle(const Particle &arg) : Pt(arg.Pt),P(arg.P),Theta(arg.Theta),Eta(arg.Eta),Phi(arg.Phi),E(arg.E),Ih(arg.Ih),corrMass(arg.corrMass),uncorrMass(arg.uncorrMass),Beta(arg.Beta),TOF(arg.TOF),PDG(arg.PDG),Event(arg.Event) {}

  //Overload the assignment operator
  void operator=(const Particle &aPart){
    Pt = aPart.Pt;
    P = aPart.P;
    Theta = aPart.Theta;
    Eta = aPart.Eta;
    Phi = aPart.Phi;
    E = aPart.E;
    Ih = aPart.Ih;
    corrMass = aPart.corrMass;
    uncorrMass = aPart.uncorrMass;
    Beta = aPart.Beta;
    TOF = aPart.TOF;    
  }
  
};


// struct Match{
//   Particle *genParticle;
//   Particle *recoParticle;
//   float dR, dEta, dPhi;
//   int event;

//   Match() : genParticle(NULL), recoParticle(NULL), dR(-1), dEta(-1), dPhi(-1) {}
//   Match( const Particle &aGen, const Particle &aReco, const float &a_dR, const float &a_dEta, const float &a_dPhi, const int &a_event ){
//     genParticle = &aGen;
//     recoParticle = &aReco;
//     dR = a_dR;
//     dEta = a_dEta;
//     dPhi = a_dPhi;
//     event = a_event;
//   }
// };

bool isMatch( const Particle &aPart1, const Particle &aPart2 ){
  //Refine this later when you look at the distributions for dR
  const float LIMIT = 0.0005;

  float dR = TMath::Sqrt(TMath::Power(aPart1.Eta-aPart2.Eta,2) + TMath::Power(aPart1.Phi-aPart2.Phi,2));
  if (dR <= LIMIT) return true;
  else return false;
}


typedef map< Key, map< int, vector< Particle > > > Events;


/*
  Function to insert available ROOT MC Sample files into a string array.
  Make the first entry in the array a string containing the number of 
  entries in the array, INCLUDING this entry.
*/

static string *ListOfFiles(const char &CC, const char &CT){
  FILE *inStream;
  char charnFiles[100];
  char charFiles[100000];
  stringstream streamnFiles;
  
  //First, count how many files we will be reading
  string command = fmt::format("ls /home/austin/HSCP_Study/HSCP_MC_Files/CC_{}_CT_{}/Histos_mchamp*.root -l | wc -l",CC,CT);
  if(!(inStream = popen(command.c_str(), "r"))){
    exit(0);
  }
    
  while(fgets(charnFiles, sizeof(charnFiles), inStream)!=NULL){
    streamnFiles << charnFiles;
  }
  pclose(inStream);

  int nFiles = atoi(streamnFiles.str().c_str());
  //cout << "nFiles form function: " << nFiles << endl;
  string *listOfFiles = new string[nFiles+1];

  //Now we have the number of files, and can dynamicall allocate a string array for it.
  //PLUS the first entry containing the number of entries in the array.
  
  listOfFiles[0] = streamnFiles.str();
  
  //Now we need to read in the list of files
  //Will do something similar to the above, but using this method to get only the file names
  stringstream testStream;
  string command2 = fmt::format("find /home/austin/HSCP_Study/HSCP_MC_Files/CC_{}_CT_{}/Histos_mchamp*.root -printf \"%f\n\"",CC,CT);
  if(!(inStream = popen(command2.c_str(),"r"))){
    exit(0);
  }
  while(fgets(charFiles, sizeof(charFiles), inStream)!=NULL){
    testStream << string(charFiles);
  }
  pclose(inStream);
  
  string line;
  int fNum = 1;
  while(getline(testStream,line)){
    listOfFiles[fNum] = line;
    fNum++;
  }
  
  return listOfFiles;
}

/*Function to parse a file name into mass and charge */
//Create a struct to allow returning both mass and charge
//as parsed from the file name
struct NameDat{
  double *charge;
  double *mass;
  string *fileNames;
  map<double,int> chargeCounts;
  map<double,int> massCounts;
};

//Actual function that will parse the file name and return a NameDat struct
static NameDat *FileNameParser( string *Names, const int nFiles ){
  NameDat *outDat = new NameDat; //Struct that will contain the parsed mass and charge
  outDat->charge = new double[nFiles];
  outDat->mass = new double[nFiles];
  outDat->fileNames = Names+1;
  
  //the one comes from the fact that the file
  //name array has the number of files as the
  //first character
  for(int iFile = 0; iFile < nFiles; iFile++){
    string aName = outDat->fileNames[iFile];
    string chunks[4]; //string array that will contain the chunks
    //Loop through each character of the file name, increasing iCh.
    //Increase iChunk at each chunk ('_' character)
    //Store each chunk in a string array (chunks)
    //Clean up the chunks of interest, and return the data in the struct
    size_t found = aName.find_first_of("_");
    
    chunks[0] = aName.substr(0,found);
    int chunkPos = 1;
    while( found != string::npos ){
      chunks[chunkPos] = aName.substr(found+1,aName.find_first_of('_',found+1)-found-1);
      found = aName.find_first_of('_',found+1);
      chunkPos++;
    }
    //Charge is in chunk 1, mass in chunk index 3
    //Clean the mass and charge chunks.
    string charge = chunks[1].substr(chunks[1].find_first_not_of("mchamp"),string::npos);
    string mass = chunks[3].substr(0,chunks[3].find_first_of(".root"));

    //Convert the strings to floats/doubles
    outDat->charge[iFile] = atof(charge.c_str())/3.0;
    outDat->mass[iFile] = atof(mass.c_str());
  }
  
  map<double, int> massCounts;
  map<double, int> chargeCounts;
  
  //Count the masses
  for( int iMass = 0; iMass < nFiles; iMass++ ){
    if( massCounts.count(outDat->mass[iMass]) != 0 ){
      massCounts.find(outDat->mass[iMass])->second += 1;
    }
    else{
      massCounts.emplace( outDat->mass[iMass], 1 );
    }
  }

  //Count the charges
  for( int iCharge = 0; iCharge < nFiles; iCharge++ ){
    if( chargeCounts.count(outDat->charge[iCharge]*3) != 0 ){
      chargeCounts.find(outDat->charge[iCharge]*3)->second += 1;
    }
    else{
      chargeCounts.emplace( outDat->charge[iCharge]*3, 1 );
    }
  }

  outDat->chargeCounts = chargeCounts;
  outDat->massCounts = massCounts;
  
  return outDat;
}


//For each distribution type, go through all files (masses and charges)
static void AllocateDistributions( distMap &argDists, const double *charges, const double *masses, const int &numFiles, const vector<string> &distNames, map<string,distProp> &distProps, const vector<string> &types){
  /*"beta","energy",
    "eta",
    "gamma",
    "Ih",
    "genP",
    "genPt",
    "recoPCorr",
    "recoPUncorr",
    "recoPtCorr",
    "recoPtUncorr",
    "phi",
    "theta",
    "TOF",
    "genMass",
    "recoMassCorr",
    "recoMassUncorr",
    "charge",
  */
  bool isGen = false;
  for( const auto &iType : types){
    for (const auto &iNames : distNames) {
      
      double *lowerLim = &distProps[iNames].axisLimits.first;
      double *upperLim = &distProps[iNames].axisLimits.second;
      double *nBins = &distProps[iNames].nBins;
      string xAxis(distProps[iNames].axisTitles.first);
      string yAxis(distProps[iNames].axisTitles.second);
    
      for (int iFile=0; iFile < numFiles; iFile++) {
        Key entryKey (iNames,iType,(int)(3* charges[iFile]),(int)masses[iFile]);
      
        string canv_name = entryKey.name + "_canv";
        string dist_name = entryKey.name + "_dist";

        yAxis = fmt::format(yAxis.c_str(), (*upperLim-*lowerLim) / *nBins);
        //TCanvas *tempCanv = new TCanvas(canv_name.c_str(),canv_name.c_str(),500,500);
        TCanvas *tempCanv = NULL;
        TH1F *tempHist = new TH1F(dist_name.c_str(),dist_name.c_str(),*nBins,*lowerLim,*upperLim);

        tempHist->GetXaxis()->SetTitle(xAxis.c_str());
        tempHist->GetYaxis()->SetTitle(yAxis.c_str());
      
        argDists.emplace(entryKey,distObjects(tempCanv,tempHist));
      }
    }
  }
}


static void AllocateOutfile( TFile *&aFile, const vector<string> &distNames ){
  aFile->mkdir("Reco");
  aFile->mkdir("Gen");
  aFile->mkdir("Bkg");
  for (const auto &iNames : distNames) {
    aFile->cd("Reco");
    gDirectory->mkdir(iNames.c_str());
    aFile->cd("Gen");
    gDirectory->mkdir(iNames.c_str());
    aFile->cd("Bkg");
    gDirectory->mkdir(iNames.c_str());
  }
  aFile->cd("/");
}

void SetRootStyle(){
  //Useful Style Tips for ROOT
  /*
    http://www.nbi.dk/~petersen/Teaching/Stat2014/PythonRootIntro/ROOT_TipsAndTricks.pdf
  */
  TStyle *aStyle = new TStyle("aStyle","Austin's Root Style");
  aStyle->SetPalette(1,0); //Get rid of terrible default color scheme
  aStyle->SetOptStat(0);
  aStyle->SetOptTitle(0);
  aStyle->SetOptDate(0);
  aStyle->SetLabelSize(0.03,"xyz"); //Axis value font size
  aStyle->SetTitleSize(0.035,"xyz"); //Axis title font size
  aStyle->SetLabelFont(22,"xyz");
  aStyle->SetTitleOffset(1.2,"y");
  //Default canvas options
  aStyle->SetCanvasDefW(500);
  aStyle->SetCanvasDefH(500);
  aStyle->SetCanvasColor(0);
  aStyle->SetCanvasBorderMode(0);
  aStyle->SetPadBottomMargin(0.1);
  aStyle->SetPadTopMargin(0.1);
  aStyle->SetPadLeftMargin(0.1);
  aStyle->SetPadRightMargin(0.1);
  aStyle->SetPadGridX(0);
  aStyle->SetPadGridY(0);
  aStyle->SetPadTickX(1);
  aStyle->SetPadTickY(1);
  aStyle->SetFrameBorderMode(0);
  aStyle->SetPaperSize(20,24); //US Letter
  gROOT->SetStyle("aStyle");
  gROOT->ForceStyle();
}


void WriteStandardDistributions( TFile *&outFile, const vector<string> &types, const distMap &distList){
  
  bool isGen;
  double norm = 1;

  TH1F *currDist;
  
  for (const auto &iTypes : types) {
    if (iTypes == "Gen")
      isGen = true;
    else
      isGen = false;
    for (const auto &distIt : distList){
      string dir = fmt::format("/{}/{}",iTypes,distIt.first.dist);
      if (isGen) {
        if (distIt.first.type != "Gen") continue;
      }
      else {
        if (distIt.first.type != "Reco") continue;
      }
      outFile->cd(dir.c_str());
      currDist = distIt.second.distribution;
      if (distIt.first.dist.find("Mass") == string::npos)//No mass normalization
        currDist->Scale(norm/currDist->Integral("width"));
      //currDist->Draw("P");
      currDist->Write();
    }
  }
}


void WriteGenFracTrackVsBeta(TFile *&outFile, const distMap &distList, const map<double,int> &massCounts, const map<string,distProp> &distProps, const char &CC, const char &CT){
  //Only want to make the plot with unity charges (key for charge = 3)
  const array<int,12> colorList = {1,2,3,4,6,41,34,1,2,3,4};
  const array<int,12> markerList = {2,3,4,5,25,26,27,20,21,22,33,34};

  const double CHARGE = 1.0; //Desired constant charge for the plots, in units of actual charge (not e/3)

  double per = (distProps[string("beta")].axisLimits.second -
                distProps[string("beta")].axisLimits.first)/
    distProps[string("beta")].nBins;
  
  THStack *LowMassOutStack = new THStack("LowMassfracPart_beta","");
  TLegend *LowMassLegend = new TLegend(0.1,0.65,0.5,0.85);
  //LowMassLegend->SetTextFont(72);
  LowMassLegend->SetTextSize(0.025);
  LowMassLegend->SetFillColor(0);
  TCanvas *LowMassCanvas = new TCanvas("low","low",2000,2000);
                                       
  THStack *HighMassOutStack = new THStack("HighMassfracPart_beta","");
  TLegend *HighMassLegend = new TLegend(0.6,0.65,0.8,0.85);
  //LowMassLegend->SetTextFont(72);
  HighMassLegend->SetTextSize(0.025);
  HighMassLegend->SetFillColor(0);
  TCanvas *HighMassCanvas = new TCanvas("high","high",2000,2000);
  
  TH1F *currDist = NULL;
  const double NORM = 1;

  //Make a vector of TGraph objects whose constituents will compose the TMultiGraph that is saved to disc
  vector<TH1F*> LowMassDists;
  vector<TH1F*> HighMassDists;
  vector<string> LowMassNames;
  vector<string> HighMassNames;
    
  //Loop through all available masses
  auto iMass = massCounts.begin();
  auto iColor = colorList.begin();
  auto iMarker = markerList.begin();
  for (; iMass != massCounts.end(); ++iMass){
    bool isLowMass;
    if (iMass->first <=800)
      isLowMass = true;
    else
      isLowMass = false;
    
    //Check to see if this mass is available with the desired charge
    Key betaKey = Key (string("beta"), string("Gen"), (int)(3 * CHARGE), (int)iMass->first);

    //map.count returns zero if the item is not found. want to skip masses that are not with desired charge
    if (distList.count(betaKey) == 0){
      continue;
    }

    currDist = distList.find(betaKey)->second.distribution;
    currDist->Scale(NORM/currDist->Integral());
    currDist->SetMarkerStyle(*iMarker);
    //currDist->SetFillColorAlpha(*iColor,0.75);
    //currDist->SetLineColorAlpha(*iColor,0.75);
    currDist->SetLineColor(*iColor);
    currDist->SetLineWidth(3.0);
    if( isLowMass ){
      LowMassLegend->AddEntry(currDist,(fmt::format("Q: {} M: {} GeV/",(int)(3*CHARGE),(int)iMass->first)+string("c^{2}")).c_str());
      LowMassOutStack->Add(currDist);
    }
    else{
      HighMassLegend->AddEntry(currDist,(fmt::format("Q: {} M: {} GeV/",(int)(3*CHARGE),(int)iMass->first)+string("c^{2}")).c_str());
      HighMassOutStack->Add(currDist);
    }

    ++iColor;
    ++iMarker;
  }//File loop

  LowMassCanvas->cd();
  LowMassOutStack->Draw("NOSTACK");
  LowMassOutStack->GetYaxis()->SetTitle(fmt::format("Fraction of Tracks/{}",per).c_str());
  LowMassOutStack->GetXaxis()->SetTitle("Gen #beta");
  gPad->Update();
  LowMassLegend->Draw("F");
  gPad->Update();
  outFile->cd();
  LowMassCanvas->Write();
  string LMName = fmt::format("LMGenFracTrackVsBeta_CC_{}_CT_{}.pdf",CC,CT);
  LowMassCanvas->Print(LMName.c_str(),"pdf");

  HighMassCanvas->cd();
  HighMassOutStack->Draw("NOSTACK");
  HighMassOutStack->GetYaxis()->SetTitle(fmt::format("Fraction of Tracks/{}",per).c_str());
  HighMassOutStack->GetYaxis()->SetRangeUser(0.0, 0.15);
  HighMassOutStack->GetXaxis()->SetTitle("Gen #beta");
  gPad->Update();
  HighMassLegend->Draw("F");
  gPad->Update();
  outFile->cd();
  HighMassCanvas->Write();
  string HMName = fmt::format("HMGenFracTrackVsBeta_CC_{}_CT_{}.pdf",CC,CT);
  HighMassCanvas->Print(HMName.c_str(),"pdf");
  
}


void WriteIhVsP(TFile *&outFile, const distMap &distList, const map<double,int> &chargeCounts, const Events &events, const char &CC, const char &CT){
  const array<int,12> colorList = {1,2,3,4,6,41,34,46,12,8,14,5};
  const array<int,12> markerList = {2,3,4,5,25,26,27,20,21,22,33,34};

  const int MASS = 900; //Desired mass to vary charges over

  TMultiGraph *recoOutGraph_KnownQ = new TMultiGraph("multrecoknownq","");
  TMultiGraph *recoOutGraph_UnknownQ = new TMultiGraph("multrecounknownq","");
  TMultiGraph *genOutGraph = new TMultiGraph("multgen","");
  
  TGraph* currGraph;
  TLegend *recoDistLegend_KnownQ = new TLegend(0.5,0.55,0.80,0.85);
  TLegend *recoDistLegend_UnknownQ = new TLegend(0.5,0.55,0.80,0.85);
  TLegend *genDistLegend = new TLegend(0.5,0.55,0.80,0.85);
  
  recoDistLegend_KnownQ->SetTextSize(0.025);
  recoDistLegend_KnownQ->SetFillColor(0);
  recoDistLegend_UnknownQ->SetTextSize(0.025);
  recoDistLegend_UnknownQ->SetFillColor(0);
  genDistLegend->SetTextSize(0.025);
  genDistLegend->SetFillColor(0);
 
  vector<string> chargeNames;

  //Check to see if this mass is available with the desired charge
  auto iCharge = chargeCounts.begin();
  auto iColor = colorList.begin();
  for( ; iCharge!=chargeCounts.end(); iCharge++) {
    Key corr_pKey = Key (string("recoPCorr"), string("Reco"), (int)(iCharge->first), MASS);
    Key uncorr_pKey = Key (string("recoPUncorr"), string("Reco"), (int)(iCharge->first), MASS);
    Key reco_ihKey = Key (string("Ih"), string("Reco"), (int)(iCharge->first), MASS);
  
    //Reco will generally have fewer points than gen, so let's check this one
    unsigned int NPOINTS = distList.find(corr_pKey)->second.data.size();
    if ( NPOINTS < 100 ){
      cout << fmt::format("Bad file at mchamp{}_M_{}",iCharge->first,MASS) << endl;
      continue;
    }
    //map.count returns zero if the item is not found. want to skip masses that are not with desired charge  
    if (distList.count(corr_pKey) == 0){
      continue;
    }    
   
    chargeNames.push_back(fmt::format("Q={}",iCharge->first)+string("#frac{e}{3}"));

    /*--------------------------------------------------------------------*/
    //Fill the reco graph, known Q
    currGraph = new TGraph(distList.find(corr_pKey)->second.data.size());
    currGraph->SetMarkerStyle(20);
    currGraph->SetMarkerColorAlpha(*iColor,0.8);
    currGraph->SetMarkerSize(0.3);
    
    //Now we want to loop through all momenta and Ih for this charge and mass
    //Have the data stored in distList
    auto iP = distList.find(corr_pKey)->second.data.begin();
    auto iIh = distList.find(reco_ihKey)->second.data.begin();
    unsigned int i = 0;
    //Correcting momentum here. Assuming correct charge isknown
    //Therefore need to multiply reco by Q in order to get corrected
    for( ; iP != distList.find(corr_pKey)->second.data.end(); ++iP, ++iIh ){
      currGraph->SetPoint(i,(*iP) * (iCharge->first/3.0),*iIh);
      ++i;
    }
    recoOutGraph_KnownQ->Add(currGraph);
    recoDistLegend_KnownQ->AddEntry(currGraph, (fmt::format("Q: {} M: {}",iCharge->first,MASS)+string(" GeV/c^{2}")).c_str(), "P");
    /*--------------------------------------------------------------------*/

    /*--------------------------------------------------------------------*/
    //Fill the reco graph, unknown Q
    currGraph = new TGraph(distList.find(uncorr_pKey)->second.data.size());
    currGraph->SetMarkerStyle(20);
    currGraph->SetMarkerColorAlpha(*iColor,0.8);
    currGraph->SetMarkerSize(0.3);
    
    //Now we want to loop through all momenta and Ih for this charge and mass
    //Have the data stored in distList
    iP = distList.find(uncorr_pKey)->second.data.begin();
    iIh = distList.find(reco_ihKey)->second.data.begin();
    i = 0;
    //Correcting momentum here. Assuming correct charge isknown
    //Therefore need to multiply reco by Q in order to get corrected
    for( ; iP != distList.find(uncorr_pKey)->second.data.end(); ++iP, ++iIh ){
      currGraph->SetPoint(i,(*iP),*iIh);
      ++i;
    }
    recoOutGraph_UnknownQ->Add(currGraph);
    recoDistLegend_UnknownQ->AddEntry(currGraph, (fmt::format("Q: {} M: {}",iCharge->first,MASS)+string(" GeV/c^{2}")).c_str(), "P");
    /*--------------------------------------------------------------------*/
    
    /*--------------------------------------------------------------------*/
    //Now we can compare the above two to the version below, where we use gen (absolutely correct) P, and the matched reco Ih.
    //Will have as many points as we do reco ihs.
    currGraph = new TGraph(distList.find(reco_ihKey)->second.data.size());
    currGraph->SetMarkerStyle(20);
    currGraph->SetMarkerColorAlpha(*iColor,0.8);
    currGraph->SetMarkerSize(0.3);
    Key recoKey (string("particle"), "Reco", (int)(iCharge->first), MASS);
    Key genKey (string("particle"), "Gen", (int)(iCharge->first), MASS);
    i=0;
    //Need to loop the gen and reco particles, looking for matches
    for( const auto &iGenEvents : events[genKey] ){
      //Get the vector of all gen particles from this event
      auto genParticles = iGenEvents.second;
      //Loop over the gen particles from this event
      for( const auto &iGenParticle : genParticles ){
        //Skip all non gen HSCP from event
        if (TMath::Abs(iGenParticle.PDG)!=17) continue;
        for( const auto &iRecoParticle : events[recoKey][iGenEvents.first] ){
          //Only fill the dist if we have a match between the particles
          if (!isMatch(iGenParticle,iRecoParticle)) continue;
          currGraph->SetPoint(i,iGenParticle.P,iRecoParticle.Ih);
          i++;
        }//reco
      }//gen
    }//events
    genOutGraph->Add(currGraph);
    genDistLegend->AddEntry(currGraph, (fmt::format("Q: {} M: {}",iCharge->first,MASS)+string(" GeV/c^{2}")).c_str(), "P");
    /*--------------------------------------------------------------------*/
    
    ++iColor;    
  }
  
  //Now all distributions have been added to their own TH2F's
  //They have also been added to the THStack.
  TCanvas *genOutCanvas = new TCanvas("gen_IhVsP","gen_IhVsP",500,500);
  TCanvas *recoOutCanvas_KnownQ = new TCanvas("reco_IhVsP_KnownQ","reco_IhVsP_KnownQ",500,500);
  TCanvas *recoOutCanvas_UnknownQ = new TCanvas("reco_IhVsP_UnknownQ","reco_IhVsP_UnknownQ",500,500);
  
  //Write the reco level, unknown Q plot
  recoOutCanvas_UnknownQ->cd();
  recoOutGraph_UnknownQ->Draw("ap");
  recoOutGraph_UnknownQ->SetTitle("reco;P_{r} [GeV/c];I_{h} [MeV/cm]");
  gPad->Update();
  recoDistLegend_UnknownQ->Draw();
  gPad->Update();
  outFile->cd();
  recoOutCanvas_UnknownQ->Write();
  string IhVsPName = fmt::format("IhVsP_CC_{}_CT_{}_RecoUnknownQ.pdf",CC,CT);
  recoOutCanvas_UnknownQ->Print(IhVsPName.c_str(),"pdf");

  //Write the reco level, known Q plot
  recoOutCanvas_KnownQ->cd();
  recoOutGraph_KnownQ->Draw("ap");
  recoOutGraph_KnownQ->SetTitle("reco;P_{r} [GeV/c];I_{h} [MeV/cm]");
  gPad->Update();
  recoDistLegend_KnownQ->Draw();
  gPad->Update();
  outFile->cd();
  recoOutCanvas_KnownQ->Write();
  IhVsPName = fmt::format("IhVsP_CC_{}_CT_{}_RecoKnownQ.pdf",CC,CT);
  recoOutCanvas_KnownQ->Print(IhVsPName.c_str(),"pdf");

  //Write the matched gen P plot
  genOutCanvas->cd();
  genOutGraph->Draw("ap");
  genOutGraph->SetTitle("Gen;P_{g} [GeV/c];I_{h} [MeV/cm]");
  gPad->Update();
  genDistLegend->Draw();
  gPad->Update();
  outFile->cd();
  genOutCanvas->Write();
  IhVsPName = fmt::format("IhVsP_CC_{}_CT_{}_GenP.pdf",CC,CT);
  genOutCanvas->Print(IhVsPName.c_str(),"pdf");
}


void WritePrecoVsPgen(TFile *&outFile, const distMap &distList, const map<double,int> &chargeCounts, const Events &events, const char &CC, const char &CT){
  const array<int,12> colorList = {1,2,3,4,6,41,34,46,12,8,14,5};
  const array<int,12> markerList = {2,3,4,5,25,26,27,20,21,22,33,34};
  const int MASS = 900; //Desired mass to vary charges over
  THStack *outDist = new THStack("PrVsPg","");
  THStack *outDistScaled = new THStack("PrVsPg_Scaled","");
  
  TF1 *currLine = NULL;
  TH2F *currDist, *currDistScaled= NULL;
  TLegend *distLegend = new TLegend(0.65,0.65,0.95,0.95);
  distLegend->SetTextSize(0.025);
  distLegend->SetFillColor(0);
  
  vector<string> chargeNames;
  vector<TF1*> lines;
  
  //Check to see if this mass is available with the desired charge
  auto iCharge = chargeCounts.begin();
  auto iColor = colorList.begin();
  for( ; iCharge!=chargeCounts.end(); iCharge++) {
    Key recoKey (string("particle"), "Reco", (int)(iCharge->first), MASS);
    Key genKey (string("particle"), "Gen", (int)(iCharge->first), MASS);
    //iCharge->fisrt is in units of e/3
    //map.count returns zero if the item is not found. want to skip masses that are not with desired charge. Also, only want files with > 200 events  
    if (events.count(genKey) == 0 || events[recoKey].size() < 200){
      continue;
    }    
    
    chargeNames.push_back(fmt::format("Q = {}",iCharge->first));
    
    string name = fmt::format("P{}",iCharge->first);
    currDist = new TH2F(name.c_str(),name.c_str(),100,0,2000,100,0,2000);
    currDist->SetFillColorAlpha(*iColor,0.75);
    currDistScaled = currDist->Clone();
    string function = fmt::format("(1/{})*x",(iCharge->first/3.0));
    currLine = new TF1(name.c_str(),function.c_str(),0,2000);
    currLine->SetLineColor(*iColor);
    currLine->SetLineStyle(2);
    currLine->SetLineWidth(2);
    lines.push_back( currLine );
    //We want to make sure that we are only taking the ratio of reco to
    //gen momentum for HSCP's only, and on top of that, only those HSCP's
    //that are in the same event. So, for each Gen HSCP, we need to
    //loop over ALL reco HSCP's that are in the same event.
    //Get the Gen key particles
    //typedef map< Key, map< int, vector< Particle > > > Events;
    
    for( const auto &iGenEvents : events[genKey] ){
      //Get the vector of all gen particles from this event
      auto genParticles = iGenEvents.second;
      //Loop over the gen particles from this event
      for( const auto &iGenParticle : genParticles ){
        //Skip all non gen HSCP from event
        if (TMath::Abs(iGenParticle.PDG)!=17) continue;
        for( const auto &iRecoParticle : events[recoKey][iGenEvents.first] ){
          if( !isMatch(iGenParticle,iRecoParticle) ) continue;
          currDist->Fill(iGenParticle.Pt,iRecoParticle.Pt);
          currDistScaled->Fill(iGenParticle.Pt,(iCharge->first/3.0)*iRecoParticle.Pt);
        }//reco
      }//gen
    }//events
   
    outDist->Add(currDist);
    outDistScaled->Add(currDistScaled);
    distLegend->AddEntry(currDist, (fmt::format("Q: {} M: {} GeV/",iCharge->first,MASS)+string("c^{2}")).c_str());
    ++iColor;
  }//file loop

  TCanvas *outCanvas = new TCanvas("Poutcanv","Poutcanv",2000,2000);
  TCanvas *outCanvasScaled = new TCanvas("PoutcanvS","PoutcanvS",2000,2000);
  outCanvas->cd();
  outDist->Draw("BOX1");
  outDist->GetXaxis()->SetTitle("P_{t}^{g} [GeV/c]");
  outDist->GetYaxis()->SetTitle("P_{t}^{r} [GeV/c]");
  outDist->GetYaxis()->SetTitleOffset(1.4);
  outDist->GetXaxis()->SetRangeUser(0,2000);
  outDist->GetYaxis()->SetRangeUser(0,2000);
  outDist->SetTitle("HSCP Reco Momentum Vs Gen Momentum");
  gPad->Update();
  //gPad->Update();
  //Draw all of the lines
  for( const auto &iLine : lines ){
    iLine->Draw("SAME");
  }
  gPad->Update();
  distLegend->Draw("F");
  //distLegend->SetTextAlign(22);
  gPad->Update();
  outFile->cd();
  outCanvas->Write();
  string PrVsPgName = fmt::format("PrVsPg_CC_{}_CT_{}.pdf",CC,CT);
  outCanvas->Print(PrVsPgName.c_str(),"pdf");


  outCanvasScaled->cd();
  outDistScaled->Draw("BOX1");
  outDistScaled->GetXaxis()->SetTitle("P_{t}^{g} [GeV/c]");
  outDistScaled->GetYaxis()->SetTitle("Q #times P_{t}^{r} [GeV/c]");
  outDistScaled->GetYaxis()->SetTitleOffset(1.4);
  outDistScaled->GetXaxis()->SetRangeUser(0,2000);
  outDistScaled->GetYaxis()->SetRangeUser(0,2000);
  outDistScaled->SetTitle("HSCP Scaled Reco Momentum Vs Gen Momentum");
  gPad->Update();
  currLine = new TF1("45","x",0,2000);
  currLine->SetLineColor(kYellow);
  currLine->SetLineStyle(2);
  currLine->SetLineWidth(2);
  currLine->Draw("SAME");
  gPad->Update();
  distLegend->Draw("F");
  distLegend->SetTextAlign(22);
  gPad->Update();
  outFile->cd();
  outCanvasScaled->Write();
  string PrVsPgScaledName = fmt::format("PrVsPgScaled_CC_{}_CT_{}.pdf",CC,CT);
  outCanvasScaled->Print(PrVsPgScaledName.c_str(),"pdf");
}

void WriteBrecoVsBgen(TFile *&outFile, const distMap &distList, const map<double,int> &chargeCounts, const map<double,int> &massCounts, const Events &events, const char &CC, const char &CT){
  //Want a mass scan and a charge scan
  const Int_t NRGBs = 5;
  const Int_t NCont = 255;
 
  Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
  Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
  Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
  Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
  TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
  gStyle->SetNumberContours(NCont);
  
  const array<int,12> colorList = {1,2,3,4,6,41,34,46,12,8,14,5};
  const array<int,12> markerList = {2,3,4,5,25,26,27,20,21,22,33,34};
  
  gStyle->SetPalette(54);
  const int MASS = 300; //Desired mass to vary charges over
  const int CHARGE = 3; //Desired charge (in e/3) to vary masses over

  vector<TCanvas*> constMassOutCanvases;
  vector<TCanvas*> constChargeOutCanvases;
  vector<TH2F*> constMassOutDists;
  vector<TH2F*> constChargeOutDists;
  vector<string> chargeNames;
  vector<string> massNames;

  TH2F *currBetaDist = NULL;
  
  auto iCharge = chargeCounts.begin();
  auto iColor = colorList.begin();

  /*------------------------------------------------------------------------*/
  //Conduct the charge scan
  for( ; iCharge!=chargeCounts.end(); iCharge++) {
    Key recoKey_Q (string("particle"), "Reco", (int)(iCharge->first), MASS);
    Key genKey_Q (string("particle"), "Gen", (int)(iCharge->first), MASS);

    if (events.count(genKey_Q) == 0 || events[recoKey_Q].size() < 200){
      continue;
    }    
    
    chargeNames.push_back(fmt::format("{}",iCharge->first));
    
    string name = fmt::format("Beta_Q_{}_M_{}_M",iCharge->first,MASS);
    currBetaDist = new TH2F(name.c_str(),name.c_str(),50,0,1.0,50,0,1.5);
    currBetaDist->SetFillColorAlpha(*iColor,0.75);
    
    for( const auto &iGenEvents : events[genKey_Q] ){
      //Get the vector of all gen particles from this event
      auto genParticles = iGenEvents.second;
      //Loop over the gen particles from this event
      for( const auto &iGenParticle : genParticles ){
        //Skip all non gen HSCP from event
        if (TMath::Abs(iGenParticle.PDG)!=17) continue;
        for( const auto &iRecoParticle : events[recoKey_Q][iGenEvents.first] ){
          if( !isMatch(iGenParticle,iRecoParticle) ) continue;
          currBetaDist->Fill(iGenParticle.Beta,iRecoParticle.Beta);
        }//reco
      }//gen
    }//events

    constMassOutDists.push_back(currBetaDist);
    constMassOutCanvases.push_back(new TCanvas(fmt::format("BrVsBg_Q_{}_MCONST",iCharge->first).c_str(),"",2000,2000)); 
    ++iColor;
  }//charge loop
  /*------------------------------------------------------------------------*/
  iColor = colorList.begin();
  auto iMass = massCounts.begin();
  /*------------------------------------------------------------------------*/
  //Conduct the mass scan
  for( ; iMass!=massCounts.end(); iMass++) {
    Key recoKey_M (string("particle"), "Reco", CHARGE, (int)(iMass->first));
    Key genKey_M (string("particle"), "Gen", CHARGE, (int)(iMass->first));

    if (events.count(genKey_M) == 0 || events[recoKey_M].size() < 200){
      continue;
    }    
    
    massNames.push_back(fmt::format("{}",iMass->first));
    
    string name = fmt::format("Beta_Q_{}_M_{}_Q",CHARGE,iMass->first);
    currBetaDist = new TH2F(name.c_str(),name.c_str(),50,0,1.0,50,0,1.5);
    currBetaDist->SetFillColorAlpha(*iColor,0.75);
    
    for( const auto &iGenEvents : events[genKey_M] ){
      //Get the vector of all gen particles from this event
      auto genParticles = iGenEvents.second;
      //Loop over the gen particles from this event
      for( const auto &iGenParticle : genParticles ){
        //Skip all non gen HSCP from event
        if (TMath::Abs(iGenParticle.PDG)!=17) continue;
        for( const auto &iRecoParticle : events[recoKey_M][iGenEvents.first] ){
          currBetaDist->Fill(iGenParticle.Beta,iRecoParticle.Beta);
        }//reco
      }//gen
    }//events

    constChargeOutDists.push_back(currBetaDist);
    constChargeOutCanvases.push_back(new TCanvas(fmt::format("BrVsBg_M_{}_QCONST",iMass->first).c_str(),"",2000,2000)); 
    ++iColor;
  }//mass loop
  /*------------------------------------------------------------------------*/
  
  /*------------------------------------------------------------------------*/
  //Write the constant mass distributions
  auto iDists = constMassOutDists.begin();
  auto iCanvases = constMassOutCanvases.begin();
  auto iChargeName = chargeNames.begin();
  TF1 *currLine;
  currLine = new TF1("45","x",0,1);
  currLine->SetLineColor(kBlack);
  currLine->SetLineStyle(2);
  currLine->SetLineWidth(2);
  for(; iDists != constMassOutDists.end(); ++iDists,++iCanvases,++iChargeName){
    string name = fmt::format("BrVsBg_Q_{}_CC_{}_CT_{}.pdf",(*iChargeName),CC,CT);
    (*iCanvases)->cd();
    (*iDists)->Scale(1.0/(*iDists)->Integral());
    gStyle->SetOptTitle(true);
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.15);
    gPad->SetTopMargin(0.10);
    gPad->SetBottomMargin(0.10);
    (*iDists)->Draw("colz");
    (*iDists)->GetYaxis()->SetTitle("#beta_{r}");
    (*iDists)->GetXaxis()->SetTitle("#beta_{g}");
    (*iDists)->SetTitle((fmt::format("Q: {} M: {} GeV/",(*iChargeName),MASS)+string("c^{2}")).c_str());
    (*iDists)->GetYaxis()->SetTitleOffset(1.4);
    (*iDists)->GetXaxis()->SetRangeUser(0,1.0);
    (*iDists)->GetYaxis()->SetRangeUser(0,1.5);
    currLine->Draw("SAME");
    outFile->cd();
    (*iCanvases)->Write(); 
    (*iCanvases)->Print(name.c_str(),"pdf");
  }
  /*------------------------------------------------------------------------*/

  
  /*------------------------------------------------------------------------*/
  //Write the constant charge distributions
  iDists = constChargeOutDists.begin();
  iCanvases = constChargeOutCanvases.begin();
  auto iMassName = massNames.begin();
  currLine = new TF1("45","x",0,1);
  currLine->SetLineColor(kBlack);
  currLine->SetLineStyle(2);
  currLine->SetLineWidth(2);
  for(; iDists != constChargeOutDists.end(); ++iDists,++iCanvases,++iMassName){
    string name = fmt::format("BrVsBg_M_{}_CC_{}_CT_{}.pdf",(*iMassName),CC,CT);
    (*iCanvases)->cd();
    (*iDists)->Scale(1.0/(*iDists)->Integral());
    gStyle->SetOptTitle(true);
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.15);
    gPad->SetTopMargin(0.10);
    gPad->SetBottomMargin(0.10);
    (*iDists)->Draw("colz");
    (*iDists)->GetYaxis()->SetTitle("#beta_{r}");
    (*iDists)->GetXaxis()->SetTitle("#beta_{g}");
    (*iDists)->SetTitle((fmt::format("Q: {} M: {} GeV/",CHARGE,(*iMassName))+string("c^{2}")).c_str());
    (*iDists)->GetYaxis()->SetTitleOffset(1.4);
    (*iDists)->GetXaxis()->SetRangeUser(0,1.0);
    (*iDists)->GetYaxis()->SetRangeUser(0,1.5);
    currLine->Draw("SAME");
    outFile->cd();
    (*iCanvases)->Write(); 
    (*iCanvases)->Print(name.c_str(),"pdf");
  }
  /*------------------------------------------------------------------------*/
}


void WriteMrecoVsMgen(TFile *&outFile, const distMap &distList, const map<double,int> &chargeCounts, const map<double,int> &massCounts, const Events &events, const char &CC, const char &CT){
  const array<int,12> qcolorList = {1,2,3,4,6,41,34,46,12,8,14,5};
  const array<int,12> mcolorList = {1,2,3,4,6,41,34,46,12,8,14,5};
  const array<int,12> markerList = {2,3,4,5,25,26,27,20,21,22,33,34};
  
  vector<string> Names;
  vector<TCanvas*> outCorrCanvs, outUncorrCanvs;
  vector<TCanvas*> outCorrCanvs2, outUncorrCanvs2;
  vector<TH1F*> outCorrDists, outUncorrDists;
  vector<TH2F*> outCorrHistos, outUncorrHistos;
  TH1F *currCorrDist, *currUncorrDist;
  TH2F *currCorrHisto, *currUncorrHisto;

  const int MASS = 300;
  const int CHARGE = 3;

  auto iMColor = mcolorList.begin();
  auto iQColor = qcolorList.begin();
  
  //We now want to make stacked plots for constant charge with varying mass & constant mass with varying charge for both corrected and not
  THStack *corrMassScanStack = new THStack("corrMassScanS","");
  THStack *uncorrMassScanStack = new THStack("uncorrMassScanS","");
  THStack *corrChargeScanStack = new THStack("corrChargeScanS","");
  THStack *uncorrChargeScanStack = new THStack("uncorrChargeScanS","");

  TCanvas *corrMassScanCanv = new TCanvas("corrMassScanC","corrMassScanC",2000,2000);
  TCanvas *uncorrMassScanCanv = new TCanvas("uncorrMassScanC","uncorrMassScanC",2000,2000);
  TCanvas *corrChargeScanCanv = new TCanvas("corrChargeScanC","corrChargeScanC",2000,2000);
  TCanvas *uncorrChargeScanCanv = new TCanvas("uncorrChargeScanC","uncorrChargeScanC",2000,2000);

  TLegend *corrMassScanLegend = new TLegend(0.6,0.7,0.85,0.9);
  TLegend *uncorrMassScanLegend = new TLegend(0.6,0.7,0.85,0.9);
  TLegend *corrChargeScanLegend = new TLegend(0.6,0.7,0.85,0.9);
  TLegend *uncorrChargeScanLegend = new TLegend(0.6,0.7,0.85,0.9);

  //Use the particles, divide by the mass in the loop, make new distributions.
  //Take the means of these TH1Fs, plot on a multigraph, but keep these as well.
  for(auto iMass = massCounts.begin() ; iMass!=massCounts.end(); iMass++) {
    auto iChargeColor = qcolorList.begin();
    for(auto iCharge = chargeCounts.begin(); iCharge!=chargeCounts.end(); iCharge++){
      Key recoKey (string("particle"), "Reco", iCharge->first, iMass->first);
      Key genKey (string("particle"), "Gen", iCharge->first, iMass->first);

      if (events.count(genKey) == 0 || events[recoKey].size() < 200){
        continue;
      }    
      string corrName = fmt::format("MrDivMg_Q_{}_M_{}_Corr",iMass->first,iCharge->first);
      string corrName2 = fmt::format("MrVsMg_Q_{}_M_{}_Corr",iMass->first,iCharge->first);
      string uncorrName = fmt::format("MrDivMg_Q_{}_M_{}_Uncorr",iMass->first,iCharge->first);
      string uncorrName2 = fmt::format("MrVsMg_Q_{}_M_{}_Uncorr",iMass->first,iCharge->first);
      
      currCorrDist = new TH1F(corrName.c_str(),corrName.c_str(),20,0,6);
      currCorrHisto = new TH2F(corrName2.c_str(),corrName2.c_str(),100,0,2500,100,0,2500);
      currUncorrDist = new TH1F(uncorrName.c_str(),uncorrName.c_str(),20,0,6);
      currUncorrHisto = new TH2F(uncorrName2.c_str(),uncorrName2.c_str(),100,0,2500,100,0,2500);
      
      for( const auto &iGenEvents : events[genKey] ){
        //Get the vector of all gen particles from this event
        auto genParticles = iGenEvents.second;
        //Loop over the gen particles from this event
        for( const auto &iGenParticle : genParticles ){
          //Skip all non gen HSCP from event
          if (TMath::Abs(iGenParticle.PDG)!=17) continue;
          //Now for each of these reco particles we want to compare the gen particles
          for( const auto &iRecoParticle : events[recoKey][iGenEvents.first] ){
            if ( !isMatch(iGenParticle,iRecoParticle) ) continue;
            currCorrDist->Fill( iRecoParticle.corrMass / iGenParticle.corrMass );
            currCorrHisto->Fill( iGenParticle.corrMass, iRecoParticle.corrMass );
            currUncorrDist->Fill( iRecoParticle.uncorrMass / iMass->first );
            currUncorrHisto->Fill( iGenParticle.corrMass, iRecoParticle.uncorrMass );
          }//reco
        }//gen
      }//events

      outCorrDists.push_back(currCorrDist);
      outCorrHistos.push_back(currCorrHisto);
      outUncorrDists.push_back(currUncorrDist);
      outUncorrHistos.push_back(currUncorrHisto);
      
      Names.push_back( fmt::format("M_{}_Q_{}",iMass->first,iCharge->first).c_str() );
      
      outCorrCanvs.push_back(new TCanvas(fmt::format("M_{}_Q_{}_MrDivMgCorr",iMass->first,iCharge->first).c_str(),"",2000,2000));
      outCorrCanvs2.push_back(new TCanvas(fmt::format("M_{}_Q_{}_MrVsMgCorr",iMass->first,iCharge->first).c_str(),"",2000,2000));
      outUncorrCanvs.push_back(new TCanvas(fmt::format("M_{}_Q_{}_MrDivMgUncorr",iMass->first,iCharge->first).c_str(),"",2000,2000));
      outUncorrCanvs2.push_back(new TCanvas(fmt::format("M_{}_Q_{}_MrVsMgUncorr",iMass->first,iCharge->first).c_str(),"",2000,2000));

      if( iCharge->first == CHARGE ){
        currCorrDist->SetLineColor(*iMColor);
        corrMassScanStack->Add(currCorrDist);
        currCorrDist->Scale(1.0/currCorrDist->Integral());
        corrMassScanLegend->AddEntry(currCorrDist,fmt::format("M_{}_Q_{}",iMass->first,CHARGE).c_str());

        currUncorrDist->SetLineColor(*iMColor);
        uncorrMassScanStack->Add(currUncorrDist);
        currUncorrDist->Scale(1.0/currUncorrDist->Integral());
        uncorrMassScanLegend->AddEntry(currUncorrDist,fmt::format("M_{}_Q_{}",iMass->first,CHARGE).c_str());
                
        iMColor++;
      }
      if( iMass->first == MASS ){
        currCorrDist->SetLineColor(*iQColor);
        //Normalize before addin
        corrChargeScanStack->Add(currCorrDist);
        currCorrDist->Scale(1.0/currCorrDist->Integral());
        corrChargeScanLegend->AddEntry(currCorrDist,fmt::format("M_{}_Q_{}",MASS,iCharge->first).c_str());

        currUncorrDist->SetLineColor(*iQColor);
        //Normalize before adding
        uncorrChargeScanStack->Add(currUncorrDist);
        currUncorrDist->Scale(1.0/currUncorrDist->Integral());
        uncorrChargeScanLegend->AddEntry(currUncorrDist,fmt::format("M_{}_Q_{}",MASS,iCharge->first).c_str());

        iQColor++;
      }
    }//charge
  }//mass

  //Write the distributions to file. Should have as many corr dists as uncorr dists
  auto iCorrDists = outCorrDists.begin();
  auto iCorrHistos = outCorrHistos.begin();
  auto iCorrCanvs = outCorrCanvs.begin();
  auto iCorrCanvs2 = outCorrCanvs2.begin();
  
  auto iUncorrDists = outUncorrDists.begin();
  auto iUncorrHistos = outUncorrHistos.begin();
  auto iUncorrCanvs = outUncorrCanvs.begin();
  auto iUncorrCanvs2 = outUncorrCanvs2.begin();
  
  auto iNames = Names.begin();

  //Write the individual distributions to the TFile, as well as the TH2s
  string fName;
  for( ; iCorrDists != outCorrDists.end(); iCorrDists++, iCorrHistos++, iUncorrDists++, iUncorrHistos++, iCorrCanvs++, iCorrCanvs2++, iUncorrCanvs++, iUncorrCanvs2++, iNames++ ){
    (*iCorrCanvs)->cd();
    (*iCorrDists)->Draw();
    (*iCorrDists)->GetXaxis()->SetTitle("m_{r}^{c}/m_{g}");
    (*iCorrDists)->GetYaxis()->SetTitle("Events");
    (*iCorrDists)->SetTitle(iNames->c_str());
    outFile->cd();
    fName = string("MrDivMg_") + *iNames + fmt::format("_CC_{}_CT_{}_",CC,CT) +string("Corr.pdf");
    (*iCorrCanvs)->Write();

    (*iUncorrCanvs)->cd();
    (*iUncorrDists)->Draw();
    (*iUncorrDists)->GetXaxis()->SetTitle("m_{r}^{u}/m_{g}");
    (*iUncorrDists)->GetYaxis()->SetTitle("Events");
    (*iUncorrDists)->SetTitle(iNames->c_str());
    outFile->cd();
    fName = string("MrDivMg_") + *iNames + fmt::format("_CC_{}_CT_{}_",CC,CT) + string("Uncorr.pdf");
    (*iUncorrCanvs)->Write();

    (*iCorrCanvs2)->cd();
    (*iCorrHistos)->Draw("colz");
    (*iCorrHistos)->GetXaxis()->SetTitle("m_{g}");
    (*iCorrHistos)->GetXaxis()->SetTitle("m_{r}^{c}");
    (*iCorrHistos)->SetTitle(iNames->c_str());
    outFile->cd();
    (*iCorrHistos)->Write();
    fName = string("MrVsMg_") + *iNames + fmt::format("_CC_{}_CT_{}_",CC,CT) +string("Corr.pdf");
    (*iCorrCanvs2)->Print(fName.c_str(),"pdf");
    
    (*iUncorrCanvs2)->cd();
    (*iUncorrHistos)->Draw("colz");
    (*iUncorrHistos)->GetXaxis()->SetTitle("m_{g}");
    (*iUncorrHistos)->GetXaxis()->SetTitle("m_{r}^{c}");
    (*iUncorrHistos)->SetTitle(iNames->c_str());
    outFile->cd();
    (*iUncorrHistos)->Write();
    fName = string("MrVsMg_") + *iNames + fmt::format("_CC_{}_CT_{}_",CC,CT) +string("Uncorr.pdf");
    (*iUncorrCanvs2)->Print(fName.c_str(),"pdf");
    
  }

  //Now write the stacks
  corrMassScanCanv->cd();
  corrMassScanStack->Draw();
  gPad->Update();
  corrMassScanStack->GetXaxis()->SetTitle("m_{r}^{c}/m_{g}");
  corrMassScanStack->GetYaxis()->SetTitle("Events");
  corrMassScanStack->SetTitle("Corrected Mass");
  gPad->Update();
  corrMassScanLegend->Draw("SAME");
  gPad->Update();
  fName = fmt::format("MrDivMg_MScan_Q_{}_CC_{}_CT_{}_Corr.pdf",CHARGE,CC,CT);
  outFile->cd();
  corrMassScanCanv->Print(fName.c_str(),"pdf");

  uncorrMassScanCanv->cd();
  uncorrMassScanStack->Draw();
  gPad->Update();
  uncorrMassScanStack->GetXaxis()->SetTitle("m_{r}^{u}/m_{g}");
  uncorrMassScanStack->GetYaxis()->SetTitle("Events");
  uncorrMassScanStack->SetTitle("Uncorrected Mass");
  gPad->Update();
  uncorrMassScanLegend->Draw("SAME");
  gPad->Update();
  fName = fmt::format("MrDivMg_MScan_Q_{}_CC_{}_CT_{}_Uncorr.pdf",CHARGE,CC,CT);
  outFile->cd();
  uncorrMassScanCanv->Print(fName.c_str(),"pdf");

  corrChargeScanCanv->cd();
  corrChargeScanStack->Draw();
  gPad->Update();
  corrChargeScanStack->GetXaxis()->SetTitle("m_{r}^{c}/m_{g}");
  corrChargeScanStack->GetYaxis()->SetTitle("Events/0.2");
  corrChargeScanStack->SetTitle("Corrected Mass");
  gPad->Update();
  corrChargeScanLegend->Draw("SAME");
  gPad->Update();
  fName = fmt::format("MrDivMg_QScan_M_{}_CC_{}_CT_{}_Corr.pdf",MASS,CC,CT);
  outFile->cd();
  corrChargeScanCanv->Print(fName.c_str(),"pdf");

  uncorrChargeScanCanv->cd();
  uncorrChargeScanStack->Draw();
  gPad->Update();
  uncorrChargeScanStack->GetXaxis()->SetTitle("m_{r}^{u}/m_{g}");
  uncorrChargeScanStack->GetYaxis()->SetTitle("Events/0.2");
  uncorrChargeScanStack->SetTitle("Uncorrected Mass");
  gPad->Update();
  uncorrChargeScanLegend->Draw("SAME");
  gPad->Update();
  fName = fmt::format("MrDivMg_MScan_M_{}_CC_{}_CT_{}_Uncorr.pdf",MASS,CC,CT);
  outFile->cd();
  uncorrChargeScanCanv->Print(fName.c_str(),"pdf");
  
  
}

void WriteBgenVsEtagen(TFile *&outFile, const distMap &distList, const map<double,int> &massCounts, const Events &events, const char &CC, const char &CT){
  const Int_t NRGBs = 5;
  const Int_t NCont = 255;
 
  Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
  Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
  Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
  Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
  TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
  gStyle->SetNumberContours(NCont);
  
  const array<int,12> colorList = {1,2,3,4,6,41,34,46,12,8,14,5};
  const array<int,12> markerList = {2,3,4,5,25,26,27,20,21,22,33,34};
  
  gStyle->SetPalette(54);
  const int CHARGE = 3; //Desired charge to vary masses over (in units of e/3)
  
  vector<TH2F*> outDists;
  vector<TCanvas*> outCanvases;
  vector<string> massNames;

  TH2F *currDist = NULL;
  
  //Check to see if this mass is available with the desired charge
  auto iMass = massCounts.begin();
  auto iColor = colorList.begin();
  for( ; iMass!=massCounts.end(); iMass++) {
    Key recoKey (string("particle"), "Reco", CHARGE, iMass->first);
    Key genKey (string("particle"), "Gen", CHARGE, iMass->first);

    if (events.count(genKey) == 0 || events[recoKey].size() < 200){
      continue;
    }    
    
    massNames.push_back(fmt::format("{}",iMass->first));
    string name = fmt::format("BetaVsEta_Q_{}_M_{}",iMass->first,CHARGE);
    currDist = new TH2F(name.c_str(),name.c_str(),20,-4,4,20,10,1);
    currDist->SetFillColorAlpha(*iColor,0.75);
    
    for( const auto &iGenEvents : events[genKey] ){
      //Get the vector of all gen particles from this event
      auto genParticles = iGenEvents.second;
      //Loop over the gen particles from this event
      for( const auto &iGenParticle : genParticles ){
        //Skip all non gen HSCP from event
        if (TMath::Abs(iGenParticle.PDG)!=17) continue;
        currDist->Fill(iGenParticle.Eta,iGenParticle.P/iGenParticle.E);
      }//gen
    }//events

    outDists.push_back(currDist);
    outCanvases.push_back(new TCanvas(fmt::format("{}",iMass->first).c_str(),"",2000,2000));
    ++iColor;
  }//file loop

 
  auto iDists = outDists.begin();
  auto iCanvases = outCanvases.begin();
  auto iMassName = massNames.begin();
  
  for(; iDists != outDists.end(); ++iDists,++iCanvases,++iMassName){
    string name = fmt::format("BgVsEtag_M_{}_CC_{}_CT_{}.pdf",(*iMassName),CC,CT);
    (*iCanvases)->cd();
    (*iDists)->Scale(1.0/(*iDists)->Integral());
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.15);
    gPad->SetTopMargin(0.10);
    gPad->SetBottomMargin(0.10);
    (*iDists)->Draw("colz");
    (*iDists)->GetYaxis()->SetTitle("#beta_{g}");
    (*iDists)->GetXaxis()->SetTitle("#eta_{g}");
    (*iDists)->SetTitle((fmt::format("Q: {} M: {} GeV",CHARGE,(*iMassName))+string("/c^{2}")).c_str());
    (*iDists)->GetYaxis()->SetTitleOffset(1.4);
    (*iDists)->GetXaxis()->SetRangeUser(-4,4);
    (*iDists)->GetYaxis()->SetRangeUser(0,1);
  
    outFile->cd();
    (*iCanvases)->Write(); 
    (*iCanvases)->Print(name.c_str(),"pdf");
  }
}



/*Main*/
int main(int argc, char **argv){
  TApplication theApp("App",0,0);

  //Parse the command line arguments.
  //argv[1] is T or F for Cluster Cleaning
  //argv[2] it T or F for Cross Talk Inversion
  char CC = toupper(argv[1][0]);
  char CT = toupper(argv[1][2]);

  SetRootStyle();
  gROOT->SetBatch(kTRUE); //Don't draw things when created.
  
  string *fileList = ListOfFiles(CC,CT);
  
  //Total number of files to be processed
  int numFiles = atoi(fileList[0].c_str());

  NameDat *fileNameData = FileNameParser( fileList, numFiles );

  //Arrays of the charges and masses associated with each file
  double *masses = fileNameData->mass;
  double *charges = fileNameData->charge; //These are the ACTUAL charges

  //List of the names of the files to be opened
  fileList = fileNameData->fileNames; //Can index from 0 now

  // //Number of each charge and each mass
  map<double,int> chargeCounts = fileNameData->chargeCounts;
  map<double,int> massCounts = fileNameData->massCounts;

  distMap distList;

  Events events;
  
  vector<string> distNames = {
    "beta",
    "energy",
    "eta",
    "gamma",
    "Ih",
    "genP",
    "genPt",
    "recoPCorr",
    "recoPUncorr",
    "recoPtCorr",
    "recoPtUncorr",
    "phi",
    "theta",
    "TOF",
    "genMass",
    "recoMassCorr",
    "recoMassUncorr",
    "charge"
  };
  
  //Right now, sets the same limits and bins for both reco and gen particles.
  map<string,distProp> distProps = {
    {"beta",{limits(0,1.5),50,axes("#beta","events/{}")}},
    {"energy",{limits(0,2000),200,axes("E [GeV]","events/{} GeV")}},
    {"eta",{limits(-5,5),50,axes("#eta", "events/{}")}},
    {"gamma",{limits(0,100.0),200,axes("#gamma", "events/{}")}},
    {"Ih",{limits(0,50),25,axes("I_{h} [MeV/cm]", "events/{} MeV/cm")}},
    {"genP",{limits(0,2600),200,axes("P [GeV/c]", "events/{} GeV/c")}},
    {"genPt",{limits(0,2600),200,axes("P#_t [GeV/c]", "events/{} GeV/c")}},
    {"recoPCorr",{limits(0,2600),200,axes("P [GeV/c]", "events/{} GeV/c")}},
    {"recoPUncorr",{limits(0,2600),200,axes("P [GeV/c]", "events/{} GeV/c")}},
    {"recoPtCorr",{limits(0,2600),200,axes("P [GeV/c]", "events/{} GeV/c")}},
    {"recoPtUncorr",{limits(0,2600),200,axes("P [GeV/c]", "events/{} GeV/c")}},
    {"phi",{limits(-200,200),100,axes("#phi [deg]", "events/{} deg")}},
    {"theta",{limits(0,180),180,axes("#theta [deg]", "events/{} deg")}},
    {"TOF",{limits(0.0,4.0),20,axes("t [arb. units]", "events/{}")}},
    {"genMass",{limits(0,2600),200,axes("mass [GeV/c#^2]", "events/{} GeV/c#^2")}},
    {"recoMassCorr",{limits(0,2600),200,axes("mass [GeV/c#^2]", "events/{} GeV/c#^2")}},
    {"recoMassUncorr",{limits(0,2600),200,axes("mass [GeV/c#^2]", "events/{} GeV/c#^2")}},
    {"charge",{limits(0,40),40,axes("charge [e/3]", "events/{} e/3")}}
  };

  vector<string> types = {"Gen","Reco"};
  
  AllocateDistributions(distList,charges,masses,numFiles,distNames,distProps,types);

  //Pointer for the tree being read
  TTree *tree;
 
  unsigned int reco_event;
  int gen_event;
  unsigned int Hscp;
  float E;
  float reco_Pt;
  double gen_Pt;
  float P;
  float I;
  float Ih; //This is the energy deposition that we are interested in 
  float reco_Mass;
  double gen_Mass;
  float dZ, dXY, dR;
  float reco_eta, reco_phi;
  float recoMassCorr, recoMassUncorr;
  double gen_eta, gen_phi;
  float TOF;
  float theta, thetaDeg;
  float beta, gamma;
  bool hasMuon;
  int genCharge; //Charge given from gen
  int PDG;

  float eta, Pt, Mass, phi;
  unsigned int event;

  bool isGen;
  //Loop over gen/reco
  for (const auto &iTypes : types) {
    //Loop over the files
    if (iTypes == "Gen")
      isGen = true;
    else isGen = false;
    
    for( int iFile = 0; iFile < numFiles; iFile++ ){
      string file = fmt::format("/home/austin/HSCP_Study/HSCP_MC_Files/CC_{}_CT_{}/",CC,CT);
      file += fileList[iFile];

      TFile *datFile = new TFile(file.c_str(), "READ"); //Open the current file

      double *charge = &charges[iFile];
      double *mass = &masses[iFile];
      string type = iTypes;
      
      //Create a string for the appropriate file directory
      ///storage/6/work/askeeters/HSCPStudy/HSCP_Root_Files
      string fileDir;
      if (type == "Reco")
        fileDir = string(fmt::format("mchamp{}_M_{}/HscpCandidates",3*charges[iFile],masses[iFile]));
      else 
        fileDir = string(fmt::format("mchamp{}_M_{}/GENinfo",3*charges[iFile],masses[iFile]));
      
      //Set the appropriate branches for the MC particles
      tree = (TTree*)datFile->Get(fileDir.c_str());
      
      if (!isGen){
        tree->SetBranchAddress("I",&I);
        tree->SetBranchAddress("Ih",&Ih);
        tree->SetBranchAddress("TOF",&TOF);
        tree->SetBranchAddress("dZ",&dZ);
        tree->SetBranchAddress("dXY",&dXY);
        tree->SetBranchAddress("dR",&dR);
        tree->SetBranchAddress("hasMuon",&hasMuon);
        tree->SetBranchAddress("Hscp", &Hscp);
        tree->SetBranchAddress("Event", (&reco_event));
        tree->SetBranchAddress("Pt",(&reco_Pt));
        tree->SetBranchAddress("Mass",(&reco_Mass));
        tree->SetBranchAddress("eta",(&reco_eta));
        tree->SetBranchAddress("phi",(&reco_phi));
      }
      else{
        tree->SetBranchAddress("charge", &genCharge);
        tree->SetBranchAddress("Event", (&gen_event));
        tree->SetBranchAddress("Pt",(&gen_Pt));
        tree->SetBranchAddress("Mass",(&gen_Mass));
        tree->SetBranchAddress("eta",(&gen_eta));
        tree->SetBranchAddress("phi",(&gen_phi));
        tree->SetBranchAddress("PDG",&PDG);
      }

        
      for(int iEvt = 0; iEvt < tree->GetEntries(); iEvt++){
        tree->GetEntry(iEvt);
        if(isGen){
          eta = gen_eta;
          event = gen_event;
          Pt = gen_Pt;
          Mass = gen_Mass;
          eta = gen_eta;
          phi = gen_phi;
        }
        else{
          eta = (float)reco_eta;
          event = (unsigned int)reco_event;
          Pt = (float)reco_Pt;
          Mass = (float)reco_Mass;
          eta = (float)reco_eta;
          phi = (float)reco_phi;
        }

        //We only want the HSCPs for gen. So continue the loop if we do not have one.
        if(isGen)
          if(TMath::Abs(PDG) != 17) continue;
                
        //Calculate theta
        theta = 2*TMath::ATan( TMath::Exp( -1 * (eta) ) );
        thetaDeg = theta * (180.0/TMath::Pi());

        //Fill the theta distribution
        Key thetaKey (string("theta"), type, (int)(3* *charge), (int)*mass);
        distList[thetaKey].distribution->Fill(thetaDeg);
        distList.find(thetaKey)->second.data.push_back(thetaDeg);
                
        //fill the eta distribution
        Key etaKey (string("eta"), type, (int)(3* *charge), (int)*mass);
        distList[etaKey].distribution->Fill((eta));
        distList.find(etaKey)->second.data.push_back(eta);

        //fill the phi distribution
        Key phiKey (string("phi"), type, (int)(3* *charge), (int)*mass);
        phi *= (180.0/TMath::Pi());
        distList[phiKey].distribution->Fill( phi );
        distList.find(phiKey)->second.data.push_back(phi);

        //Calculate Momentum from transverse and theta (rad)
        P = (Pt) / TMath::Sin(theta);
   
        //Fill momentum distributions
        if( !isGen ){
          Key corrPKey (string("recoPCorr"), type, (int)(3* *charge), (int)*mass);
          distList[corrPKey].distribution->Fill( P * (*charge) );
          distList.find(corrPKey)->second.data.push_back( P * (*charge) );

          Key uncorrPKey (string("recoPUncorr"), type, (int)(3* *charge), (int)*mass);
          distList[uncorrPKey].distribution->Fill( P );
          distList.find(uncorrPKey)->second.data.push_back( P );
        }
        else{
          Key PKey (string("genP"), type, (int)(3* *charge), (int)*mass);
          distList[PKey].distribution->Fill( P );
          distList.find(PKey)->second.data.push_back( P );
        }
      
        //Fill the transverse momentum distribution
        if( !isGen ){
          Key corrPtKey (string("recoPtCorr"), type, (int)(3* *charge), (int)*mass);
          distList[corrPtKey].distribution->Fill( Pt * (*charge) );
          distList.find(corrPtKey)->second.data.push_back( Pt * (*charge) );

          Key uncorrPtKey (string("recoPtUncorr"), type, (int)(3* *charge), (int)*mass);
          distList[uncorrPtKey].distribution->Fill( Pt );
          distList.find(uncorrPtKey)->second.data.push_back( Pt );
        }
        else{
          Key PtKey (string("genPt"), type, (int)(3* *charge), (int)*mass);
          distList[PtKey].distribution->Fill( Pt );
          distList.find(PtKey)->second.data.push_back( Pt );
        }
        
        //Calculate the relativistic energy
        //Calculate it correctly (known Q) for HSCP each time.
        //No separate variable for incorrect E needed.
        //Fill the energy distribution
        Key enKey (string("energy"), type, (int)(3* *charge), (int)*mass);
        if( !isGen ){
          E = TMath::Sqrt( (*charge)*(*charge)*P*P + TMath::Power((float)*mass,2) );
        }
        else{
          E = TMath::Sqrt( (P*P) + ((Mass) * (Mass)) );
        }
        distList[enKey].distribution->Fill(E);
        distList.find(enKey)->second.data.push_back(E);

        //Fill the Ih distribution
        Key ihKey (string("Ih"), type, (int)(3* *charge), (int)*mass);
        distList[ihKey].distribution->Fill( Ih );
        distList.find(ihKey)->second.data.push_back(Ih);
       
        //Fill the mass distribution
        //If we have a gen particle, fill with the known mass.
        //If we have a reco particle, use the formula that relies on Ih
        if( !isGen ){
          Key corrMassKey (string("recoMassCorr"), type, (int)(3* *charge), (int)*mass);
          Key uncorrMassKey (string("recoMassUncorr"), type, (int)(3* *charge), (int)*mass);
          recoMassCorr = (*charge) * P * TMath::Sqrt( (Ih - 2.557) / 2.579 );
          recoMassUncorr = P * TMath::Sqrt( (Ih - 2.557) / 2.579 );
          //recoMassCorr = Mass * (*charge);
          //recoMassUncorr = Mass;
          
          distList[corrMassKey].distribution->Fill( recoMassCorr );
          distList.find(corrMassKey)->second.data.push_back( recoMassCorr );

          distList[uncorrMassKey].distribution->Fill( recoMassUncorr );
          distList.find(uncorrMassKey)->second.data.push_back( recoMassUncorr );
        }
        else{
          Key MassKey (string("genMass"), type, (int)(3* *charge), (int)*mass);
          distList[MassKey].distribution->Fill(Mass);
          distList.find(MassKey)->second.data.push_back(Mass);
        }
        
        //Fill TOF distribution
        Key tofKey (string("TOF"), type, (int)(3* *charge), (int)*mass);
        //If we have a gen particle, we estimate TOF by 1/B.
        //If we have a reco particle, we estimate B by 1/TOF
        if( !isGen ){
          beta = 1/TOF;
          distList[tofKey].distribution->Fill( TOF );
          distList.find(tofKey)->second.data.push_back( TOF );
        }
        else{
          beta = P / E;
          TOF = 1/beta;
          distList[tofKey].distribution->Fill( TOF );
          distList.find(tofKey)->second.data.push_back( TOF );
        }
        
        //Calculate gamma
        gamma = 1.0 / TMath::Sqrt( 1 - beta*beta );

        //Fill the gamma distribution
        Key gammaKey (string("gamma"), type, (int)(3* *charge), (int)*mass);
        distList[gammaKey].distribution->Fill(gamma);
        distList.find(gammaKey)->second.data.push_back(gamma);

        //Fill beta distribution
        Key betaKey (string("beta"), type, (int)(3* *charge), (int)*mass);
        distList[betaKey].distribution->Fill( beta );
        distList.find(betaKey)->second.data.push_back(beta);
        
        //Create a particle to contain the information.
        Key particleKey (string("particle"), type, (int)(3* *charge), (int)*mass);
        
        //Note: here, P and Pt are NOT corrected for the reco particles, so we're okay.
        //Note: we have a corrected and uncorrected mass for each particle. Maybe do the same with P and Pt?
        //Would have to change the P Vs P plot(s).
        //Always giving the correct particle energy though (corrected mass and p)
        if( !isGen ){
          events[particleKey][(int)event].push_back(Particle(Pt,P,theta,eta,phi,E,Ih,recoMassCorr,recoMassUncorr,beta,TOF,PDG,event));
        }
        else{
          events[particleKey][(int)event].push_back(Particle(Pt,P,theta,eta,phi,E,Ih,Mass,Mass,beta,TOF,PDG,event));
        }
          
      }//end event loop
      datFile->Close();
    }//end file loop
  }//end types loop
  //Loop through all available distributions, writing them to the disc.

  //Data output file
  string outFileName = fmt::format("HSCP_MC_Analysis_CC_{}_CT_{}.root",CC,CT);
  TFile *outFile = new TFile(outFileName.c_str(),"RECREATE");
  //Set up the outfile to have a folder for each data type
  AllocateOutfile(outFile, distNames);
  
  //Loop through each standard distribution (eta,beta,gamma,etc)
  WriteStandardDistributions( outFile, types, distList );

  // //Create stacked beta distributions
  WriteGenFracTrackVsBeta(outFile, distList, massCounts, distProps, CC, CT);

  // //Create Ih vs P Distribution for gen, reco,
  WriteIhVsP( outFile, distList, chargeCounts, events, CC, CT);

  // //Create the PReco/Q vs PGen plot
  WritePrecoVsPgen( outFile, distList, chargeCounts, events, CC, CT);

  // Create plots of reco mass vs gen mass. One assuming correct Q, on
  WriteMrecoVsMgen( outFile, distList, chargeCounts, massCounts, events, CC, CT);

  
  // //Same as above but for beta
  WriteBrecoVsBgen( outFile, distList, chargeCounts, massCounts, events, CC, CT);
  
  // //Beta Versus Eta
  WriteBgenVsEtagen( outFile, distList, massCounts, events, CC, CT);

  
  cout << "Finished writing to file" << endl;
  outFile->Close();
  theApp.Run();
  return 0;
}
