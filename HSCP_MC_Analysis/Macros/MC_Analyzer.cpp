#include <cstddef> //size_t
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <string>
#include <cmath>
#include <map>
#include <unordered_map>
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


/*
  Function to insert available ROOT MC Sample files into a string array.
  Make the first entry in the array a string containing the number of 
  entries in the array, INCLUDING this entry.
*/

static string *ListOfFiles(){
  FILE *inStream;
  char charnFiles[10];
  char charFiles[10000];
  stringstream streamnFiles;
  
  //First, count how many files we will be reading
  if(!(inStream = popen("ls /home/austin/HSCP_Study/HSCP_MC_Files/Histos_mchamp*.root -l | wc -l", "r"))){
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
  //PLUS the first entry containing the number of entries in the arary.
  
  listOfFiles[0] = streamnFiles.str();
  
  //Now we need to read in the list of files
  //Will do something similar to the above, but using this method to get only the file names
  stringstream testStream;
  if(!(inStream = popen("find /home/austin/HSCP_Study/HSCP_MC_Files/Histos_mchamp*.root -printf \"%f\n\"","r"))){
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
  //const string type = "Reco";
  bool isGen = false;
  for( const auto &iType : types){
    for (const auto &iNames : distNames) {
      if (isGen){
        if (iNames == "I" || iNames == "Ih" || iNames == "dZ" || iNames == "dXY" || iNames == "dR" || iNames == "hasMuon")
          continue;
      }
      else{
        if (iNames == "charge") continue;
      }
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
      currDist->Scale(norm/currDist->Integral("width"));
      //currDist->Draw("P");
      currDist->Write();
    }
  }

}


void WriteGenFracTrackVsBeta(TFile *&outFile, const distMap &distList, const map<double,int> &massCounts, const map<string,distProp> &distProps){
  //Have to integrate between two beta values in order to get the fraction of particles
  //that lie between those values (assuming that the primitive distribution is normalized to unity)
  //Only want to make the plot with unity charges (key for charge = 3)
  const array<int,12> colorList = {1,2,3,4,6,41,34,46,12,8,14,5};
  const array<int,12> markerList = {2,3,4,5,25,26,27,20,21,22,33,34};

  const double CHARGE = 1.0; //Desired constant charge for the plots, in units of actual charge (not e/3)

  int nBins = distProps.find(string("beta"))->second.nBins;
  double betaMin = distProps.find(string("beta"))->second.axisLimits.first; 
  double betaMax = distProps.find(string("beta"))->second.axisLimits.second;
  double BETASTEP = (betaMax-betaMin)/nBins; //change to math limits and step in main
  
  unsigned int nPoints = TMath::Ceil((betaMax - betaMin) / BETASTEP);
  
  TMultiGraph *LowMassOutGraph = new TMultiGraph("LowMassfracPart_beta","LowMassfracPart_beta");
  TLegend *LowMassLegend = new TLegend(0.2,0.65,0.5,0.85);
  //LowMassLegend->SetTextFont(72);
  LowMassLegend->SetTextSize(0.025);
  LowMassLegend->SetFillColor(0);
  TCanvas *LowMassCanvas = new TCanvas("low","low",500,500);
                                       
  TMultiGraph *HighMassOutGraph = new TMultiGraph("HighMassfracPart_beta","HighMassfracPart_beta");
  TLegend *HighMassLegend = new TLegend(0.5,0.75,0.78,0.90);
  //HighMassLegend->SetTextFont(72);
  HighMassLegend->SetTextSize(0.025);
  HighMassLegend->SetFillColor(0);
  TCanvas *HighMassCanvas = new TCanvas("high","high",500,500);
  
  TH1F *currDist = NULL;
  TGraph *currGraph = NULL;
  const double NORM = 1;

  //Make a vector of TGraph objects whose constituents will compose the TMultiGraph that is saved to disc
  vector<TGraph*> LowMassGraphs;
  vector<TGraph*> HighMassGraphs;
  vector<string> LowMassNames;
  vector<string> HighMassNames;
  //These are the steps that will be taken during the beta integrations, and for the distance
  //between the points in the resulting distribution
  //Fill an array of beta values based on the step size and range of the beta distributions
  
  vector<double> xPoints; //Vectors automatically order increasing from min
  for (double iP = betaMin; iP <= betaMax+BETASTEP; iP+=BETASTEP) {
    xPoints.push_back(iP);
  }

  vector<double> LowMassYPoints;
  vector<double> HighMassYPoints;
  
  TAxis *xAxis = NULL;
  int bMin, bMax;
  
  //Loop through all available masses
  auto iMass = massCounts.begin();
  auto iColor = colorList.begin();
  auto iMarker = markerList.begin();
  for (; iMass != massCounts.end(); ++iMass,++iColor,++iMarker){
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
    xAxis = currDist->GetXaxis();

    if( isLowMass ){
      LowMassNames.push_back( fmt::format("Q={}e M={} GeV/",CHARGE,(int)iMass->first)+string("c^{2}"));
    }
    else{
      HighMassNames.push_back( fmt::format("Q={}e M={} GeV/",CHARGE,(int)iMass->first)+string("c^{2}"));
    }
    
    //Calculate the fraction of tracks at each point. For each point, the marker represents the fraction
    //of tracks between the previous point and the current point
    if (isLowMass)
      LowMassYPoints.push_back(0.0);
    else
      HighMassYPoints.push_back(0.0);
    
    for (auto iPoint = xPoints.begin(); iPoint != xPoints.end(); ++iPoint) {
      if (iPoint == xPoints.begin()) continue;
      double Prev = *prev(iPoint,1);
      bMax = xAxis->FindBin(*iPoint);

      if( isLowMass ) 
        LowMassYPoints.push_back(currDist->GetBinContent(bMax));
      else
        HighMassYPoints.push_back(currDist->GetBinContent(bMax));
    }
   
    xPoints.shrink_to_fit();

    if (isLowMass)
      LowMassYPoints.shrink_to_fit();
    else
      HighMassYPoints.shrink_to_fit();

    //Now have all x and y points for the current mass distribution
    //Need to convert these vectors to TVectors for entry into
    //TGraphs
    TVectorD root_xPoints(xPoints.size(),&xPoints[0]);

    TVectorD *root_LowMassYPoints;
    TVectorD *root_HighMassYPoints;

    if (isLowMass){
      root_LowMassYPoints = new TVectorD(LowMassYPoints.size(),&LowMassYPoints[0]);
      currGraph = new TGraph(root_xPoints,*root_LowMassYPoints);
    }
    else{
      root_HighMassYPoints = new TVectorD(HighMassYPoints.size(),&HighMassYPoints[0]);
      currGraph = new TGraph(root_xPoints,*root_HighMassYPoints);
    }
    currGraph->SetMarkerColorAlpha(*iColor,0.8);
    currGraph->SetLineColor(*iColor);
    currGraph->SetLineWidth(2);
    currGraph->SetFillColorAlpha(*iColor,0.8);
    currGraph->SetMarkerStyle(*iMarker);

    if (isLowMass){
      LowMassGraphs.push_back(currGraph);
      LowMassYPoints.clear();
    }
    else{
      HighMassGraphs.push_back(currGraph);
      HighMassYPoints.clear();
    }
    
  }
  //Loop through the constituents, adding them to the multi
  auto iDistLow = LowMassGraphs.begin();
  auto iLNames = LowMassNames.begin();
  for (; iDistLow != LowMassGraphs.end(); ++iDistLow, ++iLNames){
    LowMassOutGraph->Add(*iDistLow,"AFC");
    LowMassLegend->AddEntry(*iDistLow, iLNames->c_str(), "PL");
  }

  auto iDistHigh = HighMassGraphs.begin();
  auto iHNames = HighMassNames.begin();
  for (; iDistHigh != HighMassGraphs.end(); ++iDistHigh, ++iHNames){
    HighMassOutGraph->Add(*iDistHigh,"AFC");
    cout << iHNames->c_str() << endl;
    HighMassLegend->AddEntry(*iDistHigh, iHNames->c_str(), "PL");
  }


  HighMassCanvas->cd();
  HighMassOutGraph->Draw("ACP");
  gPad->Update();
  HighMassOutGraph->GetXaxis()->SetTitle("#beta");
  HighMassOutGraph->GetYaxis()->SetTitle(fmt::format("Fraction of Particles/{}",BETASTEP).c_str());
  HighMassOutGraph->SetMaximum(4.0);
  gPad->Update();
  HighMassLegend->Draw();

  LowMassCanvas->cd();
  LowMassOutGraph->Draw("ACP");
  gPad->Update();
  LowMassOutGraph->GetXaxis()->SetTitle("#beta");
  LowMassOutGraph->GetYaxis()->SetTitle(fmt::format("Fraction of Particles/{}",BETASTEP).c_str());
  LowMassOutGraph->SetMaximum(8.0);
  gPad->Update();
  LowMassLegend->Draw();


  outFile->cd();
  LowMassCanvas->Write();
  LowMassCanvas->Print("LowMassPartFrac.pdf","pdf");
  HighMassCanvas->Write();
  HighMassCanvas->Print("HighMassPartFrac.pdf","pdf");
}


void WriteIhVsP(TFile *&outFile, const distMap &distList, const map<double,int> &chargeCounts, const map<string,distProp> &distProps){
  const array<int,12> colorList = {1,2,3,4,6,41,34,46,12,8,14,5};
  const array<int,12> markerList = {2,3,4,5,25,26,27,20,21,22,33,34};

  const int MASS = 900; //Desired mass to vary charges over

  //TH2F *outDist; //= new TH2F("IhVsP","IhVsP",200,0,2000,70,0,35);
  TMultiGraph *outGraph = new TMultiGraph("mult","");
  TGraph* currGraph;
  TLegend *distLegend = new TLegend(0.5,0.75,0.80,0.93);
  distLegend->SetTextSize(0.025);
  distLegend->SetFillColor(0);
  
  vector<string> chargeNames;

  //Check to see if this mass is available with the desired charge
  auto iCharge = chargeCounts.begin();
  auto iColor = colorList.begin();
  for( ; iCharge!=chargeCounts.end(); iCharge++) {
    Key pKey = Key (string("momentum"), string("Reco"), (int)(iCharge->first), MASS);
    Key ihKey = Key (string("Ih"), string("Reco"), (int)(iCharge->first), MASS);
    unsigned int NPOINTS = distList.find(pKey)->second.data.size();
    if ( NPOINTS < 20 ){
      cout << fmt::format("Bad file at mchamp{}_M_{}",iCharge->first,MASS) << endl;
      continue;
    }
    cout << NPOINTS << endl;
    //map.count returns zero if the item is not found. want to skip masses that are not with desired charge  
    if (distList.count(pKey) == 0){
      continue;
    }    
   
    chargeNames.push_back(fmt::format("Q={}",iCharge->first)+string("#frac{e}{3}"));
    
    string name = fmt::format("{}",iCharge->first);
    currGraph = new TGraph(distList.find(pKey)->second.data.size());
    currGraph->SetMarkerStyle(20);
    currGraph->SetMarkerColorAlpha(*iColor,0.8);
    currGraph->SetMarkerSize(0.3);
    
    //Now we want to loop through all momenta and Ih for this charge and mass
    //Have the data stored in distList
    auto iP = distList.find(pKey)->second.data.begin();
    auto iIh = distList.find(ihKey)->second.data.begin();
    unsigned int i = 0;
    for( ; iP != distList.find(pKey)->second.data.end(); ++iP, ++iIh ){
      currGraph->SetPoint(i,*iP,*iIh);
      ++i;
    }
    outGraph->Add(currGraph);
    distLegend->AddEntry(currGraph, (fmt::format("Q={}",iCharge->first)+string("e/3")).c_str(), "P");
    
    ++iColor;
    
  }
  //Now all distributions have been added to their own TH2F's
  //They have also been added to the THStack.
  TCanvas *outCanvas = new TCanvas("IhVsP","IhVsP",500,500);
  outCanvas->cd();
  outGraph->Draw("ap");
  gPad->Update();
  distLegend->Draw();
  gPad->Update();
  outFile->cd();
  outCanvas->Write();
  outCanvas->Print("IhVsP.pdf","pdf");
}

//NEED TO MATCH PARTICLES FOR THESE
//Reco Pt(scaled by 1/Q) div Gen Pt vs Q for a given mass
//Reco beta div gen beta vs Q for a given mass

/*Main*/
int main(int argc, char **argv){
  TApplication theApp("App",0,0);

  // gROOT->SetStyle("Pub");
  SetRootStyle();
  //gROOT->SetBatch(kTRUE); //Don't draw things when created.

  
  string *fileList = ListOfFiles(); 

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
  vector<string> distNames = {
    "beta",
    "energy",
    "eta",
    "gamma",
    "Ih",
    "momentum",
    "phi",
    "trans_momentum",
    "theta",
    "TOF",
    "TOFdivTRel",
    "mass",
    "charge"
  };

  //Right now, sets the same limits and bins for both reco and gen particles.
  map<string,distProp> distProps = {
    {"beta",{limits(0.0,1.0),20,axes("#beta","events/{}")}},
    {"energy",{limits(0,2200),200,axes("E [GeV]","events/{} GeV")}},
    {"eta",{limits(-3,3),30,axes("#eta", "events/{}")}},
    {"gamma",{limits(0,100.0),200,axes("#gamma", "events/{}")}},
    {"Ih",{limits(0,35),70,axes("I_{h} [MeV/cm]", "events/{} MeV/cm")}},
    {"momentum",{limits(0.0,2000),200,axes("P [GeV/c]", "events/{} GeV/c")}},
    {"phi",{limits(-200,200),100,axes("#phi [deg]", "events/{} deg")}},
    {"trans_momentum",{limits(0,2000),200,axes("P_{t} [GeV/c]", "events/{} GeV/c")}},
    {"theta",{limits(0,180),180,axes("#theta [deg]", "events/{} deg")}},
    {"TOF",{limits(0,7),70,axes("t [s]", "events/{} s")}},
    {"TOFdivTRel",{limits(0.9,1.0),200,axes("#frac{t}{t_{r}}", "events/{}")}},
    {"mass",{limits(0,1000),100,axes("mass [GeV/c^{2}]", "events/{} GeV/c^2")}},
    {"charge",{limits(0,40),40,axes("charge [e/3]", "events/{} e/3")}}
  };

  vector<string> types = {"Gen","Reco"};
  
  AllocateDistributions(distList,charges,masses,numFiles,distNames,distProps, types);

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
  float TOF; 
  float reco_Mass;
  double gen_Mass;
  float dZ, dXY, dR;
  float reco_eta, reco_phi;
  double gen_eta, gen_phi;
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
      string file = "/home/austin/HSCP_Study/HSCP_MC_Files/";
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

        //Calculate theta
        theta = 2*TMath::ATan( TMath::Exp( -1 * (eta) ) );
        thetaDeg = theta * (180.0/TMath::Pi());
        
        //Fill the theta distribution
        Key thetaKey (string("theta"), type, (int)(3* *charge), (int)*mass);
        distList[thetaKey].distribution->Fill(thetaDeg);
        distList[thetaKey].data.push_back(thetaDeg);
        
        //fill the eta distribution
        Key etaKey (string("eta"), type, (int)(3* *charge), (int)*mass);
        distList[etaKey].distribution->Fill((eta));
        distList[etaKey].data.push_back(eta);

        //fill the phi distribution
        Key phiKey (string("phi"), type, (int)(3* *charge), (int)*mass);
        phi *= (180.0/TMath::Pi());
        distList[phiKey].distribution->Fill( phi );
        distList[phiKey].data.push_back(phi);
      
        //Calculate Momentum from transverse and theta (rad)
        P = (Pt) / TMath::Sin(theta);
        //Fill momentum distribution
        Key pKey (string("momentum"), type, (int)(3* *charge), (int)*mass);
        distList[pKey].distribution->Fill(P);
        distList[pKey].data.push_back(P);

        //Fill the transverse momentum distribution
        Key ptKey (string("trans_momentum"), type, (int)(3* *charge), (int)*mass);
        distList[ptKey].distribution->Fill( Pt );
        distList[ptKey].data.push_back( Pt );
      
        //Calculate the relativistic energy
        E = TMath::Sqrt( P*P + TMath::Power(Mass,2) );
        //Fill the energy distribution
        Key enKey (string("energy"), type, (int)(3* *charge), (int)*mass);
        distList[enKey].distribution->Fill(E);
        distList[enKey].data.push_back( E );

        //Calculate gamma
        gamma = 1.0 / TMath::Sqrt( 1 - beta*beta );

        //Fill the gamma distribution
        Key gammaKey (string("gamma"), type, (int)(3* *charge), (int)*mass);
        distList[gammaKey].distribution->Fill(gamma);
        distList[gammaKey].data.push_back( gamma );

        //Fill the Ih distribution
        Key ihKey (string("Ih"), type, (int)(3* *charge), (int)*mass);
        distList[ihKey].distribution->Fill( Ih );
        distList[ihKey].data.push_back( Ih );

        //Fill the tof key
        Key tofKey (string("TOF"), type, (int)(3* *charge), (int)*mass);
        distList[tofKey].distribution->Fill( TOF );
        distList[tofKey].data.push_back( TOF );

        //Fill the mass distribution
        Key massKey (string("mass"), type, (int)(3* *charge), (int)*mass);
        distList[massKey].distribution->Fill( Mass );
        distList[massKey].data.push_back( Mass );
        
        //Calculate beta
        beta = P / E;
        //Fill beta distribution
        Key betaKey (string("beta"), type, (int)(3* *charge), (int)*mass);
        //Only want the gen particle HSCP's 
        if( !isGen || (TMath::Abs(PDG) != 17) )
          continue;
        else{
          distList[betaKey].distribution->Fill( beta );
          distList[betaKey].data.push_back( beta );
        }

        
          
      }//end event loop
      datFile->Close();
    }//end file loop
  }//end types loop
  //Loop through all available distributions, writing them to the disc.

  //Data output file
  TFile *outFile = new TFile("HSCP_MC_Analysis.root","RECREATE");
  //Set up the outfile to have a folder for each data type
  AllocateOutfile(outFile, distNames);
  
  //Loop through each standard distribution (eta,beta,gamma,etc)
  WriteStandardDistributions( outFile, types, distList );

  WriteGenFracTrackVsBeta(outFile, distList, massCounts, distProps);
  //Create Ih vs P Distribution for gen, reco,


  //WriteIhVsP(TFile *&outFile, const distMap &distList, const map<double,int> &chargeCounts, const map<string,distProp> &distProps){
  WriteIhVsP( outFile, distList, chargeCounts, distProps );

  
  cout << "Finished writing to file" << endl;
  outFile->Close();
  theApp.Run();
  return 0;
}
