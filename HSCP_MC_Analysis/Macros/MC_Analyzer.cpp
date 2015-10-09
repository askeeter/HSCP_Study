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
  float Pt,P,Theta,Eta,Phi,E,Ih,Mass;
  int Event,PDG;
  Particle() : Pt(0),P(0),Theta(0),Eta(0),Phi(0),E(0),Ih(0),Mass(0),Event(-1),PDG(-999) {}
  Particle(const float Pt, const float P, const float Theta, const float Eta, const float Phi, const float E, const float Ih, const float Mass, const int PDG, const int Event) : Pt(Pt),P(P),Theta(Theta),Eta(Eta),Phi(Phi),E(E),Ih(Ih),Mass(Mass),PDG(PDG),Event(Event) {}
  Particle(const Particle &arg) : Pt(arg.Pt),P(arg.P),Theta(arg.Theta),Eta(arg.Eta),Phi(arg.Phi),E(arg.E),Ih(arg.Ih),Mass(arg.Mass),PDG(arg.PDG),Event(arg.Event) {}
};

typedef map< Key, map< int, vector< Particle > > > Events;

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
  TLegend *HighMassLegend = new TLegend(0.5,0.65,0.9,0.85);
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
    currDist->SetFillColorAlpha(*iColor,0.75);
    
    
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
  LowMassOutStack->Draw("C");
  LowMassOutStack->GetYaxis()->SetTitle(fmt::format("Fraction of Tracks/{}",per).c_str());
  LowMassOutStack->GetXaxis()->SetTitle("Gen #beta");
  gPad->Update();
  LowMassLegend->Draw("F");
  gPad->Update();
  outFile->cd();
  LowMassCanvas->Write();
  LowMassCanvas->Print("LMGenFracTracVsBeta.pdf","pdf");

  HighMassCanvas->cd();
  HighMassOutStack->Draw("C");
  HighMassOutStack->GetYaxis()->SetTitle(fmt::format("Fraction of Tracks/{}",per).c_str());
  HighMassOutStack->GetXaxis()->SetTitle("Gen #beta");
  gPad->Update();
  HighMassLegend->Draw("F");
  gPad->Update();
  outFile->cd();
  HighMassCanvas->Write();
  HighMassCanvas->Print("HMGenFracTracVsBeta.pdf","pdf");
  
}


void WriteIhVsP(TFile *&outFile, const distMap &distList, const map<double,int> &chargeCounts){
  const array<int,12> colorList = {1,2,3,4,6,41,34,46,12,8,14,5};
  const array<int,12> markerList = {2,3,4,5,25,26,27,20,21,22,33,34};

  const int MASS = 900; //Desired mass to vary charges over

  //TH2F *outDist; //= new TH2F("IhVsP","IhVsP",200,0,2000,70,0,35);
  TMultiGraph *outGraph = new TMultiGraph("mult","");
  TGraph* currGraph;
  TLegend *distLegend = new TLegend(0.5,0.55,0.80,0.85);
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
    if ( NPOINTS < 100 ){
      cout << fmt::format("Bad file at mchamp{}_M_{}",iCharge->first,MASS) << endl;
      continue;
    }
    //map.count returns zero if the item is not found. want to skip masses that are not with desired charge  
    if (distList.count(pKey) == 0){
      continue;
    }    
   
    chargeNames.push_back(fmt::format("Q={}",iCharge->first)+string("#frac{e}{3}"));
    
    string name = fmt::format("Ih{}",iCharge->first);
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
    distLegend->AddEntry(currGraph, (fmt::format("Q: {} M: {}",iCharge->first,MASS)+string(" GeV/c^{2}")).c_str(), "P");
    
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

void WritePrecoVsPgen(TFile *&outFile, const distMap &distList, const map<double,int> &chargeCounts, const Events &events){
  const array<int,12> colorList = {1,2,3,4,6,41,34,46,12,8,14,5};
  const array<int,12> markerList = {2,3,4,5,25,26,27,20,21,22,33,34};
  const int MASS = 900; //Desired mass to vary charges over
  THStack *outDist = new THStack("PrVsPg","");
  THStack *outDistScaled = new THStack("PrVsPg_Scaled","");
  
  TF1 *currLine = 0;
  TH2F *currDist, *currDistScaled= 0;
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
  outDist->GetXaxis()->SetTitle("P_{t g} [GeV/c]");
  outDist->GetYaxis()->SetTitle("P_{t r} [GeV/c]");
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
  outCanvas->Print("PrVsPg.pdf","pdf");


  outCanvasScaled->cd();
  outDistScaled->Draw("BOX1");
  outDistScaled->GetXaxis()->SetTitle("P_{t g} [GeV/c]");
  outDistScaled->GetYaxis()->SetTitle("Q #times P_{t r} [GeV/c]");
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
  outCanvasScaled->Print("PrVsPg_Scaled.pdf","pdf");
}

void WriteBrecoVsBgen(TFile *&outFile, const distMap &distList, const map<double,int> &chargeCounts, const Events &events){

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
  
  vector<TH2F*> outDists;
  vector<TCanvas*> outCanvases;
  vector<string> chargeNames;

  
  TH2F *currBetaDist = 0;
  
  //Check to see if this mass is available with the desired charge
  auto iCharge = chargeCounts.begin();
  auto iColor = colorList.begin();
  for( ; iCharge!=chargeCounts.end(); iCharge++) {
    Key recoKey (string("particle"), "Reco", (int)(iCharge->first), MASS);
    Key genKey (string("particle"), "Gen", (int)(iCharge->first), MASS);
    //iCharge->fisrt is in units of e/3
    //map.count returns zero if the item is not found. want to skip masses that are not with desired charge. Also, only want files with > 200 HSCPs
    if (events.count(genKey) == 0 || events[recoKey].size() < 200){
      continue;
    }    
    
    chargeNames.push_back(fmt::format("{}",iCharge->first));
    
    string name = fmt::format("Beta_Q_{}_M_{}",iCharge->first,MASS);
    currBetaDist = new TH2F(name.c_str(),name.c_str(),50,0,1.0,45,0,1.5);
    currBetaDist->SetFillColorAlpha(*iColor,0.75);
    
    for( const auto &iGenEvents : events[genKey] ){
      //Get the vector of all gen particles from this event
      auto genParticles = iGenEvents.second;
      //Loop over the gen particles from this event
      for( const auto &iGenParticle : genParticles ){
        //Skip all non gen HSCP from event
        if (TMath::Abs(iGenParticle.PDG)!=17) continue;
        for( const auto &iRecoParticle : events[recoKey][iGenEvents.first] ){
          currBetaDist->Fill(iGenParticle.P/iGenParticle.E,iRecoParticle.P/iRecoParticle.E);
        }//reco
      }//gen
    }//events

    outDists.push_back(currBetaDist);
    outCanvases.push_back(new TCanvas(fmt::format("{}",iCharge->first).c_str(),"",2000,2000));
    ++iColor;
  }//file loop

  auto iDists = outDists.begin();
  auto iCanvases = outCanvases.begin();
  auto iChargeName = chargeNames.begin();
  TF1 *currLine;
  currLine = new TF1("45","x",0,1);
  currLine->SetLineColor(kBlack);
  currLine->SetLineStyle(2);
  currLine->SetLineWidth(2);
  for(; iDists != outDists.end(); ++iDists,++iCanvases,++iChargeName){
    string name = fmt::format("BrVsBg_{}.pdf",(*iChargeName));
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
}


void WriteBgenVsEtagen(TFile *&outFile, const distMap &distList, const map<double,int> &massCounts, const Events &events){

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

  TH2F *currDist = 0;
  
  //Check to see if this mass is available with the desired charge
  auto iMass = massCounts.begin();
  auto iColor = colorList.begin();
  for( ; iMass!=massCounts.end(); iMass++) {
    Key recoKey (string("particle"), "Reco", CHARGE, iMass->first);
    Key genKey (string("particle"), "Gen", CHARGE, iMass->first);
    //iCharge->fisrt is in units of e/3
    //map.count returns zero if the item is not found. want to skip masses that are not with desired charge. Also, only want files with > 200 HSCPs
    if (events.count(genKey) == 0 || events[recoKey].size() < 200){
      continue;
    }    
    
    massNames.push_back(fmt::format("{}",iMass->first));
    string name = fmt::format("BetaVsEta_Q_{}_M_{}",iMass->first,CHARGE);
    currDist = new TH2F(name.c_str(),name.c_str(),40,-4,4,50,10,1);
    currDist->SetFillColorAlpha(*iColor,0.75);
    
    for( const auto &iGenEvents : events[genKey] ){
      //Get the vector of all gen particles from this event
      auto genParticles = iGenEvents.second;
      //Loop over the gen particles from this event
      for( const auto &iGenParticle : genParticles ){
        //Skip all non gen HSCP from event
        if (TMath::Abs(iGenParticle.PDG)!=17) continue;
        for( const auto &iRecoParticle : events[recoKey][iGenEvents.first] ){
          currDist->Fill(iGenParticle.Eta,iGenParticle.P/iGenParticle.E);
        }//reco
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
    string name = fmt::format("BgVsEtag_{}.pdf",(*iMassName));
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

  Events events;
  
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
    "mass",
    "charge",
  };

  //Right now, sets the same limits and bins for both reco and gen particles.
  map<string,distProp> distProps = {
    {"beta",{limits(0.0,1.5),45,axes("#beta","events/{}")}},
    {"energy",{limits(0,2200),200,axes("E [GeV]","events/{} GeV")}},
    {"eta",{limits(-3,3),30,axes("#eta", "events/{}")}},
    {"gamma",{limits(0,100.0),200,axes("#gamma", "events/{}")}},
    {"Ih",{limits(0,35),70,axes("I_{h} [MeV/cm]", "events/{} MeV/cm")}},
    {"momentum",{limits(0.0,2000),200,axes("P [GeV/c]", "events/{} GeV/c")}},
    {"phi",{limits(-200,200),100,axes("#phi [deg]", "events/{} deg")}},
    {"trans_momentum",{limits(0,2000),200,axes("P_{t} [GeV/c]", "events/{} GeV/c")}},
    {"theta",{limits(0,180),180,axes("#theta [deg]", "events/{} deg")}},
    {"TOF",{limits(0.0,4.0),20,axes("t [arb. units]", "events/{}")}},
    {"mass",{limits(0,1000),100,axes("mass [GeV/c^{2}]", "events/{} GeV/c^2")}},
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

      vector<Particle> partVec;
      
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
        //Fill momentum distribution
        Key pKey (string("momentum"), type, (int)(3* *charge), (int)*mass);
        distList[pKey].distribution->Fill(P);
        distList.find(pKey)->second.data.push_back(P);

        //Fill the transverse momentum distribution
        Key ptKey (string("trans_momentum"), type, (int)(3* *charge), (int)*mass);
        distList[ptKey].distribution->Fill( Pt );
        distList.find(ptKey)->second.data.push_back(Pt);
      
        //Calculate the relativistic energy
        E = TMath::Sqrt( P*P + TMath::Power(Mass,2) );
        //Fill the energy distribution
        Key enKey (string("energy"), type, (int)(3* *charge), (int)*mass);
        distList[enKey].distribution->Fill(E);
        distList.find(enKey)->second.data.push_back(E);

        //Calculate gamma
        gamma = 1.0 / TMath::Sqrt( 1 - beta*beta );

        //Fill the gamma distribution
        Key gammaKey (string("gamma"), type, (int)(3* *charge), (int)*mass);
        distList[gammaKey].distribution->Fill(gamma);
        distList.find(gammaKey)->second.data.push_back(gamma);

        //Fill the Ih distribution
        Key ihKey (string("Ih"), type, (int)(3* *charge), (int)*mass);
        distList[ihKey].distribution->Fill( Ih );
        distList.find(ihKey)->second.data.push_back(Ih);

       
        //Fill the mass distribution
        Key massKey (string("mass"), type, (int)(3* *charge), (int)*mass);
        distList[massKey].distribution->Fill( Mass );
        distList.find(massKey)->second.data.push_back(Mass);

        
        //Calculate beta
        beta = P / E;

        //Create a particle to contain the information.
        //Pt(0),P(0),Theta(0),Eta(0),Phi(0),E(0),Ih(0),Mass(0),Event(-1),PDG(-999) {}
        Key particleKey (string("particle"), type, (int)(3* *charge), (int)*mass);
        events[particleKey][(int)event].push_back(Particle(Pt,P,theta,eta,phi,E,Ih,Mass,PDG,event));
        //Fill beta distribution
        Key betaKey (string("beta"), type, (int)(3* *charge), (int)*mass);
        Key tofKey (string("TOF"), type, (int)(3* *charge), (int)*mass);
        distList[tofKey].distribution->Fill( 1.0 / beta );
        distList.find(tofKey)->second.data.push_back( 1.0 / beta );
        //Only want the gen particle HSCP's 
        if( !isGen || (TMath::Abs(PDG) != 17) ){
          continue;
        }
        else{
          //Estimating the time of flight by 1/B
          distList[betaKey].distribution->Fill( beta );
          distList.find(betaKey)->second.data.push_back(beta);
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
  WriteIhVsP( outFile, distList, chargeCounts );

  //Create the PReco/Q vs PGen plot
  //void WritePrecoVsPgen(TFile *&outFile, const distMap &distList, const map<double,int> &chargeCounts, const Events &events){
  WritePrecoVsPgen( outFile, distList, chargeCounts, events );

  //Same as above but for beta
  WriteBrecoVsBgen( outFile, distList, chargeCounts, events );

  //Beta Versus Eta
  WriteBgenVsEtagen( outFile, distList, massCounts, events );
  
  cout << "Finished writing to file" << endl;
  outFile->Close();
  theApp.Run();
  return 0;
}
