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
  distObjects() : canvas(NULL), distribution(NULL) {}
  distObjects(const distObjects &arg) : canvas(arg.canvas), distribution(arg.distribution) {}
  distObjects(TCanvas *&aCanv, TH1F *&aDist) : canvas((TCanvas*)aCanv->Clone()), distribution((TH1F*)aDist->Clone()) {}
};

typedef map<Key,distObjects> distMap;
typedef pair<string,string> axes;
typedef pair<double,double> limits;

struct distProp{
  limits axisLimits;
  double nBins;
  axes axisTitles;
};

//typedef map< Key, pair<TCanvas, TH1F> > distMap;

// TCanvas *canv_Preco_Vs_Pgen[nCHARGE][nMASS];
// TCanvas *canv_Tofreco_Vs_Tofgen[nCHARGE][nMASS];
// TCanvas *canv_Ihreco_Vs_Ihgen[nCHARGE][nMASS];

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
  if(!(inStream = popen("ls /home/austin/HSCP_Study/HSCP_MC_Files/*.root -l | wc -l", "r"))){
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
  if(!(inStream = popen("find /home/austin/HSCP_Study/HSCP_MC_Files/*.root -printf \"%f\n\"","r"))){
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
    if( !massCounts.insert( make_pair(outDat->mass[iMass], 1 ) ).second ){
      //Element is alread present. Need to increment the count
      massCounts[outDat->mass[iMass]] += 1;
    }
  }

  //Count the charges
  for( int iCharge = 0; iCharge < nFiles; iCharge++ ){
    if( !chargeCounts.insert( make_pair(outDat->charge[iCharge], 1 ) ).second ){
      //Element is alread present. Need to increment the count
      chargeCounts[outDat->charge[iCharge]] += 1;
    }
  }

  outDat->chargeCounts = chargeCounts;
  outDat->massCounts = massCounts;
  return outDat;
}


//For each distribution type, go through all files (masses and charges)
static void AllocateDistributions( distMap &argDists, const double *charges, const double *masses, const int &numFiles, const vector<string> &distNames, map<string,distProp> &distProps){
  const string type = "Reco";
  for (auto iNames = distNames.begin(); iNames != distNames.end(); ++iNames) {
    double *lowerLim = &distProps[*iNames].axisLimits.first;
    double *upperLim = &distProps[*iNames].axisLimits.second;
    double *nBins = &distProps[*iNames].nBins;
    string xAxis(distProps[*iNames].axisTitles.first);
    string yAxis(distProps[*iNames].axisTitles.second);
    
    for (int iFile=0; iFile < numFiles; iFile++) {
      Key entryKey (*iNames,type,(int)(3* charges[iFile]),(int)masses[iFile]);
      
      string canv_name = entryKey.name + "_canv";
      string dist_name = entryKey.name + "_dist";

      yAxis = fmt::format(yAxis.c_str(), (*upperLim-*lowerLim) / *nBins);
      TCanvas *tempCanv = new TCanvas(canv_name.c_str(),canv_name.c_str(),500,500);
      TH1F *tempHist = new TH1F(dist_name.c_str(),dist_name.c_str(),*nBins,*lowerLim,*upperLim);

      tempHist->GetXaxis()->SetTitle(xAxis.c_str());
      tempHist->GetYaxis()->SetTitle(yAxis.c_str());

      //tempHist->SetLineWidth(2); //Better visibility
      
      argDists.emplace(entryKey,distObjects(tempCanv,tempHist));
    }
  }
}

static void AllocateOutfile( TFile *&aFile, const vector<string> &distNames ){
  aFile->mkdir("Reco");
  aFile->mkdir("Gen");
  for (const auto &iNames : distNames) {
    aFile->cd("Reco");
    gDirectory->mkdir(iNames.c_str());
    aFile->cd("Gen");
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

/*Main*/
int main(int argc, char **argv){
  //TApplication theApp("App",&argc,argv);

  // gROOT->SetStyle("Pub");
  SetRootStyle();
  gROOT->SetBatch(kTRUE); //Don't draw things when created.

  
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
  map<double,int> massCounts = fileNameData->chargeCounts;

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
    "mass"
  };
  
  map<string,distProp> distProps = {
    {"beta",{limits(0.9,1.0),200,axes("#beta","events/{}")}},
    {"energy",{limits(0,2200),200,axes("E [GeV]","events/{} GeV")}},
    {"eta",{limits(0.9,1.0),200,axes("#eta", "events/{}")}},
    {"gamma",{limits(0,100.0),200,axes("#gamma", "events/{}")}},
    {"Ih",{limits(0.9,1.0),200,axes("I_{h} [MeV/cm]", "events/{} MeV/cm")}},
    {"momentum",{limits(0.0,2000),200,axes("P [GeV/c]", "events/{} GeV/c")}},
    {"phi",{limits(0.9,1.0),200,axes("#phi [deg]", "events/{} deg")}},
    {"trans_momentum",{limits(0.9,1.0),200,axes("P_{t} [GeV/c]", "events/{} GeV/c")}},
    {"theta",{limits(0,180),180,axes("#theta [deg]", "events/{} deg")}},
    {"TOF",{limits(0.9,1.0),200,axes("t [s]", "events/{} s")}},
    {"TOFdivTRel",{limits(0.9,1.0),200,axes("#frac{t}{t_{r}}", "events/{}")}},
    {"mass",{limits(0,1000),100,axes("mass [GeV/c^{2}]", "events/{} GeV/c^2")}}
  };

  
  AllocateDistributions(distList,charges,masses,numFiles,distNames,distProps);
  
  
  //Pointer for the tree being read
  TTree *tree;
  unsigned int event;
  unsigned int Hscp;
  float E;
  float Pt, P;
  float I;
  float Ih; //This is the energy deposition that we are interested in 
  float TOF; 
  float Mass; //Reco mass. Not correct. Assumes unit charge
  float dZ, dXY, dR;
  float eta, phi, theta;
  float beta, gamma;
  bool hasMuon;

  //Loop over the files
  for( int iFile = 0; iFile < numFiles; iFile++ ){
    string file = "/home/austin/HSCP_Study/HSCP_MC_Files/";
    file += fileList[iFile];
    TFile *datFile = new TFile(file.c_str(), "READ"); //Open the current file

    double *charge = &charges[iFile];
    double *mass = &masses[iFile];
    string type = "Reco";
    //Create a string for the appropriate file directory
    ///storage/6/work/askeeters/HSCPStudy/HSCP_Root_Files
    string mcDir(fmt::format("mchamp{}_M_{}/HscpCandidates",3*charges[iFile],masses[iFile]));

    //Set the appropriate branches for the MC particles
    tree = (TTree*)datFile->Get(mcDir.c_str());
    
    tree->SetBranchAddress("Event", &event);
    tree->SetBranchAddress("Hscp", &Hscp);
    tree->SetBranchAddress("Pt",&Pt);
    tree->SetBranchAddress("I",&I);
    tree->SetBranchAddress("Ih",&Ih);
    tree->SetBranchAddress("TOF",&TOF);
    tree->SetBranchAddress("Mass",&Mass);
    tree->SetBranchAddress("dZ",&dZ);
    tree->SetBranchAddress("dXY",&dXY);
    tree->SetBranchAddress("dR",&dR);
    tree->SetBranchAddress("eta",&eta);
    tree->SetBranchAddress("phi",&phi);
    tree->SetBranchAddress("hasMuon",&hasMuon);
    
    //Loop over the events
    //cout << fmt::format("{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}", "event", "Pt(GeV)", "Ih(MeV/cm)", "Mass(GeV)", "eta", "phi", "theta(deg)","E") << endl;
    for(int iEvt = 0; iEvt < tree->GetEntries(); iEvt++){
      tree->GetEntry(iEvt);
      
      //Calculate theta
      theta = 2*TMath::ATan( TMath::Exp( -1 * eta ) );
      //Fill the theta distribution
      Key thetaKey (string("theta"), type, (int)(3* *charge), (int)*mass);
      distList[thetaKey].distribution->Fill(theta * (180.0/TMath::Pi()));
      
      //Calculate Momentum from transverse and theta (rad)
      P = Pt / TMath::Sin(theta);
      //Fill momentum distribution
      Key pKey (string("momentum"), type, (int)(3* *charge), (int)*mass);
      distList[pKey].distribution->Fill(P);
      
      //Calculate the relativistic energy
      E = TMath::Sqrt( P*P + Mass*Mass );
      //Fill the energy distribution
      Key enKey (string("energy"), type, (int)(3* *charge), (int)*mass);
      distList[enKey].distribution->Fill(E);

      //Calculate beta
      beta = P / E;
      //Fill beta distribution
      Key betaKey (string("beta"), type, (int)(3* *charge), (int)*mass);
      distList[betaKey].distribution->Fill(beta);
      

      //Calculate gamma
      gamma = 1.0 / TMath::Sqrt( 1 - beta*beta );

      //Fill the gamma distribution
      Key gammaKey (string("gamma"), type, (int)(3* *charge), (int)*mass);
      distList[gammaKey].distribution->Fill(gamma);

      
    }//end event loop
    datFile->Close();
  }//end file loop
  
  //Loop through all available distributions, writing them to the disc.
  TCanvas *currCanv;
  TH1F *currDist;
  double norm = 1;

  //Data output file
  TFile *outFile = new TFile("HSCP_MC_Analysis.root","RECREATE");
  //Set up the outfile to have a folder for each data type
  AllocateOutfile(outFile, distNames);
  
  //Loop through each distribution (eta,beta,gamma,etc)
  outFile->cd("Reco");
  for(const auto &distIt : distList){
    string dir = string("Reco/") + distIt.first.dist;
    outFile->cd(dir.c_str());
    currDist = distIt.second.distribution;
    currDist->Scale(norm/currDist->Integral("width"));
    currDist->Write();
    outFile->cd("/");
  }

  // for(const auto &distIt : GammaDist){
  //   currDist = distIt.second.second;
  //   currDist->Scale(norm/currDist->Integral("width"));
  //   currDist->Write();
  // }
  
  cout << "Finished writing to file" << endl;
  outFile->Close();
  //theApp.Run();
  return 0;
}
