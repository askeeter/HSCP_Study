
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
#include "format.h"
#include "TROOT.h" //for gROOT
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TH1.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TMathBase.h"
#include "TMath.h"
#include "TNetFile.h"
#include "TAuthenticate.h"
#include "TObject.h"
#include "TApplication.h"

using namespace std;


/*Declare canvases for each plot*/
const unsigned int RECO = 0;
const unsigned int GEN = 1;
const unsigned int BKG = 2;
const unsigned int nTYPES = 3;

struct distObjects{
  TCanvas *canvas;
  TH1F *distribution;
  distObjects() : canvas(NULL), distribution(NULL) {}
  distObjects(const distObjects &arg) : canvas(arg.canvas), distribution(arg.distribution) {}
  distObjects(TCanvas *&aCanv, TH1F *&aDist) : canvas(aCanv), distribution(aDist) {}
};

class Key{
public:
  Key(int charge,int mass,int ID){
    this->charge = charge;
    this->mass = mass;
    this->ID = ID;
  }
  int charge;
  int mass;
  int ID;
  
  bool operator<(const Key &k) const{
    return (this->ID < k.ID);
  }
  bool operator==(const Key &rhs){
    return this->charge == rhs.charge && this->mass == rhs.mass && this->ID == rhs.ID;
  }
  
};

//Create a map of canvases and histoscorresponding to each mass and charge
//typedef map< Key, distObjects > distMap;
typedef map< Key, pair<TCanvas*, TH1F*> > distMap;

// distMap EtaDist;
// distMap EDist;

// distMap IhDist;
// distMap PDist;
// distMap PhiDist; //azimuthal angle
// distMap PtDist;
// distMap ThetaDist; //polar angle
// distMap TofDist;
// distMap TofTrelRatioDist;

//Graph (Scatter) Canvases (dep vs independent)
//TCanvas *canv_Ih_Vs_P[nTYPESS];
//one for given mass with varied charge colors
//one for given charge with varied masses
//for the following three distributions
// TCanvas *canv_Preco_Vs_Pgen[nCHARGE][nMASS];
// TCanvas *canv_Tofreco_Vs_Tofgen[nCHARGE][nMASS];
// TCanvas *canv_Ihreco_Vs_Ihgen[nCHARGE][nMASS];

/*Declare the corresponding distributions (TH1D)*/

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

/*Main*/
int main(int argc, char **argv){
  TApplication theApp("App",&argc,argv);
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
  map<double,int> massCounts = fileNameData->chargeCounts;
  
  // //Add the appropriate canvases to the map and initiate them
  //TCanvas*, TH1F*
  distMap BetaDist;
  distMap GammaDist;//1/sqrt(a-beta^2)
  
  TCanvas tempCanv;
  TH1F tempDist;
  
  for( int iFile = 0; iFile < numFiles; iFile++ ){
    double *charge = &charges[iFile];
    double *mass = &masses[iFile];
    string name = fmt::format("beta_{}_{}",3 * *charge,*mass);
    Key entryKey ( (int)(3* *charge), (int)*mass, iFile );
    
    BetaDist.emplace(entryKey, pair<TCanvas*,TH1F*>(new TCanvas(name.c_str(),name.c_str(),500,500), new TH1F(name.c_str(),name.c_str(),200,0,2)));

    name = fmt::format("#gamma_{}_{}",3 * *charge,*mass);
    GammaDist.emplace(entryKey, pair<TCanvas*,TH1F*>(new TCanvas(name.c_str(),name.c_str(),500,500), new TH1F(name.c_str(),name.c_str(),200,0,2)));
  }

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
    Key entryKey ( (int)(3* *charge), (int)*mass, iFile );
    //Create a string for the appropriate file directory
    ///storage/6/work/askeeters/HSCPStudy/HSCP_Root_Files
    string mcDir(fmt::format("mchamp{}_M_{}/HscpCandidates",3*charges[iFile],masses[iFile]));

    //Key accessKey ((int)(3 * *charge),(int)*mass, iFile);
    //cout << &BetaDist[accessKey].second << endl;
    
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

      //Calculate Momentum from transverse and theta (rad)
      P = Pt / TMath::Sin(theta);
      
      //Calculate the relativistic energy
      E = TMath::Sqrt( P*P + masses[iFile] );

      //Calculate beta
      beta = P / E;
      //Fill beta distribution
      auto it = BetaDist.find(entryKey);
      it->second.second->Fill(beta);

      //Calculate gamma
      gamma = 1.0 / TMath::Sqrt( 1 - beta*beta );
      //Fill the gamma distribution
      it = GammaDist.find(entryKey);
      it->second.second->Fill(gamma);
      
    }//end event loop
  }//end file loop

  //Loop through all available distributions, writing them to the disc.
  TCanvas *currCanv;
  TH1F *currDist;
  double norm = 1;
  // for(const auto &it : BetaDist){
  //   currCanv = it.second.first;
  //   currDist = it.second.second;
  //   currDist->Scale(norm/currDist->Integral("width"));
  //   currCanv->cd();
  //   currDist->Draw();
  //   currCanv->Update();
  //   currCanv->Draw();
  // }

  for(const auto &it : GammaDist){
    currCanv = it.second.first;
    currDist = it.second.second;
    currDist->Scale(norm/currDist->Integral("width"));
    currCanv->cd();
    currDist->Draw();
    currCanv->Update();
    currCanv->Draw();
  }
  
  theApp.Run();
  return 0;
}
