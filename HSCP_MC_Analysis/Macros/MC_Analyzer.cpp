
#include <cstddef> //size_t
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <string>
#include <cmath>
#include <map>
#include <array>
#include "format.h"
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TH1.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TMathBase.h"
#include "TMath.h"

using namespace std;

/*Declare canvases for each plot*/
const unsigned int RECO = 0;
const unsigned int GEN = 1;
const unsigned int BKG = 2;
const unsigned int nTYPES = 3;

//Create a map of canvases corresponding to each mass and charge
//A map with masses as key which has a value of a map of charges
//with charges as key and canvas as value
// [gen/reco][charge][mass]
map< string,map< double,map< double ,TCanvas > > > canv_BetaDist;//Pc/E
map< double,map<double,TCanvas> > canv_EtaDist;
map< double,map<double,TCanvas> > canv_GammaDist;//1/sqrt(a-beta^2)
map< double,map<double,TCanvas> > canv_IhDist;
map< double,map<double,TCanvas> > canv_PDist;
map< double,map<double,TCanvas> > canv_PhiDist; //azimuthal angle
map< double,map<double,TCanvas> > canv_PtDist;
map< double,map<double,TCanvas> > canv_ThetaDist; //polar angle
map< double,map<double,TCanvas> > canv_TofDist;
map< double,map<double,TCanvas> > canv_TofTrelRatioDist;

//Graph (Scatter) Canvases (dep vs independent)
TCanvas *canv_Ih_Vs_P[nTYPES];
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
  if(!(inStream = popen("ls /storage/6/work/askeeters/HSCPStudy/HSCP_Root_Files/Histos_mchamp*.root -l | wc -l", "r"))){
    exit(0);
  }
    
  while(fgets(charnFiles, sizeof(charnFiles), inStream)!=NULL){
    streamnFiles << charnFiles;
  }
  pclose(inStream);

  int nFiles = atoi(streamnFiles.str().c_str());
  string *listOfFiles = new string[nFiles+1];

  //Now we have the number of files, and can dynamicall allocate a string array for it.
  //PLUS the first entry containing the number of entries in the arary.
  
  listOfFiles[0] = streamnFiles.str();
  
  //Now we need to read in the list of files
  //Will do something similar to the above, but using this method to get only the file names
  stringstream testStream;
  if(!(inStream = popen("find /storage/6/work/askeeters/HSCPStudy/HSCP_Root_Files/*.root -printf \"%f\n\"","r"))){
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
  outDat->fileNames = &Names[1];
  
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
int main(void){
  string *fileList = ListOfFiles(); 

  //Total number of files to be processed
  int numFiles = atoi(fileList[0].c_str());

  NameDat *fileNameData = FileNameParser( fileList, numFiles );

  //Arrays of the charges and masses associated with each file
  double *masses = fileNameData->mass;
  double *charges = fileNameData->charge;

  //List of the names of the files to be opened
  fileList = fileNameData->fileNames; //Can index from 0 now
 
  //Number of each charge and each mass
  map<double,int> *chargeCounts = &fileNameData->chargeCounts;
  map<double,int> *massCounts = &fileNameData->chargeCounts;


  //Pointer for the current file
  TFile *datFile;

  //Pointer for the tree being read
  TTree *tree;
 
  
  int event;
  int Hscp;
  float Pt, P;
  float I;
  float Ih; //This is the energy deposition that we are interested in 
  float TOF; 
  float Mass; //Reco mass. Not correct. Assumes unit charge
  float dZ, dXY, dR;
  float eta, phi, theta;
  bool hasMuon;

  //Loop over the files
  for( int iFile = 0; iFile < numFiles; iFile++ ){
    string file = "/storage/6/work/askeeters/HSCPStudy/HSCP_Root_Files/";
    file += fileList[iFile];
    datFile = new TFile(file.c_str()); //Open the current file

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
    cout << file << endl;
    cout << fmt::format("{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}", "event", "Pt(GeV)", "Ih(MeV/cm)", "Mass(GeV)", "eta", "phi", "theta(deg)") << endl;
    for(int iEvt = 0; iEvt < tree->GetEntriesFast(); iEvt++){
      tree->GetEntry(iEvt);
      
      //Calculate theta
      theta = 2*TMath::ATan( TMath::Exp( -1 * eta ) );

      //Calculate Momentum from transverse and theta (rad)
      P = Pt / TMath::Sin(theta);
      
      cout << fmt::format("{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}", event, Pt, Ih, Mass, eta, phi, theta*180.0/TMath::Pi())<< endl;
    }

    datFile->Close();
    
  }//end file loop
  return 0;
}
