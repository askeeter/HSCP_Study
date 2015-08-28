
#include <cstddef> //size_t
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <string>
#include <cmath>
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TH1.h"
#include "TGraph.h"
#include "TCanvas.h"

using namespace std;

/*Declare canvases for each plot*/
const unsigned int RECO = 0;
const unsigned int GEN = 1;
const unsigned int BKG = 2;
const unsigned int nTYPES = 3;
//use static constants when calculating the number of masses and 
//number of charges. Declare the corresponding constants here

//Declare distribution (histogram) Canvases
TCanvas *canv_BetaDist[nTYPES];
TCanvas *canv_EtaDist[nTYPES];
TCanvas *canv_GammaDist[nTYPES];
TCanvas *canv_IhDist[nTYPES];
TCanvas *canv_PDist[nTYPES];
TCanvas *canv_PhiDist[nTYPES]; //azimuthal angle
TCanvas *canv_PtDist[nTYPES];
TCanvas *canv_ThetaDist[nTYPES]; //polar angle
TCanvas *canv_TofDist[nTYPES];
TCanvas *canv_TofTrelRatioDist[nTYPES];

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
  
  return outDat;
}


/*Main*/
int main(void){
  //TFile *outFile = new TFile("MC_Analysis.root","RECREATE");
  //outFile->Close();
  string *fileList = ListOfFiles();
  int numFiles = atoi(fileList[0].c_str());
  NameDat *fileNameData = FileNameParser( fileList, numFiles );
  
  
  return 0;
}
