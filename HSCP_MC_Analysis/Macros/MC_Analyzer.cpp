
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <string>
#include <cmath>
#include "TROOT.h"
#include "TFile.h"
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
  cout << "There are " << nFiles << " files to be processed" << endl << flush;
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
  double charge;
  double mass;
};

//Actual function that will parse the file name and return a NameDat struct
NameDat FileNameParser( const string &aName ){
  NameDat outDat; //Struct that will contain the parsed mass and charge
  string chunks[4]; //string array that will contain the chunks
  //Loop through each character of the file name, increasing iCh.
  //Increase iChunk at each chunk ('_' character)
  //Store each chunk in a string array (chunks)
  //Clean up the chunks of interest, and return the data in the struct
    
  for (int iCh = 0, iChunk = 0; iCh < aName.length(); iCh++){
    //Assign the current character to a variable
    char currChar = aName[iCh];
    stringstream chunk; //stringstream object to contain pieced chunk
      
    if ( currChar == '_' ){
      chunks[iChunk] = chunk.str();
      chunk.str(string(""));
      iChunk++;
    }
    else
      chunk << currChar;
  }

  //Same process as the above, but to clean the mass chunk
  stringstream massStream;
  for (int iCh = 0; iCh < chunks[3].length(); iCh++){
    if ( chunks[3].at(iCh) == '.' )
      break;
    else
      massStream << iCh;
  }

  //Convert the strings to floats/doubles
  outDat.charge = (double)atof(string(chunks[1],6).c_str());
  outDat.mass = (double)atof(massStream.str().c_str());
  return outDat;
}

/*Function to return a string array of all of the available MC files*/


/*Main*/
int main(void){
  TFile *outFile = new TFile("MC_Analysis.root","RECREATE");
  outFile->Close();
  
  string *fileList = ListOfFiles();
  int numFiles = atoi(fileList[0].c_str());
  for(int i=1; i<numFiles; i++){
    cout << fileList[i] << endl;
  }
  return 0;
}
