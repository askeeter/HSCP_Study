
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
const unsigned int nTYPES = 2;

//Declare distribution (histogram) Canvases
TCanvas *canv_BetaDist[nTYPES];
TCanvas *canv_EtaDist[nTYPES];
TCanvas *canv_GammaDist[nTYPES];
TCanvas *canv_IhDist[nTYPES];
TCanvas *canv_PDist[nTYPES];
TCanvas *canv_PhiDist[nTYPES];
TCanvas *canv_PtDist[nTYPES];
TCanvas *canv_TofDist[nTYPES];


//Graph (Scatter) Canvases (dep vs independent)

/*Declare the corresponding distributions (TH1D)*/


/*Function to parse a file name into mass and charge */

//Create a struct to allow returning both mass and charge
//as parsed from the file name
struct NameDat{
  double charge = 0.0;
  double mass = 0.0;
};

//Actual function that will parse the file name and return a NameDat struct
NameDat FileNameParser( const string &aName ){
  NameDat outDat; //Struct that will contain the parsed mass and charge

  //Loop through each character of the file name, increasing iCh.
  //Increase iChunk at each chunk ('_' character)
  //Store each chunk in a string array (chunks)
  //Clean up the chunks of interest, and return the data in the struct
  for (int iCh = 0, int iChunk = 0; iCh < aName.length(); iCh++){
    //Assign the current character to a variable
    char currChar = aName[iCh];
    string chunks[4]; //string array that will contain the chunks
    stringstream chunk; //stringstream object to contain pieced chunk
    
    if ( iCh == '_' ){
      chunks[j] = chunk.str();
      chunk = stringstream(string(""));
      iChunk++;
    }
    else
      chunk << currChar;
  }

  //Same process as the above, but to clean the mass chunk
  stringstream massStream;
  for (int iCh = 0; iCh < chunks[3].length(); iCh++){
    if ( chunks[3][iCh] == '.' )
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
  system("!source");
  system("cd /storage/6/work/askeeters/HSCPStudy/HSCP_Root_Files/");
  TFile *outFile = new TFile(fname.c_str(),"RECREATE");
  outFile.close();
  return outFile;
  return 0;
}
