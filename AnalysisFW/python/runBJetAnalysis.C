#include "BJetAnalysis.h"
#include "BJetAnalysis.C"
#include "TSystem.h"
#include "TCanvas.h"

void runBJetAnalysis (TString foldername, int mctrue, TString rangeOfpt)
{
   //if (fChain == 0) return;
   //TString _dataPath = filename;
   bool    _ismc = mctrue;
  
   // Read the tree from an input filename = /fullpath/file.root
   //TFile* infile  = new TFile( foldername, "read" );
   //TTree* myInput  = (TTree*) infile -> Get( "ak5ak7/OpenDataTree");
   
  // Create a TChain to read several output_*.root trees in the foldername folder.
  TChain* myInput = new TChain("ak5ak7/OpenDataTree", "");
  TString a = gSystem->GetFromPipe("ls " + foldername + "output_*.root | wc -l"); 
  int b = atoi( a );
  for (int nfile=0; nfile < b; nfile++)
   {
    //if (nfile == 1) continue; // just QCDPt300to470
    TString numb; 
    numb += nfile;
    myInput->Add(foldername + "output_" + numb + ".root"); 
   }

  // call the analysis functions
   BJetAnalysis BJetAnalysis(myInput); 
   BJetAnalysis.Loop(foldername, _ismc, rangeOfpt); 
 }

# ifndef __CINT__
int main(int argc, char ** argv)
{
 printf("\n  ismc = 0 -> data; ismc = 1 -> mc  \n");
 if (argc < 4)
  { 
   printf("\n ./runBJetAnalysis <filename> <bool ismc> <samplePt>\n");
   return -1;
  }
 else if (argc > 4 ) 
  {
   printf("\n too many arguments: \n");
   printf("\n      ./runBJetAnalysis <filename> <bool ismc> <samplePt>\n");
   return -1;
  }
  else if (argc == 4)  runBJetAnalysis(argv[1], atoi(argv[2]), argv[3]);
  return 0;
}
#endif
