#include "BJetAnalysis.C"
#include "TSystem.h"
#include "TChain.h"

void runBJetAnalysis (TString foldername, int mctrue, TString rangeOfpt)
{
   //if (fChain == 0) return;
   //TString _dataPath = filename;
   bool    _ismc = mctrue;
  
   // Read the tree from an input filename
  // TFile* infile  = new TFile( filename, "read" );
  // TTree* mytree  = (TTree*) infile -> Get( "ak5ak7/OpenDataTree");
  // Int_t nentries = (Int_t) mytree->GetEntries();
  // // call the analysis functions
  // BJetAnalysis BJetAnalysis(mytree);  // remeber change BJetAnalysis.C , BJetAnalysis.h and Jet
  // BJetAnalysis.Loop(filename, _ismc, rangeOfpt); 
   
  // Create a TChain to read several output_*.root trees in the foldername folder.
  TChain* mychain("ak5ak7/OpenDataTree");
  TString a = gSystem->GetFromPipe("ls " + foldername + "/output_*.root | wc -l"); 
  int b = atoi( a );
  for (int nfile=0; nfile < b; nfile++)
   {
    TString numb; 
    numb += nfile;
    mychain->Add(foldername + "output_" + numb + ".root"); 
   }
   Int_t nentries = (Int_t) mychain->GetEntries();

   BJetAnalysis BJetAnalysis(mychain); 
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
