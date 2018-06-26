#include "BJetAnalysis.C"

void runBJetAnalysis (TString filename, bool mctrue, TString rangeOfpt)
{
   //if (fChain == 0) return;
   //TString _dataPath = filename;
   //bool    _ismc = mctrue;
   //cout << " _ismc " << _ismc << endl;
   //TString _ptRange = rangeOfpt;
  // open the input file
   TFile* infile  = new TFile( filename, "read" );
  // get the tree
   TTree* mytree  = (TTree*) infile -> Get( "ak5ak7/OpenDataTree");
   Int_t nentries = (Int_t) mytree->GetEntries();

   BJetAnalysis BJetAnalysis(mytree); 
   BJetAnalysis.Loop(filename, mctrue, _ptRange); 
   //BJetAnalysis.Loop(_dataPath, _ismc, _ptRange); 
 }

# ifndef __CINT__
int main(int argc, char ** argv)
{
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
  else if (argc == 4)  runBJetAnalysis(argv[1], argv[2], argv[3]);
  return 0;
}
#endif
