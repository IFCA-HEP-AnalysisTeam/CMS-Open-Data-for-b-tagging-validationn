#include "TCanvas.h"
#include "TLegend.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2.h"
#include "TGraph.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TGraphAsymmErrors.h"
#include "TMath.h"
#include "TPaletteAxis.h"
#include "TROOT.h"
#include "TString.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "TAxis.h"
#include "TPaletteAxis.h"
#include "THStack.h"
#include "TAttFill.h"
#include "TColor.h"
#include "TString.h"
#include "TCut.h"
#include <fstream>
#include <iostream>




void Analyzer( )
{
TString filename = "OpenDataTree_mc.root";

/////////////////////////////////////////////////////////////////////////////////
//                                                                             //
//   This code opens OpenDataTree_mc.root input file. Extracts the branch      //
//   information (variables information). Stablishes a selection criteria.     //
//   Creates Histo_BranchName.root output file. After the selection            //
//   requirements, it stores one histogram per cut.                             //
//                                                                             //
/////////////////////////////////////////////////////////////////////////////////

// Name of branches to read
const int nvar = 5;
TString bname [nvar]= {  "jet_pt", 
                         "jet_CSV", 
                         "jet_JBP",
                         "jet_TCHP",
                         "nSVertex"
                      };
// Name of the output file
TString ofname [nvar]= { "Histo_" + bname[0],  
                         "Histo_" + bname[1],
                         "Histo_" + bname[2],
                         "Histo_" + bname[3],
                         "Histo_" + bname[4]
                       };
// Set the #of bins, x-axix limits in the histo
TString hdimension[nvar] = { "(250, 0, 250)",
                             "(50, 0,   1)",
			     "(50, 0,   8)",
                             "(50,  0,   30)",
                             "(  5, 0,   5)"
                           };

std::cout << "1" << std::endl; 

// Open  input file
TFile* ifile  = new TFile( filename, "read" );
// Get the tree
TTree* mytree = (TTree*) ifile -> Get( "ak5ak7/OpenDataTree");

// Flavour selection
const int ncut = 5;  
TCut Flavour [ncut] = { "PartonF == 5", 
                        "PartonF == 4", 
                        "PartonF != 5 && PartonF != 4",
                        "PartonF == 5 && nBHadrons == 2",
                        ""};
                     
// Trigger selection
TCut triggerCut0 = "triggernames == \"jt30\"";
TCut triggerCut1 = "triggers";

 
// MC Weigth
float  cross_section = 784265000; // fb
int    Ngen = 100000;
float  Lumi = 2.33; // /fb
float eventW = (cross_section * Lumi )/ Ngen;

// Name of histograms to save
TString hnameMC [ncut]= { "b_quark",
                          "c_quark",
		          "lgluon",
                          "b_gsplitting",
                          "allFlavours"
                        };


//TCanvas* c1 = new TCanvas ("c1", "c1");
std::cout << "2" << std::endl;

// Create the out put file 
// Make the loop over nvar
for (int i = 0; i < nvar; i++) { 
  
   TFile* ofile;
    
   if (filename.Contains("mc")){
     ofile  = new TFile( ofname [i]+ "_MC.root", "recreate" );

     // Make the loop over ncut 
     for (int j = 0; j < ncut; j++){
       TH1F* myhisto = new TH1F();  
       mytree -> Draw(bname[i] + ">>" + hnameMC[j] + hdimension[i], ( triggerCut0 && triggerCut1 && Flavour[j] && "jet_pt >30"));
       myhisto = (TH1F*) gDirectory -> Get(hnameMC[j]);
       // Move the over flow
       int nbins = myhisto -> GetNbinsX ();     
       myhisto -> SetBinContent(nbins, myhisto -> GetBinContent(nbins+1) +  myhisto -> GetBinContent(nbins));   
       // Save the histogram 
       myhisto -> Write();
      }
    }

   else  if ( filename.Contains("data")){
     ofile  = new TFile( ofname [i]+ "_Data.root", "recreate" );
     TH1F* myhisto = new TH1F();  
     mytree -> Draw(bname[i] + ">>" + "data" + hdimension[i], triggerCut0 && triggerCut1 && "jet_pt > 30");
     myhisto = (TH1F*) gDirectory -> Get("data");
     // Move the over flow
     int nbins = myhisto -> GetNbinsX ();
     myhisto -> SetBinContent(nbins, myhisto -> GetBinContent(nbins+1) +  myhisto -> GetBinContent(nbins));
     // Save the histogram
     myhisto -> Write();
    }

    ofile -> Close(); 
}
std::cout << "3" << std::endl; 
//c1-> Destructor();
ifile -> Close();
std::cout << "4" << std::endl; 
}


 
