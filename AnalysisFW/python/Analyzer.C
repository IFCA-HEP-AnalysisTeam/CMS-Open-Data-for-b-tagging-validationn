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




void Analyzer(TString filename, TString ptRange, bool _ismc)
//void Analyzer()
{
//TString filename = "/eos/user/b/bchazinq/QCDPt30to50/tuples_MC30-50.root";
//TString ptRange = "_pthat30to50";
//TString filename = "/eos/user/b/bchazinq/Data/tuples_Data.root";
//TString filename = "OpenDataTree_mc.root";
//bool _ismc   = true;

/////////////////////////////////////////////////////////////////////////////////
//                                                                             //
//   This code opens OpenDataTree_mc.root input file. Extracts the branch      //
//   information (variables information). Stablishes a selection criteria.     //
//   Creates Histo_BranchName_ptRange.root output file. After the selection            //
//   requirements, it stores one histogram per cut.                             //
//                                                                             //
/////////////////////////////////////////////////////////////////////////////////

// Open  input file
// ----------------
TFile* infile  = new TFile( filename, "read" );
// Get the tree
TTree* mytree = (TTree*) infile -> Get( "ak5ak7/OpenDataTree");

// Name of branches to read
// ------------------------
const int nvar = 3;
TString bname [nvar]= { //"pthat"
                         "jet_pt", 
                         "jet_eta", 
                         "jet_phi",
                         
                      };
// Name of the output file
TString hname [nvar]= { "Histo_" + bname[0] + ptRange,  
                        "Histo_" + bname[1] + ptRange,
                        "Histo_" + bname[2] + ptRange,
                        //"Histo_" + bname[3] + ptRange,
                       };
// Set the #of bins, x-axix limits in the histo
TString hdimension[nvar] = { //"(130, 0, 650)"
                             "(100, 0, 250)",
                             "( 50,-5,   5)",
			     "( 50,-4,   4)",
                           };



// Flavour selection
// -----------------
const int ncut = 1;  
TCut Flavour [ncut] = {// "fabs(PartonF) == 5", 
                       // "fabs(PartonF) == 4", 
                       // "fabs(PartonF) != 5 && fabs(PartonF) != 4",
                       // "fabs(PartonF) == 5 && nBHadrons == 2",
                        "1>0"};
                     
// Name of histograms to save
TString hnameMC [ncut]= {// "b_quark",
                         // "c_quark",
		         // "lgluon",
                         // "b_gsplitting",
                          "allFlavours"
                        };
// cut selection
// -----------------
TCut triggerName = "triggernames == \"jt30\"";
//TCut triggerName = "triggernames == \"jt30\"";
TCut triggerPass = "triggers";
TCut minimumPt   = "jet_pt >60";
//TCut minimumPt   = "jet_pt >30";
TCut maximumEta  = "jet_eta <2.4";
TCut minimumEta  = "jet_eta >-2.4";

TCut tightID     = "jet_tightID";
//TCut flavour     = "fabs(PartonF) == 5 || fabs(PartonF) == 4 || fabs(PartonF) != 5 && fabs(PartonF) != 4 || fabs(PartonF) == 5 && nBHadrons == 2";
//TCut flavourb    = "fabs(PartonF) == 5";
//TCut flavourc    = "fabs(PartonF) == 4";
//TCut flavourlg   = "fabs(PartonF) != 5 && fabs(PartonF) != 4";
//TCut flavourb_gsplit = "fabs(PartonF) == 5 && nBHadrons == 2";
//TCut mainCut     = triggerName && triggerPass && minimumPt && (flavourb || flavourc || flavourlg || flavourb_gsplit) ; 
//TCut mainCut     = triggerName && triggerPass && minimumPt && (flavourb || flavourc || flavourlg || flavourb_gsplit) ; 
TCut mainCut = triggerName && triggerPass && minimumPt && minimumEta && maximumEta; 

//Values to normalize the mc
//---------------------------
//float eventw = (mcweight * lumi )/ ngen;
//float mcweight;
//luminosity value for 2011 runA JSON file 
//float  lumi = 2.33; // /fb
//number of generated events
TString ngen;
if (filename.Contains("15to30"))   ngen = "9978850"; //http://opendata.cern.ch/record/1366
if (filename.Contains("30to50"))   ngen = "5837856"; //http://opendata.cern.ch/record/1539
if (filename.Contains("50to80"))   ngen = "5766430"; //http://opendata.cern.ch/record/1555
if (filename.Contains("80to120"))  ngen = "5867864"; //http://opendata.cern.ch/record/1562
if (filename.Contains("120to170")) ngen = "5963264"; //http://opendata.cern.ch/record/1348 
if (filename.Contains("170to300")) ngen = "5975592"; //http://opendata.cern.ch/record/1369 
if (filename.Contains("300to470")) ngen = "5975016"; //http://opendata.cern.ch/record/1469 
if (filename.Contains("470to600")) ngen = "3967154"; //http://opendata.cern.ch/record/1553 
//For data: http://opendata.cern.ch/record/21

// Create the out put file
// ----------------------- 
// Make the loop over nvar
for (int i = 0; i < nvar; i++) { 
  
   TFile* outfile;
   if (_ismc)  
   {
     outfile  = new TFile( hname [i]+ "_MC.root", "recreate" );

     // Make the loop over ncut 
     for (int j = 0; j < ncut; j++){
       TH1D* myhisto = new TH1D();  
       //mytree -> Draw( bname[i] +">>" + hnameMC[j] + hdimension[i], "mcweight *2.33/" + ngen); // To check the pthat
       mytree -> Draw( bname[i] +">>" + hnameMC[j] + hdimension[i], (triggerName+triggerPass+minimumPt+maximumEta+minimumEta+tightID) * ("mcweight *2.33/" + ngen));
       //mytree -> Draw( bname[i] +">>" + hnameMC[j] + hdimension[i], (triggerCut0 && triggerCut1 && Flavour[j] && "jet_pt >30") * ("mcweight *2.33/" + ngen));
       //std::cout<< "the total cut is:   " <<"("<< triggerCut0 << " && " <<  triggerCut1 << " && " << Flavour[j] << " && " << "jet_pt >30) * (mcweight *2.33/" << ngen << ")" <<std::endl;
       myhisto = (TH1D*) gDirectory -> Get(hnameMC[j]);
       std::cout<< bname[i] << " Integral: " << myhisto -> Integral() << std::endl;   
       // Move the over flow
       int nbins = myhisto -> GetNbinsX ();     
       myhisto -> SetBinContent(nbins, myhisto -> GetBinContent(nbins+1) +  myhisto -> GetBinContent(nbins)); 
       // Save the histogram 
       myhisto -> Write();
       }
      
    }
   else
    {
     outfile  = new TFile( hname [i]+ "_Data.root", "recreate" );
     TH1D* myhisto = new TH1D();  
     mytree -> Draw(bname[i] + ">>" + "data" + hdimension[i], triggerName+triggerPass+minimumPt+maximumEta+minimumEta+tightID);
     ////mytree -> Draw(bname[i] + ">>" + "data" + hdimension[i], triggerCut0 && triggerCut1 && "jet_pt > 30");
     myhisto = (TH1D*) gDirectory -> Get("data");
     // Move the over flow
     int nbins = myhisto -> GetNbinsX ();
     myhisto -> SetBinContent(nbins, myhisto -> GetBinContent(nbins+1) +  myhisto -> GetBinContent(nbins));
     // Save the histogram
     myhisto -> Write();
    }

    outfile -> Close(); 
}
//c1-> Destructor();
infile -> Close();
}


void runAnalyzer(bool _ismc)
{
 // bool _ismc; 
  if (_ismc) 
   {
   /*TString filename0 = "/eos/user/b/bchazinq/QCDPt15to30/tuples_MC15-30.root";
    TString ptRange0 = "_pthat15to30";
    TString filename1 = "/eos/user/b/bchazinq/QCDPt30to50/tuples_MC30-50.root";
    TString ptRange1 = "_pthat30to50";
    TString filename2 = "/eos/user/b/bchazinq/QCDPt50to80/tuples_MC50-80.root";
    TString ptRange2 = "_pthat50to80";
    TString filename3 = "/eos/user/b/bchazinq/QCDPt80to120/tuples_MC80-120.root";
    TString ptRange3 = "_pthat80to120";
    TString filename4 = "/eos/user/b/bchazinq/QCDPt120to170/tuples_MC120-170.root";
    TString ptRange4 = "_pthat120to170";
    */TString filename5 = "/eos/user/b/bchazinq/QCDPt170to300/tuples_MC170-300.root";
    TString ptRange5 = "_pthat170to300";
    TString filename6 = "/eos/user/b/bchazinq/QCDPt300to470/tuples_MC300-470.root";
    TString ptRange6 = "_pthat300to470";
    TString filename7 = "/eos/user/b/bchazinq/QCDPt470to600/tuples_MC470-600.root";
    TString ptRange7 = "_pthat470to600";
   /* Analyzer (filename0,ptRange0,_ismc); 
    Analyzer (filename1,ptRange1,_ismc); 
    Analyzer (filename2,ptRange2,_ismc); 
    Analyzer (filename3,ptRange3,_ismc); 
    Analyzer (filename4,ptRange4,_ismc); 
    */Analyzer (filename5,ptRange5,_ismc); 
    Analyzer (filename6,ptRange6,_ismc);
    Analyzer (filename7,ptRange7,_ismc);
   }else{ 
   // _ismc   = false;
    TString filenameData = "/eos/user/b/bchazinq/Data/tuples_Data.root";
    TString ptRange = "";
    Analyzer (filenameData,ptRange,_ismc);
   }
}
 
