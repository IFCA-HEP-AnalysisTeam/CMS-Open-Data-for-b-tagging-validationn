#ifndef BJetAnalysis_h
#define BJetAnalysis_h
//#include "JetAnalysisBase.h"
//#include "DataOutput2Class.h"
//#include "TreeClass_QCD15to30.h"
#include "ChainClass.h"
#include "Constants.h"
#include <TTree.h>
#include <TH1.h>
#include <TMath.h>
#include <TChain.h>
#include <fstream>
#include <iostream>


class BJetAnalysis : public ChainClass
//class BJetAnalysis : public TreeClass_QCD15to30
//class BJetAnalysis : public JetAnalysisBase
//class BJetAnalysis : public DataOutput2Class 
{
 public:

 BJetAnalysis (TChain* chain);
 //BJetAnalysis (TTree* tree);

 void  BeginJob         (TString filename,
                         bool _ismc);


 void  Loop             (TString _dataPath, 
                         bool _ismc, 
                         TString _ptRange);


 void  DefineHistograms (int iflv, 
                         int sflv);


 void  SaveHistograms   (TH1F* myHistogram, 
                         TString hRootName, 
                         TString ptRange);


 void   PrintProgress   (Long64_t counter,
                         Long64_t total);


 void  ResetHistograms  (TH1F* rootHisto);

 
 float xOverflow        (TH1F* rootHisto);


 float xUnderflow       (TH1F* rootHisto);


 void  EndJob           ();



 // Global variables
 // --------------------------------------------------------------------------------------------
  Int_t nentries;
  TString filename; 
  TString ptRange;
  bool ismc;  
  float eventw;     
  float ngen; 
      
 
 // TH1 histograms to plot
 // --------------------------------------------------------------------------------------------

  // minimun dR between ak5CaloJets and ak5PFJets. The current cut is dRmin < 0.3
  TH1F* dRmin[nflavour];
 
 // jet variables
  TH1F* jetPt [nflavour];
  TH1F* jetEta[nflavour];
  TH1F* jetPhi[nflavour];
 
 // selected tracks variables
  TH1F* IP3D                [nflavour]; 
  TH1F* IP3Dsignif          [nflavour]; 
  TH1F* avgTrackMultiplicity[nflavour];

 // ordinary tracks
  TH1F* nrPixelHits      [nflavour]; 
 // TH1F* nrTrackerHits    [nflavour];
  TH1F* trackPt          [nflavour];
  TH1F* distanceToJetAxis[nflavour]; 

 // secondary vertex
  TH1F* flight3Dsignif[nflavour]; 
  TH1F* massSV        [nflavour];
  TH1F* nrSV          [nflavour];

 // b-discriminants
  TH1F* TCHE[nflavour]; 
  TH1F* TCHP[nflavour];
  TH1F* JP  [nflavour];  
  TH1F* JBP [nflavour];
  TH1F* CSV [nflavour];   
};
# endif   
