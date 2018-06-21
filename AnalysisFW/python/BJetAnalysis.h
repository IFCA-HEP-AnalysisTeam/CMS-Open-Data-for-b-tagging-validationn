#ifndef BJetAnalysis_h
#define BJetAnalysis_h

#include "TH1F.h"


class BJetAnalysis : JetAnalysisBase
{
 public:

 //BJetAnalysis (TTree* tree = 0, TString systematic = "nominal");


 void  BeginJob         ();
 
 void  EndJob           ();

 void  Loop             (TString fileName);

 void  DefineHistograms (TH1F     rootHisto, 
                         TString  hname,
                         int      nbins, 
                         float    xmin, 
                         float    xmax);


 void  ResetHistograms  (TH1F rootHisto);

 
 float Overflow         (TH1F rootHisto);


 float Underflow        (TH1F rootHisto);


 void  FillHistograms   (TH1F rootHisto, varyyyyyyy)



 // Global variables
 // --------------------------------------------------------------------------------------------
  Long64_t nentries = 0;
  bool ismc;  
  TString filename; 
  float lumi
  unsigned int nvariables; 
  unsigned int ncuts;

  float lumi = 2.33 // 2011 legacy runA
  Int_t ngen; 
  float eventw;     
      
  float ptRange;

// Flavour selection
// -----------------  
  enum {
         allflavour,
         b_quark,
         c_quark,
         lgluon,
         b_gsplitting,
         nflavour,
        }

 const TString sflavour [nflavour+1] = {
         "allflavours",
         "b_quark",
         "c_quark",
         "lgluon",
         "b_gsplitting",
         "data"
        } 


 // TH1 histograms to plot
 // --------------------------------------------------------------------------------------------
 
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
   
