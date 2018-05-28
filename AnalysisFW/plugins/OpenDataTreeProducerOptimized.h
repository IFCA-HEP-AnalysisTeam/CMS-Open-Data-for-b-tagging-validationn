#ifndef OpenDataTreeProducerOptimized_h
#define OpenDataTreeProducerOptimized_h

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "JetMETCorrections/Objects/interface/JetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/JetCorrFactors.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
// Discriminator test
#include "DataFormats/BTauReco/interface/JetTag.h"
#include "DataFormats/BTauReco/interface/JetTagInfo.h"
//////////////

using namespace edm;
using namespace reco;
using namespace std;
using namespace trigger;

class OpenDataTreeProducerOptimized : public edm::EDAnalyzer 
{
  public:

    explicit OpenDataTreeProducerOptimized(edm::ParameterSet const& cfg);
    virtual void beginJob();
    virtual void beginRun(edm::Run const &, edm::EventSetup const& iSetup);
    virtual void analyze(edm::Event const& evt, edm::EventSetup const& iSetup);
    virtual void endRun(edm::Run const &, edm::EventSetup const& iSetup);
    virtual void endJob();
    virtual ~OpenDataTreeProducerOptimized();


  private:  

    // Function to help sort the jet wrt. pT
    static bool cmp_patjets(const pat::Jet &pj1, const pat::Jet &pj2) {
        return pj1.pt() > pj2.pt();
    }

    //---- Configurable parameters --------  
    bool            mIsMCarlo;
    bool            mUseGenInfo;
    bool            mPrintTriggerMenu;
    int             mMinNPFJets;
    double          mMinPFPt, mMinGenPt, mMaxY, mMinJJMass;
    int             mGoodVtxNdof;
    double          mGoodVtxZ; 
    edm::InputTag   mPFak5JetsName;
    edm::InputTag   mPFak7JetsName;

    // Test flavour 
    edm::InputTag jetFlavourInfos_;
    // Test SV
    edm::InputTag secondaryVertexTagInfos_;   
    // Test Tracks IP
    edm::InputTag impactParameterTagInfos_; 
    edm::InputTag m_ipassoc;

    // ---- PF Jet input tags ----- //
    edm::InputTag   mGenJetsName;
    edm::InputTag   mSrcPFRho;
    edm::InputTag   mPFMET; 
    edm::InputTag   mOfflineVertices;
    
    //---- Trigger----------------------
    std::string                 processName_;
    std::vector<std::string>    triggerNames_;
    std::vector<unsigned int>   triggerIndex_;
    edm::InputTag               triggerResultsTag_;
    HLTConfigProvider           hltConfig_;
   
   
 
    // Output variables
    edm::Service<TFileService>  fs;
    TTree                       *mTree;

    
    //---- TTree variables --------
    
    static const UInt_t kMaxNjet = 64;
    static const UInt_t kMaxNtrg = 32;



    // PF jets
    UInt_t njet;
    Float_t jet_pt[kMaxNjet];
    Float_t jet_eta[kMaxNjet];
    Float_t jet_phi[kMaxNjet];
    Float_t jet_E[kMaxNjet];
    Bool_t jet_tightID[kMaxNjet];
    Float_t jet_area[kMaxNjet];
    Float_t jet_jes[kMaxNjet];
    Int_t jet_igen[kMaxNjet];

    // tracks
    UInt_t  trackIndex;                      // number of tracks in the jet
    Int_t   ntracks         [kMaxNjet];  // track multiplicity in the jet
    Float_t track_pt        [kMaxNjet];  // pt of tracks in the jet
    Int_t track_nValidPixelHits   [kMaxNjet]; 
    Int_t track_nValidTrackerHits [kMaxNjet]; 
    Float_t track_IPz         [kMaxNjet];  // longitudinal Impact parameter
    Float_t track_IP2D        [kMaxNjet];
    Float_t track_IP2Dsig     [kMaxNjet];
    Float_t track_IP3D        [kMaxNjet];
    Float_t track_IP3Dsig     [kMaxNjet];
    Float_t track_distToJetaxis [kMaxNjet];
    // ..... chi2 normalized
    // ..... decay length

    // loose WP for commisionnig
    ////////////////////////////
    Bool_t jet_looseID[kMaxNjet];
    ////////////////////////////

    // b discriminants
    /////////////////////////// 
    Float_t jet_CSV[kMaxNjet];
    Float_t jet_JBP[kMaxNjet];
    Float_t jet_TCHP[kMaxNjet];
    //////////////////////////

    // Test Flavour
    /////////////////////////// 
    Float_t ptF       [kMaxNjet];   
    Float_t etaF      [kMaxNjet];   
    Float_t phiF      [kMaxNjet];  
    Float_t HadronF   [kMaxNjet];   
    Float_t PartonF   [kMaxNjet];  
    Float_t nBHadrons [kMaxNjet];  
    ///////////////////////////

    // Secondary Vertex
    ///////////////////////////
    UInt_t  nSVertex;
    ///////////////////////////

    // PF jets
    UInt_t njet_ak7;
    Float_t jet_pt_ak7[kMaxNjet];
    Float_t jet_eta_ak7[kMaxNjet];
    Float_t jet_phi_ak7[kMaxNjet];
    Float_t jet_E_ak7[kMaxNjet];
    Float_t jet_area_ak7[kMaxNjet];
    Float_t jet_jes_ak7[kMaxNjet];
    Int_t ak7_to_ak5[kMaxNjet];

    // Jet composition
    Float_t chf[kMaxNjet];
   	Float_t nhf[kMaxNjet];
   	Float_t phf[kMaxNjet];
   	Float_t elf[kMaxNjet];
   	Float_t muf[kMaxNjet];
   	Float_t hf_hf[kMaxNjet];
   	Float_t hf_phf[kMaxNjet];
   	Int_t hf_hm[kMaxNjet];
   	Int_t hf_phm[kMaxNjet];
   	Int_t chm[kMaxNjet];
   	Int_t nhm[kMaxNjet];
   	Int_t phm[kMaxNjet];
   	Int_t elm[kMaxNjet];
   	Int_t mum[kMaxNjet];   
    Float_t beta[kMaxNjet];   
    Float_t bstar[kMaxNjet];
   
        // for loose WP commisioning plots
        //////////////////////////////////
        Float_t nhfJet[kMaxNjet];  
        Float_t nemfJet[kMaxNjet]; 
        Float_t chemfJet[kMaxNjet];
        Float_t chmJet[kMaxNjet];  
        //////////////////////////////////


    // Generated jets
    UInt_t ngen;
    Float_t gen_pt[kMaxNjet];
    Float_t gen_eta[kMaxNjet];
    Float_t gen_phi[kMaxNjet];
    Float_t gen_E[kMaxNjet];

    // Event identification
    UInt_t run;
    UInt_t lumi;
    ULong64_t event;

    // Test to get the N generated in MC, N processed in data
    //     mTree->Branch("nevent", nevent,"nevent/i");
    UInt_t nevent = 0; 

    // Triggers
    UInt_t ntrg;
    Bool_t triggers[kMaxNtrg];
    std::vector<std::string> triggernames;
    UInt_t prescales[kMaxNtrg];

    // MET, SuMET, rho, eventfilter
    Float_t met;
    Float_t sumet;
    Float_t rho;

    // MC variables
    Float_t pthat;
    Float_t mcweight;

    // Jet correction labels
    std::string mJetCorr_ak5;
    std::string mJetCorr_ak7;
};

#endif
