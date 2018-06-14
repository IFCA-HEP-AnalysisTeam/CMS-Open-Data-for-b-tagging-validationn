

// Forked from SMPJ Analysis Framework
// https://twiki.cern.ch/twiki/bin/viewauth/CMS/SMPJAnalysisFW
// https://github.com/cms-smpj/SMPJ/tree/v1.0

#include <typeinfo>
#include <iostream>
#include <sstream>
#include <istream>
#include <fstream>
#include <iomanip>
#include <string>
#include <cmath>
#include <functional>
#include "TTree.h"
#include <vector>
#include <cassert>
#include <TLorentzVector.h>

#include "OpenDataTreeProducerOptimized.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/Common/interface/TriggerResultsByName.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/GeometryCommonDetAlgo/interface/Measurement1D.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/Common/interface/RefToBase.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/JetReco/interface/JetCollection.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/JetReco/interface/JetExtendedAssociation.h"
#include "DataFormats/JetReco/interface/JetFloatAssociation.h"
#include "DataFormats/JetReco/interface/JetID.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/PFMETCollection.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/BTauReco/interface/JetTag.h"
#include "DataFormats/BTauReco/interface/JetTagInfo.h"
#include "DataFormats/BTauReco/interface/SecondaryVertexTagInfo.h"
#include "DataFormats/BTauReco/interface/TaggingVariable.h"  

#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenRunInfoProduct.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "SimDataFormats/JetMatching/interface/JetFlavour.h"
#include "SimDataFormats/JetMatching/interface/JetFlavourMatching.h"
#include "SimDataFormats/JetMatching/interface/MatchedPartons.h"
#include "SimDataFormats/JetMatching/interface/JetMatchedPartons.h"
#include "SimDataFormats/JetMatching/interface/JetFlavourInfo.h"
#include "SimDataFormats/JetMatching/interface/JetFlavourInfoMatching.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/JetReco/interface/JetTracksAssociation.h"
#include "DataFormats/TrackReco/interface/HitPattern.h"
#include "DataFormats/BTauReco/interface/TrackIPTagInfo.h"
#include "DataFormats/Math/interface/Vector3D.h"

#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"
#include "RecoJets/JetAssociationProducers/src/JetTracksAssociatorAtVertex.h"

// track IP info
#include "DataFormats/BTauReco/interface/TrackIPTagInfo.h"
#include "DataFormats/BTauReco/interface/TrackProbabilityTagInfo.h"
#include "RecoBTag/ImpactParameter/interface/TrackCountingComputer.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "DataFormats/GeometryVector/interface/GlobalVector.h"
#include "RecoVertex/VertexPrimitives/interface/ConvertToFromReco.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
//https://github.com/cms-sw/cmssw/blob/CMSSW_5_3_X/PhysicsTools/JetExamples/test/printJetFlavourInfo.py
//https://github.com/cms-sw/cmssw/blob/CMSSW_5_3_X/PhysicsTools/JetExamples/test/printJetFlavourInfo.cc

OpenDataTreeProducerOptimized::OpenDataTreeProducerOptimized(edm::ParameterSet const &cfg) {
  mMinPFPt           = cfg.getParameter<double>                    ("minPFPt");
  mMinJJMass         = cfg.getParameter<double>                    ("minJJMass");
  mMaxY              = cfg.getParameter<double>                    ("maxY");
  mMinNPFJets        = cfg.getParameter<int>                       ("minNPFJets");
  mPFak5JetsName     = cfg.getParameter<edm::InputTag>             ("pfak5jets");
  mPFak7JetsName     = cfg.getParameter<edm::InputTag>             ("pfak7jets");
  mCaloak5JetsName   = cfg.getParameter<edm::InputTag>             ("caloak5jets");
  mOfflineVertices   = cfg.getParameter<edm::InputTag>             ("offlineVertices");
  mGoodVtxNdof       = cfg.getParameter<double>                    ("goodVtxNdof");
  mGoodVtxZ          = cfg.getParameter<double>                    ("goodVtxZ");
  mSrcPFRho          = cfg.getParameter<edm::InputTag>             ("srcPFRho");
  mPFMET             = cfg.getParameter<edm::InputTag>             ("pfmet");
  mGenJetsName       = cfg.getUntrackedParameter<edm::InputTag>    ("genjets",edm::InputTag(""));
  mPrintTriggerMenu  = cfg.getUntrackedParameter<bool>             ("printTriggerMenu",false);
  mIsMCarlo          = cfg.getUntrackedParameter<bool>             ("isMCarlo",false);
  mUseGenInfo        = cfg.getUntrackedParameter<bool>             ("useGenInfo",false);
  mMinGenPt          = cfg.getUntrackedParameter<double>           ("minGenPt",30);
  processName_       = cfg.getParameter<std::string>               ("processName");
  triggerNames_      = cfg.getParameter<std::vector<std::string> > ("triggerNames");
  triggerResultsTag_ = cfg.getParameter<edm::InputTag>             ("triggerResults");
  mJetCorr_ak5       = cfg.getParameter<std::string>               ("jetCorr_ak5");
  mJetCorr_ak7       = cfg.getParameter<std::string>               ("jetCorr_ak7");
  // tracks IP
  m_ipassoc = cfg.getParameter<edm::InputTag>("ipassociation");

  // test flavour
  if (mIsMCarlo){mJetFlavourInfos = cfg.getParameter<edm::InputTag>("jetFlavourInfos");}
  // test SV
//  impactParameterTagInfoTrack_  = cfg.getParameter<edm::InputTag>("impactParameterTagInfoTrack");
  secondaryVertexTagInfos_  = cfg.getParameter<edm::InputTag>("secondaryVertexTagInfos"); 
}

void OpenDataTreeProducerOptimized::beginJob() {
    mTree = fs->make< TTree >("OpenDataTree", "OpenDataTree");

    // Variables of the flat tuple
    mTree->Branch("njet", &njet, "njet/i");
    mTree->Branch("jet_pt", jet_pt, "jet_pt[njet]/F");
    mTree->Branch("jet_eta", jet_eta, "jet_eta[njet]/F");
    mTree->Branch("jet_phi", jet_phi, "jet_phi[njet]/F");
    mTree->Branch("jet_E", jet_E, "jet_E[njet]/F");   
    mTree->Branch("jet_tightID", jet_tightID, "jet_tightID[njet]/O");
    mTree->Branch("jet_area", jet_area, "jet_area[njet]/F");
    mTree->Branch("jet_jes", jet_jes, "jet_jes[njet]/F");
    mTree->Branch("jet_igen", jet_igen, "jet_igen[njet]/I");
    // b discriminant
    mTree->Branch("jet_CSV", jet_CSV, "jet_CSV[njet]/F");
    mTree->Branch("jet_JP", jet_JP, "jet_JP[njet]/F");
    mTree->Branch("jet_JBP", jet_JBP, "jet_JBP[njet]/F");
    mTree->Branch("jet_TCHP", jet_TCHP, "jet_TCHP[njet]/F");
    mTree->Branch("jet_TCHE", jet_TCHE, "jet_TCHE[njet]/F");
    mTree->Branch("dRmin_matching", dRmin_matching, "dRmin_matching[njet]/F");

    // AK7 variables
    mTree->Branch("njet_ak7", &njet_ak7, "njet_ak7/i");
    mTree->Branch("jet_pt_ak7", jet_pt_ak7, "jet_pt_ak7[njet_ak7]/F");
    mTree->Branch("jet_eta_ak7", jet_eta_ak7, "jet_eta_ak7[njet_ak7]/F");
    mTree->Branch("jet_phi_ak7", jet_phi_ak7, "jet_phi_ak7[njet_ak7]/F");
    mTree->Branch("jet_E_ak7", jet_E_ak7, "jet_E_ak7[njet_ak7]/F");
    mTree->Branch("jet_area_ak7", jet_area_ak7, "jet_area_ak7[njet_ak7]/F");
    mTree->Branch("jet_jes_ak7", jet_jes_ak7, "jet_jes_ak7[njet_ak7]/F");
    mTree->Branch("ak7_to_ak5", ak7_to_ak5, "ak7_to_ak5[njet_ak7]/I");

    mTree->Branch("ngen", &ngen, "ngen/i");
    mTree->Branch("gen_pt", gen_pt, "gen_pt[ngen]/F");
    mTree->Branch("gen_eta", gen_eta, "gen_eta[ngen]/F");
    mTree->Branch("gen_phi", gen_phi, "gen_phi[ngen]/F");
    mTree->Branch("gen_E", gen_E, "gen_E[ngen]/F");

    mTree->Branch("run", &run, "run/i");
    mTree->Branch("lumi", &lumi, "lumi/i");
    mTree->Branch("event", &event, "event/l");
    mTree->Branch("ntrg", &ntrg, "ntrg/i");
    mTree->Branch("triggers", triggers, "triggers[ntrg]/O");
    mTree->Branch("triggernames", &triggernames);
    mTree->Branch("prescales", prescales, "prescales[ntrg]/i");
    mTree->Branch("met", &met, "met/F");
    mTree->Branch("sumet", &sumet, "sumet/F");
    mTree->Branch("rho", &rho, "rho/F");
    mTree->Branch("pthat", &pthat, "pthat/F");
    mTree->Branch("mcweight", &mcweight, "mcweight/F");

    mTree->Branch("chf", chf, "chf[njet]/F");   
    mTree->Branch("nhf", nhf, "nhf[njet]/F");   
    mTree->Branch("phf", phf, "phf[njet]/F");   
    mTree->Branch("elf", elf, "elf[njet]/F");   
    mTree->Branch("muf", muf, "muf[njet]/F");   
    mTree->Branch("hf_hf", hf_hf, "hf_hf[njet]/F");   
    mTree->Branch("hf_phf", hf_phf, "hf_phf[njet]/F");   
    mTree->Branch("hf_hm", hf_hm, "hf_hm[njet]/i");    
    mTree->Branch("hf_phm", hf_phm, "hf_phm[njet]/i");
    mTree->Branch("chm", chm, "chm[njet]/i");   
    mTree->Branch("nhm", nhm, "nhm[njet]/i");   
    mTree->Branch("phm", phm, "phm[njet]/i");   
    mTree->Branch("elm", elm, "elm[njet]/i");   
    mTree->Branch("mum", mum, "mum[njet]/i");
    mTree->Branch("beta", beta, "beta[njet]/F");   
    mTree->Branch("bstar", bstar, "bstar[njet]/F");
    
    //loose WP for commisionning
    mTree->Branch("nhfJet", nhfJet, "nhfJet[njet]/F");
    mTree->Branch("nemfJet", nemfJet, "nemfJet[njet]/F");
    mTree->Branch("chemfJet", chemfJet, "chemfJet[njet]/F");
    mTree->Branch("chmJet", chmJet, "chmJet[njet]/i"); 
    mTree->Branch("jet_looseID", jet_looseID, "jet_looseID[njet]/O");
 
    // Test flavour  
   // mTree->Branch("ptF",     ptF,    "ptF[njet]/F");    
   // mTree->Branch("etaF",    etaF,   "etaF[njet]/F");    
   // mTree->Branch("phiF",    phiF,   "phiF[njet]/F");    
    mTree->Branch("HadronF", HadronF,"HadronF[njet]/F");    
    mTree->Branch("PartonF", PartonF,"PartonF[njet]/F");    
    mTree->Branch("nBHadrons", nBHadrons,"nBHadrons[njet]/F");   
    // Test to get the N generated in MC, N processed in data
    mTree->Branch("nevent", &nevent,"nevent/i");   
    // Secondary Vertex
    mTree->Branch("nSVinEvent",  &nSVinEvent, "nSVinEvent/I");    
    mTree->Branch("jetSVIndex",  &jetSVIndex, "jetSVIndex[nSVinEvent]/i");    
    mTree->Branch("svmass",      svmass,    "svmass[nSVinEvent]/F");
    mTree->Branch("flight3DSignificance", flight3DSignificance, "flight3DSignificance[nSVinEvent]/F");    
    // B-tag(IPTagInfo) selected tracks (associated to ak5CaloJets)
    mTree->Branch("seltracksInEvent", &seltracksInEvent, "seltracksInEvent/i");
    mTree->Branch("jetSeltrackIndex", jetSeltrackIndex, "jetSeltrackIndex[seltracksInEvent]/I"); 
    mTree->Branch("seltrack_pt", seltrack_pt, "seltrack_pt[seltracksInEvent]/F");
    mTree->Branch("seltrack_nValidPixelHits", seltrack_nValidPixelHits, "seltrack_nValidPixelHits[seltracksInEvent]/I");
    mTree->Branch("seltrack_nValidTrackerHits", seltrack_nValidTrackerHits, "seltrack_nValidTrackerHits[seltracksInEvent]/I");
    mTree->Branch("seltrack_IP2D", seltrack_IP2D, "seltrack_IP2D [seltracksInEvent]/F");
    mTree->Branch("seltrack_IP2Dsig", seltrack_IP2Dsig, "seltrack_IP2Dsig [seltracksInEvent]/F");
    mTree->Branch("seltrack_IP3D", seltrack_IP3D, "seltrack_IP3D [seltracksInEvent]/F");
    mTree->Branch("seltrack_IP3Dsig", seltrack_IP3Dsig, "seltrack_IP3Dsig [seltracksInEvent]/F");
    mTree->Branch("seltrack_distToJetAxis", seltrack_distToJetAxis, "seltrack_distToJetAxis [seltracksInEvent]/F");
//  mTree->Branch("track_IPz", "track_IPz", "track_IPz [seltracksInEvent]/F");
  
    // tracks (associated to ak5PFJets)
    mTree->Branch("tracks_inEvent", &tracks_inEvent, "tracks_inEvent/i");
    mTree->Branch("tracks_jetIndex", tracks_jetIndex, "tracks_jetIndex[tracks_inEvent]/I");
    mTree->Branch("tracks_nValidPixelHits", tracks_nValidPixelHits, "tracks_nValidPixelHits[tracks_inEvent]/I");
    mTree->Branch("tracks_nValidTrackerHits", tracks_nValidTrackerHits, "tracks_nValidTrackerHits[tracks_inEvent]/I");
    mTree->Branch("tracks_pt", tracks_pt, "tracks_pt[tracks_inEvent]/F");
    mTree->Branch("tracks_chi2", tracks_chi2, "tracks_chi2[tracks_inEvent]/F");
    mTree->Branch("tracks_IPz", tracks_IPz, "tracks_IPz[tracks_inEvent]/F");
    mTree->Branch("tracks_distToJetAxis", tracks_distToJetAxis, "tracks_distToJetAxis[tracks_inEvent]/D");
    mTree->Branch("tracks_decayLength", tracks_decayLength, "tracks_decayLength[tracks_inEvent]/D");
}

void OpenDataTreeProducerOptimized::endJob() {
}


void OpenDataTreeProducerOptimized::beginRun(edm::Run const &iRun,
                                     edm::EventSetup const &iSetup) {

    // Mapping trigger indices 
    int iTrigger=0;//mi trigger***
    int jTrigger=0; // trigger de la coleccion***
    bool changed(true);
    if (hltConfig_.init(iRun, iSetup, processName_, changed) && changed) {

        // List of trigger names and indices 
        // are not emptied between events, must be done heckIndex re
        triggerIndex_.clear();
        triggernames.clear();

        // Iterate over all active triggers of the AOD file
        iTrigger++;
        auto name_list = hltConfig_.triggerNames();
        for (std::string name_to_search: triggerNames_) {
            // Find the version of jet trigger that is active in this run 
            
            for (std::string name_candidate: name_list) {
	   
                jTrigger++;
	        //printf(" Mi trigger %i: %s y el trigger de la coleccion %i: %s\n", iTrigger++, triggerNames_, jTrigger++, name_list);
		//printf(" Mi trigger %i: %s y el trigger de la coleccion %i: %s\n", iTrigger++, name_to_search, jTrigger++, name_candidate);
//                std::cout<< "Mi trigger  : " << iTrigger << "    name_to_search     " << name_to_search << std::endl; 
//                std::cout<< "el trigger de la coleccion  : " << jTrigger << "    name_candidate     " << name_candidate << std::endl; 

                // Match the prefix to the full name (eg. HLT_Jet30 to HLT_Jet30_v10)
                if ( name_candidate.find(name_to_search + "_v") != std::string::npos ) {
                    // Save index corresponding to the trigger
                    triggerIndex_.push_back(hltConfig_.triggerIndex(name_candidate));

                    // Save the trigger name
                    triggernames.push_back("jt" + name_to_search.substr(7, string::npos));
                    break;            
                }
            }
        }
    }

    // Retrieve cross section of the simulated process
    mcweight = 0;
    if (mIsMCarlo) {

        edm::Handle<GenRunInfoProduct> genRunInfo;
        iRun.getByLabel("generator", genRunInfo );

        // Save only the cross section, since the total number of 
        // generated events is not available in this context (!!)
        mcweight = genRunInfo->crossSection();
        std::cout << "Cross section: " <<  mcweight << std::endl;
    }
    
}


void OpenDataTreeProducerOptimized::analyze(edm::Event const &event_obj,
                                    edm::EventSetup const &iSetup) {

    // Event info
    run = event_obj.id().run();
    lumi = event_obj.luminosityBlock();
    event = event_obj.id().event();
    cout << " ##########################################################################################  " << endl; 
    std::cout << " run number: " << run <<std::endl; 
    std::cout << " lumi = lumiosityBlock: " << lumi <<std::endl; 
    cout << " event number: " << event << endl; 
    /////////////////////////////////////////////////////////////////////////////////////////////////////////  
    // Test Discriminant 
    //---------------------------- Jet CSV discriminantor -----------------------
    edm::Handle<reco::JetTagCollection> tagHandle_CSV;
    event_obj.getByLabel("combinedSecondaryVertexBJetTags", tagHandle_CSV);
    //const reco::JetTagCollection & tag_CSV = *(tagHandle_CSV.product());
    //---------------------------- Jet JBP tag discriminantor -------------------
    edm::Handle<reco::JetTagCollection> tagHandle_JBP;
    event_obj.getByLabel("jetBProbabilityBJetTags", tagHandle_JBP); 
    //const reco::JetTagCollection & tag_JBP = *(tagHandle_JBP.product());
    //---------------------------- Jet JP tag discriminantor -------------------
    edm::Handle<reco::JetTagCollection> tagHandle_JP;
    event_obj.getByLabel("jetProbabilityBJetTags", tagHandle_JP); 
    //const reco::JetTagCollection & tag_JP = *(tagHandle_JP.product());
    //---------------------------- Jet TCHP discriminator -----------------------
    edm::Handle<reco::JetTagCollection> tagHandle_TCHP;
    event_obj.getByLabel("trackCountingHighPurBJetTags", tagHandle_TCHP); 
    //const reco::JetTagCollection & tag_TCHP = *(tagHandle_TCHP.product());
    //---------------------------- Jet TCHE discriminator -----------------------
    edm::Handle<reco::JetTagCollection> tagHandle_TCHE;
    event_obj.getByLabel("trackCountingHighEffBJetTags", tagHandle_TCHE); 
    //const reco::JetTagCollection & tag_TCHE = *(tagHandle_TCHE.product());
    
    // Print out the info
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////
//    std::cout << "event " << event << std::endl;  
//    std::cout << "-----------------------------------------------------------------------" << std::endl;    
//    std::cout << "-----------------------------------------------------------------------" << std::endl;    
//
//    std::cout <<      "---------------------------- Jet CSV tag Info -------------------" << std::endl;  
//    std::cout << "tag_CSV.size()    " << tag_CSV.size() << std::endl; 
//    for (int i = 0; i != (int)tag_CSV.size(); i++)
//     {
//     std::cout << "ptCSV   " << tag_CSV[i].first -> pt() << "    etaCSV   " << tag_CSV[i].first -> eta() << "   phiCSV   " << tag_CSV[i].first -> phi() << "    disc1   "<< tag_CSV[i].second << std::endl;   
//     } 
//    std::cout <<      "---------------------------- Jet JBP tag Info -------------------" << std::endl; 
//      std::cout << "tag_JBP.size()    " << tag_JBP.size() << std::endl; 
//    for (int i = 0; i != (int)tag_JBP.size(); i++)
//    { 
//      std::cout << "ptJBP   " << tag_JBP[i].first -> pt() << "  etaJBP    " << tag_JBP[i].first -> eta() << "    phiJBP   " << tag_JBP[i].first -> phi() << "    disc1   "<< tag_JBP[i].second << std::endl;   
//    }
//  //////////////////////////////////////////////////////////////////////////////////////////////////////////// 
    /////////////////////////////////////////////////////////////////////////////////////////////////////////  
   
    // Test Flavour
    // Print out the info
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //if (mIsMCarlo){
    edm::Handle<reco::JetFlavourInfoMatchingCollection> theJetFlavourInfos;
    //event_obj.getByLabel(mJetFlavourInfos, theJetFlavourInfos );
    // for ( reco::JetFlavourInfoMatchingCollection::const_iterator j  = theJetFlavourInfos->begin(); j != theJetFlavourInfos->end(); ++j ) 
    // {
    //      std::cout << "-------------------- Jet Flavour Info --------------------" << std::endl;
    //      const reco::Jet *aJet = (*j).first.get(); 
    //      reco::JetFlavourInfo aInfo = (*j).second;
    //      // ----------------------- Hadrons -------------------------------
    //      std::cout << "                      Hadron-based flavour: " << aInfo.getHadronFlavour() << std::endl;
    //      // ----------------------- Partons -------------------------------
    //      std::cout << "                      Parton-based flavour: " << aInfo.getPartonFlavour() << std::endl;
    // }
    //}
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Triggers
    edm::Handle<edm::TriggerResults>   triggerResultsHandle_;
    event_obj.getByLabel(triggerResultsTag_, triggerResultsHandle_);

    // Sanity checks
    assert(triggerResultsHandle_.isValid() && "Error in getting TriggerResults from Event!");
    assert(triggerResultsHandle_->size() == hltConfig_.size() && "Size mismatch between triggerResultsHandle_ and hltConfig_");
    
    // Number of triggers to be saved
    ntrg = triggerIndex_.size();

    // Iterate only over the selected jet triggers
    for (unsigned itrig = 0; itrig < ntrg; itrig++) {

        // Trigger bit
        Bool_t isAccepted = triggerResultsHandle_->accept(triggerIndex_[itrig]);
        triggers[itrig] = isAccepted;

        // Trigger prescales are retrieved using the trigger name
        std::string trgName = hltConfig_.triggerName(triggerIndex_[itrig]);
        const std::pair< int, int > prescalePair(hltConfig_.prescaleValues(event_obj, iSetup, trgName));

        // Total prescale: PreL1*PreHLT 
        prescales[itrig] = prescalePair.first*prescalePair.second;   
    }    

    // Rho
    Handle< double > rho_handle;
    event_obj.getByLabel(mSrcPFRho, rho_handle);
    rho = *rho_handle;


    // Generator Info

    // Retrieve pthat and mcweight (only MC)
    pthat = 0;
    if (mIsMCarlo && mUseGenInfo) {
        Handle< GenEventInfoProduct > hEventInfo;
        event_obj.getByLabel("generator", hEventInfo);

        // Monte Carlo weight (NOT AVAILABLE FOR 2011 MC!!)
        //mcweight = hEventInfo->weight();
        
        // Pthat 
        if (hEventInfo->hasBinningValues()) {
            pthat = hEventInfo->binningValues()[0];
        }
    }

    // Generator-level jets
    ngen = 0;
    if (mIsMCarlo) {

        Handle< GenJetCollection > genjets;
        event_obj.getByLabel(mGenJetsName, genjets);
    
        // Index of the simulated jet
        int gen_index = 0; 

        for (GenJetCollection::const_iterator i_gen = genjets->begin(); i_gen != genjets->end(); i_gen++)  {

            // pT and rapidity selection
            if (i_gen->pt() > mMinGenPt && fabs(i_gen->y()) < mMaxY) {
                gen_pt[gen_index] = i_gen->pt();
                gen_eta[gen_index] = i_gen->eta();
                gen_phi[gen_index] = i_gen->phi();
                gen_E[gen_index] = i_gen->energy();
                gen_index++;
            }
        }

        // Number of generated jets in this event
        ngen = gen_index;
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////////////  
    // Secondary Vertex Info
    edm::Handle<reco::SecondaryVertexTagInfoCollection> svTagInfosHandle;
    event_obj.getByLabel(secondaryVertexTagInfos_, svTagInfosHandle);
    //event_obj.getByLabel("secondaryVertexTagInfos", svTagInfosHandle);
    const reco::SecondaryVertexTagInfoCollection & svTagInfoColl = *(svTagInfosHandle.product());
    // counter of number of secondary vertex in the event
    int nSVinEvent_index = 0; 
    /////////////////////////////////////////////////////////////////////////////////////////////////////////  
    // print the info
    /////////////////////////////////////////////////////////////////////////////////////////////////////////  
//    for(reco::SecondaryVertexTagInfoCollection::const_iterator iter = svTagInfoColl.begin(); iter != svTagInfoColl.end(); ++iter) 
//    {
//      // if there are reconstructed vertices in this jet
//      nSVertex = iter->nVertices();
//      cout << " a jet with   pt = " << iter->jet()->pt() << " eta = " << iter->jet()->eta() << " phi = " << iter->jet()->phi() << '\n';       
//      cout << " contains " << iter->nVertices() << " secondary vertices " << endl;
//      if(iter->nVertices() > 0 )
//       {      
//        for (unsigned int vtx = 0; vtx < iter->nVertices(); ++vtx)
//         {
//          std::cout<<"with" << endl;  
//          std::cout<<"        3D flight distance of " << iter->flightDistance(vtx,false).value() << " cm"<<std::endl;
//          std::cout<<"        3D flight significance of " << iter->flightDistance(vtx,false).significance() << " cm"<<std::endl;
//          std::cout<<"        2D flight distance of " << iter->flightDistance(vtx,true).value() << " cm"<<std::endl;
//          std::cout<<"        2D flight significance of " << iter->flightDistance(vtx,true).significance() << " cm"<<std::endl;
//          std::cout<<"        invariant mass of     " << iter->secondaryVertex(vtx).p4().mass() << " GeV"<<endl;
//          const Vertex &vertex = iter->secondaryVertex(vtx); 
//          std::cout<<"        invariant mass 2nd way of     " << vertex.p4().mass() << " GeV"<<endl;
//         }
//       }
//    }
//    /////////////////////////////////////////////////////////////////////////////////////////////////////////  
   
   

    // Transient track for IP calculation
    edm::ESHandle<TransientTrackBuilder> builder;
    iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", builder);
    
    
    // Vertex Info
    Handle<reco::VertexCollection> recVtxs;
    event_obj.getByLabel(mOfflineVertices, recVtxs);
    // use first pv of the collection
    //   Vertex dummy;
    //   const Vertex *pv = &dummy
    const reco::Vertex *pv;

    // PF AK5 Jets

    edm::Handle<reco::PFJetCollection> ak5_handle;
    event_obj.getByLabel(mPFak5JetsName, ak5_handle);
    const JetCorrector* corrector_ak5 = JetCorrector::getJetCorrector(mJetCorr_ak5, iSetup);


    // Jet Track Association (JTA)
    edm::Handle <reco::TrackCollection> tracks_h;
    event_obj.getByLabel ("generalTracks", tracks_h);
    std::auto_ptr<reco::JetTracksAssociation::Container> tracksInJets (new reco::JetTracksAssociation::Container (reco::JetRefBaseProd(ak5_handle)));
    // format inputs
    std::vector <edm::RefToBase<reco::Jet> > allJets;
    allJets.reserve (ak5_handle->size());
    for (unsigned i = 0; i < ak5_handle->size(); ++i)
    {
      edm::RefToBase<reco::Jet> jetRef(edm::Ref<reco::PFJetCollection>(ak5_handle, i));
      allJets.push_back(jetRef);
    }
    std::vector <reco::TrackRef> allTracks;
    allTracks.reserve(tracks_h->size());
    for (unsigned i = 0; i < tracks_h->size(); ++i) 
      allTracks.push_back (reco::TrackRef(tracks_h, i));
    // run JTA algorithm
    JetTracksAssociationDRVertex mAssociator(0.5); // passed argument: 0.5 cone size
    mAssociator.produce (&*tracksInJets, allJets, allTracks);
  
    // Index of the selected jet 
    int ak5_index = 0;

    // Jet energy correction factor
    double jec = -1.0;

    // Jets will be unsorted in pT after applying JEC,  
    // therefore store corrected jets in a new collection (map): key (double) is pT * -1 (key), 
    // value (std::pair<PFJet*, double>) is pair of original jet iterator and corresponding JEC factor
    
    int njetNoCORR =0;

    std:: map<int, int> ak5_handle_index;
    std::map<double, std::pair<reco::PFJetCollection::const_iterator, double> > sortedJets;
    int ak5jet_orig_index =0;


    //######################################################
    // My AK5CaloJets reco::Jets
    //######################################################
         edm::Handle<edm::View<reco::Jet> > myJets;
         event_obj.getByLabel(mCaloak5JetsName, myJets);
      /////////////////////////////////////////////////////////////////////////////////////////////////////////  
      // print the info
      /////////////////////////////////////////////////////////////////////////////////////////////////////////
//    std::cout <<      "---------------------------- mCaloak5JetsName -------------------" << std::endl; 
//    for (unsigned index = 0; index < myJets -> size(); ++index)
//     {
//       edm::RefToBase<reco::Jet> jetRef = myJets->refAt(index);
//    
// 
//       float disc = (*tagHandle_CSV)[myJets->refAt(index)];
//       //float disc = (*tagHandle_JBP)[jetRef];
//       cout << " disc_csv = " << disc<< endl; 
//       cout << " pt = " << (*myJets)[index].pt() << "   eta = "  << (*myJets)[index].eta() <<  "   phi  " << (*myJets)[index].phi() << "  disc2 = " << disc<< endl; 
//     }
//    
//    std::cout <<      "---------------------------- mCaloak5JetsName 2nd way ---------------------" << std::endl; 
//    for( edm::View<reco::Jet>::const_iterator jet = myJets->begin(); jet != myJets->end(); ++jet )
//     {
//      cout << "  mCaloak5JetsName     pt = " << jet->pt() << "   eta = "  << jet->eta() <<  "   phi  " << jet->phi() << endl;
//     }
//    std::cout <<      "---------------------------------------------------------------------------" << std::endl; 
      /////////////////////////////////////////////////////////////////////////////////////////////////////////
    //######################################################


    //################################################################
    // B-tagging selected tracks from TrackIPTagInfoCollection 
    //################################################################
     Handle<TrackIPTagInfoCollection> ipHandle;
     event_obj.getByLabel(m_ipassoc, ipHandle);
     const TrackIPTagInfoCollection & ip = *(ipHandle.product());
     //index for the  track in the event
     int tracks_inEvent_index = 0; 
     //index for the B-tag(IPTagInfo) selected track in the event
     int seltracksInEvent_index =0;
     //################################################################


    for (auto i_ak5jet_orig = ak5_handle->begin(); i_ak5jet_orig != ak5_handle->end(); ++i_ak5jet_orig) {
        
        /////////////////////////////////////////////////////////////////////////////////////////////
        // Para chequear la correccion 
        /////////////////////////////////////////////////////////////////////////////////////////////
      /*  printf( "\nmi jet sin corregir %i\n", njetNoCORR);
        printf("  pt = %f\n",                           i_ak5jet_orig->pt());
        printf("  eta = %f\n",                           i_ak5jet_orig->eta());
        printf("  phi = %f\n",                           i_ak5jet_orig->phi());
        printf("  chargedHadronEnergyFraction  = %f\n", i_ak5jet_orig->chargedHadronEnergyFraction());
        printf("  muonEnergyFraction  = %f\n",          i_ak5jet_orig->muonEnergyFraction());
        printf("  neutralHadronMultiplicity = %i\n",    i_ak5jet_orig->neutralHadronMultiplicity());
        printf("  chargedHadronMultiplicity = %i\n",    i_ak5jet_orig->chargedHadronMultiplicity());
        printf("  chargedMultiplicity       = %i\n",    i_ak5jet_orig->chargedMultiplicity());
        printf("  neutralMultiplicity       = %i\n",    i_ak5jet_orig->chargedHadronMultiplicity());
       */ 
        ////////////////////////////////////////////////////////////////////////////////////////////

        // take jet energy correction and get corrected pT
        jec = corrector_ak5->correction(*i_ak5jet_orig, event_obj, iSetup);
        // Multiply pT by -1 in order to have largest pT jet first (sorted in ascending order by default)
        sortedJets.insert(std::pair<double, std::pair<reco::PFJetCollection::const_iterator, double> >(-1 * i_ak5jet_orig->pt() * jec, std::pair<reco::PFJetCollection::const_iterator, double>(i_ak5jet_orig, jec)));
        // Fill the index map 
        ak5_handle_index.insert(std::make_pair (int(-1000 * i_ak5jet_orig->pt() * jec), ak5jet_orig_index));  
        
        njetNoCORR ++;
        ak5jet_orig_index ++;

     }
   
    //cout << "the number of NON corrected PF Jets in the event is: " << njetNoCORR << endl; 
    int njetCORR =0;


    // Iterate over the jets (sorted in pT) of the event
    for (auto i_ak5jet_orig = sortedJets.begin(); i_ak5jet_orig != sortedJets.end(); ++i_ak5jet_orig) {

        // Apply jet energy correction "on the fly":
        // copy original (uncorrected) jet;
        PFJet corjet = *((i_ak5jet_orig->second).first);
        // take stored JEC factor
        jec = (i_ak5jet_orig->second).second;
        // apply JEC
        corjet.scaleEnergy(jec);
        // pointer for further use
        const PFJet* i_ak5jet = &corjet;

        // jet index in the original collection
        int indexKey = -1000 * corjet.pt();
        int jetRefIndex = ak5_handle_index.at(indexKey);

        /////////////////////////////////////////////////////////////////////////////////////////////
        // Para chequear la correccion 
        /////////////////////////////////////////////////////////////////////////////////////////////
        /*printf("\nmi jet  corregido %i\n", njetCORR);
        printf("  pt = %f\n",i_ak5jet->pt());
        printf("  eta = %f\n",i_ak5jet->eta());
        printf("  phi = %f\n",i_ak5jet->phi());
        printf("  chargedHadronEnergyFraction  = %f\n", i_ak5jet->chargedHadronEnergyFraction());
        printf("  muonEnergyFraction  = %f\n",          i_ak5jet->muonEnergyFraction());
        printf("  neutralHadronMultiplicity = %i\n",    i_ak5jet->neutralHadronMultiplicity());
        printf("  chargedHadronMultiplicity = %i\n",    i_ak5jet->chargedHadronMultiplicity());
        printf("  chargedMultiplicity       = %i\n",    i_ak5jet->chargedMultiplicity());
        printf("  neutralMultiplicity       = %i\n",    i_ak5jet->chargedHadronMultiplicity());
        */njetCORR ++;
        /////////////////////////////////////////////////////////////////////////////////////////////
 
        // Skip the current iteration if jet is not selected
        if (fabs(i_ak5jet->y()) > mMaxY || 
            (i_ak5jet->pt()) < mMinPFPt) {
            continue;
        }

        // Computing beta and beta*

        // Get tracks
        reco::TrackRefVector tracks = reco::JetTracksAssociation::getValue(*tracksInJets, *((i_ak5jet_orig->second).first));

        float sumTrkPt(0.0), sumTrkPtBeta(0.0),sumTrkPtBetaStar(0.0);
        beta[ak5_index] = 0.0;
        bstar[ak5_index] = 0.0;


        // Loop over tracks of the jet
        //cout << " -------------------------- No selected tracks --------------------------------" << endl;  
        for(auto i_trk = tracks.begin(); i_trk != tracks.end(); i_trk++) {

            if (recVtxs->size() == 0) break;
            tracks_jetIndex [tracks_inEvent_index] = ak5_index;
            //cout << "The track number " << tracks_inEvent_index << " in jet number  " << ak5_index << "  has : " << endl; 

            //cout << "  pt = " << (*i_trk)->pt() << endl;              
            //cout << "  normalized chi2 = " << (*i_trk)->normalizedChi2()<< endl;              
            //cout << "  numberOfValidPixelHits = " << (*i_trk)->hitPattern().numberOfValidPixelHits() << endl;
            //cout << "  numberOfValidTrackerHits = "   << (*i_trk)->hitPattern().numberOfValidTrackerHits() << endl;    
            tracks_pt[tracks_inEvent_index] = (*i_trk)->pt();
            tracks_nValidTrackerHits[tracks_inEvent_index] = (*i_trk)->hitPattern().numberOfValidTrackerHits();
            tracks_nValidPixelHits[tracks_inEvent_index] = (*i_trk)->hitPattern().numberOfValidPixelHits();
            tracks_chi2[tracks_inEvent_index] = (*i_trk)->normalizedChi2();          

            //Extract the IP Info
            // we always use the first vertex (at the moment)
            pv = &*recVtxs->begin();
            //cout << "  transverse impact parameter 2d (xy) = " << (*i_trk)->dxy(pv->position()) << endl; 
            //cout << "  longitudinal impact parameter (z)   = " << (*i_trk)->dz(pv->position()) << endl; 
            tracks_IPz [tracks_inEvent_index] = (*i_trk)->dz(pv->position());
            

            //TrackIPTagInfo::TrackIPData trackIP; 
            const reco::TransientTrack transientTrack = builder->build(*i_trk); 
            const GlobalVector direction(i_ak5jet->px(), i_ak5jet->py(), i_ak5jet->pz()); //!!!!!!!
            Double_t distanceToJetAxis =  IPTools::jetTrackDistance(transientTrack, direction, *pv).second.value();
            //cout << "  distance to the jet = " << distanceToJetAxis << endl;
      
            tracks_distToJetAxis[tracks_inEvent_index] = distanceToJetAxis;
           
            Double_t IP2d = -999;
            bool ipPass = IPTools::signedTransverseImpactParameter(transientTrack, direction, *pv).first;
            if (ipPass) 
             {
              IP2d = IPTools::signedTransverseImpactParameter(transientTrack, direction, *pv).second.value();
              //cout << "  if (ip2d) signedimpactparameter = " << IP2d <<endl;       
              tracks_IP2D [tracks_inEvent_index] = IP2d;
             }
            Double_t decayLength= -999;
            TrajectoryStateOnSurface closest = IPTools::closestApproachToJet(transientTrack.impactPointState(), *pv, direction,transientTrack.field());
            if (closest.isValid()) 
             {
              decayLength = (closest.globalPosition() - RecoVertex::convertPos(pv->position())).mag();
              //cout << "  decay lenght = " << decayLength << endl;
              tracks_decayLength[tracks_inEvent_index] = decayLength; 
             }
 
           // Sum pT
            sumTrkPt += (*i_trk)->pt();
            
            // Loop over vertices
            for (unsigned ivtx = 0; ivtx < recVtxs->size(); ivtx++) {
                reco::Vertex vertex = (*recVtxs)[ivtx];

                // Loop over tracks associated with the vertex
                bool flagBreak = false;
                if (!(vertex.isFake()) && 
                    vertex.ndof() >= mGoodVtxNdof && 
                    fabs(vertex.z()) <= mGoodVtxZ) {
                    
                    for(auto i_vtxTrk = vertex.tracks_begin(); i_vtxTrk != vertex.tracks_end(); ++i_vtxTrk) {
                        
                        // Match the jet track to the track from the vertex
                        reco::TrackRef trkRef(i_vtxTrk->castTo<reco::TrackRef>());
                        
                        // Check for matching vertices
                        if (trkRef == (*i_trk)) {
                            if (ivtx == 0) {
                                sumTrkPtBeta += (*i_trk)->pt();
                            }
                            else {
                                sumTrkPtBetaStar += (*i_trk)->pt();
                            }
                            flagBreak = true;
                            break;
                        } 
                    } 
                    if(flagBreak)
                      break;
                } 
            }
        tracks_inEvent_index ++; 
        }
        
    
        if (sumTrkPt > 0) 
          {
            beta[ak5_index]   = sumTrkPtBeta/sumTrkPt;
            bstar[ak5_index]  = sumTrkPtBetaStar/sumTrkPt;
          } 


        // Jet composition
        // (all energy fractions have to be multiplied by the JEC factor)
        chf[ak5_index]     = i_ak5jet->chargedHadronEnergyFraction() * jec;
        nhf[ak5_index]     = (i_ak5jet->neutralHadronEnergyFraction() + i_ak5jet->HFHadronEnergyFraction()) * jec;
        phf[ak5_index]     = i_ak5jet->photonEnergyFraction() * jec;
        elf[ak5_index]     = i_ak5jet->electronEnergyFraction() * jec;
        muf[ak5_index]     = i_ak5jet->muonEnergyFraction() * jec;
        hf_hf[ak5_index]   = i_ak5jet->HFHadronEnergyFraction() * jec;
        hf_phf[ak5_index]  = i_ak5jet->HFEMEnergyFraction() * jec;
        hf_hm[ak5_index]   = i_ak5jet->HFHadronMultiplicity();
        hf_phm[ak5_index]  = i_ak5jet->HFEMMultiplicity();
        chm[ak5_index]     = i_ak5jet->chargedHadronMultiplicity();
        nhm[ak5_index]     = i_ak5jet->neutralHadronMultiplicity();
        phm[ak5_index]     = i_ak5jet->photonMultiplicity();
        elm[ak5_index]     = i_ak5jet->electronMultiplicity();
        mum[ak5_index]     = i_ak5jet->muonMultiplicity();
    
        int npr      = i_ak5jet->chargedMultiplicity() + i_ak5jet->neutralMultiplicity();

        bool isHighEta = fabs(i_ak5jet->eta()) > 2.4;
        bool isLowEta = fabs(i_ak5jet->eta()) <= 2.4 && 
                        nhf[ak5_index] < 0.9 &&
                        phf[ak5_index] < 0.9 && 
                        elf[ak5_index] < 0.99 && 
                        chf[ak5_index] > 0 && 
                        chm[ak5_index] > 0;
        bool tightID =  npr > 1 && 
                        phf[ak5_index] < 0.99 && 
                        nhf[ak5_index] < 0.99 &&
                        (isLowEta || isHighEta);

        // Jet ID requirement for btagging commisioning plots 2011 collected from twiki: 
        // https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetID#Recommendations_for_7_TeV_data_a
        nhfJet[ak5_index]   = i_ak5jet->neutralHadronEnergyFraction() * jec; 
        nemfJet[ak5_index]  = i_ak5jet->neutralEmEnergyFraction() * jec;
        chemfJet[ak5_index] = i_ak5jet->chargedEmEnergyFraction() * jec;
        chmJet[ak5_index]   = i_ak5jet->chargedMultiplicity (); 
        // Loose WP
        bool looseID = i_ak5jet->pt() > 10.0 && fabs(i_ak5jet->eta()) < 2.4 && chf[ak5_index] > 0.0 &&  nhfJet[ak5_index] < 0.99 && chmJet[ak5_index] > 0.0 && nemfJet[ak5_index] < 0.99 &&  chemfJet[ak5_index] < 0.99 && npr > 1; 
        jet_looseID[ak5_index] = looseID; 
        
 
        // Variables of the tuple
        jet_tightID[ak5_index] = tightID;
        jet_area[ak5_index] = i_ak5jet->jetArea();
        jet_jes[ak5_index] = jec; // JEC factor

        // p4 is already corrected!
        auto p4 = i_ak5jet->p4();
        jet_pt[ak5_index]   = p4.Pt();
        jet_eta[ak5_index]  = p4.Eta();
        jet_phi[ak5_index]  = p4.Phi();
        jet_E[ak5_index]    = p4.E(); 
        
        // Matching a GenJet to this PFjet
        jet_igen[ak5_index] = 0;
        if (mIsMCarlo && ngen > 0) {

            // Index of the generated jet matching this PFjet
            jet_igen[ak5_index] = -1; // is -1 if no matching jet

            // Search generated jet with minimum distance to this PFjet   
            float r2min(999);
            for (unsigned int gen_index = 0; gen_index != ngen; gen_index++) {
                double deltaR2 = reco::deltaR2( jet_eta[ak5_index], 
                                                jet_phi[ak5_index],
                                                gen_eta[gen_index], 
                                                gen_phi[gen_index]);
                if (deltaR2 < r2min) {
                    r2min = deltaR2;
                    jet_igen[ak5_index] = gen_index;
                }
            }
        }
       
    HadronF[ak5_index] = -999; PartonF[ak5_index] = -999; nBHadrons [ak5_index] = -999;
    //Test Flavour
    if (mIsMCarlo){
      
      //edm::Handle<reco::JetFlavourInfoMatchingCollection> theJetFlavourInfos;
      event_obj.getByLabel(mJetFlavourInfos, theJetFlavourInfos );
      
      cout<<"     " <<endl;       
      cout<<" ------------------------ selected jet number "<< ak5_index<<"  ----------------------------" <<endl;
      cout<<"                                 pt   = " << jet_pt[ak5_index] << endl;      
      cout<<"                                 eta  = " << jet_eta[ak5_index] << endl;      
      cout<<"                                 phi  = " << jet_phi[ak5_index] << endl;      

      //---------------------------- Associated Jet Info -------------------------------
      const reco::Jet *aJet =(*theJetFlavourInfos)[jetRefIndex].first.get();
      cout<<" ------------------------ JetFlavourInfoMatchingCollection jet  "<< ak5_index<<"  ----------------------------" <<endl;
      cout<<"                                 pt   = " << aJet->pt()<< endl;      
      cout<<"                                 eta  = " << aJet->eta() << endl;      
      cout<<"                                 phi  = " << aJet->phi() << endl;
      cout<<"     " <<endl;       
      double deltaR2 = reco::deltaR2( jet_eta[ak5_index], jet_phi[ak5_index], aJet->eta(), aJet->phi()); 
      if (sqrt(deltaR2) < 0.1)
       {
        //---------------------------- Jet Flavour Info --------------------------------
        reco::JetFlavourInfo aInfo = (*theJetFlavourInfos)[jetRefIndex].second;
        // ----------------------- Hadrons-based flavour -------------------------------
        HadronF[ak5_index] = aInfo.getHadronFlavour();
        // ----------------------- Parton-based flavour  -------------------------------
        PartonF[ak5_index] = aInfo.getPartonFlavour();
        //------------------------ # of clustered b-hadrons ----------------------------
        const reco::GenParticleRefVector & bHadrons = aInfo.getbHadrons();
        nBHadrons [ak5_index] = bHadrons.size();
        //check if the method getbHadrons() gives b hadrons that do not have other b hadrons as daughters  
        for(reco::GenParticleRefVector::const_iterator it = bHadrons.begin(); it != bHadrons.end(); ++it)
         {
           //cout << "b hadron " << endl; 
           //cout << "  pdgId = "   << (*it)->pdgId() << endl; 
           //cout << "  number of daughters = " << (*it)->numberOfDaughters() << endl; 
           unsigned int nDaughters = (*it)->numberOfDaughters();  
           int hasBHadronDaughter = 0;
           for (unsigned int d=0; d<nDaughters; ++d) 
            {
              int daughterID = abs((*it)->daughter(d)->pdgId());
              if ( (daughterID/100)%10 == 5 || (daughterID/1000)%10 == 5 ) { hasBHadronDaughter = 1; cout << " Warning !! the hadron hasBHadronDaughter  = " << hasBHadronDaughter << " nBHadrons variable can not be used for gluonsplitting tagging flavour" << endl; break; }
            } 
         }

        }else{
          cout<<" Warning !! No matching between flavourjet collection and the selected PF jet number " << ak5_index << " review lines 835-881 "<< endl; 
        }
                  
       }
    //################################################################
    // B-discriminants from JetTagCollection
    //  negative vaules means not enough information 
    //  for the tagger algorithm produces the discriminant value    
    //################################################################
      jet_CSV[ak5_index]  = -999;
      jet_JBP [ak5_index] = -999;
      jet_JP [ak5_index]  = -999;
      jet_TCHP[ak5_index] = -999; 
      jet_TCHE[ak5_index] = -999; 
    
      //matching ak5CaloJets and ak5PFJets
      float dR2min = 999; 
      unsigned indexmin = 999;
      for (unsigned index = 0; index < myJets -> size(); ++index)
      {
       float dR2 = reco::deltaR2 (jet_eta[ak5_index], jet_phi[ak5_index], (*myJets)[index].eta(), (*myJets)[index].phi());   
       if (dR2 < dR2min) 
       { 
        dR2min = dR2;
        indexmin = index;
       }
      }
      dRmin_matching[ak5_index] = sqrt(dR2min);
      if (sqrt(dR2min) < 0.1) 
      {
       //cout <<" PFindex " << ak5_index << "  caloJet " << indexmin << endl; 
       jet_CSV[ak5_index]  =  (*tagHandle_CSV)[myJets->refAt(indexmin)];
       //cout << "Jet_CSV disc =  " << jet_CSV[ak5_index]<<endl;; 
       jet_JBP[ak5_index]  =  (*tagHandle_JBP)[myJets->refAt(indexmin)];
       jet_JP[ak5_index]   =  (*tagHandle_JP)[myJets->refAt(indexmin)];
       //cout << "Jet_JP disc =  " << jet_JP[ak5_index]<<endl;
       jet_TCHP[ak5_index] =  (*tagHandle_TCHP)[myJets->refAt(indexmin)];
       //cout << "Jet_TCHP disc =  " << jet_TCHP[ak5_index]<<endl; 
       jet_TCHE[ak5_index] =  (*tagHandle_TCHE)[myJets->refAt(indexmin)];
       //cout << "Jet_TCHE disc =  " << jet_TCHE[ak5_index] << endl; 
      } else 
        {
         cout <<" Warning!! There is not matching between ak5CaloJets and ak5PFJets => No b-tagging info available for the original PFJet number =  " << jetRefIndex << endl;
        }  

    //################################################################
    // B-tagging selected tracks from TrackIPTagInfoCollection 
    //################################################################
     /*cout << " JETS FROM TrackIPTagInfoCollection " << endl; 
     cout << " " << endl; 
     cout << "Found " << ip.size() << " TagInfo" << endl;
     cout << " " << endl; 
     */
     unsigned tagindex_test = 0; 
     float dR2min_test = 999; 
     unsigned indexmin_test = 999;
     float tagindexmin = 999; 

     TrackIPTagInfoCollection::const_iterator it = ip.begin();
     for(; it != ip.end(); it++)
     {
       float dR2_test = reco::deltaR2 (jet_eta[ak5_index], jet_phi[ak5_index], it->jet()->eta(), it->jet()->phi());   
       if (dR2_test < dR2min_test) 
       { 
        dR2min_test = dR2_test;
        indexmin_test = tagindex_test;
       }
    /*cout << "TagInfoCollection Jet number: "<< TagInfoJet_index << endl;  
      cout << "  pt = "  << it->jet()->pt() << endl;
      cout << "  eta = : " << it->jet()->eta() << endl;
      cout << "  phi =  " << it->jet()->phi() << endl;
      cout << " " << endl;  
      */
      tagindex_test ++; 
     }
     if (sqrt(dR2min_test) < 0.1) 
      { 
       tagindexmin = indexmin_test; 
       TrackRefVector selTracks = (*ipHandle)[tagindexmin].selectedTracks();
       int n=selTracks.size();
       if (n > 0)
        {
         //cout << " -------------------------- Selected tracks --------------------" << endl; 
         for(int j=0;j< n;j++)
         {
          jetSeltrackIndex [seltracksInEvent_index] = ak5_index;
          //cout << " track number in the event  " << seltracksInEvent_index << " associated to the selected jet number  " << jetSeltrackIndex [seltracksInEvent_index] << " has : " << endl;  
          //cout << " pt = " << selTracks[j]->pt() << "\n";
          seltrack_pt [seltracksInEvent_index] = selTracks[j]->pt();
          // Extract the Impact paramenter info for this track
          TrackIPTagInfo::TrackIPData data = (*ipHandle)[tagindexmin].impactParameterData()[j];  
          // Hit pattern of the track
          const reco::HitPattern& p = selTracks[j]->hitPattern();
          //const reco::HitPattern& p = (*i_trk)->hitPattern();
          // Loop over the hits of the track
          //for (int i=0; i<p.numberOfHits(); i++) {
          seltrack_nValidPixelHits[seltracksInEvent_index] = p.numberOfValidPixelHits() ;
          seltrack_nValidTrackerHits [seltracksInEvent_index] = p.numberOfValidTrackerHits() ;
          //}
          //cout << " numberOfValidPixelHits = "  << p.numberOfValidPixelHits() << "\n";    
          //cout << " numberOfValidTrackerHits = " << p.numberOfValidTrackerHits() << "\n";    
          //cout << " ip3d.value = " << data.ip3d.value() << "\n";
          //cout << " ip3d.significance = " << data.ip3d.significance() << "\n";
          //cout << " distanceToJetAxis.value = " << data.distanceToJetAxis.value() << "\n";
          //cout << " distanceToJetAxis.significance = " << data.distanceToJetAxis.significance() << "\n";
          // cout << data.distanceToGhostTrack.value() << "\t";
          // cout << data.distanceToGhostTrack.significance() << "\t";
          // cout << data.closestToJetAxis << "\t";
          // cout << (data.closestToJetAxis - pv).mag() << "\t";
          // cout << data.closestToGhostTrack << "\t";
          // cout << (data.closestToGhostTrack - pv).mag() << "\t";
          //cout <<  " ip2d.value = " << data.ip2d.value() << "\n";
          //cout <<  " ip2d.significance = " << data.ip2d.significance() <<  endl;     
          seltrack_IP2D [seltracksInEvent_index] = data.ip2d.value();  
          seltrack_IP2Dsig [seltracksInEvent_index] = data.ip2d.significance();  
          seltrack_IP3D [seltracksInEvent_index] = data.ip3d.value();  
          seltrack_IP3Dsig [seltracksInEvent_index] = data.ip3d.significance();
          seltrack_distToJetAxis [seltracksInEvent_index] = data.distanceToJetAxis.value();  
          //track_IPz [trackPV_index] = ???;
          seltracksInEvent_index ++; 
         }
        } 
       }  
     if (tagindexmin != indexmin) cout << " Recorcholis !! TagCollection Jets and TrackIPTagInfoCollection does not match one to one!! " << endl; 
     
 
    //################################################################
    // Secondary vertices 
    //################################################################
     unsigned svindex_test = 0; 
     float svdR2min_test = 999; 
     unsigned svindexmin_test = 999;
     float svindexmin = 999;

     nSVinJet[ak5_index] = 0;  
     // loop over the SecondaryVertexTagInfoCollection and match the ak5CaloJets with at least one secondary vertex 
     for(reco::SecondaryVertexTagInfoCollection::const_iterator iter = svTagInfoColl.begin(); iter != svTagInfoColl.end(); ++iter) 
     {
       if(iter->nVertices() > 0 )
       {
        float svdR2_test = reco::deltaR2 (jet_eta[ak5_index], jet_phi[ak5_index], iter->jet()->eta(), iter->jet()->phi());   
        if (svdR2_test < svdR2min_test) 
        { 
         svdR2min_test = svdR2_test;
         svindexmin_test = svindex_test;
        }
       }
      svindex_test ++; 
     }
     if (sqrt(svdR2min_test) < 0.1)
     {
      svindexmin = svindexmin_test; 
      unsigned int nVertices = (*svTagInfosHandle)[svindexmin].nVertices();  
      nSVinJet[ak5_index] = nVertices;

      //cout<< "the selected jet number " << ak5_index << " with   pt = " << jet_pt[ak5_index] << "  eta = " << jet_eta[ak5_index] << "  phi = " <<jet_phi[ak5_index]<< endl; 
      //cout<< "that matches with the calo jet number " << svindexmin << " with  pt = " << (*svTagInfosHandle)[svindexmin].jet() ->pt() << " eta = " << (*svTagInfosHandle)[svindexmin].jet() ->eta() << "  phi = " << (*svTagInfosHandle)[svindexmin].jet() ->phi() <<endl;  
      //cout<<" has " << nSVinJet[ak5_index]  << " secondary vertices " << endl;   
      for (unsigned int vtx = 0; vtx < nVertices; ++vtx)
         {
          jetSVIndex[nSVinEvent_index] = ak5_index;  
          flight3DSignificance[nSVinEvent_index] = (*svTagInfosHandle)[svindexmin].flightDistance(vtx,false).significance();
          const Vertex &vertex = (*svTagInfosHandle)[svindexmin].secondaryVertex(vtx); 
          svmass[nSVinEvent_index] = vertex.p4().mass(); 
         // std::cout<<"with" << endl;  
         // std::cout<<"        3D flight distance of " << (*svTagInfosHandle)[svindexmin].flightDistance(vtx,false).value() << " cm"<<std::endl;
         // std::cout<<"        3D flight significance of " << flight3DSignificance[nSVinEvent_index] << " cm"<<std::endl;
         // std::cout<<"        2D flight distance of " << (*svTagInfosHandle)[svindexmin].flightDistance(vtx,true).value() << " cm"<<std::endl;
         // std::cout<<"        2D flight significance of " << (*svTagInfosHandle)[svindexmin].flightDistance(vtx,true).significance() << " cm"<<std::endl;
         // std::cout<<"        invariant mass of     " <<  svmass[nSVinEvent_index]  << " GeV"<<endl;
         // std::cout<<"        invariant mass 2nd way of     " << (*svTagInfosHandle)[svindexmin].secondaryVertex(vtx).p4().mass() << " GeV"<<endl;
          nSVinEvent_index ++;
         }
       }
     if (nSVinJet[ak5_index] > 0 && svindexmin != indexmin) cout << " Recorcholis !! TagCollection Jets and SecondaryVertexTagInfoCollection does not match one to one!! " << endl; 
     
    /////////////////////////////////////////////////////////////////////////////////////////////////////////  
    

     //cout << "the number of TagInfoCollection Jets in the event is: " << TagInfoJet_index << endl; 
    //################################################################
    /* cout << " " << endl;
     cout << "PFCollection Selected Jet number: "<< njetCORR-1 << endl;    
     cout << "  pt = "  << jet_pt[ak5_index]   << endl;
     cout << "  eta = " << jet_eta[ak5_index]  << endl;
     cout << "  phi = " << jet_phi[ak5_index]  << endl;
     cout << " " << endl;
    *///################################################################
    
    ak5_index++;
    }  
  
   
    // Number of selected jets in the event
    njet = ak5_index;    
    // Number of tracks in the event
    tracks_inEvent = tracks_inEvent_index; 
    // Number of B-tag(IPTagInfo) selected tracks in the event
    seltracksInEvent = seltracksInEvent_index;
    // Number of secondary vertex in the event
    nSVinEvent = nSVinEvent_index;   
    
    //cout << " " << endl; 
    //cout << "the number of CORRECTED PF Jets in the event is: " << njetCORR << endl; 
    //cout << "the number of selected PF Jets in the event is: "  << njet << endl; 
    //cout << " " << endl; 

    
    // Four leading AK7 Jets

    edm::Handle<reco::PFJetCollection> ak7_handle;
    event_obj.getByLabel(mPFak7JetsName, ak7_handle);
    const JetCorrector* corrector_ak7 = JetCorrector::getJetCorrector(mJetCorr_ak7, iSetup);

    // Index of the selected jet 
    int ak7_index = 0;

    // Jets will be unsorted in pT after applying JEC,  
    // therefore store corrected jets in a new collection (map): key (double) is pT * -1 (key), 
    // value (std::pair<PFJet*, double>) is pair of original jet iterator and corresponding JEC factor
    sortedJets.clear();
    for (auto i_ak7jet_orig = ak7_handle->begin(); i_ak7jet_orig != ak7_handle->end(); ++i_ak7jet_orig) {
        // take jet energy correction and get corrected pT
        jec = corrector_ak7->correction(*i_ak7jet_orig, event_obj, iSetup);
        // Multiply pT by -1 in order to have largest pT jet first (sorted in ascending order by default)
        sortedJets.insert(std::pair<double, std::pair<reco::PFJetCollection::const_iterator, double> >(-1 * i_ak7jet_orig->pt() * jec, std::pair<reco::PFJetCollection::const_iterator, double>(i_ak7jet_orig, jec)));
    }

    // Iterate over the jets (sorted in pT) of the event
    for (auto i_ak7jet_orig = sortedJets.begin(); i_ak7jet_orig != sortedJets.end() && ak7_index < 4; ++i_ak7jet_orig) {

        // Apply jet energy correction "on the fly":
        // copy original (uncorrected) jet;
        PFJet corjet = *((i_ak7jet_orig->second).first);
        // take stored JEC factor
        jec = (i_ak7jet_orig->second).second;
        // apply JEC
        corjet.scaleEnergy(jec);
        // pointer for further use
        const PFJet* i_ak7jet = &corjet;

        // Skip the current iteration if jet is not selected
        if (fabs(i_ak7jet->y()) > mMaxY || 
            (i_ak7jet->pt()) < mMinPFPt) {
            continue;
        }

        // Variables of the tuple
        jet_area_ak7[ak7_index] = i_ak7jet->jetArea();
        jet_jes_ak7[ak7_index] = jec; // JEC factor

        // p4 is already corrected!
        auto p4 = i_ak7jet->p4();
        jet_pt_ak7[ak7_index]   = p4.Pt();
        jet_eta_ak7[ak7_index]  = p4.Eta();
        jet_phi_ak7[ak7_index]  = p4.Phi();
        jet_E_ak7[ak7_index]    = p4.E(); 
        
        // Matching AK5 jet to this AK7 jet
        // Index of the generated jet matching this PFjet
        ak7_to_ak5[ak7_index] = -1; // -1 if no matching jet

        float r2min(999);
        for (unsigned int ak5_index = 0; ak5_index != njet; ak5_index++) {

            // Compute distance squared
            double deltaR2 = reco::deltaR2( jet_eta_ak7[ak7_index], 
                                            jet_phi_ak7[ak7_index],
                                            jet_eta[ak5_index], 
                                            jet_phi[ak5_index]);
            if (deltaR2 < r2min) {
                r2min = deltaR2;
                ak7_to_ak5[ak7_index] = ak5_index;
            }
        }
        
    ak7_index++;
    }  
    // Number of saved jets in the event
    njet_ak7 = ak7_index;    


    // MET
    Handle< PFMETCollection > met_handle;
    event_obj.getByLabel("pfMet", met_handle);

    met = (*met_handle)[0].et();
    sumet = (*met_handle)[0].sumEt();

    // Finally, fill the tree
    if (njet >= (unsigned)mMinNPFJets && 
        njet_ak7 >= (unsigned)mMinNPFJets ) {            
            mTree->Fill();
    }
}


void OpenDataTreeProducerOptimized::endRun(edm::Run const &iRun, edm::EventSetup const &iSetup) {

}

OpenDataTreeProducerOptimized::~OpenDataTreeProducerOptimized() {
}


DEFINE_FWK_MODULE(OpenDataTreeProducerOptimized);
