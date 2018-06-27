//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Jun 26 18:54:39 2018 by ROOT version 5.32/00
// from TChain ak5ak7/OpenDataTree/
//////////////////////////////////////////////////////////

#ifndef ChainClass_h
#define ChainClass_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include <vector>

// Fixed size dimensions of array or collections stored in the TTree if any.

class ChainClass {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   UInt_t          njet;
   Float_t         jet_pt[11];   //[njet]
   Float_t         jet_eta[11];   //[njet]
   Float_t         jet_phi[11];   //[njet]
   Float_t         jet_E[11];   //[njet]
   Bool_t          jet_tightID[11];   //[njet]
   Float_t         jet_area[11];   //[njet]
   Float_t         jet_jes[11];   //[njet]
   Int_t           jet_igen[11];   //[njet]
   Float_t         jet_CSV[11];   //[njet]
   Float_t         jet_JP[11];   //[njet]
   Float_t         jet_JBP[11];   //[njet]
   Float_t         jet_TCHP[11];   //[njet]
   Float_t         jet_TCHE[11];   //[njet]
   Float_t         dRmin_matching[11];   //[njet]
   UInt_t          njet_ak7;
   Float_t         jet_pt_ak7[4];   //[njet_ak7]
   Float_t         jet_eta_ak7[4];   //[njet_ak7]
   Float_t         jet_phi_ak7[4];   //[njet_ak7]
   Float_t         jet_E_ak7[4];   //[njet_ak7]
   Float_t         jet_area_ak7[4];   //[njet_ak7]
   Float_t         jet_jes_ak7[4];   //[njet_ak7]
   Int_t           ak7_to_ak5[4];   //[njet_ak7]
   UInt_t          ngen;
   Float_t         gen_pt[10];   //[ngen]
   Float_t         gen_eta[10];   //[ngen]
   Float_t         gen_phi[10];   //[ngen]
   Float_t         gen_E[10];   //[ngen]
   UInt_t          run;
   UInt_t          lumi;
   ULong64_t       event;
   UInt_t          ntrg;
   Bool_t          triggers[6];   //[ntrg]
   vector<string>  *triggernames;
   UInt_t          prescales[6];   //[ntrg]
   Float_t         met;
   Float_t         sumet;
   Float_t         rho;
   Float_t         pthat;
   Float_t         mcweight;
   Float_t         chf[11];   //[njet]
   Float_t         nhf[11];   //[njet]
   Float_t         phf[11];   //[njet]
   Float_t         elf[11];   //[njet]
   Float_t         muf[11];   //[njet]
   Float_t         hf_hf[11];   //[njet]
   Float_t         hf_phf[11];   //[njet]
   UInt_t          hf_hm[11];   //[njet]
   UInt_t          hf_phm[11];   //[njet]
   UInt_t          chm[11];   //[njet]
   UInt_t          nhm[11];   //[njet]
   UInt_t          phm[11];   //[njet]
   UInt_t          elm[11];   //[njet]
   UInt_t          mum[11];   //[njet]
   Float_t         beta[11];   //[njet]
   Float_t         bstar[11];   //[njet]
   Float_t         nhfJet[11];   //[njet]
   Float_t         nemfJet[11];   //[njet]
   Float_t         chemfJet[11];   //[njet]
   UInt_t          chmJet[11];   //[njet]
   Bool_t          jet_looseID[11];   //[njet]
   Float_t         HadronF[11];   //[njet]
   Float_t         PartonF[11];   //[njet]
   Float_t         nBHadrons[11];   //[njet]
   UInt_t          nevent;
   Int_t           nSVinEvent;
   UInt_t          jetSVIndex[6];   //[nSVinEvent]
   Float_t         svmass[6];   //[nSVinEvent]
   Float_t         flight3DSignificance[6];   //[nSVinEvent]
   UInt_t          seltracksInEvent;
   Int_t           jetSeltrackIndex[147];   //[seltracksInEvent]
   Float_t         seltrack_pt[147];   //[seltracksInEvent]
   Int_t           seltrack_nValidPixelHits[147];   //[seltracksInEvent]
   Int_t           seltrack_nValidTrackerHits[147];   //[seltracksInEvent]
   Float_t         seltrack_IP2D[147];   //[seltracksInEvent]
   Float_t         seltrack_IP2Dsig[147];   //[seltracksInEvent]
   Float_t         seltrack_IP3D[147];   //[seltracksInEvent]
   Float_t         seltrack_IP3Dsig[147];   //[seltracksInEvent]
   Float_t         seltrack_distToJetAxis[147];   //[seltracksInEvent]
   UInt_t          tracks_inEvent;
   Int_t           tracks_jetIndex[351];   //[tracks_inEvent]
   Int_t           tracks_nValidPixelHits[351];   //[tracks_inEvent]
   Int_t           tracks_nValidTrackerHits[351];   //[tracks_inEvent]
   Float_t         tracks_pt[351];   //[tracks_inEvent]
   Float_t         tracks_chi2[351];   //[tracks_inEvent]
   Float_t         tracks_IPz[351];   //[tracks_inEvent]
   Float_t         tracks_IP2D[351];   //[tracks_inEvent]
   Double_t        tracks_distToJetAxis[351];   //[tracks_inEvent]
   Double_t        tracks_decayLength[351];   //[tracks_inEvent]

   // List of branches
   TBranch        *b_njet;   //!
   TBranch        *b_jet_pt;   //!
   TBranch        *b_jet_eta;   //!
   TBranch        *b_jet_phi;   //!
   TBranch        *b_jet_E;   //!
   TBranch        *b_jet_tightID;   //!
   TBranch        *b_jet_area;   //!
   TBranch        *b_jet_jes;   //!
   TBranch        *b_jet_igen;   //!
   TBranch        *b_jet_CSV;   //!
   TBranch        *b_jet_JP;   //!
   TBranch        *b_jet_JBP;   //!
   TBranch        *b_jet_TCHP;   //!
   TBranch        *b_jet_TCHE;   //!
   TBranch        *b_dRmin_matching;   //!
   TBranch        *b_njet_ak7;   //!
   TBranch        *b_jet_pt_ak7;   //!
   TBranch        *b_jet_eta_ak7;   //!
   TBranch        *b_jet_phi_ak7;   //!
   TBranch        *b_jet_E_ak7;   //!
   TBranch        *b_jet_area_ak7;   //!
   TBranch        *b_jet_jes_ak7;   //!
   TBranch        *b_ak7_to_ak5;   //!
   TBranch        *b_ngen;   //!
   TBranch        *b_gen_pt;   //!
   TBranch        *b_gen_eta;   //!
   TBranch        *b_gen_phi;   //!
   TBranch        *b_gen_E;   //!
   TBranch        *b_run;   //!
   TBranch        *b_lumi;   //!
   TBranch        *b_event;   //!
   TBranch        *b_ntrg;   //!
   TBranch        *b_triggers;   //!
   TBranch        *b_triggernames;   //!
   TBranch        *b_prescales;   //!
   TBranch        *b_met;   //!
   TBranch        *b_sumet;   //!
   TBranch        *b_rho;   //!
   TBranch        *b_pthat;   //!
   TBranch        *b_mcweight;   //!
   TBranch        *b_chf;   //!
   TBranch        *b_nhf;   //!
   TBranch        *b_phf;   //!
   TBranch        *b_elf;   //!
   TBranch        *b_muf;   //!
   TBranch        *b_hf_hf;   //!
   TBranch        *b_hf_phf;   //!
   TBranch        *b_hf_hm;   //!
   TBranch        *b_hf_phm;   //!
   TBranch        *b_chm;   //!
   TBranch        *b_nhm;   //!
   TBranch        *b_phm;   //!
   TBranch        *b_elm;   //!
   TBranch        *b_mum;   //!
   TBranch        *b_beta;   //!
   TBranch        *b_bstar;   //!
   TBranch        *b_nhfJet;   //!
   TBranch        *b_nemfJet;   //!
   TBranch        *b_chemfJet;   //!
   TBranch        *b_chmJet;   //!
   TBranch        *b_jet_looseID;   //!
   TBranch        *b_HadronF;   //!
   TBranch        *b_PartonF;   //!
   TBranch        *b_nBHadrons;   //!
   TBranch        *b_nevent;   //!
   TBranch        *b_nSVinEvent;   //!
   TBranch        *b_jetSVIndex;   //!
   TBranch        *b_svmass;   //!
   TBranch        *b_flight3DSignificance;   //!
   TBranch        *b_seltracksInEvent;   //!
   TBranch        *b_jetSeltrackIndex;   //!
   TBranch        *b_seltrack_pt;   //!
   TBranch        *b_seltrack_nValidPixelHits;   //!
   TBranch        *b_seltrack_nValidTrackerHits;   //!
   TBranch        *b_seltrack_IP2D;   //!
   TBranch        *b_seltrack_IP2Dsig;   //!
   TBranch        *b_seltrack_IP3D;   //!
   TBranch        *b_seltrack_IP3Dsig;   //!
   TBranch        *b_seltrack_distToJetAxis;   //!
   TBranch        *b_tracks_inEvent;   //!
   TBranch        *b_tracks_jetIndex;   //!
   TBranch        *b_tracks_nValidPixelHits;   //!
   TBranch        *b_tracks_nValidTrackerHits;   //!
   TBranch        *b_tracks_pt;   //!
   TBranch        *b_tracks_chi2;   //!
   TBranch        *b_tracks_IPz;   //!
   TBranch        *b_tracks_IP2D;   //!
   TBranch        *b_tracks_distToJetAxis;   //!
   TBranch        *b_tracks_decayLength;   //!

   ChainClass(TTree *tree=0);
   virtual ~ChainClass();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef ChainClass_cxx
ChainClass::ChainClass(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {

#ifdef SINGLE_TREE
      // The following code should be used if you want this class to access
      // a single tree instead of a chain
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("Memory Directory");
      if (!f || !f->IsOpen()) {
         f = new TFile("Memory Directory");
      }
      f->GetObject("ak5ak7/OpenDataTree",tree);

#else // SINGLE_TREE

      // The following code should be used if you want this class to access a chain
      // of trees.
      TChain * chain = new TChain("ak5ak7/OpenDataTree","");
      chain->Add("/eos/user/b/bchazinq/QCDPt300to470/output_0.root/ak5ak7/OpenDataTree");
      chain->Add("/eos/user/b/bchazinq/QCDPt300to470/output_1.root/ak5ak7/OpenDataTree");
      chain->Add("/eos/user/b/bchazinq/QCDPt300to470/output_2.root/ak5ak7/OpenDataTree");
      chain->Add("/eos/user/b/bchazinq/QCDPt300to470/output_3.root/ak5ak7/OpenDataTree");
      chain->Add("/eos/user/b/bchazinq/QCDPt300to470/output_4.root/ak5ak7/OpenDataTree");
      chain->Add("/eos/user/b/bchazinq/QCDPt300to470/output_5.root/ak5ak7/OpenDataTree");
      chain->Add("/eos/user/b/bchazinq/QCDPt300to470/output_6.root/ak5ak7/OpenDataTree");
      chain->Add("/eos/user/b/bchazinq/QCDPt300to470/output_7.root/ak5ak7/OpenDataTree");
      chain->Add("/eos/user/b/bchazinq/QCDPt300to470/output_8.root/ak5ak7/OpenDataTree");
      chain->Add("/eos/user/b/bchazinq/QCDPt300to470/output_9.root/ak5ak7/OpenDataTree");
      chain->Add("/eos/user/b/bchazinq/QCDPt300to470/output_10.root/ak5ak7/OpenDataTree");
      chain->Add("/eos/user/b/bchazinq/QCDPt300to470/output_11.root/ak5ak7/OpenDataTree");
      chain->Add("/eos/user/b/bchazinq/QCDPt300to470/output_12.root/ak5ak7/OpenDataTree");
      chain->Add("/eos/user/b/bchazinq/QCDPt300to470/output_13.root/ak5ak7/OpenDataTree");
      chain->Add("/eos/user/b/bchazinq/QCDPt300to470/output_14.root/ak5ak7/OpenDataTree");
      chain->Add("/eos/user/b/bchazinq/QCDPt300to470/output_15.root/ak5ak7/OpenDataTree");
      chain->Add("/eos/user/b/bchazinq/QCDPt300to470/output_16.root/ak5ak7/OpenDataTree");
      chain->Add("/eos/user/b/bchazinq/QCDPt300to470/output_17.root/ak5ak7/OpenDataTree");
      chain->Add("/eos/user/b/bchazinq/QCDPt300to470/output_18.root/ak5ak7/OpenDataTree");
      chain->Add("/eos/user/b/bchazinq/QCDPt300to470/output_19.root/ak5ak7/OpenDataTree");
      chain->Add("/eos/user/b/bchazinq/QCDPt300to470/output_20.root/ak5ak7/OpenDataTree");
      chain->Add("/eos/user/b/bchazinq/QCDPt300to470/output_21.root/ak5ak7/OpenDataTree");
      chain->Add("/eos/user/b/bchazinq/QCDPt300to470/output_22.root/ak5ak7/OpenDataTree");
      chain->Add("/eos/user/b/bchazinq/QCDPt300to470/output_23.root/ak5ak7/OpenDataTree");
      chain->Add("/eos/user/b/bchazinq/QCDPt300to470/output_24.root/ak5ak7/OpenDataTree");
      chain->Add("/eos/user/b/bchazinq/QCDPt300to470/output_25.root/ak5ak7/OpenDataTree");
      chain->Add("/eos/user/b/bchazinq/QCDPt300to470/output_26.root/ak5ak7/OpenDataTree");
      chain->Add("/eos/user/b/bchazinq/QCDPt300to470/output_27.root/ak5ak7/OpenDataTree");
      chain->Add("/eos/user/b/bchazinq/QCDPt300to470/output_28.root/ak5ak7/OpenDataTree");
      chain->Add("/eos/user/b/bchazinq/QCDPt300to470/output_29.root/ak5ak7/OpenDataTree");
      chain->Add("/eos/user/b/bchazinq/QCDPt300to470/output_30.root/ak5ak7/OpenDataTree");
      chain->Add("/eos/user/b/bchazinq/QCDPt300to470/output_31.root/ak5ak7/OpenDataTree");
      chain->Add("/eos/user/b/bchazinq/QCDPt300to470/output_32.root/ak5ak7/OpenDataTree");
      chain->Add("/eos/user/b/bchazinq/QCDPt300to470/output_33.root/ak5ak7/OpenDataTree");
      chain->Add("/eos/user/b/bchazinq/QCDPt300to470/output_34.root/ak5ak7/OpenDataTree");
      chain->Add("/eos/user/b/bchazinq/QCDPt300to470/output_35.root/ak5ak7/OpenDataTree");
      chain->Add("/eos/user/b/bchazinq/QCDPt300to470/output_36.root/ak5ak7/OpenDataTree");
      chain->Add("/eos/user/b/bchazinq/QCDPt300to470/output_37.root/ak5ak7/OpenDataTree");
      chain->Add("/eos/user/b/bchazinq/QCDPt300to470/output_38.root/ak5ak7/OpenDataTree");
      chain->Add("/eos/user/b/bchazinq/QCDPt300to470/output_39.root/ak5ak7/OpenDataTree");
      chain->Add("/eos/user/b/bchazinq/QCDPt300to470/output_40.root/ak5ak7/OpenDataTree");
      chain->Add("/eos/user/b/bchazinq/QCDPt300to470/output_41.root/ak5ak7/OpenDataTree");
      chain->Add("/eos/user/b/bchazinq/QCDPt300to470/output_42.root/ak5ak7/OpenDataTree");
      chain->Add("/eos/user/b/bchazinq/QCDPt300to470/output_43.root/ak5ak7/OpenDataTree");
      chain->Add("/eos/user/b/bchazinq/QCDPt300to470/output_44.root/ak5ak7/OpenDataTree");
      chain->Add("/eos/user/b/bchazinq/QCDPt300to470/output_45.root/ak5ak7/OpenDataTree");
      chain->Add("/eos/user/b/bchazinq/QCDPt300to470/output_46.root/ak5ak7/OpenDataTree");
      chain->Add("/eos/user/b/bchazinq/QCDPt300to470/output_47.root/ak5ak7/OpenDataTree");
      chain->Add("/eos/user/b/bchazinq/QCDPt300to470/output_48.root/ak5ak7/OpenDataTree");
      chain->Add("/eos/user/b/bchazinq/QCDPt300to470/output_49.root/ak5ak7/OpenDataTree");
      chain->Add("/eos/user/b/bchazinq/QCDPt300to470/output_50.root/ak5ak7/OpenDataTree");
      chain->Add("/eos/user/b/bchazinq/QCDPt300to470/output_51.root/ak5ak7/OpenDataTree");
      chain->Add("/eos/user/b/bchazinq/QCDPt300to470/output_52.root/ak5ak7/OpenDataTree");
      chain->Add("/eos/user/b/bchazinq/QCDPt300to470/output_53.root/ak5ak7/OpenDataTree");
      chain->Add("/eos/user/b/bchazinq/QCDPt300to470/output_54.root/ak5ak7/OpenDataTree");
      chain->Add("/eos/user/b/bchazinq/QCDPt300to470/output_55.root/ak5ak7/OpenDataTree");
      chain->Add("/eos/user/b/bchazinq/QCDPt300to470/output_56.root/ak5ak7/OpenDataTree");
      chain->Add("/eos/user/b/bchazinq/QCDPt300to470/output_57.root/ak5ak7/OpenDataTree");
      chain->Add("/eos/user/b/bchazinq/QCDPt300to470/output_58.root/ak5ak7/OpenDataTree");
      chain->Add("/eos/user/b/bchazinq/QCDPt300to470/output_59.root/ak5ak7/OpenDataTree");
      chain->Add("/eos/user/b/bchazinq/QCDPt300to470/output_60.root/ak5ak7/OpenDataTree");
      chain->Add("/eos/user/b/bchazinq/QCDPt300to470/output_61.root/ak5ak7/OpenDataTree");
      chain->Add("/eos/user/b/bchazinq/QCDPt300to470/output_62.root/ak5ak7/OpenDataTree");
      chain->Add("/eos/user/b/bchazinq/QCDPt300to470/output_63.root/ak5ak7/OpenDataTree");
      chain->Add("/eos/user/b/bchazinq/QCDPt300to470/output_64.root/ak5ak7/OpenDataTree");
      chain->Add("/eos/user/b/bchazinq/QCDPt300to470/output_65.root/ak5ak7/OpenDataTree");
      chain->Add("/eos/user/b/bchazinq/QCDPt300to470/output_66.root/ak5ak7/OpenDataTree");
      chain->Add("/eos/user/b/bchazinq/QCDPt300to470/output_67.root/ak5ak7/OpenDataTree");
      chain->Add("/eos/user/b/bchazinq/QCDPt300to470/output_68.root/ak5ak7/OpenDataTree");
      chain->Add("/eos/user/b/bchazinq/QCDPt300to470/output_69.root/ak5ak7/OpenDataTree");
      chain->Add("/eos/user/b/bchazinq/QCDPt300to470/output_70.root/ak5ak7/OpenDataTree");
      chain->Add("/eos/user/b/bchazinq/QCDPt300to470/output_71.root/ak5ak7/OpenDataTree");
      chain->Add("/eos/user/b/bchazinq/QCDPt300to470/output_72.root/ak5ak7/OpenDataTree");
      chain->Add("/eos/user/b/bchazinq/QCDPt300to470/output_73.root/ak5ak7/OpenDataTree");
      chain->Add("/eos/user/b/bchazinq/QCDPt300to470/output_74.root/ak5ak7/OpenDataTree");
      chain->Add("/eos/user/b/bchazinq/QCDPt300to470/output_75.root/ak5ak7/OpenDataTree");
      chain->Add("/eos/user/b/bchazinq/QCDPt300to470/output_76.root/ak5ak7/OpenDataTree");
      chain->Add("/eos/user/b/bchazinq/QCDPt300to470/output_77.root/ak5ak7/OpenDataTree");
      chain->Add("/eos/user/b/bchazinq/QCDPt300to470/output_78.root/ak5ak7/OpenDataTree");
      chain->Add("/eos/user/b/bchazinq/QCDPt300to470/output_79.root/ak5ak7/OpenDataTree");
      chain->Add("/eos/user/b/bchazinq/QCDPt300to470/output_80.root/ak5ak7/OpenDataTree");
      chain->Add("/eos/user/b/bchazinq/QCDPt300to470/output_81.root/ak5ak7/OpenDataTree");
      chain->Add("/eos/user/b/bchazinq/QCDPt300to470/output_82.root/ak5ak7/OpenDataTree");
      chain->Add("/eos/user/b/bchazinq/QCDPt300to470/output_83.root/ak5ak7/OpenDataTree");
      chain->Add("/eos/user/b/bchazinq/QCDPt300to470/output_84.root/ak5ak7/OpenDataTree");
      chain->Add("/eos/user/b/bchazinq/QCDPt300to470/output_85.root/ak5ak7/OpenDataTree");
      chain->Add("/eos/user/b/bchazinq/QCDPt300to470/output_86.root/ak5ak7/OpenDataTree");
      chain->Add("/eos/user/b/bchazinq/QCDPt300to470/output_87.root/ak5ak7/OpenDataTree");
      chain->Add("/eos/user/b/bchazinq/QCDPt300to470/output_88.root/ak5ak7/OpenDataTree");
      chain->Add("/eos/user/b/bchazinq/QCDPt300to470/output_89.root/ak5ak7/OpenDataTree");
      chain->Add("/eos/user/b/bchazinq/QCDPt300to470/output_90.root/ak5ak7/OpenDataTree");
      chain->Add("/eos/user/b/bchazinq/QCDPt300to470/output_91.root/ak5ak7/OpenDataTree");
      chain->Add("/eos/user/b/bchazinq/QCDPt300to470/output_92.root/ak5ak7/OpenDataTree");
      chain->Add("/eos/user/b/bchazinq/QCDPt300to470/output_93.root/ak5ak7/OpenDataTree");
      chain->Add("/eos/user/b/bchazinq/QCDPt300to470/output_94.root/ak5ak7/OpenDataTree");
      chain->Add("/eos/user/b/bchazinq/QCDPt300to470/output_95.root/ak5ak7/OpenDataTree");
      chain->Add("/eos/user/b/bchazinq/QCDPt300to470/output_96.root/ak5ak7/OpenDataTree");
      chain->Add("/eos/user/b/bchazinq/QCDPt300to470/output_97.root/ak5ak7/OpenDataTree");
      chain->Add("/eos/user/b/bchazinq/QCDPt300to470/output_98.root/ak5ak7/OpenDataTree");
      chain->Add("/eos/user/b/bchazinq/QCDPt300to470/output_99.root/ak5ak7/OpenDataTree");
      tree = chain;
#endif // SINGLE_TREE

   }
   Init(tree);
}

ChainClass::~ChainClass()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t ChainClass::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t ChainClass::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void ChainClass::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   triggernames = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("njet", &njet, &b_njet);
   fChain->SetBranchAddress("jet_pt", jet_pt, &b_jet_pt);
   fChain->SetBranchAddress("jet_eta", jet_eta, &b_jet_eta);
   fChain->SetBranchAddress("jet_phi", jet_phi, &b_jet_phi);
   fChain->SetBranchAddress("jet_E", jet_E, &b_jet_E);
   fChain->SetBranchAddress("jet_tightID", jet_tightID, &b_jet_tightID);
   fChain->SetBranchAddress("jet_area", jet_area, &b_jet_area);
   fChain->SetBranchAddress("jet_jes", jet_jes, &b_jet_jes);
   fChain->SetBranchAddress("jet_igen", jet_igen, &b_jet_igen);
   fChain->SetBranchAddress("jet_CSV", jet_CSV, &b_jet_CSV);
   fChain->SetBranchAddress("jet_JP", jet_JP, &b_jet_JP);
   fChain->SetBranchAddress("jet_JBP", jet_JBP, &b_jet_JBP);
   fChain->SetBranchAddress("jet_TCHP", jet_TCHP, &b_jet_TCHP);
   fChain->SetBranchAddress("jet_TCHE", jet_TCHE, &b_jet_TCHE);
   fChain->SetBranchAddress("dRmin_matching", dRmin_matching, &b_dRmin_matching);
   fChain->SetBranchAddress("njet_ak7", &njet_ak7, &b_njet_ak7);
   fChain->SetBranchAddress("jet_pt_ak7", jet_pt_ak7, &b_jet_pt_ak7);
   fChain->SetBranchAddress("jet_eta_ak7", jet_eta_ak7, &b_jet_eta_ak7);
   fChain->SetBranchAddress("jet_phi_ak7", jet_phi_ak7, &b_jet_phi_ak7);
   fChain->SetBranchAddress("jet_E_ak7", jet_E_ak7, &b_jet_E_ak7);
   fChain->SetBranchAddress("jet_area_ak7", jet_area_ak7, &b_jet_area_ak7);
   fChain->SetBranchAddress("jet_jes_ak7", jet_jes_ak7, &b_jet_jes_ak7);
   fChain->SetBranchAddress("ak7_to_ak5", ak7_to_ak5, &b_ak7_to_ak5);
   fChain->SetBranchAddress("ngen", &ngen, &b_ngen);
   fChain->SetBranchAddress("gen_pt", gen_pt, &b_gen_pt);
   fChain->SetBranchAddress("gen_eta", gen_eta, &b_gen_eta);
   fChain->SetBranchAddress("gen_phi", gen_phi, &b_gen_phi);
   fChain->SetBranchAddress("gen_E", gen_E, &b_gen_E);
   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("lumi", &lumi, &b_lumi);
   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("ntrg", &ntrg, &b_ntrg);
   fChain->SetBranchAddress("triggers", triggers, &b_triggers);
   fChain->SetBranchAddress("triggernames", &triggernames, &b_triggernames);
   fChain->SetBranchAddress("prescales", prescales, &b_prescales);
   fChain->SetBranchAddress("met", &met, &b_met);
   fChain->SetBranchAddress("sumet", &sumet, &b_sumet);
   fChain->SetBranchAddress("rho", &rho, &b_rho);
   fChain->SetBranchAddress("pthat", &pthat, &b_pthat);
   fChain->SetBranchAddress("mcweight", &mcweight, &b_mcweight);
   fChain->SetBranchAddress("chf", chf, &b_chf);
   fChain->SetBranchAddress("nhf", nhf, &b_nhf);
   fChain->SetBranchAddress("phf", phf, &b_phf);
   fChain->SetBranchAddress("elf", elf, &b_elf);
   fChain->SetBranchAddress("muf", muf, &b_muf);
   fChain->SetBranchAddress("hf_hf", hf_hf, &b_hf_hf);
   fChain->SetBranchAddress("hf_phf", hf_phf, &b_hf_phf);
   fChain->SetBranchAddress("hf_hm", hf_hm, &b_hf_hm);
   fChain->SetBranchAddress("hf_phm", hf_phm, &b_hf_phm);
   fChain->SetBranchAddress("chm", chm, &b_chm);
   fChain->SetBranchAddress("nhm", nhm, &b_nhm);
   fChain->SetBranchAddress("phm", phm, &b_phm);
   fChain->SetBranchAddress("elm", elm, &b_elm);
   fChain->SetBranchAddress("mum", mum, &b_mum);
   fChain->SetBranchAddress("beta", beta, &b_beta);
   fChain->SetBranchAddress("bstar", bstar, &b_bstar);
   fChain->SetBranchAddress("nhfJet", nhfJet, &b_nhfJet);
   fChain->SetBranchAddress("nemfJet", nemfJet, &b_nemfJet);
   fChain->SetBranchAddress("chemfJet", chemfJet, &b_chemfJet);
   fChain->SetBranchAddress("chmJet", chmJet, &b_chmJet);
   fChain->SetBranchAddress("jet_looseID", jet_looseID, &b_jet_looseID);
   fChain->SetBranchAddress("HadronF", HadronF, &b_HadronF);
   fChain->SetBranchAddress("PartonF", PartonF, &b_PartonF);
   fChain->SetBranchAddress("nBHadrons", nBHadrons, &b_nBHadrons);
   fChain->SetBranchAddress("nevent", &nevent, &b_nevent);
   fChain->SetBranchAddress("nSVinEvent", &nSVinEvent, &b_nSVinEvent);
   fChain->SetBranchAddress("jetSVIndex", jetSVIndex, &b_jetSVIndex);
   fChain->SetBranchAddress("svmass", svmass, &b_svmass);
   fChain->SetBranchAddress("flight3DSignificance", flight3DSignificance, &b_flight3DSignificance);
   fChain->SetBranchAddress("seltracksInEvent", &seltracksInEvent, &b_seltracksInEvent);
   fChain->SetBranchAddress("jetSeltrackIndex", jetSeltrackIndex, &b_jetSeltrackIndex);
   fChain->SetBranchAddress("seltrack_pt", seltrack_pt, &b_seltrack_pt);
   fChain->SetBranchAddress("seltrack_nValidPixelHits", seltrack_nValidPixelHits, &b_seltrack_nValidPixelHits);
   fChain->SetBranchAddress("seltrack_nValidTrackerHits", seltrack_nValidTrackerHits, &b_seltrack_nValidTrackerHits);
   fChain->SetBranchAddress("seltrack_IP2D", seltrack_IP2D, &b_seltrack_IP2D);
   fChain->SetBranchAddress("seltrack_IP2Dsig", seltrack_IP2Dsig, &b_seltrack_IP2Dsig);
   fChain->SetBranchAddress("seltrack_IP3D", seltrack_IP3D, &b_seltrack_IP3D);
   fChain->SetBranchAddress("seltrack_IP3Dsig", seltrack_IP3Dsig, &b_seltrack_IP3Dsig);
   fChain->SetBranchAddress("seltrack_distToJetAxis", seltrack_distToJetAxis, &b_seltrack_distToJetAxis);
   fChain->SetBranchAddress("tracks_inEvent", &tracks_inEvent, &b_tracks_inEvent);
   fChain->SetBranchAddress("tracks_jetIndex", tracks_jetIndex, &b_tracks_jetIndex);
   fChain->SetBranchAddress("tracks_nValidPixelHits", tracks_nValidPixelHits, &b_tracks_nValidPixelHits);
   fChain->SetBranchAddress("tracks_nValidTrackerHits", tracks_nValidTrackerHits, &b_tracks_nValidTrackerHits);
   fChain->SetBranchAddress("tracks_pt", tracks_pt, &b_tracks_pt);
   fChain->SetBranchAddress("tracks_chi2", tracks_chi2, &b_tracks_chi2);
   fChain->SetBranchAddress("tracks_IPz", tracks_IPz, &b_tracks_IPz);
   fChain->SetBranchAddress("tracks_IP2D", tracks_IP2D, &b_tracks_IP2D);
   fChain->SetBranchAddress("tracks_distToJetAxis", tracks_distToJetAxis, &b_tracks_distToJetAxis);
   fChain->SetBranchAddress("tracks_decayLength", tracks_decayLength, &b_tracks_decayLength);
   Notify();
}

Bool_t ChainClass::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void ChainClass::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t ChainClass::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef ChainClass_cxx
