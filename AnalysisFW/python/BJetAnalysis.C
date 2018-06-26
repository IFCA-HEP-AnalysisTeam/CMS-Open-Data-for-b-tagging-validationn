#define BJetAnalysis_cxx
#include "BJetAnalysis.h"


 
// --------------------------------------------------------------------------------------
// Constructor
// --------------------------------------------------------------------------------------
BJetAnalysis::BJetAnalysis(TTree* tree) : JetAnalysisBase(tree)
{
 cout << " b-tagging validation of 2011 Legacy dataset " << endl; 
}

// --------------------------------------------------------------------------------------
// Loop
// --------------------------------------------------------------------------------------
void BJetAnalysis::Loop (TString _dataPath, bool _ismc, TString _ptRange)
{
// if (fChain == 0) return;
 filename = _dataPath;
 ismc = _ismc;
 cout << " ismc ? " << ismc << endl; 
 ptRange = _ptRange;
// open the input file
// TFile* infile  = new TFile( filename, "read" );
// get the tree
// TTree* mytree  = (TTree*) infile -> Get( "ak5ak7/OpenDataTree");
// nentries = (Int_t) mytree->GetEntries();
 BeginJob(filename, ismc);

 TH1::SetDefaultSumw2();

 
// Define histograms
// ----------------------------------------------------------------------------------
 if (ismc)
 { 
  for (int i = 0; i < nflavour; i++)
   {
     DefineHistograms(i , i);    
   }
 } else {
 DefineHistograms(0, nflavour+1);
 }

 
 // avgTrackMultiplicity  global variables
 // ----------------------------------------------------------------------------------
 // selected-track-multiplicity per jet-pt histogram
  float njet_ptCounter[nflavour][30];
  float ntracksxjet_ptCounter[nflavour][30];
  float avg[30];
 
  for (int fv = 0; fv < nflavour; fv ++)
  {
   for (int ibin = 0; ibin < 30; ibin ++)
   { 
    njet_ptCounter[fv][ibin] = 0;
    ntracksxjet_ptCounter[fv][ibin] = 0;
   }
  } 

//float njet_ptCounter[nflavour][30];
//float ntracksxjet_ptCounter[nflavour][30];
//int   ntracksxjet[njet];


// Loop over events
// ----------------------------------------------------------------------------------
 for (Long64_t jentry=0; jentry< nentries;jentry++) 
  {
   Long64_t ientry = LoadTree(jentry);
   if (ientry < 0) break;   
   fChain->GetEntry(jentry);
 
  // number of generated events (not stored in ntuples)
   if (filename.Contains("15to30"))   ngen = 9978850; //http://opendata.cern.ch/record/1366
   if (filename.Contains("30to50"))   ngen = 5837856; //http://opendata.cern.ch/record/1539
   if (filename.Contains("50to80"))   ngen = 5766430; //http://opendata.cern.ch/record/1555
   if (filename.Contains("80to120"))  ngen = 5867864; //http://opendata.cern.ch/record/1562
   if (filename.Contains("120to170")) ngen = 5963264; //http://opendata.cern.ch/record/1348 
   if (filename.Contains("170to300")) ngen = 5975592; //http://opendata.cern.ch/record/1369 
   if (filename.Contains("300to470")) ngen = 5975016; //http://opendata.cern.ch/record/1469 
   if (filename.Contains("470to600")) ngen = 3967154; //http://opendata.cern.ch/record/1553 
   //For data: http://opendata.cern.ch/record/21
  
  // set the event weigth 
   eventw = 1; // for data
   if (ismc) { eventw = 2; 
   //            eventw = (lumi*mcweight)/ngen; //for mc 
   //            cout << " mcweight = " << mcweight << endl;
   //            cout << " ngen = " << ngen << endl;
             }
   //eventw = 0.0005 / 0.1;
   cout << "eventw = " << eventw << endl;

  // set the passed trigger
  if (triggers[1] == false) continue;
  //if (triggernames != jt60 || triggers[1] == false) continue;
 
 
  // Fill histograms
  // --------------------------------------------------------------------------------------
  
  // avgTrackMultiplicity per event variables
  int   ntracksxjet[njet];

  // loop over the jets in the event
  for (int j = 0; j < njet; j++)
  {
   // baseline selection
   if (jet_pt[j] < 60) continue;
   if (jet_eta[j] > 2.4 || jet_eta[j] < -2.4) continue;
   if (jet_tightID[j] = false) continue; // check with the loose
   

   // all flavour (data and mc)
     // only jet variables
   jetPt  [0] -> Fill(jet_pt [j],  eventw);  
   jetEta [0] -> Fill(jet_eta[j],  eventw);  
   jetPhi [0] -> Fill(jet_phi[j],  eventw);
   TCHE   [0] -> Fill(jet_TCHE[j], eventw); 
   TCHP   [0] -> Fill(jet_TCHP[j], eventw); 
   JP     [0] -> Fill(jet_JP[j],   eventw); 
   JBP    [0] -> Fill(jet_JBP[j],  eventw); 
   CSV    [0] -> Fill(jet_CSV[j],  eventw); 

     // selected tracks variables
   ntracksxjet[j] = 0;
   for (int k = 0; k < seltracksInEvent; k++)
   {
    if (jetSeltrackIndex[k]==j)
    { 
     ntracksxjet[j] += 1 ;
     IP3D [0] -> Fill(seltrack_IP3D[k], eventw);
     IP3Dsignif [0] -> Fill(seltrack_IP3Dsig[k], eventw);
    }
   }
    // selected secondary vertices variables
   if (nSVinEvent== 0) nrSV [0] -> Fill(nSVinEvent, eventw);
   else if (nSVinEvent > 0)
   {
    for (int i =0; i < nSVinEvent; i++)
    { 
     if (jetSVIndex[i]==j)
     {
      massSV[0] -> Fill(svmass[i], eventw);
      flight3Dsignif[0] -> Fill(flight3DSignificance[i], eventw);
      nrSV  [0] -> Fill(1, eventw);
     }
    }
   }

    // calculate the input of the selected-track-multiplicity per jet-pt histogram
     // loop over the jet pt bins
   for (int m =60; m <= 350; m+= 10)
   { 
    //bin number
    int ibin = (m-60)/10;
    //select the jets for the current pt bin 
    if(jet_pt[j] < m+1 && jet_pt[j] >= m ) 
     { 
      // save the global variables to fill the histogram  
       njet_ptCounter[0][ibin] += 1;
       ntracksxjet_ptCounter[0][ibin] += ntracksxjet[j];
     }  
    } 
   
    // set the ordinary tracks cuts 
   bool plot_tracknrPixelHits_sel = false;
   bool plot_trackPt_sel = false;
   bool plot_distToJetAxis_sel = false;
 
   for (int p = 0; p < tracks_inEvent; p++)
   {
    if (tracks_jetIndex[p] == j)
    { 
     if (tracks_nValidTrackerHits[p] >= 8 && tracks_chi2[p] < 5 && tracks_IP2D[p] < 0.2 && tracks_IPz[p] < 17 && tracks_distToJetAxis[p] < 0.07 && tracks_pt[p] > 1 )
     {
       plot_tracknrPixelHits_sel = true;
       nrPixelHits[0] -> Fill(tracks_nValidPixelHits[p], eventw);
     }
     if (tracks_nValidTrackerHits[p] >= 8 && tracks_chi2[p] < 5 && tracks_IP2D[p] < 0.2 && tracks_IPz[p] < 17 && tracks_distToJetAxis[p] < 0.07 && tracks_nValidPixelHits[p] >= 2) 
     { 
       plot_trackPt_sel = true;
       trackPt[0] -> Fill(tracks_pt[p], eventw);
     }
     if (tracks_nValidTrackerHits[p] >= 8 && tracks_chi2[p] < 5 && tracks_IP2D[p] < 0.2 && tracks_IPz[p] < 17 && tracks_pt[p] > 1 && tracks_nValidPixelHits[p] >= 2) 
     {
      plot_distToJetAxis_sel = true;
      distanceToJetAxis[0] -> Fill (tracks_distToJetAxis[p], eventw);
     }
    }
   }

   // set the flavour
   int jetFlavour = -999; 
   if (ismc)
   {
    //float flavour = fabs(PartonF[j]);
    if (fabs(PartonF[j]) == 5) jetFlavour = 1; // "b_quark"
    if (fabs(PartonF[j]) == 4) jetFlavour = 2; // "c_quark"
    if (fabs(PartonF[j]) != 5 && fabs(PartonF[j]) != 4) jetFlavour = 3; // "lgluon"
    if (fabs(PartonF[j]) == 5 && nBHadrons[j] == 2) jetFlavour = 4;// "b_gsplitting"
    if (jetFlavour < 0 ) cout << " Warning! jetFlavour = " << jetFlavour << " this does not work! " << endl;      

    jetPt  [jetFlavour] -> Fill(jet_pt[j],   eventw); 
    jetEta [jetFlavour] -> Fill(jet_eta[j],  eventw);  
    jetPhi [jetFlavour] -> Fill(jet_phi[j],  eventw);
    TCHE   [jetFlavour] -> Fill(jet_TCHE[j], eventw); 
    TCHP   [jetFlavour] -> Fill(jet_TCHP[j], eventw); 
    JP     [jetFlavour] -> Fill(jet_JP[j],   eventw); 
    JBP    [jetFlavour] -> Fill(jet_JBP[j],  eventw); 
    CSV    [jetFlavour] -> Fill(jet_CSV[j],  eventw); 
   
    for (int k = 0; k < seltracksInEvent; k++)
    {
     if (jetSeltrackIndex[k]==j)
     { 
      IP3D [jetFlavour] -> Fill(seltrack_IP3D[k], eventw);
      IP3Dsignif [jetFlavour] -> Fill(seltrack_IP3Dsig[k], eventw);
     }
    }

    if (nSVinEvent== 0) nrSV [jetFlavour] -> Fill(nSVinEvent, eventw);
    else if (nSVinEvent > 0)
    {
     for (int i =0; i < nSVinEvent; i++)
     { 
      if (jetSVIndex[i]==j)
      {
       massSV[jetFlavour] -> Fill(svmass[i], eventw);
       flight3Dsignif[jetFlavour] -> Fill(flight3DSignificance[i], eventw);
       nrSV  [jetFlavour] -> Fill(1, eventw);
      }
     }
    }
    
    for (int p = 0; p < tracks_inEvent; p++)
   {
    if (tracks_jetIndex[p] == j)
    {
     if (plot_tracknrPixelHits_sel) nrPixelHits[jetFlavour] -> Fill(tracks_nValidPixelHits[p], eventw);
     if (plot_trackPt_sel)          trackPt[jetFlavour] -> Fill(tracks_pt[p], eventw);
     if (plot_distToJetAxis_sel )   distanceToJetAxis[jetFlavour] -> Fill (tracks_distToJetAxis[p], eventw);
    }
   }    
    
   // calculate the input of the selected-track-multiplicity per jet-pt histogram
     // loop over the jet pt bins
   for (int m =60; m <= 350; m+= 10)
   { 
    //bin number
    int ibin = (m-60)/10;
    //select the jets for the current pt bin 
    if(jet_pt[j] < m+1 && jet_pt[j] >= m ) 
     { 
      // save the global variables to fill the histogram  
       njet_ptCounter[jetFlavour][ibin] += 1;
       ntracksxjet_ptCounter[jetFlavour][ibin] += ntracksxjet[j];
     }  
   }
  }//end of setting mc flavour
 }//end of the jet loop

  // average track multiplicity histo
  for ( int s = 0; s < nflavour; s++)
  { 
    for (int ibin = 0; ibin < 30; ibin++)
    {
     avg[ibin] = 0;  
     avg[ibin] = ntracksxjet_ptCounter[s][ibin]/njet_ptCounter[s][ibin];
     avgTrackMultiplicity[s]->Fill( avg[ibin], eventw); 
    }
  }
  // Save histograms
  // --------------------------------------------------------------------------------------
    for (int fv = 0; fv < nflavour; fv ++)
    {
      if (ismc == false && fv > 0) break;  

      SaveHistograms (jetPt[fv],  "jet_pt", _ptRange);  
      SaveHistograms (jetEta[fv], "jet_eta", _ptRange);  
      SaveHistograms (jetPhi[fv], "jet_phi", _ptRange); 
      SaveHistograms (TCHE[fv], "jet_TCHE", _ptRange); 
      SaveHistograms (TCHP[fv], "jet_TCHP", _ptRange); 
      SaveHistograms (JP[fv], "jet_JP", _ptRange); 
      SaveHistograms (JBP[fv], "jet_JBP", _ptRange); 
      SaveHistograms (CSV[fv], "jet_CSV", _ptRange); 
 
      SaveHistograms (IP3D[fv],  "seltrack_IP3D", _ptRange);  
      SaveHistograms (IP3Dsignif[fv], "seltrack_IP3Dsignif", _ptRange);  
      SaveHistograms (avgTrackMultiplicity[fv], "avgTrackMultiplicity", _ptRange);  

      SaveHistograms (nrPixelHits[fv],   "tracks_nrPixelHits", _ptRange);  
      SaveHistograms (trackPt[fv],   "tracks_Pt", _ptRange);  
      SaveHistograms (distanceToJetAxis[fv],   "tracks_distanceToJetAxis", _ptRange);  

      SaveHistograms (flight3Dsignif[fv], "flight3Dsignif", _ptRange);  
      SaveHistograms (massSV[fv],         "massSV",         _ptRange);  
      SaveHistograms (nrSV[fv],           "nrSV",           _ptRange);
   }  
 }
} // End of Loop

// --------------------------------------------------------------------------------------
// Begin job
// --------------------------------------------------------------------------------------
void BJetAnalysis::BeginJob(TString filename, bool _ismc)
{
 cout<< " begin the job ..."<<endl;
 cout<< " reading the dataset"<<endl;
 cout<< " filename = "<<filename<<endl; 
 cout<< " ismc     = "<<ismc<<endl;
// cout<< " nentries = "<<nentries;
 if (ismc) cout << " the ptRange = " << ptRange << endl; 
}

// --------------------------------------------------------------------------------------
// Define Histograms
// --------------------------------------------------------------------------------------
void BJetAnalysis::DefineHistograms(int iflv, int sflv)
 {
  /* Function to define the histograms with a readable hname for plotter.C
  *   if ismc:  hname correspond to the flavour selection
  *        "b_quark", "c_quark", "lgluon", "b_gsplitting"
  *   if not ismc: hname correspond "data"
  */
 // cout << "  " << endl;  
 // cout << "  " << endl;  
 // cout << "  " << endl;  
 // cout << " defining histos " << endl;  
 // cout << " flavour " << iflv << endl;  
 // cout << " sflavour " << sflavour[sflv] << endl;  
 // cout << "  " << endl;  
 // cout << "  " << endl;  
 // cout << "  " << endl;  
  // jet variables
  jetPt [iflv] = new TH1F (sflavour[sflv], " " , 100,  0,  250);
  jetEta[iflv] = new TH1F (sflavour[sflv], " " , 50,  -5,  5);
  jetPhi[iflv] = new TH1F (sflavour[sflv], " " , 40,  -4,  4);
  // selected tracks variables
  IP3D                [iflv] = new TH1F (sflavour[sflv], " " , 100, -0.1, 0.1);  
  IP3Dsignif          [iflv] = new TH1F (sflavour[sflv], " " , 100, -35, 35);  
  avgTrackMultiplicity[iflv] = new TH1F (sflavour[sflv], " " , 30, 60, 350);
  // ordinary tracks
  nrPixelHits      [iflv] = new TH1F (sflavour[sflv], " ", 9, 0, 9); 
  trackPt          [iflv] = new TH1F (sflavour[sflv], " ", 150, 0, 15); 
  distanceToJetAxis[iflv] = new TH1F (sflavour[sflv], " ", 100, 0, 0.3); 
  // secondary vertex
  flight3Dsignif[iflv] = new TH1F (sflavour[sflv], " ", 50, 0, 80); 
  massSV        [iflv] = new TH1F (sflavour[sflv], " ", 50, 0, 8); 
  nrSV          [iflv] = new TH1F (sflavour[sflv], " ", 5, 0, 5); 
  // b-discriminants
  TCHE[iflv] = new TH1F (sflavour[sflv], " ", 50, 0, 30);  
  TCHP[iflv] = new TH1F (sflavour[sflv], " ", 50, 0, 30);
  JP  [iflv] = new TH1F (sflavour[sflv], " ", 50, 0, 2.5);
  JBP [iflv] = new TH1F (sflavour[sflv], " ", 50, 0, 8);
  CSV [iflv] = new TH1F (sflavour[sflv], " ", 50, 0, 1);
 }
 
// --------------------------------------------------------------------------------------
// Save Histograms
// --------------------------------------------------------------------------------------
void BJetAnalysis::SaveHistograms (TH1F* myHistogram, TString hRootName, TString ptRange)
{
 /*
 * This function save the histograms with a readable name for plotter.C
 *  if ismc: ptRange correspond to the pthat of the sample
 *  if not ismc: ptRange correspond to "data"
 */
 // Histo_jet_phi_pthat170to300_MC.root
 // Histo_jet_phi_Data.root  
 TString hname = " ";
 if (ismc) hname = "test_Histo_" +  hRootName + "_" + ptRange + "_MC.root"; 
 else      hname = "test_Histo_" +  hRootName + "_Data.root" ;
 TFile* OutFileName = new TFile( hname , "update"); 
 myHistogram -> Write();
 OutFileName -> Close();
}

// --------------------------------------------------------------------------------------
// Auxiliar Histograms Functions 
// --------------------------------------------------------------------------------------

float BJetAnalysis::xOverflow (TH1F* rootHisto)
 {
  /*
  * This function gives the xmax of the histogram
  */  
  float max_value = rootHisto-> GetBinLowEdge(rootHisto -> GetNbinsX()+1);
  return max_value;
 }

float BJetAnalysis::xUnderflow (TH1F* rootHisto)
  {
  /*
  * This function gives the xmin of the histogram
  */  
   float min_value = rootHisto-> GetBinLowEdge(1); 
   return min_value;
  }


void BJetAnalysis::ResetHistograms(TH1F* rootHisto)
 {
  rootHisto -> Reset();
 }
