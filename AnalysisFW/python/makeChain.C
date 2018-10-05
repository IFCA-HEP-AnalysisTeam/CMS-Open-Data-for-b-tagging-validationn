#include "TSystem.h"
#include <fstream>
#include <iostream>
#include <string>
#include "TChain.h"
#include "TH1.h"
#include "TString.h"
#include "TROOT.h"


void  makeChain (TString infoldername)
{
  ////gSystem->cd(foldername);
  ////gSystem->Init();
  ////int a = gSystem->Exec("ls " + infoldername + "/output_*.root | wc");
/*  TString a = gSystem->GetFromPipe("ls " + infoldername + "/output_*.root | wc -l");
  TString b = gSystem->pwd();
  cout << " a = " << a << endl; 
  cout << " typeid(a).name() = " << typeid(a).name() << endl;  
  cout << " typeid(TString).name() = " << typeid(TString).name() << endl;  
  cout << " b = " << b << endl;
  int c = atoi( a );
  cout << " c = " << c << endl;  
  cout << " typeid(c).name() = " << typeid(c).name() << endl;  
  cout << " typeid(int).name() = " << typeid(int).name() << endl; 

  TString d;
  d += c; 
  cout << " d = " << c << endl;  
  cout << " typeid(d).name() = " << typeid(d).name() << endl;  
 */  
  ////  It does not work 
  ////  -----------------
  ////  TString e = std::to_string(c); 
  ////  cout << " e = " << c << endl;  
  ////  cout << " typeid(e).name() = " << typeid(e).name() << endl;


  // Create a TChain to read several output_*.root trees in the infoldername folder.
  TChain mychain ("ak5ak7/OpenDataTree");
  TString a = gSystem->GetFromPipe("ls " + infoldername + "output_*.root | wc -l");
  int b = atoi( a );
  for (int nfile=0; nfile < b; nfile++)
   {
    TString numb;
    numb += nfile;
    mychain.Add(infoldername + "output_" + numb + ".root");
   }
 int nentries = (Int_t) mychain.GetEntries();
 cout << "The chain has " << nentries << " entries" << endl;
 /*Loop
  for (Long64_t jentry=0; jentry< nentries; jentry++) 
  {
   //Long64_t ientry = LoadTree(jentry); 
   //if (ientry < 0) break;
   mychain.GetEntry(jentry);
  } 
 */

 // Draw a variable in an histogram
 //TString b_variable = "jet_pt";
 //TString b_variable = "njet";
 //TString b_variable = "tracks_inEvent";
 //----------------------------------------------------------------------------
 // (no) cuts
 // =======
 // TString b_variable = "triggernames";
 // TString b_variable2 = "ntrg";
 // mychain.Scan("run:triggernames");
  //mychain.Scan("run:triggernames[1]");
  //mychain.Scan("run:triggernames:triggers:jet_pt",  "1>0");
  //mychain.Scan("run:triggernames:jet_pt",  "triggernames[1] != \"jt60\"");
 //mychain.Draw(b_variable + ">>" + b_variable, "triggers[1]&&jet_eta>-2.4&&jet_eta<2.4&&jet_pt>60&&jet_looseID");
 //mychain.Draw(b_variable + ">>" + b_variable, "1>0");
 //
 // jet index dependent variables
 // =============================
 // TString b_variable = "jet_pt"; TString h_dimension = "(100,  0,  250)";
 // TString b_variable = "jet_eta"; TString h_dimension = "(50,  -5,  5)";
 // TString b_variable = "jet_phi"; TString h_dimension = "(40,  -4,  4)";
 // mychain.Draw(b_variable + ">>" + b_variable + h_dimension , "triggers[1] && jet_eta <2.4 && jet_eta >-2.4 && jet_pt >60 && jet_tightID");  
 //
 // Selected trackindex  dependent variables
 // ========================================
 // TString b_variable = "seltrack_IP3D"; TString h_dimension = "(100, -0.1, 0.1)";
 // mychain.Draw(b_variable + ">>" + b_variable + h_dimension , "triggers[1] && jet_eta[jetSeltrackIndex] <2.4 && jet_eta[jetSeltrackIndex] >-2.4 && jet_pt[jetSeltrackIndex] >60 && jet_tightID[jetSeltrackIndex]");
 
 // Secondary Vertex index dependent variables
 // ==========================================
 //TString b_variable = "nSVinEvent"; TString h_dimension = "(5, 0, 5)";
 //TString b_variable = "flight3DSignificance"; TString h_dimension = "(50, 0, 80)";
 //TString b_variable = "svmass"; TString h_dimension = "(50, 0, 8)";
 //TString b_variable = "jetSVIndex"; TString h_dimension = "(8, 0, 8)";
 //mychain.Draw(b_variable + ">>" + b_variable + h_dimension , "(triggers[1] && jet_eta[jetSVIndex] <2.4 && jet_eta[jetSVIndex] >-2.4 && jet_pt[jetSVIndex] >60 && jet_tightID[jetSVIndex])*(53122368*2.33/5837856)");
 //mychain.Draw(b_variable + ">>" + b_variable, "(triggers[1] && jet_eta[jetSVIndex] <2.4 && jet_eta[jetSVIndex] >-2.4 && jet_pt[jetSVIndex] >60 && jet_tightID[jetSVIndex])*(53122368*2.33/5837856)");
 
 // Ordinary good tracks index  dependent variables
 // ===============================================
  //TString b_variable = "goodtracks_nValidPixelHits"; TString h_dimension = "(9, 0, 9)";
  //mychain.Draw(b_variable + ">>" + b_variable + h_dimension,"goodtracks_distToJetAxis < 0.07 && goodtracks_pt > 1 && triggers[1] && jet_eta[goodtracks_jetIndex] <2.4 && jet_eta[goodtracks_jetIndex] >-2.4 && jet_pt[goodtracks_jetIndex] >60 && jet_looseID[goodtracks_jetIndex]");

 //----------------------------------------------------------------------------
 
 /*TH1F* myHisto = new TH1F();
 myHisto = (TH1F*) gDirectory -> Get(b_variable);
 cout << " Integral = " << myHisto->Integral() << endl;
 myHisto -> Draw();
*/ 
 mychain.MakeClass("ChainClass");   
}

