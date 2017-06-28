#include "TH1F.h"
#include "TFile.h"
#include "TCanvas.h"
#include "THStack.h"
#include "TString.h"
#include "TColor.h"
#include "TAxis.h"
#include "TGaxis.h"
//#include "TStyle.h"
#include "TAttMarker.h"
#include <fstream>
#include <iostream>

//TCanvas* Plotter()
void Plotter()
{

 
   bool UnityNorma = false; 
   bool DataNorma  = false;

  // Number of variables to plot
  const int nplot = 5;
  // Name of variables to plot
  TString vname [nplot] = { "jet_pt",
                            "jet_CSV",
                            "jet_JBP",
                            "jet_TCHP",
                            "nSVertex" 
                         };


  // Name of saved MC histograms to plot from each file
  const int nhistoMC = 4;
  TString hnameMC [nhistoMC]= {  "b_quark",
                                 "c_quark",
                                 "lgluon",
                                 "b_gsplitting"
                           };
  
  TString hnameData = "data";

  // Loop over nplot
  for (int i = 0; i < nplot; i++) { 
     
     TFile* ifileData = new TFile ("Histo_" + vname[i] + "_Data.root", "read");
     TH1F* myHistoData = (TH1F*) ifileData -> Get (hnameData);

     TFile* ifileMC   = new TFile ("Histo_" + vname[i] + "_MC.root", "read");
     TH1F* myHistoMC[nhistoMC];
     
     // Integral
     float integralDat = myHistoData -> Integral();
     std::cout << "dat" << integralDat << std::endl;
     TH1F* allMC = (TH1F*) ifileMC -> Get ("allFlavours");
     float integral = allMC -> Integral();
     std::cout << "mc" << integral << std::endl;
     
     // Define the canvas
     TCanvas* currentCanvas = new TCanvas ("current_Histo_" + vname[i], "Histo_" + vname[i], 10,10,700,500);

     // Loop over nhistoMC 
     for (int j = 0; j < nhistoMC; j++){
         myHistoMC[j] = (TH1F*) ifileMC -> Get (hnameMC[j]);
         
         // Normalize to unity
         //if (UnityNorma && !DataNorma)  myHistoMC[j] -> Scale(1/integral);
         //if (DataNorma) 
         myHistoMC[j] -> Scale(integralDat/integral);
         //if (DataNorma && !UnityNorma) myHistoMC[j] -> Scale(1/integralDat);
     }
    
     // Normalize to unity 
     //if (UnityNorma && !DataNorma) myHistoData -> Scale(1/integralDat);
    
     // Set cosmetics 
     myHistoData -> SetMarkerStyle (kFullCircle);
     myHistoData -> SetMarkerColor  (kBlack);
     myHistoData -> SetMarkerSize(0.6);
 
     myHistoMC[0] -> SetFillColor(kRed);
     myHistoMC[1] -> SetFillColor(kGreen+2);
     myHistoMC[2] -> SetFillColor(kBlue+2);
     myHistoMC[3] -> SetFillColor(kCyan+1);
   
     // Stack the MC histograms
     THStack* st1 = new THStack(vname[i],vname[i]);
     st1 -> Add(myHistoMC[0]);
     st1 -> Add(myHistoMC[3]);
     st1 -> Add(myHistoMC[1]);
     st1 -> Add(myHistoMC[2]);

     currentCanvas -> cd();
     //Set log Y scale
     gPad-> SetLogy();
     // Draw
     st1 -> Draw("hist");
     myHistoData -> Draw("ep,same");
     // Save the canvas
     currentCanvas -> Print(vname[i] + ".png");
  }
}     




