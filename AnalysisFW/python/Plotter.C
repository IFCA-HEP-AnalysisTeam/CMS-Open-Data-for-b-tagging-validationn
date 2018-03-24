#include "TH1F.h"
#include "TFile.h"
#include "TCanvas.h"
#include "THStack.h"
#include "TString.h"
#include "TColor.h"
#include "TAxis.h"
#include "TGaxis.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TAttMarker.h"
#include "TLatex.h"
#include <fstream>
#include <iostream>
#include "TPad.h"
#include "TFrame.h" 



// Use paperStyle of AnalysisCMS
// Change the hard code as lumi value
//-----------------------------------------------------------------------------
// DrawLatex 
//------------------------------------------------------------------------------
void DrawLatex(Font_t      tfont,
               Float_t     x,
	       Float_t     y,
	       Float_t     tsize,
	       Short_t     align,
	       const char* text,
	       Bool_t      setndc=true)
{
  TLatex* tl = new TLatex(x, y, text);

  tl->SetNDC      (setndc);
  tl->SetTextAlign( align);
  tl->SetTextFont ( tfont);
  tl->SetTextSize ( tsize);

  tl->Draw("same");
}

//-----------------------------------------------------------------------------
// Draw 
//------------------------------------------------------------------------------
void Plotter()
{

 
// bool UnityNorm = false; 
   bool DataNorm  = true;
   bool adjustXlimits = false;
   // set the xvalue from which you want to draw your histo
   float lowerXlimit = 0; float upperXlimit = 0; //This is applied on FindFirstBinAbove, FindLastBinAbove functions which return the first/last bin above lowerXlimit/upperXlimit  
  // Number of variables to plot
  const int nplot = 4;
  // Name of variables to plot
  TString vname [nplot] = { "pthat"//,"jet_pt", "jet_phi", "jet_eta"
                            //"jet_CSV",
                            //"jet_JBP",
                            //"jet_TCHP",
                            //"nSVertex" 
                          };

  TString xname [nplot] = {"#bf{pthat}"//,"#bf{jet p_{T}}", "#bf{jet #phi}", "#bf{jet #eta}" 
                            //,"CSV discriminator", "JBP discriminator", "TCHP discriminator", "nSVertex"
                           };
  TString units [nplot] = {"#bf{[GeV/c]}"//, "#bf{[GeV/c]}","#bf{[rad]}",     ""
                            //,                  "",                  "",                   "",         ""
                           }; 
  // Name of saved MC histograms to plot from each file
  const int nhistoMC = 1;
  TString hnameMC [nhistoMC]= {"allFlavours"
                                // "b_quark",
                                // "c_quark",
                                // "lgluon",
                                // "b_gsplitting"
                           };
  
  TString hnameData = "data";

  // Loop over nplot
  for (int i = 0; i < 1/*nplot*/; i++) { 
     
     TFile* infileData = new TFile ("Histo_" + vname[i] + "_Data.root", "read");
     TH1D* myHistoData = (TH1D*) infileData -> Get (hnameData);

     TFile* infileMC   = new TFile ("Histo_" + vname[i] + "_MC.root", "read");
     TH1D* myHistoMC[nhistoMC];
     
     // Integral
     float integralDat = myHistoData -> Integral();
     std::cout << "dat" << integralDat << std::endl;
     TH1D* allMC = (TH1D*) infileMC -> Get ("allFlavours");
     float integral = allMC -> Integral();
     std::cout << "mc" << integral << std::endl;
   
     // Loop over nhistoMC 
     for (int j = 0; j < nhistoMC; j++){
         myHistoMC[j] = (TH1D*) infileMC -> Get (hnameMC[j]);
         
         // Normalize to unity
         //if (UnityNorm && !DataNorm)  myHistoMC[j] -> Scale(1/integral);
         if (DataNorm) myHistoMC[j] -> Scale(integralDat/integral);
         //if (DataNorm && !UnityNorm) myHistoMC[j] -> Scale(1/integralDat);
     }
    
     // Normalize to unity 
     //if (UnityNorma && !DataNorma) myHistoData -> Scale(1/integralDat);

     //Create the directory to save the plots
     gSystem->mkdir("plots/", kTRUE);
     
     // Define the canvas
     TCanvas* currentcanvas = NULL;
     TPad* pad1 = NULL;
     TPad* pad2 = NULL;
     
     TCanvas* currentCanvas = new TCanvas ("Histo_" + vname[i], "", 900,600);
     //currentCanvas -> SetTitle("public 2011 run A    2.33/fb (7 TeV)");
     pad1 = new TPad();
     pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
     pad1->SetTopMargin   (0.08);
     pad1->SetBottomMargin(0.02);
     pad1->SetLeftMargin(0.15); 
     pad1->Draw();

     pad2 = new TPad("pad2", "pad2", 0, 0.0, 1, 0.31);
     pad2->SetTopMargin   (0.08);//0.09, 0.08
     pad2->SetBottomMargin(0.35);
     pad2->SetLeftMargin(0.15); 
     pad2->Draw();


    //--------------------------------------------------------------------------
    // pad1 : variable
    //--------------------------------------------------------------------------

     pad1->cd();
     gStyle->SetOptStat("");
     pad1->SetLogy();

     //axis title
     TAxis* xaxis = (TAxis*)myHistoData->GetXaxis();
     xaxis->SetTitleSize(.10);
     xaxis->SetTitle(xname[i] +units[i]);
     xaxis->SetLabelSize(0.08);

     // Adjust xaxis and yaxis

     //set the minimum
     float min = myHistoData->GetMinimum() > 1.e-1 ? myHistoData->GetMinimum(): 1.e-1; 
     float max = myHistoData->GetMaximum()*10.;     
     for (int j =1; j<=nhistoMC; j++) 
      { 
       if (myHistoMC[j]->GetMinimum()> 1.e-1 && myHistoMC[j]->GetMinimum() < myHistoData->GetMinimum()) min = myHistoMC[j]->GetMinimum();
       if (myHistoMC[j]->GetMaximum() > max) max = myHistoMC[j]->GetMaximum();
      }
     //Set the range user
     float binmin; float binmax;
     if(adjustXlimits)
     {
      binmin = (myHistoData->FindFirstBinAbove(lowerXlimit) < allMC->FindFirstBinAbove(lowerXlimit)) ? myHistoData->FindFirstBinAbove(lowerXlimit) : allMC->FindFirstBinAbove(lowerXlimit);
      binmax = (myHistoData->FindLastBinAbove(upperXlimit) > allMC->FindLastBinAbove(upperXlimit)) ? myHistoData->FindLastBinAbove(upperXlimit) : allMC->FindLastBinAbove(upperXlimit);
      xaxis->SetRangeUser(myHistoData->GetBinLowEdge(binmin), myHistoData->GetBinLowEdge(binmax)+myHistoData->GetBinWidth(binmax));
     }    
     // Set cosmetics 
     myHistoData -> SetMarkerStyle (kFullCircle);
     myHistoData -> SetMarkerColor  (kBlack);
     myHistoData -> SetMarkerSize(0.6);
           
     myHistoMC[0] -> SetFillColor(kRed);
     //myHistoMC[1] -> SetFillColor(kGreen+2);
     //myHistoMC[2] -> SetFillColor(kCyan);
     //myHistoMC[3] -> SetFillColor(kCyan+1);
     
     std::cout << "2" <<std::endl;   
     // Stack the MC histograms
     THStack* st1 = new THStack(vname[i], "");
     st1 -> Add(myHistoMC[0]);
     //st1 -> Add(myHistoMC[3]);
     //st1 -> Add(myHistoMC[1]);
     //st1 -> Add(myHistoMC[2]);

     //currentCanvas -> cd();
     //Set log Y scale
     //gPad-> SetLogy();
     // Draw
     st1 -> SetMinimum(min);
     st1 -> SetMaximum(max);
     
     //TAxis* yaxis = (TAxis*)st1->GetYaxis();
     //yaxis->SetTitle("#bf{entries}");
     //yaxis->SetTitleSize(0.09);
     //yaxis->SetTitleOffset(0);// 1
     //yaxis->CenterTitle();
     st1 -> Draw("hist");
     myHistoData -> Draw("ep,same");
     st1 -> GetYaxis()->SetTitle("#bf{entries}");
     st1 -> GetYaxis()->SetTitleSize(0.05);
     st1 -> GetYaxis()->SetTitleOffset(0.55);
     st1 -> GetYaxis()->CenterTitle();
     // apply the range user
     if(adjustXlimits) st1 -> GetXaxis()->SetRangeUser(myHistoData->GetBinLowEdge(binmin), myHistoData->GetBinLowEdge(binmax)+myHistoData->GetBinWidth(binmax));
     //pad1->GetFrame()->DrawClone();
     //pad1->RedrawAxis();
     pad1->Update();
     // draw header
     DrawLatex(61, 0.165, 0.945, 0.050, 11,  "CMS 2011");
     DrawLatex(52, 0.282, 0.945, 0.030, 11, "Public Data");
     DrawLatex(42, 0.895, 0.945, 0.050, 31, Form("%.2f fb^{-1} (7TeV)", 2.3));
     std::cout << "3" <<std::endl;   
     
     //--------------------------------------------------------------------------
     // pad2 : Data/MC 
     //--------------------------------------------------------------------------
     pad2->cd();
     TH1D* ratio = (TH1D*)myHistoData->Clone("ratio");
     ratio->SetTitle("");
     float ymin =0.4, ymax =1.6;
     for (Int_t ibin=1; ibin<=ratio->GetNbinsX(); ibin++) 
      {
       float ratioVal = myHistoMC[0]->GetBinContent(ibin) != 0 ? (myHistoData->GetBinContent(ibin)/myHistoMC[0]->GetBinContent(ibin)) : -999;
       float ratioErr = myHistoData->GetBinContent(ibin)  != 0 ? (myHistoData->GetBinError(ibin)/myHistoMC[0]->GetBinContent(ibin))    :    0;
      // if (myHistoMC[0]->GetBinContent(ibin) != 0) 
      //    ratioVal = myHistoData->GetBinContent(ibin)/myHistoMC[0]->GetBinContent(ibin);
       std::cout << "variable=  " << vname[i] << "  bin= " << ibin << "  data=  " << myHistoData->GetBinContent(ibin) 
       << "  mc=  "<< myHistoMC[0]->GetBinContent(ibin) << "  ratio  =  " << ratioVal <<std::endl; 
       ratio -> SetBinContent(ibin, ratioVal);
       ratio -> SetBinError(ibin, ratioErr);
      }
     // cosmetics
     ratio->SetMarkerColor(1);
     ratio->SetMarkerSize(0.8);
     ratio->SetMarkerStyle(kFullSquare);
     ratio->GetYaxis()->SetRangeUser(ymin, ymax);
     // ratio axix title
     TAxis* yaxis2 = (TAxis*)ratio->GetYaxis();
     yaxis2->SetLabelSize(.08);
     yaxis2->SetTitle("#bf{Data/MC}");
     yaxis2->SetTitleSize(.09);
     yaxis2->SetTitleOffset(0.3);
     yaxis2->CenterTitle();     
     // draw
     ratio->Draw("hist ep");
     currentCanvas ->Update();
     pad2->Update();
   
      
     // draw line in ratio=1
     float xlow; float xhigh; float rbins;  
     if(adjustXlimits)
     {
     // apply the range user
      xlow  = ratio->GetBinLowEdge(binmin);
      xhigh = ratio->GetBinLowEdge(binmax) + ratio->GetBinWidth(binmax); 
     }else{
     // default range
     rbins = ratio->GetNbinsX(); 
     xlow  = ratio->GetBinLowEdge(1);
     xhigh = ratio->GetBinLowEdge(rbins) + ratio->GetBinWidth(rbins); 
     }
     // Set the line  
     TLine *line=new TLine(xlow,1.0,xhigh,1.);
     line->SetLineColor(1);
     line->SetLineWidth(2);
     line->SetLineStyle(kDotted);
     line->Draw("same");
     currentCanvas->Modified();
     currentCanvas->Update();
     // Save the canvas
     currentCanvas -> Print("plots/"+ vname[i] + ".png");
  }
}     




