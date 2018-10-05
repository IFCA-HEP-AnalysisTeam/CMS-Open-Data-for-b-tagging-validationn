#include "TH1F.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TLegend.h"
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

//-----------------------------------------------------
// Set Legend 
//-----------------------------------------------------

TLegend *myLegend ( )
{

 TLegend *leg1;
 leg1 = new TLegend(0.75, 0.9, 0.9, 0.7);
 leg1->SetFillColor(kWhite); leg1->SetBorderSize(0.);
 leg1->SetTextColor(1); leg1->SetTextSize(0.020);
 leg1->SetTextFont(50);
 //leg1->SetHeader(svar[vv] + "  " + lHeader[R]); //leg1_fit->SetMargin(0.2); 
 leg1->AddEntry((TObject*)0, " ", "");
                  
  return leg1; 
}


//-----------------------------------------------------
// MoveOverflows
//-----------------------------------------------------
//
//   For all histogram types: nbins, xlow, xup
//
//     bin = 0;       underflow bin
//     bin = 1;       first bin with low-edge xlow INCLUDED
//     bin = nbins;   last bin with upper-edge xup EXCLUDED
//     bin = nbins+1; overflow bin
//
//-----------------------------------------------------
void MoveOverflows(TH1F* hist, Float_t xmin, Float_t xmax)
{
  int nentries = hist->GetEntries();
  int nbins    = hist->GetNbinsX();
  
  TAxis* xaxis = (TAxis*)hist->GetXaxis();
  // Underflow
  //---------------------------------------------------------------------------- 
   if (xmin != -999)
    {
      Int_t   firstBin = -1;
      Float_t firstVal = 0;
      Float_t firstErr = 0;
      
      for (Int_t i=0; i<=nbins+1; i++)
	{
	  if (xaxis->GetBinLowEdge(i) < xmin)
	    {
	      firstVal += hist->GetBinContent(i);
	      firstErr += (hist->GetBinError(i)*hist->GetBinError(i));
	      hist->SetBinContent(i, 0);
	      hist->SetBinError  (i, 0);
	    }
	  else if (firstBin == -1)
	    {
	      firstVal += hist->GetBinContent(i);
	      firstErr += (hist->GetBinError(i)*hist->GetBinError(i));
	      firstBin = i;
	    }
	}

      firstErr = sqrt(firstErr);
  
      hist->SetBinContent(firstBin, firstVal);
      hist->SetBinError  (firstBin, firstErr);
    }

   // Overflow
   //----------------------------------------------------------------------------
   if (xmax != -999)
    {
      Int_t   lastBin = -1;
      Float_t lastVal = 0;
      Float_t lastErr = 0;
      
      for (Int_t i=nbins+1; i>=0; i--)
	{
	  Float_t lowEdge = xaxis->GetBinLowEdge(i);
      
	  if (lowEdge >= xmax)
	    {
	      lastVal += hist->GetBinContent(i);
	      lastErr += (hist->GetBinError(i)*hist->GetBinError(i));
	      hist->SetBinContent(i, 0);
	      hist->SetBinError  (i, 0);
	    }
	  else if (lastBin == -1)
	    {
	      lastVal += hist->GetBinContent(i);
	      lastErr += (hist->GetBinError(i)*hist->GetBinError(i));
	      lastBin = i;
	    }
	}

      lastErr = sqrt(lastErr);
  
      hist->SetBinContent(lastBin, lastVal);
      hist->SetBinError  (lastBin, lastErr);
    }
 
 hist->SetEntries(nentries);

} 

// ----------------------------------------------------------------------------
// Auxiliar -> to remove 
// -------------------------------------------------------------------------
void AuxiliarPlot (TString inFolder, TString outFolder)
{
  TFile* npV_MC       = new TFile( inFolder + "Histo_nPV_MC.root", "read");
  TFile* npV_Data     = new TFile( inFolder + "Histo_nPV_TotalData.root", "read");
  TFile* npV_rewMC = new TFile( inFolder + "Histo_nPV_puRew_MC.root", "read");
  TH1F* nPVdata    =  (TH1F*)npV_Data->Get("data");
  float nPVdata_Integral = nPVdata->Integral(); 
  TH1F* nPVmc      =  (TH1F*)npV_MC->Get("allflavours");
  float nPVmc_Integral = nPVmc ->Integral();
  TH1F* nPVmc_rew  =  (TH1F*)npV_rewMC->Get("allflavours");
  float nPVmc_rewIntegral = nPVmc_rew ->Integral();
  //Normallize MC to data
  nPVmc->Scale(nPVdata_Integral/nPVmc_Integral);
  nPVmc_rew->Scale(nPVdata_Integral/nPVmc_rewIntegral);

  // Define the canvas
  TCanvas* currentCanvas = NULL;
  TPad* pad1 = NULL;
  TPad* pad2 = NULL;

  currentCanvas = new TCanvas ("nPV_Histo_" , "", 900,600);
  TLegend* _leg = myLegend ();
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

  pad1->cd();
  gStyle->SetOptStat("");

  //axis title
  TAxis* xaxis = (TAxis*)nPVdata->GetXaxis();
  xaxis->SetTitleSize(.10);
  xaxis->SetTitle("number of PV");
  xaxis->SetLabelSize(0.08);

  _leg -> AddEntry(nPVdata,   "data", "p" );
  _leg -> AddEntry(nPVmc,     "allmc before reweighting", "p");
  _leg -> AddEntry(nPVmc_rew, "allmc after reweigthing", "p");
  //Cosmetics
  nPVdata -> SetMarkerStyle (kFullCircle);
  nPVdata -> SetMarkerColor  (kBlack);
  nPVdata -> SetMarkerSize(0.8);
  nPVmc -> SetMarkerStyle (kFullCircle);
  nPVmc -> SetMarkerColor  (kGreen + 1);
  nPVmc -> SetMarkerSize(0.8);
  nPVmc_rew -> SetMarkerStyle (kFullCircle);
  nPVmc_rew -> SetMarkerColor  (kRed + 1);
  nPVmc_rew -> SetMarkerSize(0.8);

  //draw
  nPVdata->Draw();
  nPVmc->Draw("same");
  nPVmc_rew->Draw("same");
 
  nPVdata->GetYaxis()->SetTitle("#bf{entries}");
  nPVdata -> GetYaxis()->SetTitleSize(0.05);
  nPVdata -> GetYaxis()->SetTitleOffset(0.55);
  nPVdata -> GetYaxis()->CenterTitle();

  // draw header
  DrawLatex(61, 0.165, 0.945, 0.050, 11,  "CMS 2011");
  DrawLatex(52, 0.282, 0.945, 0.030, 11, "Public Data");
  DrawLatex(42, 0.895, 0.945, 0.050, 31, Form("%.2f fb^{-1} (7TeV)", 2.33));
  _leg ->Draw();

  // pad2 : Data/MC 
  pad2->cd();
  TH1F* ratio = (TH1F*)nPVdata->Clone("ratio");
  TH1F* rewratio = (TH1F*)nPVdata->Clone("rewratio");

  ratio->SetTitle("");
  float rymin =0.50, rymax =1.50;
  for (Int_t ibin=1; ibin<=ratio->GetNbinsX(); ibin++) 
   {
    float ratioVal = nPVmc->GetBinContent(ibin) != 0 ? (nPVdata->GetBinContent(ibin)/nPVmc->GetBinContent(ibin)) : -999;
    float ratioErr = nPVdata->GetBinContent(ibin) != 0 ? (nPVdata->GetBinError(ibin)/nPVmc->GetBinContent(ibin)) : 0;
    ratio -> SetBinContent(ibin, ratioVal);
    ratio -> SetBinError(ibin, ratioErr);

    float rewratioVal = nPVmc_rew->GetBinContent(ibin) != 0 ? (nPVdata->GetBinContent(ibin)/nPVmc_rew->GetBinContent(ibin)) : -999;
    float rewratioErr = nPVdata->GetBinContent(ibin) != 0 ? (nPVdata->GetBinError(ibin)/nPVmc_rew->GetBinContent(ibin)) : 0;
    rewratio -> SetBinContent(ibin, rewratioVal);
    rewratio -> SetBinError(ibin, rewratioErr);
   }
  // cosmetics 
  ratio->SetMarkerColor(1);
  ratio->SetMarkerSize(0.8);
  ratio->SetMarkerStyle(kFullSquare);
  ratio->SetMarkerColor(kGreen + 1);
  rewratio->SetMarkerStyle(kFullCircle);
  rewratio->SetMarkerColor(kRed + 1);
  ratio->GetYaxis()->SetRangeUser(rymin, rymax);

  // ratio axis title
  TAxis* yaxis2 = (TAxis*)ratio->GetYaxis();
  yaxis2->SetLabelSize(.08);
  yaxis2->SetTitle("#bf{Data/MC}");
  yaxis2->SetTitleSize(.09);
  yaxis2->SetTitleOffset(0.3);
  yaxis2->CenterTitle();

  //draw
  ratio->Draw("hist ep");
  rewratio->Draw("same, hist ep");
  currentCanvas ->Update();
  pad2->Update();


  // draw line in ratio=1
  float xlow; float xhigh; float rbins;  
  /*if(adjustXlimits)
  {
  // apply the range user
   xlow  = ratio->GetBinLowEdge(binmin);
   xhigh = ratio->GetBinLowEdge(binmax) + ratio->GetBinWidth(binmax); 
  }else{
  // default range
  rbins = ratio->GetNbinsX(); 
  xlow  = ratio->GetBinLowEdge(1);
  xhigh = ratio->GetBinLowEdge(rbins) + ratio->GetBinWidth(rbins); 
  }*/
  // default range
  rbins = ratio->GetNbinsX(); 
  xlow  = ratio->GetBinLowEdge(1);
  xhigh = ratio->GetBinLowEdge(rbins) + ratio->GetBinWidth(rbins);
  // Set the line  
  TLine *line=new TLine(xlow,1.0,xhigh,1.);
  line->SetLineColor(1);
  line->SetLineWidth(2);
  line->SetLineStyle(kDotted);
  line->Draw("same");
  currentCanvas->Modified();
  currentCanvas->Update();

  // Save the canvas
  currentCanvas -> Print(outFolder+"/nPV_nPVrew" + ".png");
  currentCanvas -> Print(outFolder+"/nPV_nPVrew" + ".jpg");
  currentCanvas -> Print(outFolder+"/nPV_nPVrew" + ".pdf");

}
//-----------------------------------------------------------------------------
// Draw 
//------------------------------------------------------------------------------
void Plotter(TString inFolder, TString outFolder)
{
  
  //float lumi= 2.33; // (/fb) 2011 legacy runA 
  float efflumifb = 0.000287; // (0.000287/fb)  HlT_jet60v* of 2011 legacy runA 
  float efflumipb = 0.287; // (0.287/pb)  HlT_jet60v* of 2011 legacy runA 
  TString sefflumipb = 0.287;  
 
  bool effLumiNorm  = true;
  bool dataNorm  = false;
  bool moveoverflow = true;
  
  bool setLinY = false; //massSV is a linear plot 
  bool adjustXlimits = false;
  // set the xvalue from which you want to draw your histo
  float lowerXlimit = 0; float upperXlimit = 0; //This is applied on FindFirstBinAbove, FindLastBinAbove functions which return the first/last bin above lowerXlimit/upperXlimit  
  //bool UnityNorm = false; 
  
// Number of variables to plot
  const int nplot =18;
  //const int nplot =20;

  // Name of variables to plot
  TString vname [nplot] = { //"nPV_puRew"           ,  "nPV"                     ,
                            "jet_pt"              ,  "jet_phi"                 ,  "jet_eta"           , 
                            "seltrack_IP3Dsignif" ,  "seltrack_IP3D"           , 
                            "nrSV"                ,  "flight3Dsignif"          ,  "massSV"            ,
                            "jet_CSV"             ,  "jet_JBP"                 ,  "jet_JP"            , "jet_TCHP" , "jet_TCHE" ,
                            "tracks_Pt"           ,  "tracks_distanceToJetAxis",  "tracks_nrPixelHits",            
                            "avgTrackMultiplicity",  "dRmin_matching" };

  TString xname [nplot] = {//"#bf{nPV after PUreweigth}",   "#bf{nPV}"                   ,   
                           "#bf{jet p_{T}}"           ,   "#bf{jet #phi}"              , "#bf{jet #eta}"         , 
                           "#bf{3D IP significance}"  ,   "#bf{3D IP}"                 , 
                           "#bf{nr.of SV}"            ,   "#bf{3D flight significance}", "#bf{SV mass}"          ,  
                           "#bf{CSV discriminator}"   ,   "#bf{JBP discriminator}"     , "#bf{JP discriminator}" , "#bf{TCHP discriminator}", "#bf{TCHE discriminator}",
                           "#bf{track p_{T}}"         ,   "#bf{distance to jet axis}"  , "#bf{nr. of pixel hits}",
                           "#bf{jet p_{T}}"           ,    "dRminak5(CalovsPF)"};

  TString units [nplot] = { //""            ,    ""          ,
                            "#bf{[GeV/c]}",    "#bf{[rad]}",          ""       ,        
                              ""          ,    "#bf{[cm]}" , 
                              ""          ,       ""       ,     "#bf{[GeV/c]}",
                              ""          ,       ""       ,          ""       ,    ""  ,   "" ,
                            "#bf{[GeV/c]}",    "#bf{[cm]}" ,          ""       , 
                            "#bf{[GeV/c]}"        ""};
 
  // Name of saved MC histograms to plot from each file
  const int nhistoMC = 5;
  TString hnameMC [nhistoMC]= {  "b_quark",
                                 "c_quark",
                                 "lgluon",
                                 "b_gsplitting",
                                 "allflavours"
                               };
  
  TString hnameData = "data";

  // Loop over nplot
  for (int i = 0; i < nplot ; i++) { 
     
     TFile* infileData = new TFile (inFolder + "/Histo_" + vname[i] + "_TotalData.root", "read");
     //TFile* infileData = new TFile ("testHistoTightID_9july/test_Histo_" + vname[i] + "_Data_CHECKING1.root", "read");
     TH1F* myHistoData = (TH1F*) infileData -> Get (hnameData);

     TFile* infileMC   = new TFile (inFolder + "/Histo_" + vname[i] + "_MC.root", "read");
     //TFile* infileMC   = new TFile ("testHistoTightID_9july/test_Histo_" + vname[i] + "_MC.root", "read");
     TH1F* myHistoMC[nhistoMC];
     
     // Integral
     float integralDat = myHistoData -> Integral();
     cout<< " " <<endl; 
     std::cout << "Integral data :" << integralDat << std::endl;
     TH1F* allMC = (TH1F*) infileMC -> Get ("allflavours");
     float integral = allMC -> Integral();
     std::cout << "Integral mc" << integral << std::endl; 
     cout<< " " <<endl; 
   
     // Loop over nhistoMC 
     for (int j = 0; j < nhistoMC; j++)
     {
         myHistoMC[j] = (TH1F*) infileMC -> Get (hnameMC[j]);
         // Not normalize
         if (vname[i] == "avgTrackMultiplicity" || vname[i] == "dRmin_matching" ) { dataNorm = false; effLumiNorm = true;}

         if ((dataNorm== true || effLumiNorm == true) && (vname[i] == "avgTrackMultiplicity" || vname[i] == "dRmin_matching" )) cout << " Warning !! The figure " << vname[i] << " is reescaled !!" <<endl; //for debuging !!!! 
     cout<< " " <<endl; 
            //if (vname[i] == "avgTrackMultiplicity" || vname[i] == "dRmin_matching" ) cout << " Warning !! The figure " << vname[i] << " is reescaled !!" <<endl; //for debuging !!!!
            if (dataNorm) myHistoMC[j] -> Scale(integralDat/integral);
            // we need to apply a correction factor of 1000 to match the mc crosssection (pb) with lumi /fb
            // if (effLumiNorm) myHistoMC[j] -> Scale(1000*efflumifb/lumi); //for plots without PUreweigthing applied
            if (effLumiNorm) myHistoMC[j] -> Scale(efflumipb);
     }
     if (dataNorm) cout << "The figure " << vname[i] <<" has the MC normalized to Data " << endl; 
     else if (effLumiNorm) cout << "The figure "<< vname[i] <<" has the MC is normalized to this Lumi (/pb): " << efflumipb <<endl;
     else cout << "The figure "<< vname[i] <<" has the MC NOT normalized. "  <<endl; 
     cout<< " " <<endl; 
     // Normalize to unity 
     //if (UnityNorma && !dataNorma) myHistoData -> Scale(1/integralDat);

     //Create the directory to save the plots
     gSystem->mkdir(outFolder + "/", kTRUE);
     
     // Define the canvas
     TCanvas* currentCanvas = NULL;
     TPad* pad1 = NULL;
     TPad* pad2 = NULL;

     if (vname[i] == "massSV" || vname[i] == "avgTrackMultiplicity") setLinY = true;
     
     currentCanvas = new TCanvas ("test_Histo_" + vname[i], "", 900,600);
     TLegend* _leg = myLegend (); 
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
     if (setLinY == false) pad1->SetLogy();

     //axis title
     TAxis* xaxis = (TAxis*)myHistoData->GetXaxis();
     xaxis->SetTitleSize(.10);
     xaxis->SetTitle(xname[i] +units[i]);
     xaxis->SetLabelSize(0.08);

     // Adjust xaxis and yaxis

     //set the maximun and minimum on the y axis
     float ymin = myHistoData->GetMinimum() > 1.e-1 ? myHistoData->GetMinimum(): 1.e-1; 
     float ymax = myHistoData->GetMaximum()*10.;     
     for (int j =0; j<nhistoMC; j++) 
      { 
       if (myHistoMC[j]->GetMinimum()> 1.e-1 && myHistoMC[j]->GetMinimum() < myHistoData->GetMinimum()) ymin = myHistoMC[j]->GetMinimum();
       if (myHistoMC[j]->GetMaximum() > ymax) ymax = myHistoMC[j]->GetMaximum();
      }
     if (setLinY) ymax = myHistoData->GetMaximum() + 100.;       
     if (setLinY && vname[i] == "massSV") ymax = myHistoData->GetMaximum() + 5000.;      
     if (setLinY && vname[i] == "avgTrackMultiplicity") { ymax = myHistoData->GetMaximum() + 10; ymin = 0;}  
     setLinY = false;     
     //Set the range user on the x axis
     int nxbins = myHistoData->GetNbinsX();
     float xmin = myHistoData->GetBinLowEdge(1); 
     float xmax = myHistoData->GetBinLowEdge(nxbins+1);
     float binmin; float binmax;
     if(adjustXlimits)
     {
      binmin = (myHistoData->FindFirstBinAbove(lowerXlimit) < allMC->FindFirstBinAbove(lowerXlimit)) ? myHistoData->FindFirstBinAbove(lowerXlimit) : allMC->FindFirstBinAbove(lowerXlimit);
      binmax = (myHistoData->FindLastBinAbove(upperXlimit) > allMC->FindLastBinAbove(upperXlimit)) ? myHistoData->FindLastBinAbove(upperXlimit) : allMC->FindLastBinAbove(upperXlimit);
      xaxis->SetRangeUser(myHistoData->GetBinLowEdge(binmin), myHistoData->GetBinLowEdge(binmax)+myHistoData->GetBinWidth(binmax));
      xmin = myHistoData->GetBinLowEdge(binmin);
      xmax = myHistoData->GetBinLowEdge(binmax)+myHistoData->GetBinWidth(binmax);      
     }    
     // Set cosmetics 
     myHistoData -> SetMarkerStyle (kFullCircle);
     myHistoData -> SetMarkerColor  (kBlack);
     myHistoData -> SetMarkerSize(0.8);
     
    if (vname[i] == "avgTrackMultiplicity" || vname[i] == "dRmin_matching")
    {
     myHistoMC[0] -> SetLineColor(kRed);
     myHistoMC[1] -> SetLineColor(kGreen+1);
     myHistoMC[2] -> SetLineColor(kBlue+2);
     myHistoMC[3] -> SetLineColor(kAzure+6);
 
    }
    else
    {      
     myHistoMC[0] -> SetFillColor(kRed);
     myHistoMC[1] -> SetFillColor(kGreen+1);
     myHistoMC[2] -> SetFillColor(kBlue+2);
     myHistoMC[3] -> SetFillColor(kAzure+6);
    }
     if (moveoverflow)
     {
       MoveOverflows(myHistoData, xmin, xmax);     
       MoveOverflows(myHistoMC[0], xmin, xmax);     
       MoveOverflows(myHistoMC[1], xmin, xmax);     
       MoveOverflows(myHistoMC[2], xmin, xmax);     
       MoveOverflows(myHistoMC[3], xmin, xmax);     
     }

     // Set legend
     //_leg -> SetHeader("effective lumi HLT_Jet60 " + sefflumipb);  
     _leg -> AddEntry(myHistoData, "data" , "p");
     _leg -> AddEntry(myHistoMC[0], "b quark" , "f");
     _leg -> AddEntry(myHistoMC[3], "b from gluon splitting" , "f");
     _leg -> AddEntry(myHistoMC[1], "c quark" , "f");
     _leg -> AddEntry(myHistoMC[2], "uds quark or gluon" , "f");
 
     // Stack the MC histograms
     THStack* st1 = new THStack(vname[i], "");
     st1 -> Add(myHistoMC[0]);
     st1 -> Add(myHistoMC[3]);
     st1 -> Add(myHistoMC[1]);
     st1 -> Add(myHistoMC[2]);

     //currentCanvas -> cd();
     //Set log Y scale
     //gPad-> SetLogy();
     // Draw
     st1 -> SetMinimum(ymin);
     st1 -> SetMaximum(ymax);
     
     //TAxis* yaxis = (TAxis*)st1->GetYaxis();
     //yaxis->SetTitle("#bf{entries}");
     //yaxis->SetTitleSize(0.09);
     //yaxis->SetTitleOffset(0);// 1
     //yaxis->CenterTitle();
     if (vname[i] == "avgTrackMultiplicity" || vname[i] == "dRmin_matching") 
     {
     st1 -> Draw("nostack,hist");
     }
     else
     {
       st1 -> Draw("hist");
     }
     myHistoData -> Draw("ep,same");
     st1 -> GetYaxis()->SetTitle("#bf{entries}");
     if (vname[i] == "avgTrackMultiplicity") st1 -> GetYaxis()->SetTitle("#bf{track multiplicity}");
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
     DrawLatex(42, 0.895, 0.945, 0.050, 31, Form("%.2f fb^{-1} (7TeV)", 2.33));
     _leg ->Draw(); 
     //--------------------------------------------------------------------------
     // pad2 : Data/MC ratio 
     //--------------------------------------------------------------------------
     pad2->cd();
     TH1F* ratio = (TH1F*)myHistoData->Clone("ratio");
     ratio->SetTitle("");
     float rymin =0.50, rymax =1.50;
     for (Int_t ibin=1; ibin<=ratio->GetNbinsX(); ibin++) 
      {
       float ratioVal = myHistoMC[4]->GetBinContent(ibin) != 0 ? (myHistoData->GetBinContent(ibin)/myHistoMC[4]->GetBinContent(ibin)) : -999;
       float ratioErr = myHistoData->GetBinContent(ibin)  != 0 ? (myHistoData->GetBinError(ibin)/myHistoMC[4]->GetBinContent(ibin))    :    0;
      // if (myHistoMC[0]->GetBinContent(ibin) != 0) 
      //    ratioVal = myHistoData->GetBinContent(ibin)/myHistoMC[0]->GetBinContent(ibin);
      //  std::cout << "variable=  " << vname[i] << "  bin= " << ibin << "  data=  " << myHistoData->GetBinContent(ibin) << "  mc=  "<< myHistoMC[4]->GetBinContent(ibin) << "  ratio  =  " << ratioVal <<std::endl; 
       ratio -> SetBinContent(ibin, ratioVal);
       ratio -> SetBinError(ibin, ratioErr);
      }
     // cosmetics
     ratio->SetMarkerColor(1);
     ratio->SetMarkerSize(0.8);
     ratio->SetMarkerStyle(kFullSquare);
     ratio->GetYaxis()->SetRangeUser(rymin, rymax);
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
     currentCanvas -> Print(outFolder+"/"+ vname[i] + ".png");
     currentCanvas -> Print(outFolder+"/"+ vname[i] + ".jpg");
     currentCanvas -> Print(outFolder+"/"+ vname[i] + ".pdf");
  }
}     




