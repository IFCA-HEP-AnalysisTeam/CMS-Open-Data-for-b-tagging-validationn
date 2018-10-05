#include "TH1D.h"
#include "TFile.h"
#include <fstream>
#include <iostream>


void PUforMCreweigth() 
{
  /*
  * Calculate the PU weigths distribution for MC reweighting
  */


  //--------------------------
  //2011 MC normalized profile
  //--------------------------
  const int nbins = 35; 
  TH1D* mcPUprofile = new TH1D("mcPUprofile", "", nbins, 0, nbins);
  // https://github.com/cms-sw/cmssw/blob/CMSSW_5_3_X/SimGeneral/MixingModule/python/mix_2011_FinalDist_OOTPU_cfi.py#L22-L58
  double yvalues[nbins] = {1.30976e-05,
                           0.000148266,
                           0.00226073,
                           0.030543,
                           0.0868303,
                           0.120295,
                           0.124687,
                           0.110419,
                           0.0945742,
                           0.0837875,
                           0.0774277,
                           0.0740595,
                           0.0676844,
                           0.0551203,
                           0.0378357,
                           0.0210203,
                           0.00918262,
                           0.00309786,
                           0.000808509,
                           0.000168568,
                           3.02344e-05,
                           5.16455e-06,
                           8.83185e-07,
                           1.43975e-07,
                           2.07228e-08,
                           2.51393e-09,
                           2.52072e-10,
                           2.07328e-11,
                           1.39369e-12,
                           7.63843e-14,
                           3.4069e-15,
                           1.23492e-16,
                           3.63059e-18,
                           8.53277e-20,
                           1.33668e-22}; 

   for (int x=1; x<=nbins; x++)
   {
     mcPUprofile -> SetBinContent(x,yvalues[x-1]);
     cout << "mcPUprofile: x = " << x << " ; y = " << yvalues[x-1] << endl;
   }
   // std::cout << "mcPUprofile Integral = " << mcPUprofile->Integral()<< std::endl;
   // normalize to unity
   mcPUprofile->Scale(1/mcPUprofile->Integral());
   cout << "               the number of bins = " << mcPUprofile->GetNbinsX() <<endl;

  //-----------------------------
  //2011 Data normalized profile 
  //-----------------------------
  TH1D* dataPUprofile = new TH1D("dataPUprofile", "", 35, 0, 35);
  //PileUp distribution calculated for a HLTpath & RunA 2011
  //TFile* truePUdata = new TFile("runAOpenDataPU_hdata.root", "read");
  TFile* truePUdata = new TFile("HLT_Jet60_RunA_PU_hdata.root", "read");
  TH1D*  dummy = (TH1D*) truePUdata->Get("pileup"); 
  // pileup histogram has 100 bins
  float dataOverflow = 0;  
  for (int x=1; x<101; x++)
  {
    if (x<36)  dataPUprofile -> SetBinContent(x, dummy->GetBinContent(x));
    // set the overflow in the last bin
    else if (x >= 36) dataOverflow =+ dummy->GetBinContent(x);    
  }
  dataPUprofile -> SetBinContent(36, dataOverflow); 
  cout << "  dataPUprofile histo has the number of bins = " << dataPUprofile->GetNbinsX() <<endl;
  //std::cout << "dataPUprofile Integral = " << dataPUprofile->Integral()<< std::endl;
   dataPUprofile->Scale(1/dataPUprofile->Integral());
  
  //-------------------------
  //2011 PU weights histogram 
  //-------------------------
  TH1D* truePUweigth = (TH1D*) dataPUprofile->Clone("truePUweigth");
  truePUweigth-> Divide(mcPUprofile);
  // save 
  TFile* truePUweigths = new TFile("truePUweigths.root", "recreate");
  mcPUprofile->Write();
  dataPUprofile->Write();
  truePUweigth->Write();
  truePUweigths->Close();
  // cosmetics 
  dataPUprofile->SetLineColor(kBlue);
  mcPUprofile->SetLineColor(kRed);
  truePUweigth->SetLineColor(kBlack);
  // draw
  truePUweigth-> Draw();
  dataPUprofile-> Draw("same");
  //dataPUprofile-> Draw();
  mcPUprofile->Draw("same");
 
}


void PUdata()
{
  TH1D* dataPUprofile       = new TH1D("dataPUprofile", "", 35, 0, 35);
  TH1D* runAdataPUprofile   = new TH1D("runAdataPUprofile", "", 35, 0, 35);
  TH1D* hltPUprofile        = new TH1D("hltdataPUprofile", "", 35, 0, 35);
  TH1D* runAhltPUprofile    = new TH1D("runAhltdataPUprofile", "", 35, 0, 35);
  
  TFile* truePUdata         = new TFile("PU_hdata.root", "read");
  TFile* runAtruePUdata     = new TFile("runAOpenDataPU_hdata.root", "read");
  TFile* hlttruePUdata      = new TFile("HLT_Jet60_PU_hdata.root", "read");
  TFile* runAhlttruePUdata  = new TFile("HLT_Jet60_RunA_PU_hdata.root", "read");
  
  TH1D*  dummy         = (TH1D*) truePUdata->Get("pileup"); 
  TH1D*  runAdummy     = (TH1D*) runAtruePUdata->Get("pileup"); 
  TH1D*  hltdummy      = (TH1D*) hlttruePUdata->Get("pileup"); 
  TH1D*  runAhltdummy  = (TH1D*) runAhlttruePUdata->Get("pileup"); 
  
  // pileup histogram has 100 bins
  float dataOverflow = 0;
  float runAdataOverflow = 0;
  float hltdataOverflow = 0;
  float runAhltdataOverflow = 0;
  for (int x=1; x<101; x++)
  {
    if (x<36) 
     { 
      dataPUprofile     -> SetBinContent(x, dummy->GetBinContent(x)); 
      runAdataPUprofile -> SetBinContent(x, runAdummy->GetBinContent(x)); 
      hltPUprofile      -> SetBinContent(x, hltdummy->GetBinContent(x));
      runAhltPUprofile  -> SetBinContent(x, runAhltdummy->GetBinContent(x));
     }
    // set the overflow in the last bin
    else if (x >= 36) 
     {
      dataOverflow =+ dummy->GetBinContent(x); 
      runAdataOverflow =+ runAdummy->GetBinContent(x);
      hltdataOverflow =+ hltdummy->GetBinContent(x); 
      runAhltdataOverflow =+ runAhltdummy->GetBinContent(x); 
     }
  }
  dataPUprofile     -> SetBinContent(36, dataOverflow);
  runAdataPUprofile -> SetBinContent(36, runAdataOverflow);
  hltPUprofile      -> SetBinContent(36, hltdataOverflow);
  runAhltPUprofile  -> SetBinContent(36, runAhltdataOverflow);
  
  cout << "  dataPUprofile histo has the number of bins = "         << dataPUprofile->GetNbinsX() <<endl;
  cout << "  runAdataPUprofile histo has the number of bins = "     << runAdataPUprofile->GetNbinsX() <<endl;
  cout << "  hltdataPUprofile histo has the number of bins = "      << hltPUprofile->GetNbinsX() <<endl;
  cout << "  runAhltdataPUprofile histo has the number of bins = "  << runAhltPUprofile->GetNbinsX() <<endl;
  
  dataPUprofile->Scale(1/dataPUprofile->Integral());
  runAdataPUprofile->Scale(1/runAdataPUprofile->Integral());
  hltPUprofile->Scale(1/hltPUprofile->Integral());
  runAhltPUprofile->Scale(1/runAhltPUprofile->Integral());

  dataPUprofile     ->SetLineColor(kOrange);
  runAdataPUprofile ->SetLineColor(kOrange+5);
  hltPUprofile      ->SetLineColor(kBlue);
  runAhltPUprofile  ->SetLineColor(kBlack);
 
  runAhltPUprofile->Draw();
  runAdataPUprofile->Draw("same");   
  dataPUprofile->Draw("same");
  hltPUprofile->Draw("same");
}
