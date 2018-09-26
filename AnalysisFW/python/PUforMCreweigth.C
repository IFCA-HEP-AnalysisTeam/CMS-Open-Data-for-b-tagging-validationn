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
  int nbins = 35; 
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

   for (int x=1; x<nbins+1; x++)
   {
     mcPUprofile -> SetBinContent(x,yvalues[x-1]);
   }
   // std::cout << "mcPUprofile Integral = " << mcPUprofile->Integral()<< std::endl;
   mcPUprofile->Scale(1/mcPUprofile->Integral());


  //-----------------------------
  //2011 Data normalized profile 
  //-----------------------------
  TH1D* dataPUprofile = new TH1D("dataPUprofile", "", 35, 0, 35);
  // https://twiki.cern.ch/twiki/bin/viewauth/CMS/PileupJSONFileforData#Calculating_Your_Pileup_Distribu 
  // pileupCalc.py -i Cert_160404-180252_7TeV_ReRecoNov08_Collisions11_JSON.txt --inputLumiJSON  pileup_2011_JSON_pixelLumi.txt --calcMode true --minBiasXsec 68000 --maxPileupBin 100 --numPileupBins 100 PU_hdata.root 

  TFile* truePUdata = new TFile("PU_hdata.root", "read");
  TH1D*  dummy = (TH1D*) truePUdata->Get("pileup"); 

  for (int x=1; x<101; x++)
  {
    if (x<36)  dataPUprofile -> SetBinContent(x, dummy->GetBinContent(x));
    else if (x >= 36) dataPUprofile -> SetBinContent(36, dummy->GetBinContent(x)); // !!!!!! anadir suma
  }
 
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
  // draw 
  dataPUprofile->SetLineColor(kBlue);
  mcPUprofile->SetLineColor(kRed);
  truePUweigth->SetLineColor(kBlack);
  //truePUweigth-> Draw();
  //dataPUprofile-> Draw("same");
  dataPUprofile-> Draw();
  mcPUprofile->Draw("same");
 
}
 
