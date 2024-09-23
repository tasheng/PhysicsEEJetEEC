#ifndef _EffCorrFactor_H_
#define _EffCorrFactor_H_
#include <iostream>
#include <vector>
using namespace std;

// root includes
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"

#include "TMath.h"

class EffCorrFactor
{
public:
   TFile *_effInf=NULL;
   TH1D *_heff_1D=NULL;
   TH2D *_heff_2D=NULL;

   EffCorrFactor() {}; // default dummy constructor
   EffCorrFactor(string EffCorrFactorPath, string effArgName);
   ~EffCorrFactor();

   void init(string EffCorrFactorPath, string effArgName);
   void write(TFile& OutputFile, string effArgName, TH1D* _h_num, TH1D* _h_den);
   void write(TFile& OutputFile, string effArgName, TH2D* _h_num, TH2D* _h_den);

   // * apply the correction factor on the exact histogram binning
   //   would just be a multiplication operation
   // * this will return the error code (1) if the efficiency factor is null
   int applyEffCorrOnHisto(TH1D* h_1D_bfCorr, TH1D* h_1D_afCorr);
   int applyEffCorrOnHisto(TH2D* h_2D_bfCorr, TH2D* h_2D_afCorr);

   // * apply the correction factor on an entry-by-entry basis
   // * this will return eff=-999 if the efficiency factor is null
   Float_t efficiency(Float_t argBin);
   Float_t efficiency(Float_t argBin, Float_t normEEBin);
   Int_t counter;

};

EffCorrFactor::EffCorrFactor(string EffCorrFactorPath, string effArgName)
{
   init(EffCorrFactorPath, effArgName);
}

EffCorrFactor::~EffCorrFactor(){if(_effInf != NULL) delete _effInf;}

void EffCorrFactor::init(string EffCorrFactorPath, string effArgName)
{
   // std::string fullPath = std::getenv("ProjectBase");
   // if(fullPath.size() == 0)
   // {
   //     std::cout << "WARNING from " << __FILE__ << ", ENVIRONMENT Variable ProjectBase NOT SET" << std::endl;
   //     // std::cout << "Please \'source setStudyMultEnv.sh\' from top directory to fix" << std::endl;
   //     return;
   // }

   _effInf   = new TFile(EffCorrFactorPath.c_str(),"READ");
   _heff_1D  = (TH1D*) _effInf->Get(Form("eff_%s", effArgName.c_str()));
   _heff_2D  = (TH2D*) _effInf->Get(Form("eff_%s_normEE", effArgName.c_str()));
   counter   = 0;
   return;
}

void EffCorrFactor::write(TFile& OutputFile, string effArgName, TH1D* _h_num, TH1D* _h_den)
{
   OutputFile.cd();
   _h_num->SetName(Form("num_%s", effArgName.c_str()));
   _h_den->SetName(Form("den_%s", effArgName.c_str()));
   _heff_1D = (TH1D*) _h_num->Clone(Form("eff_%s", effArgName.c_str()));
   _heff_1D->Divide(_h_den);
   _h_num->Write();
   _h_den->Write();
   _heff_1D->Write();
}

void EffCorrFactor::write(TFile& OutputFile, string effArgName, TH2D* _h_num, TH2D* _h_den)
{
   OutputFile.cd();
   _h_num->SetName(Form("num_%s_normEE", effArgName.c_str()));
   _h_den->SetName(Form("den_%s_normEE", effArgName.c_str()));
   _heff_2D = (TH2D*) _h_num->Clone(Form("eff_%s_normEE", effArgName.c_str()));
   _heff_2D->Divide(_h_den);
   _h_num->Write();
   _h_den->Write();
   _heff_2D->Write();
}

int EffCorrFactor::applyEffCorrOnHisto(TH1D* h_1D_bfCorr, TH1D* h_1D_afCorr)
{
   if (!_heff_1D) 
   {
      printf("[Error] EffCorrFactor::applyEffCorrOnHisto couldn't find _heff_1D.\n");
      _effInf->ls();
      return 1;
   }

   for (Int_t ix=0;ix<=h_1D_bfCorr->GetNbinsX();ix++)
   { 
      // [Warning] A check-point after HP2024
      //           checking the bin boundaries, to see if it is safe to do the matrix multiplication
      if (h_1D_bfCorr->GetXaxis()->GetBinLowEdge(ix)!=_heff_1D->GetXaxis()->GetBinLowEdge(ix) ||
          h_1D_bfCorr->GetXaxis()->GetBinLowEdge(ix)!=h_1D_afCorr->GetXaxis()->GetBinLowEdge(ix))
      {
         printf("[Error] EffCorrFactor::applyEffCorrOnHisto is doing a matrix multiplication.\n" \
                "This requires the uncorrected, corrected distributions and the efficiency factor has the same bin configuration!\n");
         return 1;
      }

      double x = h_1D_bfCorr->GetXaxis()->GetBinCenter(ix);

      int hBinIdx    = h_1D_bfCorr->FindBin(x);
      int effBinIdx  = _heff_1D->FindBin(x);
      if(_heff_1D->GetBinContent(effBinIdx)>0)
      {
         double n_afCorr   = h_1D_bfCorr->GetBinContent(hBinIdx)/_heff_1D->GetBinContent(effBinIdx);
         double errrel_num = (h_1D_bfCorr->GetBinError(hBinIdx)>0 && h_1D_bfCorr->GetBinContent(hBinIdx)>0)? h_1D_bfCorr->GetBinError(hBinIdx)/h_1D_bfCorr->GetBinContent(hBinIdx): 0;
         double errrel_den = _heff_1D->GetBinError(effBinIdx)/_heff_1D->GetBinContent(effBinIdx);
         
         double errrel_n_afCorr = TMath::Sqrt(errrel_num*errrel_num+errrel_den*errrel_den);

         h_1D_afCorr->SetBinContent(hBinIdx,n_afCorr);
         h_1D_afCorr->SetBinError  (hBinIdx,n_afCorr*errrel_n_afCorr);
      }
      else
      {
         h_1D_afCorr->SetBinContent(hBinIdx,0);
         h_1D_afCorr->SetBinError(hBinIdx,0);
      }
      // if (ix<=2)
      // {
      //    printf("[ix=%d,x=%.3f] h_1D_bfCorr=%.3f+/-%.3f, _heff_1D=%.3f+/-%.3f, h_1D_afCorr=%.3f+/-%.3f\n", 
      //             ix, x,
      //             h_1D_bfCorr->GetBinContent(hBinIdx), h_1D_bfCorr->GetBinError(hBinIdx),
      //             _heff_1D->GetBinContent(effBinIdx), _heff_1D->GetBinError(effBinIdx),
      //             h_1D_afCorr->GetBinContent(hBinIdx), h_1D_afCorr->GetBinError(hBinIdx));
      // }
   }
   return 0;
}

int EffCorrFactor::applyEffCorrOnHisto(TH2D* h_2D_bfCorr, TH2D* h_2D_afCorr)
{
   if (!_heff_2D) 
   {
      printf("[Error] EffCorrFactor::applyEffCorrOnHisto couldn't find _heff_2D.\n");
      _effInf->ls();
      return 1;
   }
   for (Int_t ix=0;ix<=h_2D_bfCorr->GetNbinsX();ix++)
   { 
      // [Warning] A check-point after HP2024
      //           checking the bin boundaries, to see if it is safe to do the matrix multiplication
      if (h_2D_bfCorr->GetXaxis()->GetBinLowEdge(ix)!=_heff_2D->GetXaxis()->GetBinLowEdge(ix) ||
          h_2D_bfCorr->GetXaxis()->GetBinLowEdge(ix)!=h_2D_afCorr->GetXaxis()->GetBinLowEdge(ix))
      {
         printf("[Error] EffCorrFactor::applyEffCorrOnHisto is doing a matrix multiplication.\n" \
                "This requires the uncorrected, corrected distributions and the efficiency factor has the same bin configuration!\n");
         return 1;
      }

      for (Int_t iy=0;iy<=h_2D_bfCorr->GetNbinsY();iy++)
      {
         // [Warning] A check-point after HP2024
         //           checking the bin boundaries, to see if it is safe to do the matrix multiplication
         if (h_2D_bfCorr->GetYaxis()->GetBinLowEdge(iy)!=_heff_2D->GetYaxis()->GetBinLowEdge(iy) ||
             h_2D_bfCorr->GetYaxis()->GetBinLowEdge(iy)!=h_2D_afCorr->GetYaxis()->GetBinLowEdge(iy))
         {
            printf("[Error] EffCorrFactor::applyEffCorrOnHisto is doing a matrix multiplication.\n" \
                   "This requires the uncorrected, corrected distributions and the efficiency factor has the same bin configuration!\n");
            return 1;
         }

         double x = h_2D_bfCorr->GetXaxis()->GetBinCenter(ix);
         double y = h_2D_bfCorr->GetYaxis()->GetBinCenter(iy);

         int hBinIdx    = h_2D_bfCorr->FindBin(x,y);
         int effBinIdx  = _heff_2D->FindBin(x,y);
         if(_heff_2D->GetBinContent(effBinIdx)>0)
         {
            double n_afCorr   = h_2D_bfCorr->GetBinContent(hBinIdx)/_heff_2D->GetBinContent(effBinIdx);
            double errrel_num = (h_2D_bfCorr->GetBinError(hBinIdx)>0 && h_2D_bfCorr->GetBinContent(hBinIdx)>0)? h_2D_bfCorr->GetBinError(hBinIdx)/h_2D_bfCorr->GetBinContent(hBinIdx): 0;
            double errrel_den = _heff_2D->GetBinError(effBinIdx)/_heff_2D->GetBinContent(effBinIdx);
            
            double errrel_n_afCorr = TMath::Sqrt(errrel_num*errrel_num+errrel_den*errrel_den);

            h_2D_afCorr->SetBinContent(hBinIdx,n_afCorr);
            h_2D_afCorr->SetBinError  (hBinIdx,n_afCorr*errrel_n_afCorr);
         }
         else
         {
            h_2D_afCorr->SetBinContent(hBinIdx,0);
            h_2D_afCorr->SetBinError(hBinIdx,0);
         }
         // if (ix<=2 and iy <=2)
         // {
         //    printf("[ix=%d,x=%.3f][iy=%d,y=%.3f] h_2D_bfCorr=%.3f+/-%.3f, _heff_2D=%.3f+/-%.3f, h_2D_afCorr=%.3f+/-%.3f\n", 
         //             ix, x, iy, y,
         //             h_2D_bfCorr->GetBinContent(hBinIdx), h_2D_bfCorr->GetBinError(hBinIdx),
         //             _heff_2D->GetBinContent(effBinIdx), _heff_2D->GetBinError(effBinIdx),
         //             h_2D_afCorr->GetBinContent(hBinIdx), h_2D_afCorr->GetBinError(hBinIdx));
         // }
      }
   }
   return 0;
}

Float_t EffCorrFactor::efficiency(Float_t argBin)
{
   if (!_heff_1D) 
   {
      printf("[Error] EffCorrFactor::efficiency couldn't find _heff_1D.\n");
      _effInf->ls();
     return -999.;
   }
   Float_t e = _heff_1D->GetBinContent(_heff_1D->FindBin(argBin));
   if(e < 0.00000000001)
   {
      if(counter < 10)
      {
         std::cout << "!!!Error on efficiency correction! Zero efficiency!!! " << _heff_1D->GetName() << "=" << argBin << std::endl << std::endl;
         if(counter == 9) std::cout << " !!!Greater than ten calls of this error... TERMINATING OUTPUT, PLEASE FIX" << std::endl;
         ++counter;
      }
    e = 1;
   }  
  return e;
}

Float_t EffCorrFactor::efficiency(Float_t argBin, Float_t normEEBin)
{
   if (!_heff_2D) 
   {
      printf("[Error] EffCorrFactor::efficiency couldn't find _heff_2D.\n");
      _effInf->ls();
     return -999.;
   }
   Float_t e = _heff_2D->GetBinContent(_heff_2D->FindBin(argBin, normEEBin));
   if(e < 0.00000000001)
   {
      if(counter < 10)
      {
         std::cout << "!!!Error on efficiency correction! Zero efficiency!!! " << _heff_1D->GetName() << "=" << argBin << std::endl << std::endl;
         if(counter == 9) std::cout << " !!!Greater than ten calls of this error... TERMINATING OUTPUT, PLEASE FIX" << std::endl;
         ++counter;
      }
      e = 1;
   }  
  return e;
}


#endif
