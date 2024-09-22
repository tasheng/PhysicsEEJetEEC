#ifndef _EffCorrFactor_H_
#define _EffCorrFactor_H_
#include <iostream>
#include <vector>
using namespace std;

// root includes
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"

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
   _h_num->SetName(Form("num_%s", effArgName.c_str()));
   _h_den->SetName(Form("den_%s", effArgName.c_str()));
   _heff_2D = (TH2D*) _h_num->Clone(Form("eff_%s_normEE", effArgName.c_str()));
   _heff_2D->Divide(_h_den);
   _h_num->Write();
   _h_den->Write();
   _heff_2D->Write();
}

Float_t EffCorrFactor::efficiency(Float_t argBin)
{
  Float_t e = _heff_1D->GetBinContent(_heff_1D->FindBin(argBin));
  if(e < 0.00000000001){
    if(counter < 10){
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
  Float_t e = _heff_2D->GetBinContent(_heff_2D->FindBin(argBin, normEEBin));
  if(e < 0.00000000001){
    if(counter < 10){
      std::cout << "!!!Error on efficiency correction! Zero efficiency!!! " << _heff_1D->GetName() << "=" << argBin << std::endl << std::endl;
      if(counter == 9) std::cout << " !!!Greater than ten calls of this error... TERMINATING OUTPUT, PLEASE FIX" << std::endl;
      ++counter;
    }
    
    e = 1;
  }  
  return e;
}


#endif
