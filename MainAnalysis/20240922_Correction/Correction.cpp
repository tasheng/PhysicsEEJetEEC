#include <iostream>
#include <vector>
#include <map>
using namespace std;

// root includes
#include "TTree.h"
#include "TChain.h"
#include "TFile.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"
#include "TStyle.h"
#include "TH1D.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TGaxis.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TMath.h"

#include "Messenger.h"
#include "CommandLine.h"
#include "ProgressBar.h"
#include "TauHelperFunctions3.h"
#include "SetStyle.h"
#include "EffCorrFactor.h"



int main(int argc, char *argv[]);
int FindBin(double Value, int NBins, double Bins[]); 
void MakeCanvasZ(vector<TH1D > Histograms, TGraphErrors DataSyst, vector<string> Labels, string Output, string X, string Y, double WorldMin, double WorldMax, bool DoRatio, bool LogX); 
void MakeCanvas(vector<TH1D> Histograms, TGraphErrors DataSyst, vector<string> Labels, string Output, string X, string Y, double WorldMin, double WorldMax, bool DoRatio, bool LogX); 
void SetPad(TPad &P); 
void DivideByBin(TH1D &H, double Bins[]); 

// can we put this in the CommonCode/ ?
void FillChain(TChain &chain, const vector<string> &files) {
  for (auto file : files) {
    chain.Add(file.c_str());
  }
}

// JC: the projection is copied here just so a quick check of the projected 1D histogram is accessible
void projection(TH2D* h_2D, TH1D* h_1D)
{
   h_1D->Reset(); 
   for (int i = 1; i <= h_2D->GetNbinsX(); ++i) {
      double weight = 0;
      double error = 0; 
      for (int j = 1; j <= h_2D->GetNbinsY(); ++j) {
         double binContent = h_2D->GetBinContent(i, j);
         double binError= h_2D->GetBinError(i,j);
         double binCenter = h_2D->GetYaxis()->GetBinCenter(j);
         weight += binContent*((binCenter));
         error += pow(binError*binCenter, 2);;

      }
      h_1D->SetBinContent(i, weight);
      h_1D->SetBinError(i, sqrt(error));
   }
}

int main(int argc, char *argv[])
{
   CommandLine CL(argc, argv);

   SetThesisStyle();
   static vector<int> Colors = GetCVDColors6();

   string InputDataPath          = CL.Get("InputData");
   string HistoName              = CL.Get("HistoName", "Bayesian_Unfoldediter3_Z");
   string OutputFileName         = CL.Get("Output", "CorrectedData.root");

   TFile InputData(InputDataPath.c_str(), "READ");
   TH2D HDataBfCorr( *((TH2D*) InputData.Get(HistoName.c_str())) );
   TH2D HDataAfCorr( *((TH2D*) HDataBfCorr.Clone(Form("%s_afCorr", HistoName.c_str()))) );

   // HDataBfCorr.Print("all");

   EffCorrFactor matchingEffCorrFactor;
   matchingEffCorrFactor.init("../../Unfolding/20240328_Unfolding/matchingScheme2/MatchingEff.root", "z");
   matchingEffCorrFactor.applyEffCorrOnHisto(&HDataBfCorr, &HDataAfCorr);

   EffCorrFactor EvtSelEffCorrFactor;
   EvtSelEffCorrFactor.init("../../EventSelectionEfficiency/20240922_evtSelEffCorr/EvtSelEff.root", "z");
   EvtSelEffCorrFactor.applyEffCorrOnHisto(&HDataBfCorr, &HDataAfCorr);

   TFile Output(OutputFileName.c_str(), "RECREATE");
   Output.cd();
   HDataBfCorr.Write();
   HDataAfCorr.Write();

   // JC: the projection is copied here just so a quick check of the projected 1D histogram is accessible
   TH1D* HDataBfCorr1D = (TH1D*) HDataBfCorr.ProjectionX();
   TH1D* HDataAfCorr1D = (TH1D*) HDataAfCorr.ProjectionX();
   projection(&HDataBfCorr, HDataBfCorr1D);
   projection(&HDataAfCorr, HDataAfCorr1D);
   Output.cd();
   HDataBfCorr1D->Write();
   HDataAfCorr1D->Write();

   TH1D* HDataBfCorr1D_unfoldBinCorr = (TH1D*) HDataBfCorr1D->Clone(Form("%s_unfoldBinCorr", HDataBfCorr1D->GetName()));
   TH1D* HDataAfCorr1D_unfoldBinCorr = (TH1D*) HDataAfCorr1D->Clone(Form("%s_unfoldBinCorr", HDataAfCorr1D->GetName()));
   EffCorrFactor UnfoldingBinCorrFactor;
   UnfoldingBinCorrFactor.init("../../Unfolding/20240923_UnfoldingBinningCorrection/UnfoldingBinCorr_with_z.root", "z");
   UnfoldingBinCorrFactor.applyEffCorrOnHisto(HDataBfCorr1D, HDataBfCorr1D_unfoldBinCorr);
   UnfoldingBinCorrFactor.applyEffCorrOnHisto(HDataAfCorr1D, HDataAfCorr1D_unfoldBinCorr);
   Output.cd();
   HDataBfCorr1D_unfoldBinCorr->Write();
   HDataAfCorr1D_unfoldBinCorr->Write();

   Output.Close();

}
