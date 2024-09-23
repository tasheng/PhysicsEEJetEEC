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

int FindBin(double Value, int NBins, double Bins[]); 
void MakeCanvasZ(vector<TH1D>& Histograms, vector<string> Labels, string Output, string X, string Y, double WorldMin, double WorldMax, bool DoRatio, bool LogX); 
void SetPad(TPad &P); 
void DivideByBin(TH1D &H, double Bins[]); 

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
   // For tgenBefore comparison
   string InputMCPath            = CL.Get("InputMC", "../../Samples/ALEPHMC/LEP1MC1994_recons_aftercut-001.root");
   string GenBeforeTreeName      = CL.Get("GenBefore", "tgenBefore");

   TFile InputData(InputDataPath.c_str(), "READ");
   TH2D HDataBfCorr( *((TH2D*) InputData.Get(HistoName.c_str())) );
   TH2D HDataAfCorr( *((TH2D*) HDataBfCorr.Clone(Form("%s_afCorr", HistoName.c_str()))) );

   // HDataBfCorr.Print("all");

   int applyEffCorrOnHistoErrorStatus = 0;
   EffCorrFactor matchingEffCorrFactor;
   matchingEffCorrFactor.init("../../Unfolding/20240328_Unfolding/matchingScheme2/MatchingEff.root", "z");
   applyEffCorrOnHistoErrorStatus += matchingEffCorrFactor.applyEffCorrOnHisto(&HDataBfCorr, &HDataAfCorr);

   EffCorrFactor EvtSelEffCorrFactor;
   EvtSelEffCorrFactor.init("../../EventSelectionEfficiency/20240922_evtSelEffCorr/EvtSelEff.root", "z");
   applyEffCorrOnHistoErrorStatus += EvtSelEffCorrFactor.applyEffCorrOnHisto(&HDataBfCorr, &HDataAfCorr);

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
   applyEffCorrOnHistoErrorStatus += UnfoldingBinCorrFactor.applyEffCorrOnHisto(HDataBfCorr1D, HDataBfCorr1D_unfoldBinCorr);
   applyEffCorrOnHistoErrorStatus += UnfoldingBinCorrFactor.applyEffCorrOnHisto(HDataAfCorr1D, HDataAfCorr1D_unfoldBinCorr);
   Output.cd();
   HDataBfCorr1D_unfoldBinCorr->Write();
   HDataAfCorr1D_unfoldBinCorr->Write();

   //------------------------------------
   // plotting
   //------------------------------------
   // theta binning
   const int BinCount = 100;
   double Bins[2*BinCount+1];
   double BinMin = 0.002;
   double BinMax = M_PI / 2;

   // z binning
   double zBins[2*BinCount+1];
   double zBinMin = (1- cos(0.002))/2;
   double zBinMax = 0.5;

   for(int i = 0; i <= BinCount; i++){
      // theta double log binning
      Bins[i] = exp(log(BinMin) + (log(BinMax) - log(BinMin)) / BinCount * i);
      Bins[2*BinCount-i] = BinMax * 2 - exp(log(BinMin) + (log(BinMax) - log(BinMin)) / BinCount * i);

      // z double log binning
      zBins[i] = exp(log(zBinMin) + (log(zBinMax) - log(zBinMin)) / BinCount * i);
      zBins[2*BinCount-i] = zBinMax * 2 - exp(log(zBinMin) + (log(zBinMax) - log(zBinMin)) / BinCount * i);
    
   }

   // [Warning] JC: this is a quick hack to get the normalization, we would need to parse that normalization value
  TString fnamesmeared = "/data/hbossi/PhysicsEEJetEEC/Unfolding/20240922_UnfoldingThetaZ/UnfoldingInputData_09232024.root";
  TFile *inputsmeared =TFile::Open(fnamesmeared);
  TTree *smeared=(TTree*)inputsmeared->Get("UnmatchedPairTree");
  Int_t nEv=smeared->GetEntries();

   HDataBfCorr1D->Scale(1.0/nEv); 
   HDataAfCorr1D->Scale(1.0/nEv);
   HDataBfCorr1D_unfoldBinCorr->Scale(1.0/nEv);
   HDataAfCorr1D_unfoldBinCorr->Scale(1.0/nEv);

   //------------------------------------
   // getting the EEC(z) at the gen level before the event selection
   //------------------------------------
   TFile InputMC(InputMCPath.c_str(), "READ");
   TH1D HzMCGenBeforeRef("HzMCGenBeforeRef", "HzMCGenBeforeRef", 2 * BinCount, 0, 2 * BinCount); 

   double TotalE = 91.1876;
   ParticleTreeMessenger MGenBefore(InputMC, GenBeforeTreeName); 
   int EntryCountBefore = MGenBefore.GetEntries();
   for(int iE = 0; iE < EntryCountBefore; iE++)
   {
      MGenBefore.GetEntry(iE);
      // fill the four vector
      vector<FourVector> PGenBefore;
      for(int i = 0; i < MGenBefore.nParticle; i++){
        // charged particle selection 
       if(MGenBefore.charge[i] == 0) continue;
       if(MGenBefore.highPurity[i] == false) continue;
         PGenBefore.push_back(MGenBefore.P[i]);
      } // end loop over the particles 

      // now calculate and fill the EECs
      for(int i = 0; i < PGenBefore.size(); i++){
        for(int j = i+1; j < PGenBefore.size();j++){
            FourVector Gen1 = PGenBefore.at(i);
            FourVector Gen2 = PGenBefore.at(j);

            // get the proper bins
            int BinThetaGen  = FindBin(GetAngle(Gen1,Gen2), 2 * BinCount, Bins);
            // int BinEnergyGen = FindBin(Gen1[0]*Gen2[0]/(TotalE*TotalE), EnergyBinCount, EnergyBins);
            double zGen = (1-cos(GetAngle(Gen1,Gen2)))/2; 
            int BinZGen = FindBin(zGen, 2*BinCount, zBins); 

            // calculate the EEC
            double EEC =  Gen1[0]*Gen2[0]/(TotalE*TotalE); 
            
            // fill the histograms
            HzMCGenBeforeRef.Fill(BinZGen, EEC); 
         }
      }
   } // end loop over the number of events
   // EEC is per-event so scale by the event number
   HzMCGenBeforeRef.Scale(1.0/EntryCountBefore);

   // divide by the bin width
   DivideByBin(*HDataBfCorr1D, zBins);
   DivideByBin(*HDataAfCorr1D, zBins);
   DivideByBin(*HDataBfCorr1D_unfoldBinCorr, zBins);
   DivideByBin(*HDataAfCorr1D_unfoldBinCorr, zBins);
   DivideByBin(HzMCGenBeforeRef, zBins);

   // set the style for the plots
   HDataBfCorr1D->SetMarkerColor(Colors[5]);
   HDataAfCorr1D->SetMarkerColor(Colors[3]);
   HDataBfCorr1D_unfoldBinCorr->SetMarkerColor(Colors[4]);
   HDataAfCorr1D_unfoldBinCorr->SetMarkerColor(Colors[2]);
   HzMCGenBeforeRef.SetMarkerColor(Colors[0]);
   HDataBfCorr1D->SetLineColor(Colors[5]);
   HDataAfCorr1D->SetLineColor(Colors[3]);
   HDataBfCorr1D_unfoldBinCorr->SetLineColor(Colors[4]);
   HDataAfCorr1D_unfoldBinCorr->SetLineColor(Colors[2]);
   HzMCGenBeforeRef.SetLineColor(Colors[0]);
   HDataBfCorr1D->SetMarkerStyle(20);
   HDataAfCorr1D->SetMarkerStyle(20);
   HDataBfCorr1D_unfoldBinCorr->SetMarkerStyle(20);
   HDataAfCorr1D_unfoldBinCorr->SetMarkerStyle(20);
   HzMCGenBeforeRef.SetMarkerStyle(20);
   HDataBfCorr1D->SetLineWidth(2);
   HDataAfCorr1D->SetLineWidth(2);
   HDataBfCorr1D_unfoldBinCorr->SetLineWidth(2);
   HDataAfCorr1D_unfoldBinCorr->SetLineWidth(2);
   HzMCGenBeforeRef.SetLineWidth(2);

   std::vector<TH1D> hists = {HzMCGenBeforeRef, *HDataBfCorr1D, *HDataAfCorr1D, *HDataBfCorr1D_unfoldBinCorr, *HDataAfCorr1D_unfoldBinCorr}; 

   system("mkdir -p plot/");
   MakeCanvasZ(hists, { "MC Gen before sel.",
                        "Before eff. corr.'s", "After eff. corr.'s", 
                        "Before eff. corr.'s + bin reweighted", "After eff. corr.'s + bin reweighted"}, 
               Form("plot/CorrectedData_z"), "#it{z} = (1- cos(#theta))/2", "#frac{1}{#it{N}_{event}}#frac{d(Sum E_{i}E_{j}/E^{2})}{d#it{z}}",1e-3,1e3, true, true);


   Output.Close();

   if (applyEffCorrOnHistoErrorStatus>0)
   {
      printf("[Error] Something wrong with applyEffCorrOnHisto.\n");
      return 1;
   } else return 0;
}

int FindBin(double Value, int NBins, double Bins[])
{
   for(int i = 0; i < NBins; i++)
      if(Value < Bins[i])
         return i - 1;
   return NBins;
}


void MakeCanvasZ(vector<TH1D>& Histograms, vector<string> Labels, string Output,
   string X, string Y, double WorldMin, double WorldMax, bool DoRatio, bool LogX){
   

   int NLine = Histograms.size();
   int N = Histograms[0].GetNbinsX();

   double MarginL = 180;
   double MarginR = 90;
   double MarginB = 120;
   double MarginT = 90;

   double WorldXMin = LogX ? 0 : 0;
   double WorldXMax = LogX ? N : 1;
   
   double PadWidth = 1200;
   double PadHeight = DoRatio ? 640 : 640 + 240;
   double PadRHeight = DoRatio ? 240 : 0.001;

   double CanvasWidth = MarginL + PadWidth + MarginR;
   double CanvasHeight = MarginT + PadHeight + PadRHeight + MarginB;

   MarginL = MarginL / CanvasWidth;
   MarginR = MarginR / CanvasWidth;
   MarginT = MarginT / CanvasHeight;
   MarginB = MarginB / CanvasHeight;

   PadWidth   = PadWidth / CanvasWidth;
   PadHeight  = PadHeight / CanvasHeight;
   PadRHeight = PadRHeight / CanvasHeight;

   TCanvas Canvas("Canvas", "", CanvasWidth, CanvasHeight);
   Canvas.cd(); 

   TPad Pad("Pad", "", MarginL, MarginB + PadRHeight, MarginL + PadWidth, MarginB + PadHeight + PadRHeight);
   Pad.SetLogy();
   SetPad(Pad);
   
   TPad PadR("PadR", "", MarginL, MarginB, MarginL + PadWidth, MarginB + PadRHeight);
   if(DoRatio)
      SetPad(PadR);

   Pad.cd();

   TH2D HWorld("HWorld", "", N, WorldXMin, WorldXMax, 100, WorldMin, WorldMax);
   HWorld.SetStats(0);
   HWorld.GetXaxis()->SetTickLength(0);
   HWorld.GetXaxis()->SetLabelSize(0);

   HWorld.Draw("axis");
   for(TH1D H : Histograms){
      TH1D *HClone = (TH1D *)H.Clone();
      HClone->Draw("exp same");
   }

   TGraph G;
   G.SetPoint(0, LogX ? N / 2 : 1 / 2, 0);
   G.SetPoint(1, LogX ? N / 2 : 1/ 2, 10000);
   G.SetLineStyle(kDashed);
   G.SetLineColor(kGray);
   G.SetLineWidth(1);
   G.Draw("l");

   if(DoRatio)
      PadR.cd();

   double WorldRMin = 0.5;
   double WorldRMax = 1.5;
   
   TH2D HWorldR("HWorldR", "", N, WorldXMin, WorldXMax, 100, WorldRMin, WorldRMax);
   TGraph G2;
   
   if(DoRatio)
   {
      HWorldR.SetStats(0);
      HWorldR.GetXaxis()->SetTickLength(0);
      HWorldR.GetXaxis()->SetLabelSize(0);
      HWorldR.GetYaxis()->SetNdivisions(505);

      HWorldR.Draw("axis");
      for(int i = 1; i < NLine; i++)
      {
         TH1D *H = (TH1D *)Histograms[i].Clone();
         H->Divide(&Histograms[0]);
         H->Draw("same");
      }

      G.Draw("l");

      G2.SetPoint(0, 0, 1);
      G2.SetPoint(1, 99999, 1);
      G2.Draw("l");
   }
   
   double BinMin    = (1- cos(0.002))/2;
   double BinMiddle = 0.5;
   double BinMax    = 1 - BinMin;

   Canvas.cd();
   std::cout << "Here 1" << std::endl;
   int nDiv = 505;
   std::cout << "Bin Min " << BinMin << " Bin Middle " << BinMiddle << std::endl;
   TGaxis X1(MarginL, MarginB, MarginL + PadWidth / 2, MarginB, BinMin, BinMiddle, nDiv, "GS");
   TGaxis X2(MarginL + PadWidth, MarginB, MarginL + PadWidth / 2, MarginB, BinMin, BinMiddle, nDiv, "-GS");
   TGaxis X3(MarginL, MarginB + PadRHeight, MarginL + PadWidth / 2, MarginB + PadRHeight, BinMin, BinMiddle, nDiv, "+-GS");
   TGaxis X4(MarginL + PadWidth, MarginB + PadRHeight, MarginL + PadWidth / 2, MarginB + PadRHeight, BinMin, BinMiddle, nDiv, "+-GS");

   std::cout << "Here 2" << std::endl;
   TGaxis Y1(MarginL, MarginB, MarginL, MarginB + PadRHeight, WorldRMin, WorldRMax, 505, "");
   TGaxis Y2(MarginL, MarginB + PadRHeight, MarginL, MarginB + PadRHeight + PadHeight, WorldMin, WorldMax, 510, "G");

   std::cout << "Here 3" << std::endl;
   TGaxis XL1(MarginL, MarginB, MarginL + PadWidth, MarginB,  (1- cos(0.002))/2, 1, 210, "S");
   TGaxis XL2(MarginL, MarginB + PadRHeight, MarginL + PadWidth, MarginB + PadRHeight,  (1- cos(0.002))/2, 1, 210, "+-S");

   Y1.SetLabelFont(42);
   Y2.SetLabelFont(42);
   XL1.SetLabelFont(42);
   XL2.SetLabelFont(42);

   X1.SetLabelSize(0);
   X2.SetLabelSize(0);
   X3.SetLabelSize(0);
   X4.SetLabelSize(0);
   // XL1.SetLabelSize(0);
   XL2.SetLabelSize(0);

   X1.SetTickSize(0.06);
   X2.SetTickSize(0.06);
   X3.SetTickSize(0.06);
   X4.SetTickSize(0.06);
   XL1.SetTickSize(0.03);
   XL2.SetTickSize(0.03);

   if(LogX == true)
   {
      X1.Draw();
      X2.Draw();
      if(DoRatio) X3.Draw();
      if(DoRatio) X4.Draw();
   }
   if(LogX == false)
   {
      XL1.Draw();
      if(DoRatio)
         XL2.Draw();
   }
   if(DoRatio)
      Y1.Draw();
   Y2.Draw();

   TLatex Latex;
   Latex.SetNDC();
   Latex.SetTextFont(42);
   Latex.SetTextSize(0.035);
   Latex.SetTextAlign(23);
   if(LogX) Latex.DrawLatex(MarginL + (1 - MarginR - MarginL) * 0.02, MarginB - 0.01, "10^{-6} ");
   if(LogX) Latex.DrawLatex(MarginL + (1 - MarginR - MarginL) * 0.180, MarginB - 0.01, "10^{-4} ");
   if(LogX) Latex.DrawLatex(MarginL + (1 - MarginR - MarginL) * 0.365, MarginB - 0.01, "10^{-2} ");
   if(LogX) Latex.DrawLatex(MarginL + (1 - MarginR - MarginL) * 0.500, MarginB - 0.01, "1/2");
   if(LogX) Latex.DrawLatex(MarginL + (1 - MarginR - MarginL) * 0.650, MarginB - 0.01, "1 - 10^{-2}");
   if(LogX) Latex.DrawLatex(MarginL + (1 - MarginR - MarginL) * 0.823, MarginB - 0.01, "1 - 10^{-4}");
   if(LogX) Latex.DrawLatex(MarginL + (1 - MarginR - MarginL) * 0.995, MarginB - 0.01, "1 - 10^{-6}");

   Latex.SetTextAlign(12);
   Latex.SetTextAngle(270);
   Latex.SetTextColor(kGray);
   Latex.DrawLatex(MarginL + (1 - MarginR - MarginL) * 0.5 + 0.0175, 1 - MarginT - 0.015, "#it{z} = 1/2");

   Latex.SetTextAlign(22);
   Latex.SetTextAngle(0);
   Latex.SetTextColor(kBlack);
   Latex.DrawLatex(MarginL + PadWidth * 0.9, MarginB * 0.4, X.c_str());

   Latex.SetTextAlign(22);
   Latex.SetTextAngle(90);
   Latex.SetTextColor(kBlack);
   if(DoRatio)
      Latex.DrawLatex(MarginL * 0.3, MarginB + PadRHeight * 0.5, "Data / MC Gen before");
   Latex.DrawLatex(MarginL * 0.3, MarginB + PadRHeight + PadHeight * 0.5, Y.c_str());

   Latex.SetTextAlign(11);
   Latex.SetTextAngle(0);
   Latex.DrawLatex(MarginL, MarginB + PadRHeight + PadHeight + 0.012, "ALEPH e^{+}e^{-}, #sqrt{s} = 91.2 GeV, Work-in-progress");

   Latex.SetTextAlign(11);
   Latex.SetTextAngle(0);
   Latex.SetTextColor(19);
   Latex.SetTextSize(0.02);
   Latex.DrawLatex(0.01, 0.01, "Work-in-progress, 2024 August 7th HB");

   TLegend Legend(0.15, 0.90, 0.35, 0.90 - 0.04 * min(NLine, 4));
   Legend.SetTextFont(42);
   Legend.SetTextSize(0.035);
   Legend.SetFillStyle(0);
   Legend.SetBorderSize(0);
   for(int i = 0; i < NLine && i < 4; i++)
      Legend.AddEntry(&Histograms[i], Labels[i].c_str(), "pl");
   Legend.Draw();

   TLegend Legend2(0.55, 0.90, 0.8, 0.90 - 0.04 * (NLine - 4));
   Legend2.SetTextFont(42);
   Legend2.SetTextSize(0.035);
   Legend2.SetFillStyle(0);
   Legend2.SetBorderSize(0);
   if(NLine >= 4)
   {
      for(int i = 4; i < NLine; i++)
         Legend2.AddEntry(&Histograms[i], Labels[i].c_str(), "pl");
      Legend2.Draw();
   }

   Canvas.SaveAs((Output + ".pdf").c_str());
}

void SetPad(TPad &P){
   P.SetLeftMargin(0);
   P.SetTopMargin(0);
   P.SetRightMargin(0);
   P.SetBottomMargin(0);
   P.SetTickx();
   P.SetTicky();
   P.Draw();
}



void DivideByBin(TH1D &H, double Bins[])
{
   int N = H.GetNbinsX();
   for(int i = 1; i <= N; i++)
   {
      double L = Bins[i-1];
      double R = Bins[i];
      H.SetBinContent(i, H.GetBinContent(i) / (R - L));
      H.SetBinError(i, H.GetBinError(i) / (R - L));
   }
}
