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
#include "TGraphAsymmErrors.h"
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
void MakeCanvasZ(vector<TH1D > Histograms, TGraphErrors DataSyst, vector<string> Labels, string Output, string X, string Y, double WorldMin, double WorldMax, bool DoRatio, bool LogX);
void SetPad(TPad &P); 
void DivideByBin(TH1D &H, double Bins[]); 
TGraphAsymmErrors getTheoryPlot();
TGraphAsymmErrors getOPALPlot();  
void MakeCanvasZTheory(vector<TH1D > Histograms, TGraphErrors DataSyst,  TGraphAsymmErrors TheorySyst, vector<string> Labels, string Output, string X, string Y, double WorldMin, double WorldMax, bool DoRatio, bool LogX); 


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
   matchingEffCorrFactor.init("/data/janicechen/PhysicsEEJetEEC/Unfolding/20240328_Unfolding/matchingScheme2/MatchingEff.root", "z");
   applyEffCorrOnHistoErrorStatus += matchingEffCorrFactor.applyEffCorrOnHisto(&HDataBfCorr, &HDataAfCorr);

   EffCorrFactor EvtSelEffCorrFactor;
   EvtSelEffCorrFactor.init("/data/janicechen/PhysicsEEJetEEC/EventSelectionEfficiency/20240922_evtSelEffCorr/EvtSelEff.root", "z");
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
   UnfoldingBinCorrFactor.init("/data/janicechen/PhysicsEEJetEEC/Unfolding/20240923_UnfoldingBinningCorrection/UnfoldingBinCorr_with_z.root", "z");
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

   std::cout << "zBins 17: " << zBins[17] << " zbins 183: " << zBins[183] << std::endl;

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

   TFile* sysFile = TFile::Open("SystematicsE2C_20240923_191800.root");
   TH2D* h2SysFromFile = (TH2D*)sysFile->Get("Systematics_Z_Total"); 
   TH2D* h2Sys = (TH2D*) HDataBfCorr.Clone("h2Sys");
   TH2D* h2SysPlus =(TH2D*) HDataBfCorr.Clone("h2SysPlus");  


   h2Sys->Add(h2SysFromFile, -1); 
   h2SysPlus->Add(h2SysFromFile); 

   TH2D HSysAfCorr( *((TH2D*) h2Sys->Clone("h2Sys_afCorr") ));
   TH2D HSysAfCorrPlus( *((TH2D*) h2SysPlus->Clone("h2SysPlus_afCorr") ));

   matchingEffCorrFactor.applyEffCorrOnHisto(h2Sys, &HSysAfCorr);   
   matchingEffCorrFactor.applyEffCorrOnHisto(h2SysPlus, &HSysAfCorrPlus);

   EvtSelEffCorrFactor.applyEffCorrOnHisto(h2Sys, &HSysAfCorr);
   EvtSelEffCorrFactor.applyEffCorrOnHisto(h2SysPlus, &HSysAfCorrPlus);


   // JC: the projection is copied here just so a quick check of the projected 1D histogram is accessible
   TH1D* HSysBfCorr1D = (TH1D*) h2Sys->ProjectionX();
   TH1D* HSysBfCorr1DPlus = (TH1D*) h2SysPlus->ProjectionX();
   TH1D* HSysAfCorr1D = (TH1D*) HSysAfCorr.ProjectionX();
   TH1D* HSysAfCorr1DPlus = (TH1D*) HSysAfCorrPlus.ProjectionX();

   projection(h2Sys, HSysBfCorr1D);
   projection(h2SysPlus, HSysBfCorr1DPlus);

   projection(&HSysAfCorr, HSysAfCorr1D);
   projection(&HSysAfCorrPlus, HSysAfCorr1DPlus);


   TH1D* HSysBfCorr1D_unfoldBinCorr = (TH1D*) HSysBfCorr1D->Clone(Form("%s_unfoldBinCorr", HSysBfCorr1D->GetName()));
   TH1D* HSysAfCorr1D_unfoldBinCorr = (TH1D*) HSysAfCorr1D->Clone(Form("%s_unfoldBinCorr", HSysAfCorr1D->GetName()));
   TH1D* HSysBfCorr1DPlus_unfoldBinCorr = (TH1D*) HSysBfCorr1DPlus->Clone(Form("%s_unfoldBinCorr", HSysBfCorr1DPlus->GetName()));
   TH1D* HSysAfCorr1DPlus_unfoldBinCorr = (TH1D*) HSysAfCorr1DPlus->Clone(Form("%s_unfoldBinCorr", HSysAfCorr1DPlus->GetName()));
   UnfoldingBinCorrFactor.applyEffCorrOnHisto(HSysBfCorr1D, HSysBfCorr1D_unfoldBinCorr);
   UnfoldingBinCorrFactor.applyEffCorrOnHisto(HSysAfCorr1D, HSysAfCorr1D_unfoldBinCorr);
   UnfoldingBinCorrFactor.applyEffCorrOnHisto(HSysBfCorr1DPlus, HSysBfCorr1DPlus_unfoldBinCorr);
   UnfoldingBinCorrFactor.applyEffCorrOnHisto(HSysAfCorr1DPlus, HSysAfCorr1DPlus_unfoldBinCorr);
   HSysAfCorr1D_unfoldBinCorr->Scale(1.0/nEv);
   DivideByBin(*HSysAfCorr1D_unfoldBinCorr, zBins);
   HSysAfCorr1DPlus_unfoldBinCorr->Scale(1.0/nEv);
   DivideByBin(*HSysAfCorr1DPlus_unfoldBinCorr, zBins);


   TGraphErrors GzDataSyst( HDataAfCorr1D_unfoldBinCorr ); GzDataSyst.SetName("HzDataSyst");

   for(int i = 1; i <= HDataAfCorr1D_unfoldBinCorr->GetNbinsX(); i++)
   {
      int iGraph = i-1;
      double err = HSysAfCorr1DPlus_unfoldBinCorr->GetBinContent(i) - HDataAfCorr1D_unfoldBinCorr->GetBinContent(i); 
      GzDataSyst.SetPoint(iGraph, HDataAfCorr1D_unfoldBinCorr->GetBinCenter(i), HDataAfCorr1D_unfoldBinCorr->GetBinContent(i));
      GzDataSyst.SetPointError(iGraph, HDataAfCorr1D_unfoldBinCorr->GetBinWidth(i)/2, err);
   }
   GzDataSyst.SetLineWidth(0);
   GzDataSyst.SetFillStyle(1001);
   GzDataSyst.SetFillColorAlpha(Colors[2], 0.3);

   std::vector<TH1D> hists = {HzMCGenBeforeRef, *HDataAfCorr1D_unfoldBinCorr}; 
   std::vector<TH1D> histsTheory = {*HDataAfCorr1D_unfoldBinCorr}; 


   system("mkdir -p plot/");
   MakeCanvasZ(hists,GzDataSyst, { "Archived MC","Fully Corrected Data"}, Form("plot/CorrectedData_z"), "#it{z} = (1- cos(#theta))/2", "#frac{1}{#it{N}_{event}}#frac{d(Sum E_{i}E_{j}/E^{2})}{d#it{z}}",1e-3,1e3, true, true);
   TGraphAsymmErrors theory = getOPALPlot(); 
   MakeCanvasZTheory(histsTheory,GzDataSyst, theory, {"Fully Corrected Data"}, Form("plot/CorrectedData_z_opal"), "#it{z} = (1- cos(#theta))/2", "#frac{1}{#it{N}_{event}}#frac{d(Sum E_{i}E_{j}/E^{2})}{d#it{z}}",1e-3,1e3, false, true);


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

/*
* Make the canvas for the results as a function of z. Note that these results include a cutoff at 1e-6. 
* Additionally, this is intended for MC comparisons. 
*/
void MakeCanvasZ(vector<TH1D > Histograms, TGraphErrors DataSyst, vector<string> Labels, string Output, string X, string Y, double WorldMin, double WorldMax, bool DoRatio, bool LogX)
{
   int NLine = Histograms.size();
   int N = Histograms[0].GetNbinsX();

   double MarginL = 180;
   double MarginR = 90;
   double MarginB = 120;
   double MarginT = 90;

   double WorldXMin = LogX ? 17 : 0;
   double WorldXMax = LogX ? 183: 1;
   
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

   const int BinCount = 100;


   // z binning
   double zBins[2*BinCount+1];
   double zBinMin = (1- cos(0.002))/2;
   double zBinMax = 0.5;

   for(int i = 0; i <= BinCount; i++){
      // z double log binning
      zBins[i] = exp(log(zBinMin) + (log(zBinMax) - log(zBinMin)) / BinCount * i);
      zBins[2*BinCount-i] = zBinMax * 2 - exp(log(zBinMin) + (log(zBinMax) - log(zBinMin)) / BinCount * i);
   }

   TH2D HWorld("HWorld", "", N, WorldXMin, WorldXMax, 100, WorldMin, WorldMax);
   HWorld.SetStats(0);
   HWorld.GetXaxis()->SetTickLength(0);
   HWorld.GetXaxis()->SetLabelSize(0);

   HWorld.Draw("axis");
   for(TH1D H : Histograms){
      TH1D *HClone = (TH1D *)H.Clone();
      HClone->Draw("exp same");
   }
   DataSyst.DrawClone("2 same");


   TGraph G;
   G.SetPoint(0, LogX ? N / 2 : 1 / 2, 0);
   G.SetPoint(1, LogX ? N / 2 : 1/ 2, 10000);
   G.SetLineStyle(kDashed);
   G.SetLineColor(kGray);
   G.SetLineWidth(1);
   G.Draw("l");

   if(DoRatio)
      PadR.cd();

   double WorldRMin = 0;
   double WorldRMax = 1.9999;
   
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

      for(int i = 1; i <= Histograms[0].GetNbinsX(); i++)
      {
         int iGraph = i-1;
         DataSyst.SetPoint(iGraph, 
                                 DataSyst.GetPointX(iGraph), 
                                 DataSyst.GetPointY(iGraph)/Histograms[0].GetBinContent(i));
         DataSyst.SetPointError( iGraph, 
                                 DataSyst.GetErrorX(iGraph), 
                                 DataSyst.GetErrorY(iGraph)/Histograms[0].GetBinContent(i));
      }
      DataSyst.DrawClone("2 same");

      G.Draw("l");

      G2.SetPoint(0, 0, 1);
      G2.SetPoint(1, 99999, 1);
      G2.Draw("l");
   }
   
   double BinMin    = zBins[17];//(1- cos(0.002))/2;
   double BinMiddle = 0.5;
   double BinMax    = 1 - BinMin;

   Canvas.cd();
   int nDiv = 505;
   TGaxis X1(MarginL, MarginB, MarginL + PadWidth / 2, MarginB, BinMin, BinMiddle, nDiv, "GS");
   TGaxis X2(MarginL + PadWidth, MarginB, MarginL + PadWidth / 2, MarginB, BinMin, BinMiddle, nDiv, "-GS");
   TGaxis X3(MarginL, MarginB + PadRHeight, MarginL + PadWidth / 2, MarginB + PadRHeight, BinMin, BinMiddle, nDiv, "+-GS");
   TGaxis X4(MarginL + PadWidth, MarginB + PadRHeight, MarginL + PadWidth / 2, MarginB + PadRHeight, BinMin, BinMiddle, nDiv, "+-GS");

   TGaxis Y1(MarginL, MarginB, MarginL, MarginB + PadRHeight, WorldRMin, WorldRMax, 505, "");
   TGaxis Y2(MarginL, MarginB + PadRHeight, MarginL, MarginB + PadRHeight + PadHeight, WorldMin, WorldMax, 510, "G");

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
   // if(LogX) Latex.DrawLatex(MarginL + (1 - MarginR - MarginL) * 0.02, MarginB - 0.01, "10^{-6} ");
   if(LogX) Latex.DrawLatex(MarginL + (1 - MarginR - MarginL) * 0.1, MarginB - 0.01, "10^{-4} ");
   if(LogX) Latex.DrawLatex(MarginL + (1 - MarginR - MarginL) * 0.32, MarginB - 0.01, "10^{-2} ");
   if(LogX) Latex.DrawLatex(MarginL + (1 - MarginR - MarginL) * 0.500, MarginB - 0.01, "1/2");
   if(LogX) Latex.DrawLatex(MarginL + (1 - MarginR - MarginL) * 0.68, MarginB - 0.01, "1 - 10^{-2}");
   if(LogX) Latex.DrawLatex(MarginL + (1 - MarginR - MarginL) * 0.9, MarginB - 0.01, "1 - 10^{-4}");
   // if(LogX) Latex.DrawLatex(MarginL + (1 - MarginR - MarginL) * 0.995, MarginB - 0.01, "1 - 10^{-6}");

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
      Latex.DrawLatex(MarginL * 0.3, MarginB + PadRHeight * 0.5, "Data/MC");
   Latex.DrawLatex(MarginL * 0.3, MarginB + PadRHeight + PadHeight * 0.5, Y.c_str());

   Latex.SetTextAlign(11);
   Latex.SetTextAngle(0);
   Latex.DrawLatex(MarginL, MarginB + PadRHeight + PadHeight + 0.012, "ALEPH e^{+}e^{-}, #sqrt{s} = 91.2 GeV, Preliminary");

   Latex.SetTextAlign(11);
   Latex.SetTextAngle(0);
   Latex.SetTextColor(19);
   Latex.SetTextSize(0.02);
   Latex.DrawLatex(0.01, 0.01, "2024 Hard Probes Preliminary");

   TLegend Legend(0.15, 0.90, 0.35, 0.90 - 0.04 * min(NLine, 4));
   Legend.SetTextFont(42);
   Legend.SetTextSize(0.035);
   Legend.SetFillStyle(0);
   Legend.SetBorderSize(0);
   for(int i = 0; i < NLine && i < 4; i++)
   {
      if (Labels[i]=="Data")
      {
         Histograms[i].SetFillStyle(DataSyst.GetFillStyle());
         Histograms[i].SetFillColor(DataSyst.GetFillColor());
      }
      Legend.AddEntry(&Histograms[i], Labels[i].c_str(), 
                      (Labels[i]=="Data")? "plf": "pl");
   }
   Legend.Draw();

   TLegend Legend2(0.55, 0.90, 0.8, 0.90 - 0.04 * (NLine - 4));
   Legend2.SetTextFont(42);
   Legend2.SetTextSize(0.035);
   Legend2.SetFillStyle(0);
   Legend2.SetBorderSize(0);
   if(NLine >= 4)
   {
      for(int i = 4; i < NLine; i++)
         Legend2.AddEntry(&Histograms[i], Labels[i].c_str(), 
                      (Labels[i]=="Data")? "plf": "pl");
      Legend2.Draw();
   }

   Canvas.SaveAs((Output + ".pdf").c_str());
}

// make the theory comparison canvas
void MakeCanvasZTheory(vector<TH1D > Histograms, TGraphErrors DataSyst,  TGraphAsymmErrors TheorySyst, vector<string> Labels, string Output, string X, string Y, double WorldMin, double WorldMax, bool DoRatio, bool LogX)
{
   int NLine = Histograms.size();
   int N = Histograms[0].GetNbinsX();

   double MarginL = 180;
   double MarginR = 90;
   double MarginB = 120;
   double MarginT = 90;

   double WorldXMin = LogX ? 17 : 0;
   double WorldXMax = LogX ? 183: 1;
   
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

   const int BinCount = 100;


   // z binning
   double zBins[2*BinCount+1];
   double zBinMin = (1- cos(0.002))/2;
   double zBinMax = 0.5;

   for(int i = 0; i <= BinCount; i++){
      // z double log binning
      zBins[i] = exp(log(zBinMin) + (log(zBinMax) - log(zBinMin)) / BinCount * i);
      zBins[2*BinCount-i] = zBinMax * 2 - exp(log(zBinMin) + (log(zBinMax) - log(zBinMin)) / BinCount * i);
   }

   TH2D HWorld("HWorld", "", N, WorldXMin, WorldXMax, 100, WorldMin, WorldMax);
   HWorld.SetStats(0);
   HWorld.GetXaxis()->SetTickLength(0);
   HWorld.GetXaxis()->SetLabelSize(0);

   HWorld.Draw("axis");
   TheorySyst.DrawClone("3 p same");

   for(TH1D H : Histograms){
      TH1D *HClone = (TH1D *)H.Clone();
      HClone->Draw("exp same");
   }
   DataSyst.DrawClone("2 same");



   TGraph G;
   G.SetPoint(0, LogX ? N / 2 : 1 / 2, 0);
   G.SetPoint(1, LogX ? N / 2 : 1/ 2, 10000);
   G.SetLineStyle(kDashed);
   G.SetLineColor(kGray);
   G.SetLineWidth(1);
   G.Draw("l");

   if(DoRatio)
      PadR.cd();

   double WorldRMin = 0;
   double WorldRMax = 1.9999;
   
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

      for(int i = 1; i <= Histograms[0].GetNbinsX(); i++)
      {
         int iGraph = i-1;
         TheorySyst.SetPoint(iGraph, 
                                 TheorySyst.GetPointX(iGraph), 
                                 TheorySyst.GetPointY(iGraph)/Histograms[0].GetBinContent(i));
         TheorySyst.SetPointError( iGraph, 
                                 TheorySyst.GetErrorXlow(iGraph), TheorySyst.GetErrorXhigh(iGraph),
                                 TheorySyst.GetErrorYlow(iGraph)/Histograms[0].GetBinContent(i), TheorySyst.GetErrorYhigh(iGraph)/Histograms[0].GetBinContent(i));
      }
      TheorySyst.DrawClone("2 same");

      G.Draw("l");

      G2.SetPoint(0, 0, 1);
      G2.SetPoint(1, 99999, 1);
      G2.Draw("l");
   }
   
   double BinMin    = zBins[17];//(1- cos(0.002))/2;
   double BinMiddle = 0.5;
   double BinMax    = 1 - BinMin;

   Canvas.cd();
   int nDiv = 505;
   TGaxis X1(MarginL, MarginB, MarginL + PadWidth / 2, MarginB, BinMin, BinMiddle, nDiv, "GS");
   TGaxis X2(MarginL + PadWidth, MarginB, MarginL + PadWidth / 2, MarginB, BinMin, BinMiddle, nDiv, "-GS");
   TGaxis X3(MarginL, MarginB + PadRHeight, MarginL + PadWidth / 2, MarginB + PadRHeight, BinMin, BinMiddle, nDiv, "+-GS");
   TGaxis X4(MarginL + PadWidth, MarginB + PadRHeight, MarginL + PadWidth / 2, MarginB + PadRHeight, BinMin, BinMiddle, nDiv, "+-GS");

   TGaxis Y1(MarginL, MarginB, MarginL, MarginB + PadRHeight, WorldRMin, WorldRMax, 505, "");
   TGaxis Y2(MarginL, MarginB + PadRHeight, MarginL, MarginB + PadRHeight + PadHeight, WorldMin, WorldMax, 510, "G");

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
   // if(LogX) Latex.DrawLatex(MarginL + (1 - MarginR - MarginL) * 0.02, MarginB - 0.01, "10^{-6} ");
   if(LogX) Latex.DrawLatex(MarginL + (1 - MarginR - MarginL) * 0.1, MarginB - 0.01, "10^{-4} ");
   if(LogX) Latex.DrawLatex(MarginL + (1 - MarginR - MarginL) * 0.32, MarginB - 0.01, "10^{-2} ");
   if(LogX) Latex.DrawLatex(MarginL + (1 - MarginR - MarginL) * 0.500, MarginB - 0.01, "1/2");
   if(LogX) Latex.DrawLatex(MarginL + (1 - MarginR - MarginL) * 0.68, MarginB - 0.01, "1 - 10^{-2}");
   if(LogX) Latex.DrawLatex(MarginL + (1 - MarginR - MarginL) * 0.9, MarginB - 0.01, "1 - 10^{-4}");
   // if(LogX) Latex.DrawLatex(MarginL + (1 - MarginR - MarginL) * 0.995, MarginB - 0.01, "1 - 10^{-6}");

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
      Latex.DrawLatex(MarginL * 0.3, MarginB + PadRHeight * 0.5, "Data/MC");
   Latex.DrawLatex(MarginL * 0.3, MarginB + PadRHeight + PadHeight * 0.5, Y.c_str());

   Latex.SetTextAlign(11);
   Latex.SetTextAngle(0);
   Latex.DrawLatex(MarginL, MarginB + PadRHeight + PadHeight + 0.012, "ALEPH e^{+}e^{-}, #sqrt{s} = 91.2 GeV, Preliminary");

   Latex.SetTextAlign(11);
   Latex.SetTextAngle(0);
   Latex.SetTextColor(19);
   Latex.SetTextSize(0.02);
   Latex.DrawLatex(0.01, 0.01, "2024 Hard Probes Preliminary");

   TLegend Legend(0.15, 0.90, 0.35, 0.90 - 0.04 * min(3, 4));
   Legend.SetTextFont(42);
   Legend.SetTextSize(0.03);
   Legend.SetFillStyle(0);
   Legend.SetBorderSize(0);
   for(int i = 0; i < NLine && i < 4; i++)
   {
      if (Labels[i]=="Data")
      {
         Histograms[i].SetFillStyle(DataSyst.GetFillStyle());
         Histograms[i].SetFillColor(DataSyst.GetFillColor());
      }
      Legend.AddEntry(&Histograms[i], Labels[i].c_str(), 
                      (Labels[i]=="Data")? "plf": "pl");
   }
   Legend.AddEntry(&TheorySyst, "OPAL Data", "pl");
   // Legend.AddEntry((TObject*)0, "(NNLL Collinear + NNNLL Sudakov)", ""); 
   Legend.Draw();

   TLegend Legend2(0.55, 0.90, 0.8, 0.90 - 0.04 * (NLine - 4));
   Legend2.SetTextFont(42);
   Legend2.SetTextSize(0.03);
   Legend2.SetFillStyle(0);
   Legend2.SetBorderSize(0);
   if(NLine >= 4)
   {
      for(int i = 4; i < NLine; i++)
         Legend2.AddEntry(&Histograms[i], Labels[i].c_str(), 
                      (Labels[i]=="Data")? "plf": "pl");
      Legend2.Draw();
   }

   Canvas.SaveAs((Output + ".pdf").c_str());
}

// create the settings for the tpad 
void SetPad(TPad &P){
   P.SetLeftMargin(0);
   P.SetTopMargin(0);
   P.SetRightMargin(0);
   P.SetBottomMargin(0);
   P.SetTickx();
   P.SetTicky();
   P.Draw();
}


// normalize by the bin width
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

// takes in theory points from Ian and returns TGraphAsymmErrors of the points
TGraphAsymmErrors getTheoryPlot(){
      // points for the central values of the theory curves
      double centralvals[200][2] =  {{1.14022*10e-6, 10.6999}, {1.30011*10e-6, 10.6999}, {1.48241*10e-6,
      10.6999}, {1.69028*10e-6, 10.6998}, {1.9273*10e-6,
      10.6998}, {2.19755*10e-6, 10.6997}, {2.50569*10e-6,
      10.6997}, {2.85705*10e-6, 10.6996}, {3.25767*10e-6,
      10.6995}, {3.71447*10e-6, 10.6994}, {4.23532*10e-6,
      10.6993}, {4.82921*10e-6, 10.6992}, {5.50638*10e-6,
      10.699}, {6.2785*10e-6, 10.6988}, {7.15888*10e-6,
      10.6985}, {8.16272*10e-6, 10.6982}, {9.30732*10e-6,
      10.6978}, {0.0000106124, 10.6973}, {0.0000121005,
      10.6967}, {0.0000137973, 10.6959}, {0.000015732,
      10.6951}, {0.000017938, 10.694}, {0.0000204533,
      10.6927}, {0.0000233213, 10.6911}, {0.0000265915,
      10.6892}, {0.0000303202, 10.6868}, {0.0000345718,
      10.684}, {0.0000394195, 10.6805}, {0.0000449471,
      10.6762}, {0.0000512496, 10.6711}, {0.000058436,
      10.6648}, {0.0000666301, 10.6572}, {0.0000759731,
      10.6479}, {0.0000866263, 10.6366}, {0.0000987733,
      10.6229}, {0.000112624, 10.6063}, {0.000128416,
      10.5862}, {0.000146423, 10.5617}, {0.000166955,
      10.5321}, {0.000190365, 10.4963}, {0.000217059,
      10.453}, {0.000247496, 10.4007}, {0.0002822, 10.3378}, {0.000321771,
         10.2623}, {0.00036689, 10.1717}, {0.000418337,
      10.0637}, {0.000476997, 9.93512}, {0.000543883,
      9.78298}, {0.000620148, 9.60392}, {0.000707107,
      9.39457}, {0.000806259, 9.15168}, {0.000919315,
      8.87238}, {0.00104822, 8.55452}, {0.00119521, 8.19696}, {0.0013628,
      7.80003}, {0.0015539, 7.36573}, {0.00177179, 6.89812}, {0.00202024,
      6.40316}, {0.00230352, 5.88873}, {0.00262653, 5.36402}, {0.00299483,
         4.83905}, {0.00341477, 4.32383}, {0.0038936, 3.82763}, {0.00443957,
         3.35839}, {0.0050621, 2.92221}, {0.00577192, 2.52321}, {0.00658127,
         2.16354}, {0.00750412, 1.84357}, {0.00855636,
      1.56227}, {0.00975616, 1.3175}, {0.0111242, 1.10644}, {0.0126841,
      0.925849}, {0.0144627, 0.772364}, {0.0164906, 0.641358}, {0.018803,
      0.531175}, {0.0214396, 0.448236}, {0.0244459, 0.386252}, {0.0278738,
         0.331888}, {0.0317824, 0.28539}, {0.036239, 0.245715}, {0.0413205,
      0.211811}, {0.0471146, 0.182846}, {0.0537211, 0.158078}, {0.061254,
      0.136891}, {0.0698433, 0.118765}, {0.0796369, 0.103247}, {0.0908038,
         0.0899689}, {0.103537, 0.0783422}, {0.118055,
      0.0681593}, {0.134609, 0.0594393}, {0.153484, 0.0521564}, {0.175006,
         0.0462156}, {0.199546, 0.0413229}, {0.227526,
      0.0372613}, {0.259431, 0.033881}, {0.295809, 0.0311329}, {0.337288,
      0.0289907}, {0.384583, 0.0274659}, {0.438511, 0.0266299}, {0.5,
      0.0266588}, {0.561489, 0.0266287}, {0.615417, 0.0285988}, {0.662712,
         0.0314269}, {0.704191, 0.0350517}, {0.740569,
      0.0395662}, {0.772474, 0.0450902}, {0.800454, 0.0517807}, {0.824994,
         0.0598343}, {0.846516, 0.0694895}, {0.865391,
      0.0810323}, {0.881945, 0.0948032}, {0.896463, 0.111205}, {0.909196,
      0.130712}, {0.920363, 0.153883}, {0.930157, 0.181372}, {0.938746,
      0.213943}, {0.946279, 0.252494}, {0.952885, 0.298203}, {0.958679,
      0.351671}, {0.963761, 0.409324}, {0.968218, 0.477187}, {0.972126,
      0.554457}, {0.975554, 0.642385}, {0.97856, 0.742408}, {0.981197,
      0.855938}, {0.983509, 0.984361}, {0.985537, 1.12917}, {0.987316,
      1.29184}, {0.988876, 1.47371}, {0.990244, 1.67615}, {0.991444,
      1.90046}, {0.992496, 2.14754}, {0.993419, 2.41836}, {0.994228,
      2.71303}, {0.994938, 3.03188}, {0.99556, 3.3737}, {0.996106,
      3.73788}, {0.996585, 4.12231}, {0.997005, 4.5241}, {0.997373,
      4.93954}, {0.997696, 5.36533}, {0.99798, 5.7979}, {0.998228,
      6.22944}, {0.998446, 6.65777}, {0.998637, 7.07661}, {0.998805,
      7.48292}, {0.998952, 7.87025}, {0.999081, 8.23594}, {0.999194,
      8.57615}, {0.999293, 8.88846}, {0.99938, 9.17203}, {0.999456,
      9.42418}, {0.999523, 9.64688}, {0.999582, 9.83995}, {0.999633,
      10.0011}, {0.999678, 10.1357}, {0.999718, 10.2459}, {0.999753,
      10.3323}, {0.999783, 10.3962}, {0.99981, 10.4435}, {0.999833,
      10.4743}, {0.999854, 10.4931}, {0.999872, 10.501}, {0.999887,
      10.5007}, {0.999901, 10.4941}, {0.999913, 10.4829}, {0.999924,
      10.4675}, {0.999933, 10.451}, {0.999942, 10.4304}, {0.999949,
      10.4113}, {0.999955, 10.3926}, {0.999961, 10.3713}, {0.999965,
      10.3556}, {0.99997, 10.3341}, {0.999973, 10.32}, {0.999977,
      10.2996}, {0.99998, 10.2829}, {0.999982, 10.2711}, {0.999984,
      10.2585}, {0.999986, 10.2451}, {0.999988, 10.2308}, {0.999989,
      10.2232}, {0.999991, 10.2071}, {0.999992, 10.1985}, {0.999993,
      10.1895}, {0.999994, 10.1801}, {0.999994, 10.1801}, {0.999995,
      10.1703}, {0.999996, 10.1603}, {0.999996, 10.1603}, {0.999997,
      10.1508}, {0.999997, 10.1508}, {0.999997, 10.1508}, {0.999998,
      10.1443}, {0.999998, 10.1443}, {0.999998, 10.1443}, {0.999999,
      10.1568}, {0.999999, 10.1568}, {0.999999, 10.1568}};

      // points for the higher values of the error bands
      double err_high[200][2] = {{1.14022*10e-6, 15.6}, {1.30011*10e-6, 15.6}, {1.48241*10e-6,
      15.6}, {1.69028*10e-6, 15.6}, {1.9273*10e-6, 15.6}, {2.19755*10e-6,
      15.5999}, {2.50569*10e-6, 15.5999}, {2.85705*10e-6,
      15.5999}, {3.25767*10e-6, 15.5999}, {3.71447*10e-6,
      15.5999}, {4.23532*10e-6, 15.5998}, {4.82921*10e-6,
      15.5998}, {5.50638*10e-6, 15.5997}, {6.2785*10e-6,
      15.5997}, {7.15888*10e-6, 15.5996}, {8.16272*10e-6,
      15.5995}, {9.30732*10e-6, 15.5993}, {0.0000106124,
      15.5992}, {0.0000121005, 15.599}, {0.0000137973,
      15.5987}, {0.000015732, 15.5984}, {0.000017938,
      15.5979}, {0.0000204533, 15.5974}, {0.0000233213,
      15.5967}, {0.0000265915, 15.5959}, {0.0000303202,
      15.5948}, {0.0000345718, 15.5935}, {0.0000394195,
      15.5918}, {0.0000449471, 15.5897}, {0.0000512496,
      15.587}, {0.000058436, 15.5837}, {0.0000666301,
      15.5795}, {0.0000759731, 15.5742}, {0.0000866263,
      15.5675}, {0.0000987733, 15.5592}, {0.000112624,
      15.5487}, {0.000128416, 15.5355}, {0.000146423,
      15.5189}, {0.000166955, 15.4981}, {0.000190365,
      15.472}, {0.000217059, 15.4393}, {0.000247496, 15.3983}, {0.0002822,
         15.3471}, {0.000321771, 15.2831}, {0.00036689,
      15.2034}, {0.000418337, 15.1043}, {0.000476997,
      14.9814}, {0.000543883, 14.8296}, {0.000620148,
      14.643}, {0.000707107, 14.4147}, {0.000806259,
      14.1374}, {0.000919315, 13.8033}, {0.00104822,
      13.4047}, {0.00119521, 12.9348}, {0.0013628, 12.3885}, {0.0015539,
      11.7633}, {0.00177179, 11.061}, {0.00202024, 10.2883}, {0.00230352,
      9.45702}, {0.00262653, 8.58438}, {0.00299483, 7.69144}, {0.00341477,
         6.80135}, {0.0038936, 5.93693}, {0.00443957, 5.11847}, {0.0050621,
      4.3619}, {0.00577192, 3.67792}, {0.00658127, 3.07189}, {0.00750412,
      2.54439}, {0.00855636, 2.09235}, {0.00975616, 1.7101}, {0.0111242,
      1.3905}, {0.0126841, 1.12578}, {0.0144627, 0.908237}, {0.0164906,
      0.727761}, {0.018803, 0.581304}, {0.0214396, 0.482044}, {0.0244459,
      0.415775}, {0.0278738, 0.356327}, {0.0317824, 0.305523}, {0.036239,
      0.262287}, {0.0413205, 0.225431}, {0.0471146, 0.194018}, {0.0537211,
         0.167316}, {0.061254, 0.14454}, {0.0698433, 0.125099}, {0.0796369,
      0.108495}, {0.0908038, 0.094326}, {0.103537, 0.0819665}, {0.118055,
      0.0711902}, {0.134609, 0.061999}, {0.153484, 0.0543494}, {0.175006,
      0.0481302}, {0.199546, 0.0430273}, {0.227526, 0.0388008}, {0.259431,
         0.0352873}, {0.295809, 0.0324343}, {0.337288, 0.030217}, {0.384583,
         0.0286472}, {0.438511, 0.0278006}, {0.5, 0.0278633}, {0.561489,
      0.0281555}, {0.615417, 0.0302384}, {0.662712, 0.0332074}, {0.704191,
         0.0370211}, {0.740569, 0.0417724}, {0.772474,
      0.0475863}, {0.800454, 0.0546263}, {0.824994, 0.0630977}, {0.846516,
         0.0732486}, {0.865391, 0.0853764}, {0.881945,
      0.0998347}, {0.896463, 0.11704}, {0.909196, 0.137484}, {0.920363,
      0.161742}, {0.930157, 0.190486}, {0.938746, 0.2245}, {0.946279,
      0.264701}, {0.952885, 0.312421}, {0.958679, 0.367707}, {0.963761,
      0.425388}, {0.968218, 0.494427}, {0.972126, 0.573641}, {0.975554,
      0.664938}, {0.97856, 0.772384}, {0.981197, 0.892013}, {0.983509,
      1.02653}, {0.985537, 1.17852}, {0.987316, 1.3496}, {0.988876,
      1.54126}, {0.990244, 1.75476}, {0.991444, 1.99078}, {0.992496,
      2.24881}, {0.993419, 2.52751}, {0.994228, 2.82376}, {0.994938,
      3.13375}, {0.99556, 3.45088}, {0.996106, 3.78308}, {0.996585,
      4.19338}, {0.997005, 4.68023}, {0.997373, 5.23372}, {0.997696,
      5.83477}, {0.99798, 6.47909}, {0.998228, 7.15381}, {0.998446,
      7.85235}, {0.998637, 8.56014}, {0.998805, 9.26649}, {0.998952,
      9.95397}, {0.999081, 10.6115}, {0.999194, 11.2265}, {0.999293,
      11.7897}, {0.99938, 12.296}, {0.999456, 12.7385}, {0.999523,
      13.1198}, {0.999582, 13.4397}, {0.999633, 13.6958}, {0.999678,
      13.8984}, {0.999718, 14.0534}, {0.999753, 14.1636}, {0.999783,
      14.2344}, {0.99981, 14.2754}, {0.999833, 14.2905}, {0.999854,
      14.2863}, {0.999872, 14.2675}, {0.999887, 14.2406}, {0.999901,
      14.206}, {0.999913, 14.1691}, {0.999924, 14.1248}, {0.999933,
      14.1202}, {0.999942, 14.1818}, {0.999949, 14.2509}, {0.999955,
      14.3081}, {0.999961, 14.3652}, {0.999965, 14.4031}, {0.99997,
      14.4502}, {0.999973, 14.4781}, {0.999977, 14.5148}, {0.99998,
      14.5418}, {0.999982, 14.5594}, {0.999984, 14.5766}, {0.999986,
      14.5932}, {0.999988, 14.6091}, {0.999989, 14.6168}, {0.999991,
      14.6312}, {0.999992, 14.6379}, {0.999993, 14.6441}, {0.999994,
      14.65}, {0.999994, 14.65}, {0.999995, 14.6553}, {0.999996,
      14.6602}, {0.999996, 14.6602}, {0.999997, 14.6651}, {0.999997,
      14.6651}, {0.999997, 14.6651}, {0.999998, 14.6723}, {0.999998,
      14.6723}, {0.999998, 14.6723}, {0.999999, 14.696}, {0.999999,
      14.696}, {0.999999, 14.696}};

         double err_low[200][2] = {{1.14022*10e-6, 5.79982}, {1.30011*10e-6, 5.79979}, {1.48241*10e-6,
         5.79974}, {1.69028*10e-6, 5.79969}, {1.9273*10e-6,
         5.79962}, {2.19755*10e-6, 5.79954}, {2.50569*10e-6,
         5.79944}, {2.85705*10e-6, 5.79932}, {3.25767*10e-6,
         5.79917}, {3.71447*10e-6, 5.799}, {4.23532*10e-6,
         5.79879}, {4.82921*10e-6, 5.79853}, {5.50638*10e-6,
         5.79822}, {6.2785*10e-6, 5.79784}, {7.15888*10e-6,
         5.79739}, {8.16272*10e-6, 5.79683}, {9.30732*10e-6,
         5.79617}, {0.0000106124, 5.79536}, {0.0000121005,
         5.79438}, {0.0000137973, 5.7932}, {0.000015732,
         5.79178}, {0.000017938, 5.79006}, {0.0000204533,
         5.78798}, {0.0000233213, 5.78548}, {0.0000265915,
         5.78245}, {0.0000303202, 5.77881}, {0.0000345718,
         5.77442}, {0.0000394195, 5.76913}, {0.0000449471,
         5.76277}, {0.0000512496, 5.75511}, {0.000058436,
         5.74591}, {0.0000666301, 5.73486}, {0.0000759731,
         5.72161}, {0.0000866263, 5.70572}, {0.0000987733,
         5.68671}, {0.000112624, 5.66398}, {0.000128416,
         5.63685}, {0.000146423, 5.60453}, {0.000166955,
         5.56611}, {0.000190365, 5.52055}, {0.000217059,
         5.46666}, {0.000247496, 5.40313}, {0.0002822,
         5.32855}, {0.000321771, 5.24138}, {0.00036689,
         5.14005}, {0.000418337, 5.02299}, {0.000476997,
         4.88881}, {0.000543883, 4.73634}, {0.000620148,
         4.56489}, {0.000707107, 4.37447}, {0.000806259,
         4.16599}, {0.000919315, 3.94149}, {0.00104822,
         3.70432}, {0.00119521, 3.45909}, {0.0013628, 3.21159}, {0.0015539,
         2.96817}, {0.00177179, 2.73519}, {0.00202024, 2.51806}, {0.00230352,
            2.32044}, {0.00262653, 2.14367}, {0.00299483,
         1.98666}, {0.00341477, 1.84631}, {0.0038936, 1.71833}, {0.00443957,
         1.59831}, {0.0050621, 1.48251}, {0.00577192, 1.36849}, {0.00658127,
         1.25518}, {0.00750412, 1.14275}, {0.00855636, 1.03218}, {0.00975616,
            0.924896}, {0.0111242, 0.822379}, {0.0126841, 0.72592}, {0.0144627,
            0.636492}, {0.0164906, 0.554956}, {0.018803, 0.481047}, {0.0214396,
            0.414429}, {0.0244459, 0.356729}, {0.0278738,
         0.307449}, {0.0317824, 0.265257}, {0.036239, 0.229142}, {0.0413205,
         0.19819}, {0.0471146, 0.171674}, {0.0537211, 0.148839}, {0.061254,
         0.129242}, {0.0698433, 0.112431}, {0.0796369,
         0.0979989}, {0.0908038, 0.0856118}, {0.103537,
         0.0747178}, {0.118055, 0.0651285}, {0.134609, 0.0568796}, {0.153484,
            0.0499634}, {0.175006, 0.044301}, {0.199546, 0.0396184}, {0.227526,
            0.0357218}, {0.259431, 0.0324747}, {0.295809,
         0.0298314}, {0.337288, 0.0277645}, {0.384583, 0.0262847}, {0.438511,
            0.0254593}, {0.5, 0.0254542}, {0.561489, 0.0251019}, {0.615417,
         0.0269592}, {0.662712, 0.0296464}, {0.704191, 0.0330823}, {0.740569,
            0.03736}, {0.772474, 0.0425942}, {0.800454, 0.0489351}, {0.824994,
         0.056571}, {0.846516, 0.0657305}, {0.865391, 0.0766881}, {0.881945,
         0.0897716}, {0.896463, 0.105369}, {0.909196, 0.12394}, {0.920363,
         0.146024}, {0.930157, 0.172258}, {0.938746, 0.203386}, {0.946279,
         0.240286}, {0.952885, 0.283985}, {0.958679, 0.335635}, {0.963761,
         0.393259}, {0.968218, 0.459947}, {0.972126, 0.535273}, {0.975554,
         0.619832}, {0.97856, 0.712431}, {0.981197, 0.819863}, {0.983509,
         0.942188}, {0.985537, 1.07983}, {0.987316, 1.23407}, {0.988876,
         1.40617}, {0.990244, 1.59755}, {0.991444, 1.81015}, {0.992496,
         2.04627}, {0.993419, 2.30921}, {0.994228, 2.6023}, {0.994938,
         2.93002}, {0.99556, 3.29652}, {0.996106, 3.69268}, {0.996585,
         4.05124}, {0.997005, 4.36798}, {0.997373, 4.64536}, {0.997696,
         4.89588}, {0.99798, 5.11671}, {0.998228, 5.30506}, {0.998446,
         5.46319}, {0.998637, 5.59308}, {0.998805, 5.69934}, {0.998952,
         5.78653}, {0.999081, 5.86035}, {0.999194, 5.92579}, {0.999293,
         5.98723}, {0.99938, 6.04808}, {0.999456, 6.10989}, {0.999523,
         6.17399}, {0.999582, 6.24024}, {0.999633, 6.30652}, {0.999678,
         6.3729}, {0.999718, 6.43849}, {0.999753, 6.50094}, {0.999783,
         6.55795}, {0.99981, 6.61152}, {0.999833, 6.65799}, {0.999854,
         6.69998}, {0.999872, 6.73445}, {0.999887, 6.76087}, {0.999901,
         6.78225}, {0.999913, 6.79668}, {0.999924, 6.81033}, {0.999933,
         6.78183}, {0.999942, 6.67909}, {0.999949, 6.57172}, {0.999955,
         6.47705}, {0.999961, 6.37749}, {0.999965, 6.30819}, {0.99997,
         6.21802}, {0.999973, 6.16187}, {0.999977, 6.08438}, {0.99998,
         6.02408}, {0.999982, 5.98275}, {0.999984, 5.94044}, {0.999986,
         5.89704}, {0.999988, 5.85246}, {0.999989, 5.82968}, {0.999991,
         5.78307}, {0.999992, 5.7592}, {0.999993, 5.73494}, {0.999994,
         5.71032}, {0.999994, 5.71032}, {0.999995, 5.68541}, {0.999996,
         5.66051}, {0.999996, 5.66051}, {0.999997, 5.63646}, {0.999997,
         5.63646}, {0.999997, 5.63646}, {0.999998, 5.61636}, {0.999998,
         5.61636}, {0.999998, 5.61636}, {0.999999, 5.61749}, {0.999999,
         5.61749}, {0.999999, 5.61749}};


   TH1D HTemp("HTemp", "HTemp", 200, 0, 200); 

   TGraphAsymmErrors graph(&HTemp); 
   

   for(int i = 1; i <= HTemp.GetNbinsX(); i++)
   {
      int iGraph = i-1;
      graph.SetPoint(iGraph, HTemp.GetBinCenter(i), centralvals[iGraph][1]);
      graph.SetPointError(iGraph, HTemp.GetBinWidth(i)/2, HTemp.GetBinWidth(i)/2, centralvals[iGraph][1] - err_low[iGraph][1], err_high[iGraph][1] - centralvals[iGraph][1]);
   }

   static vector<int> Colors = GetCVDColors6();

   graph.SetFillStyle(1001);
   graph.SetMarkerStyle(20);
   graph.SetMarkerColor(Colors[4]);
   graph.SetLineColor(Colors[4]);
   graph.SetLineWidth(2);
   graph.SetFillColorAlpha(Colors[4], 0.3);


   return graph; 

}


TGraphAsymmErrors getOPALPlot(){
   // points for the central values of the theory curves
   double centralvals[100][2] =  {{0.9, 3.24}, {2.7, 1.3}, {4.5, 1.38}, {6.3, 1.197}, {8.1, 0.974}, {9.9, 0.796}, {11.7, 0.661}, {13.5, 0.553}, {15.3, 0.461}, {17.1, 0.390}, {18.9, 0.334}, {20.7, 0.289}, {22.5, 0.255}, {24.3, 0.227}, {26.1, 0.206}, {27.9, 0.118}, {29.7, 0.174}, {31.5, 0.161}, {33.3, 0.150}, {35.1, 0.141}, {36.9, 0.133}, {38.7, 0.127}, {40.5, 0.121}, {42.3, 0.116}, {44.1, 0.110}, {45.9, 0.106}, {47.7, 0.102}, {49.5, 0.099}, {51.3, 0.096}, {53.1, 0.093}, {54.9, 0.091}, {56.7, 0.089}, {58.5, 0.087}, {60.3, 0.085}, {62.1, 0.083}, {63.9, 0.081}, {65.7, 0.081}, {71.1, 0.077}, {72.9, 0.076}, {74.7, 0.076}, {76.5, 0.075}, {78.3, 0.075}, {80.1, 0.075}, {81.9, 0.075}, {83.7, 0.074}, {85.5, 0.074}, {87.3, 0.074}, {89.1, 0.075}, {90.9, 0.076}, {92.7, 0.076}, {94.5, 0.076}, {96.3, 0.078}, {98.1, 0.078}, {99.9, 0.079}, {101.7, 0.080}, {103.5, 0.082}, {105.3, 0.083}, {107.1, 0.085}, {108.9, 0.087}, {110.7, 0.089}, {112.5, 0.091}, {114.3, 0.094}, {116.1, 0.096}, {117.9, 0.099}, {119.7, 0.102}, {121.5, 0.107}, {123.3, 0.110}, {125.1, 0.116}, {126.9, 0.121}, {128.7, 0.125}, {130.5, 0.131}, {132.3, 0.138}, {134.1, 0.146}, {135.9, 0.155}, {137.7, 0.164}, {1395, 0.174}, {141.3, 0.186}, {143.1, 0.186}, {143.1, 0.200}, {144.9, 0.213}, {146.7, 0.230}, {148.5, 0.250}, {150.3, 0.272}, {152.1, 0.299}, {153.9, 0.329}, {155.7, 0.365}, {157.5, 0.410}, {159.3, 0.457}, {161.1, 0.521}, {162.9, 0.595}, {164.7, 0.682}, {166.5, 0.783}, {168.3, 0.906}, {170.1, 1.049}, {171.9, 1.19}, {173.7, 1.31}, {175.5, 1.34}, {177.3, 1.12}, {179.1, 0.46}};
   
   // now need to translate the thetas to zs
   const int BinCount = 100;
   double zBins[2*BinCount+1]; 
   double zBinMin = (1- cos(0.002))/2; 
   double zBinMax = 0.5;

   for(int i = 0; i <= BinCount; i++){
      zBins[i] = exp(log(zBinMin) + (log(zBinMax) - log(zBinMin)) / BinCount * i);
      zBins[2*BinCount-i] = zBinMax * 2 - exp(log(zBinMin) + (log(zBinMax) - log(zBinMin)) / BinCount * i);
   }

   TH1D HTemp("HTemp", "HTemp", 200, 0, 200); 

   TGraphAsymmErrors graph(&HTemp); 
   

   for(int i = 1; i < 100; i++)
   {
      int iGraph = i-1;
      double zval = (1-cos(centralvals[iGraph][0]))/2; 
      double binWidth = (centralvals[iGraph+1][0] - centralvals[iGraph][0])/2; 
      int bin = FindBin(zval, 2*BinCount, zBins);
      std::cout << "theta: " << centralvals[iGraph][0] << " zval: " << zval << " bin " << bin << std::endl;
      double sumVal = centralvals[iGraph][1] * (2/sin(centralvals[iGraph][0])); 
      sumVal = sumVal*binWidth; 
      // now we need to divide by the theta bin width
      double binWidthZ = (zBins[bin+1] - zBins[bin])/2; 
      sumVal = sumVal/binWidthZ; 
      graph.SetPoint(iGraph, HTemp.GetBinCenter(bin), centralvals[iGraph][1]);
      //graph.SetPointError(iGraph, HTemp.GetBinWidth(i)/2, HTemp.GetBinWidth(i)/2, centralvals[iGraph][1] - err_low[iGraph][1], err_high[iGraph][1] - centralvals[iGraph][1]);
   }

   static vector<int> Colors = GetCVDColors6();

   graph.SetFillStyle(1001);
   graph.SetMarkerStyle(20);
   graph.SetMarkerColor(Colors[4]);
   graph.SetLineColor(Colors[4]);
   graph.SetLineWidth(2);
   graph.SetFillColorAlpha(Colors[4], 0.3);


   return graph; 

}


