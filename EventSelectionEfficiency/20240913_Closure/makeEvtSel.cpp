#include <iostream>
#include <vector>
#include <map>
using namespace std;

// root includes
#include "TTree.h"
#include "TFile.h"
#include "TStyle.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TGaxis.h"
#include "TLatex.h"
#include "TLegend.h"

#include "Messenger.h"
#include "CommandLine.h"
#include "Matching.h"
#include "ProgressBar.h"
#include "TauHelperFunctions3.h"
#include "SetStyle.h"



int main(int argc, char *argv[]);
int FindBin(double Value, int NBins, double Bins[]); 
void MakeCanvasZ(vector<TH1D> &Histograms, vector<string> Labels, string Output, string X, string Y, double WorldMin, double WorldMax, bool DoRatio, bool LogX); 
void MakeCanvas(vector<TH1D> &Histograms, vector<string> Labels, string Output, string X, string Y, double WorldMin, double WorldMax, bool DoRatio, bool LogX); 
void SetPad(TPad &P); 
void DivideByBin(TH1D &H, double Bins[]); 

int lowerSThetaBound = 17; 
int upperSThetaBound = 19; 

int main(int argc, char *argv[])
{
   CommandLine CL(argc, argv);

   SetThesisStyle();
   static vector<int> Colors = GetCVDColors6();


   string InputFileName     = CL.Get("Input");
   string InputDataFileName = CL.Get("InputData");
   string OutputFileName    = CL.Get("Output", "EvtSelEff.root");
   string DataTreeName      = CL.Get("Data", "t"); 
   string GenTreeName       = CL.Get("Gen", "tgen");
   string GenBeforeTreeName = CL.Get("GenBefore", "tgenBefore");
   TFile InputFile(InputFileName.c_str());
   TFile InputDataFile(InputDataFileName.c_str()); 

   
   double TotalE = 91.1876;

   ParticleTreeMessenger MGen(InputFile, GenTreeName);
   ParticleTreeMessenger MGenBefore(InputFile, GenBeforeTreeName); 
   ParticleTreeMessenger MData(InputDataFile, DataTreeName);

   //------------------------------------
   // define the binning
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

   // energy binning
   double EnergyBins[BinCount+1];
   double EnergyBinMin = 4e-6; 
   double EnergyBinMax = 0.2; 
   double logMin = std::log10(EnergyBinMin);
   double logMax = std::log10(EnergyBinMax);
   double logStep = (logMax - logMin) / (BinCount);

   for(int i = 0; i <= BinCount; i++){
      // theta double log binning
      Bins[i] = exp(log(BinMin) + (log(BinMax) - log(BinMin)) / BinCount * i);
      Bins[2*BinCount-i] = BinMax * 2 - exp(log(BinMin) + (log(BinMax) - log(BinMin)) / BinCount * i);

      // z double log binning
      zBins[i] = exp(log(zBinMin) + (log(zBinMax) - log(zBinMin)) / BinCount * i);
      zBins[2*BinCount-i] = zBinMax * 2 - exp(log(zBinMin) + (log(zBinMax) - log(zBinMin)) / BinCount * i);

      double logValue = logMin + i * logStep;
      EnergyBins[i] =  std::pow(10, logValue);//exp(log(EnergyBinMin) + (log(EnergyBinMax) - log(EnergyBinMin)) / BinCount * i);
   }

    // -------------------------------------------
    // allocate the histograms
    // -------------------------------------------

    // 2D histograms
    TH2D h2_EvtSelBefore_Theta("h2_EvtSelBefore_Theta", "h2_EvtSelBefore_Theta", 2 * BinCount, 0, 2 * BinCount,  BinCount, EnergyBins);
    TH2D h2_EvtSelBefore_Z("h2_EvtSelBefore_Z", "h2_EvtSelBefore_Z", 2 * BinCount, 0, 2 * BinCount,  BinCount, EnergyBins);
    TH2D h2_EvtSel_Theta("h2_EvtSel_Theta", "h2_EvtSel_Theta", 2 * BinCount, 0, 2 * BinCount,  BinCount, EnergyBins);
    TH2D h2_EvtSel_Z("h2_EvtSel_Z", "h2_EvtSel_Z", 2 * BinCount, 0, 2 * BinCount,  BinCount, EnergyBins);

    // 1D histograms
    TH1D h1_EvtSel_Z("h1_EvtSel_Z", "h2_EvtSel_Z", 2 * BinCount, 0, 2 * BinCount); 
    TH1D h1_EvtSelBefore_Z("h1_EvtSelBefore_Z", "h2_EvtSelBefore_Z", 2 * BinCount, 0, 2 * BinCount); 
    TH1D h1_EvtSel_Theta("h1_EvtSel_Theta", "h2_EvtSel_Theta", 2 * BinCount, 0, 2 * BinCount); 
    TH1D h1_EvtSelBefore_Theta("h1_EvtSelBefore_Theta", "h2_EvtSelBefore_Theta", 2 * BinCount, 0, 2 * BinCount); 

    // Data histograms
    TH1D h1_Data_Z("h1_Data_Z", "h2_EvtSel_Z", 2 * BinCount, 0, 2 * BinCount); 
    TH1D h1_Data_Theta("h1_Data_Theta", "h2_EvtSel_Theta", 2 * BinCount, 0, 2 * BinCount); 


   std::vector<int> lowerSThetaBounds = {7, 11, 15, 19, 21, 25}; 
   std::vector<int> upperSThetaBounds = {11, 15, 19, 21, 25, 29};
   std::vector<TH1D> vec_h1_EvtSel_Z; 
   std::vector<TH1D> vec_h1_Data_Z; 
   // need vectors to store the number of events in each sphericity bin
   std::vector<int> nEventsMC; 
   std::vector<int> nEventsData; 



   for(int j = 0; j < lowerSThetaBounds.size(); j++){
      TH1D h1_Data_Z(Form("h1_Data_Z_%d", j), Form("h1_Data_Z_%d", j), 2 * BinCount, 0, 2 * BinCount);
      vec_h1_Data_Z.push_back(h1_Data_Z); 
      TH1D h1_EvtSel_Z(Form("h1_EvtSel_Z_%d", j), Form("h1_EvtSel_Z_%d", j), 2 * BinCount, 0, 2 * BinCount);
      vec_h1_EvtSel_Z.push_back(h1_EvtSel_Z); 
      // initialize number of events to be 0
      nEventsData.push_back(0); 
      nEventsMC.push_back(0); 
   }
   
    // -------------------------------------
    // loop over the data tree 
    // -------------------------------------
   int EntryCountData = MData.GetEntries();
   for(int iE = 0; iE < EntryCountData; iE++){
      MData.GetEntry(iE);



      // fill the four vector
      vector<FourVector> PData;
      for(int i = 0; i < MData.nParticle; i++){
        // charged particle selection 
       if(MData.charge[i] == 0) continue;
       if(MData.highPurity[i] == false) continue;
         PData.push_back(MData.P[i]);
      } // end loop over the particles 


      // now calculate and fill the EECs
      for(int i = 0; i < PData.size(); i++){
        for(int j = i+1; j < PData.size();j++){
            FourVector Data1 = PData.at(i);
            FourVector Data2 = PData.at(j);
            
            // get the proper bins
            int BinThetaData  = FindBin(GetAngle(Data1,Data2), 2 * BinCount, Bins);
            int BinEnergyData = FindBin(Data1[0]*Data2[0]/(TotalE*TotalE), BinCount, EnergyBins);
            double zData = (1-cos(GetAngle(Data1, Data2)))/2; 
            int BinZData = FindBin(zData, 2*BinCount, zBins); 

            // calculate the EEC
            double EEC =  Data1[0]*Data2[0]/(TotalE*TotalE); 
            
            // fill the histograms
            h1_Data_Theta.Fill(BinThetaData, EEC); 
            h1_Data_Z.Fill(BinZData, EEC);
            // now fill the vectors
            for(int j = 0; j < lowerSThetaBounds.size(); j++){
               if(MData.STheta >= (lowerSThetaBounds.at(j)*M_PI)/36 &&  MData.STheta <= (upperSThetaBounds.at(j)*M_PI)/36){
                  vec_h1_Data_Z.at(j).Fill(BinZData, EEC); 
                  nEventsData.at(j) = nEventsData.at(j) + 1; 
               }
            }
         }
      }

    } // end loop over the number of events'


    // -------------------------------------
    // loop over the tree after event selections
    // -------------------------------------
   int EntryCount = MGen.GetEntries();
   for(int iE = 0; iE < EntryCount; iE++){
      MGen.GetEntry(iE);


      // fill the four vector
      vector<FourVector> PGen;
      for(int i = 0; i < MGen.nParticle; i++){
        // charged particle selection 
       if(MGen.charge[i] == 0) continue;
       if(MGen.highPurity[i] == false) continue;
         PGen.push_back(MGen.P[i]);
      } // end loop over the particles 


      // now calculate and fill the EECs
      for(int i = 0; i < PGen.size(); i++){
        for(int j = i+1; j < PGen.size();j++){
            FourVector Gen1 = PGen.at(i);
            FourVector Gen2 = PGen.at(j);
            
            // get the proper bins
            int BinThetaGen  = FindBin(GetAngle(Gen1,Gen2), 2 * BinCount, Bins);
            int BinEnergyGen = FindBin(Gen1[0]*Gen2[0]/(TotalE*TotalE), BinCount, EnergyBins);
            double zGen = (1-cos(GetAngle(Gen1, Gen2)))/2; 
            int BinZGen = FindBin(zGen, 2*BinCount, zBins); 

            // calculate the EEC
            double EEC =  Gen1[0]*Gen2[0]/(TotalE*TotalE); 
            
            // fill the histograms
            h2_EvtSel_Theta.Fill(BinThetaGen, BinEnergyGen, EEC); 
            h2_EvtSel_Z.Fill(BinZGen, BinEnergyGen, EEC); 
            h1_EvtSel_Z.Fill(BinZGen, EEC); 
            h1_EvtSel_Theta.Fill(BinThetaGen, EEC);

            // now fill the vectors
            for(int j = 0; j < lowerSThetaBounds.size(); j++){
               if(MGen.STheta >= (lowerSThetaBounds.at(j)*M_PI)/36 &&  MGen.STheta <= (upperSThetaBounds.at(j)*M_PI)/36){
                  vec_h1_EvtSel_Z.at(j).Fill(BinZGen, EEC); 
                  nEventsMC.at(j) = nEventsMC.at(j) + 1; 
               }
            }

         }
      }

    } // end loop over the number of events'


    // -------------------------------------
    // loop over the tree before event selections
    // -------------------------------------
   int EntryCountBefore = MGenBefore.GetEntries();
   for(int iE = 0; iE < EntryCountBefore; iE++){
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
            int BinThetaGenBefore  = FindBin(GetAngle(Gen1,Gen2), 2 * BinCount, Bins);
            int BinEnergyGenBefore = FindBin(Gen1[0]*Gen2[0]/(TotalE*TotalE), BinCount, EnergyBins);
            double zGenBefore = (1-cos(GetAngle(Gen1, Gen2)))/2; 
            int BinZGenBefore = FindBin(zGenBefore, 2*BinCount, zBins); 

            // calculate the EEC
            double EEC =  Gen1[0]*Gen2[0]/(TotalE*TotalE); 

            // fill the histograms
            h2_EvtSelBefore_Theta.Fill(BinThetaGenBefore, BinEnergyGenBefore, EEC); 
            h2_EvtSelBefore_Z.Fill(BinZGenBefore, BinEnergyGenBefore, EEC); 
            h1_EvtSelBefore_Z.Fill(BinZGenBefore, EEC); 
            h1_EvtSelBefore_Theta.Fill(BinThetaGenBefore, EEC);

         }
      }

    } // end loop over the number of events
   // EEC is per-event so scale by the event number
   h1_EvtSel_Z.Scale(1.0/EntryCount); 
   h1_EvtSelBefore_Z.Scale(1.0/EntryCountBefore);
   h1_EvtSel_Theta.Scale(1.0/EntryCount);
   h1_EvtSelBefore_Theta.Scale(1.0/EntryCountBefore);

   h1_Data_Theta.Scale(1.0/EntryCountData); 
   h1_Data_Z.Scale(1.0/EntryCountData); 

   // divide by the bin width
   DivideByBin(h1_EvtSel_Z, zBins);
   DivideByBin(h1_EvtSelBefore_Z, zBins);
   DivideByBin(h1_EvtSel_Theta, Bins);
   DivideByBin(h1_EvtSelBefore_Theta, Bins);
   DivideByBin(h1_Data_Theta, Bins); 
   DivideByBin(h1_Data_Z, zBins); 

   // set the style for the plots
   h1_EvtSel_Z.SetMarkerColor(Colors[2]);
   h1_EvtSelBefore_Z.SetMarkerColor(Colors[0]);
   h1_EvtSel_Theta.SetMarkerColor(Colors[2]);
   h1_EvtSelBefore_Theta.SetMarkerColor(Colors[0]);
   h1_Data_Theta.SetMarkerColor(Colors[1]); 
   h1_Data_Z.SetMarkerColor(Colors[1]); 
   h1_Data_Theta.SetLineColor(Colors[1]);
   h1_Data_Z.SetLineColor(Colors[1]);  
   h1_EvtSel_Z.SetLineColor(Colors[2]);
   h1_EvtSelBefore_Z.SetLineColor(Colors[0]);
   h1_EvtSel_Theta.SetLineColor(Colors[2]);
   h1_EvtSelBefore_Theta.SetLineColor(Colors[0]);
   h1_EvtSel_Z.SetMarkerStyle(20);
   h1_EvtSelBefore_Z.SetMarkerStyle(20);
   h1_EvtSel_Theta.SetMarkerStyle(20);
   h1_EvtSelBefore_Theta.SetMarkerStyle(20);
   h1_Data_Theta.SetMarkerStyle(20); 
   h1_Data_Z.SetMarkerStyle(20); 
   h1_EvtSel_Z.SetLineWidth(2);
   h1_EvtSelBefore_Z.SetLineWidth(2);
   h1_EvtSel_Theta.SetLineWidth(2);
   h1_Data_Theta.SetLineWidth(2);
   h1_Data_Z.SetLineWidth(2); 
   h1_EvtSelBefore_Theta.SetLineWidth(2);

   // construct the event selection efficiency
   TH1D *HEvtSelTheta = (TH1D*)h1_EvtSel_Theta.Clone();
   HEvtSelTheta->Divide(&h1_EvtSelBefore_Theta); 

   TH1D *HEvtSelZ= (TH1D*)h1_EvtSel_Z.Clone();
   HEvtSelZ->Divide(&h1_EvtSelBefore_Z); 

   TH1D* dataCloneZ     = (TH1D*)h1_Data_Z.Clone(); 
   TH1D* dataCloneTheta = (TH1D*)h1_Data_Theta.Clone(); 
   // correct for the event selection efficiency
   dataCloneZ->Divide(HEvtSelZ); 
   dataCloneTheta->Divide(HEvtSelTheta); 


    // Event selection is given by SelectedEvents/TotalEvents
   std::vector<TH1D> hists = {h1_EvtSelBefore_Z, *dataCloneZ}; 
   std::vector<TH1D> hists_theta = {h1_EvtSelBefore_Theta, *dataCloneTheta};

   MakeCanvasZ(hists, {"Gen Level - No Evt. Sel.", "Data w./ Evt. Sel. Corr."}, "EventSelectionWithEventScaling", "#it{z} = (1- cos(#theta))/2", "#frac{1}{#it{N}_{event}}#frac{d(Sum E_{i}E_{j}/E^{2})}{d#it{z}}",1e-3,1e3, true, true);
   MakeCanvas(hists_theta, {"Gen Level - No Evt. Sel.", "Data w./ Evt. Sel. Corr."},"EventSelectionWithEventScaling_Theta","#theta", "#frac{1}{#it{N}_{event}}#frac{d(Sum E_{i}E_{j}/E^{2})}{d#it{#theta}}",1e-3, 2e0, true, true);



   // now handle the vectors 
   std::vector<TH1D> evtCorrDataHists; 
   std::vector<std::string> labels; 
   labels.push_back("Gen Level - No Evt. Sel.");
   evtCorrDataHists.push_back(h1_EvtSelBefore_Z); 
    for(int j = 0; j < lowerSThetaBounds.size(); j++){
      vec_h1_Data_Z.at(j).SetMarkerColor(Colors[j+1]);
      vec_h1_Data_Z.at(j).Scale(1.0/nEventsData.at(j));
      vec_h1_Data_Z.at(j).SetLineColor(Colors[j+1]);
      vec_h1_Data_Z.at(j).SetMarkerStyle(20);
      vec_h1_Data_Z.at(j).SetLineWidth(2); 
      vec_h1_EvtSel_Z.at(j).Scale(1.0/nEventsMC.at(j));
      DivideByBin(vec_h1_EvtSel_Z.at(j), zBins);
      DivideByBin(vec_h1_Data_Z.at(j), zBins);


      // calculate the eventselection efficiency
      TH1D *HEvtSelZ = (TH1D*)vec_h1_EvtSel_Z.at(j).Clone(Form("mc_%d", j)); 
      HEvtSelZ->Divide(&h1_EvtSelBefore_Z); 

      TH1D *HData = (TH1D*)vec_h1_Data_Z.at(j).Clone(Form("data_%d", j)); 
      HData->Divide(HEvtSelZ);
      evtCorrDataHists.push_back(*HData);
      labels.push_back(Form("(%d#pi/36) < #it{S}_{theta} < (%d#pi/36)", lowerSThetaBounds.at(j), upperSThetaBounds.at(j)));
    }

    MakeCanvasZ(evtCorrDataHists, labels, "SThetaOverlay", "#it{z} = (1- cos(#theta))/2", "#frac{1}{#it{N}_{event}}#frac{d(Sum E_{i}E_{j}/E^{2})}{d#it{z}}",1e-3,1e5, true, true);




    // write to the output file
    TFile OutputFile(OutputFileName.c_str(), "RECREATE");
    OutputFile.cd(); 
    h2_EvtSel_Theta.Write(); 
    h2_EvtSel_Z.Write(); 
    h1_EvtSel_Z.Write(); 
    h1_EvtSelBefore_Z.Write(); 
    h1_EvtSelBefore_Theta.Write(); 
    h1_EvtSel_Theta.Write(); 




   // cleanup
   InputFile.Close(); 



    

   return 0;
}


int FindBin(double Value, int NBins, double Bins[])
{
   for(int i = 0; i < NBins; i++)
      if(Value < Bins[i])
         return i - 1;
   return NBins;
}

void MakeCanvasZ(vector<TH1D > &Histograms, vector<string> Labels, string Output,
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

   TLatex Latex;
   Latex.SetNDC();
   Latex.SetTextFont(42);
   Latex.SetTextSize(0.05);
   Latex.SetTextAlign(23);
   //Latex.DrawLatex(MarginL + PadWidth * 0.85, MarginB + 1.45*PadHeight, Form("(%d#pi/36) < #it{S}_{theta} < (%d#pi/36)", lowerSThetaBound, upperSThetaBound));
   Latex.SetTextSize(0.035);


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
      Latex.DrawLatex(MarginL * 0.3, MarginB + PadRHeight * 0.5, "Ratio");
   Latex.DrawLatex(MarginL * 0.3, MarginB + PadRHeight + PadHeight * 0.5, Y.c_str());

   Latex.SetTextAlign(11);
   Latex.SetTextAngle(0);
   Latex.DrawLatex(MarginL, MarginB + PadRHeight + PadHeight + 0.012, "ALEPH e^{+}e^{-}, #sqrt{s} = 91.2 GeV, Work-in-progress");

   Latex.SetTextAlign(11);
   Latex.SetTextAngle(0);
   Latex.SetTextColor(19);
   Latex.SetTextSize(0.02);
   Latex.DrawLatex(0.01, 0.01, "Work-in-progress, 2024 September 13th HB (20240913_Closure)");

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


void MakeCanvas(vector<TH1D> &Histograms, vector<string> Labels, string Output,
   string X, string Y, double WorldMin, double WorldMax, bool DoRatio, bool LogX)
{
   int NLine = Histograms.size();
   int N = Histograms[0].GetNbinsX();

   double MarginL = 180;
   double MarginR = 90;
   double MarginB = 120;
   double MarginT = 90;

   double WorldXMin = LogX ? 0 : 0;
   double WorldXMax = LogX ? N : M_PI;
   
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
      TH1D *H1 = (TH1D *)H.Clone();
      H1->Draw("same");
   }

   TLatex Latex;
   Latex.SetNDC();
   Latex.SetTextFont(42);
   Latex.SetTextSize(0.05);
   Latex.SetTextAlign(23);
   // Latex.DrawLatex(MarginL + PadWidth * 0.85, MarginB + 1.45*PadHeight, Form("(%d#pi/36) < #it{S}_{theta} < (%d#pi/36)", lowerSThetaBound, upperSThetaBound));
   Latex.SetTextSize(0.035);

   TGraph G;
   G.SetPoint(0, LogX ? N / 2 : M_PI / 2, 0);
   G.SetPoint(1, LogX ? N / 2 : M_PI / 2, 1000);
   G.SetLineStyle(kDashed);
   G.SetLineColor(kGray);
   G.SetLineWidth(1);
   G.Draw("l");

   if(DoRatio)
      PadR.cd();

   double WorldRMin = 0.6;
   double WorldRMax = 1.4;
   
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
   
   double BinMin    = 0.002;
   double BinMiddle = M_PI / 2;
   double BinMax    = M_PI - 0.002;

   Canvas.cd();
   TGaxis X1(MarginL, MarginB, MarginL + PadWidth / 2, MarginB, BinMin, BinMiddle, 510, "GS");
   TGaxis X2(MarginL + PadWidth, MarginB, MarginL + PadWidth / 2, MarginB, BinMin, BinMiddle, 510, "-GS");
   TGaxis X3(MarginL, MarginB + PadRHeight, MarginL + PadWidth / 2, MarginB + PadRHeight, BinMin, BinMiddle, 510, "+-GS");
   TGaxis X4(MarginL + PadWidth, MarginB + PadRHeight, MarginL + PadWidth / 2, MarginB + PadRHeight, BinMin, BinMiddle, 510, "+-GS");
   TGaxis Y1(MarginL, MarginB, MarginL, MarginB + PadRHeight, WorldRMin, WorldRMax, 505, "");
   TGaxis Y2(MarginL, MarginB + PadRHeight, MarginL, MarginB + PadRHeight + PadHeight, WorldMin, WorldMax, 510, "G");
   
   TGaxis XL1(MarginL, MarginB, MarginL + PadWidth, MarginB, 0, M_PI, 510, "S");
   TGaxis XL2(MarginL, MarginB + PadRHeight, MarginL + PadWidth, MarginB + PadRHeight, 0, M_PI, 510, "+-S");

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

   if(LogX) Latex.DrawLatex(MarginL + (1 - MarginR - MarginL) * 0.115, MarginB - 0.01, "0.01");
   if(LogX) Latex.DrawLatex(MarginL + (1 - MarginR - MarginL) * 0.290, MarginB - 0.01, "0.1");
   if(LogX) Latex.DrawLatex(MarginL + (1 - MarginR - MarginL) * 0.465, MarginB - 0.01, "1");
   if(LogX) Latex.DrawLatex(MarginL + (1 - MarginR - MarginL) * 0.535, MarginB - 0.01, "#pi - 1");
   if(LogX) Latex.DrawLatex(MarginL + (1 - MarginR - MarginL) * 0.710, MarginB - 0.01, "#pi - 0.1");
   if(LogX) Latex.DrawLatex(MarginL + (1 - MarginR - MarginL) * 0.885, MarginB - 0.01, "#pi - 0.01");

   Latex.SetTextAlign(12);
   Latex.SetTextAngle(270);
   Latex.SetTextColor(kGray);
   Latex.DrawLatex(MarginL + (1 - MarginR - MarginL) * 0.5 + 0.0175, 1 - MarginT - 0.015, "#theta_{L} = #pi/2");

   Latex.SetTextAlign(22);
   Latex.SetTextAngle(0);
   Latex.SetTextColor(kBlack);
   Latex.DrawLatex(MarginL + PadWidth * 0.5, MarginB * 0.3, X.c_str());



   Latex.SetTextAlign(22);
   Latex.SetTextAngle(90);
   Latex.SetTextColor(kBlack);
   if(DoRatio)
      Latex.DrawLatex(MarginL * 0.3, MarginB + PadRHeight * 0.5, "Ratio");
   Latex.DrawLatex(MarginL * 0.3, MarginB + PadRHeight + PadHeight * 0.5, Y.c_str());

   Latex.SetTextAlign(11);
   Latex.SetTextAngle(0);
   Latex.DrawLatex(MarginL, MarginB + PadRHeight + PadHeight + 0.012, "ALEPH e^{+}e^{-}, #sqrt{s} = 91.2 GeV, Work-in-progress");

   Latex.SetTextAlign(11);
   Latex.SetTextAngle(0);
   Latex.SetTextColor(19);
   Latex.SetTextSize(0.02);
   Latex.DrawLatex(0.01, 0.01, "Work-in-progress September 13th HB (20240913_Closure)");

   TLegend Legend(0.15, 0.90, 0.35, 0.90 - 0.06 * min(NLine, 4));
   Legend.SetTextFont(42);
   Legend.SetTextSize(0.035);
   Legend.SetFillStyle(0);
   Legend.SetBorderSize(0);
   for(int i = 0; i < NLine && i < 4; i++)
      Legend.AddEntry(&Histograms[i], Labels[i].c_str(), "pl");
   Legend.Draw();

   TLegend Legend2(0.7, 0.90, 0.9, 0.90 - 0.06 * (NLine - 4));
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