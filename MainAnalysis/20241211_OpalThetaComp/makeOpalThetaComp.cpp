#include <iostream>
#include <vector>
#include <map>
using namespace std;

#include "TTree.h"
#include "TFile.h"
#include "TStyle.h"
#include "Messenger.h"
#include "CommandLine.h"
#include "Matching.h"
#include "ProgressBar.h"
#include "TauHelperFunctions3.h"
#include "alephTrkEfficiency.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLegend.h"
#include "TGaxis.h"
#include "TPad.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TMath.h"

#include "SetStyle.h"
//#include "CommandLine.h"




int main(int argc, char *argv[]);
TH1D getOPALPlot(); 
TGraphAsymmErrors getOPALGraph(); 
void DivideByBin(TH1D &H, double Bins[]);
int FindBin(double Value, int NBins, double Bins[]); 
void MakeCanvas(vector<TH1D *> Histograms, vector<string> Labels, string Output,
   string X, string Y, double WorldMin, double WorldMax, bool DoRatio, bool LogX); 
void SetPad(TPad &P); 
double FindBinFraction(double Value, int NBins, double Bins[]); 
int main(int argc, char *argv[])
{
   CommandLine CL(argc, argv);
   string InputFileName  = CL.Get("Input");
   string GenTreeName    = CL.Get("Gen", "tgen");
   double Fraction         = CL.GetDouble("Fraction", 1.00);


   // set up the input file and the messenger
   TFile* InputFile  = new TFile(InputFileName.c_str());
   ParticleTreeMessenger* MGen = new ParticleTreeMessenger(InputFile, GenTreeName);

   // set up the binning and the histogram
   int BinCount = 50;
   double LinearBins[2*BinCount+1];
   for(int i = 0; i <= 2 * BinCount; i++) LinearBins[i] = 180.0 / (2 * BinCount) * i;
   TH1D HLinearEEC2("HLinearEEC2", ";EEC_{2};", 2 * BinCount, LinearBins);

   int BinCountDoubleLog = 100; 
   double Bins[2*BinCountDoubleLog+1];
   double BinMin = 0.002;
   double BinMax = M_PI / 2;
   for(int i = 0; i <= BinCountDoubleLog; i++)
   {
      Bins[i] = exp(log(BinMin) + (log(BinMax) - log(BinMin)) / BinCountDoubleLog * i);
      Bins[2*BinCountDoubleLog-i] = BinMax * 2 - exp(log(BinMin) + (log(BinMax) - log(BinMin)) / BinCountDoubleLog * i);
   }

   TH1D HEEC2("HEEC2", ";EEC_{2};", 2 * BinCountDoubleLog, 0, 2 * BinCountDoubleLog);


   // include the total energy 
   double TotalE = 91.1876;



  // loop over the number of events
   int EntryCount = MGen->GetEntries()*Fraction;
   ProgressBar Bar(cout, EntryCount);
   Bar.SetStyle(-1); 
   for(int iE = 0; iE < EntryCount; iE++){
      MGen->GetEntry(iE);

      if(EntryCount < 300 || (iE % (EntryCount / 250) == 0))
      {
         Bar.Update(iE);
         Bar.Print();
      }


      // fill the four vectors
      vector<FourVector> PGen;
      for(int i = 0; i < MGen->nParticle; i++){ 
        // charged particle selection 
        if(MGen->charge[i] != 0 && MGen->highPurity[i] == false) continue;
        PGen.push_back(MGen->P[i]);
      }
      //std::cout << "Pgen size " <<  PGen.size() << std::endl;
      // now create the E2C
      for(int i = 0; i < PGen.size(); i++){
         for(int j = 0; j < PGen.size();j++){
            FourVector Gen1 = PGen.at(i);
            FourVector Gen2 = PGen.at(j);
            double angleInDegrees = GetAngle(Gen1,Gen2) * (180.0/TMath::Pi()); 
            double angle = GetAngle(Gen1, Gen2); 

            // handle the linear case
            int bin = FindBin(angleInDegrees, 2*BinCount, LinearBins); 
            HLinearEEC2.Fill(angleInDegrees,Gen1[0]*Gen2[0]/(TotalE*TotalE));

            // handle the double log case
            int binDoubleLog = FindBin(angle, 2*BinCountDoubleLog, Bins); 
            HEEC2.Fill(binDoubleLog, Gen1[0]*Gen2[0]/(TotalE*TotalE)); 
         }
      } 
   }// end loop over the number of events
   static vector<int> Colors = GetCVDColors6();
   HLinearEEC2.Scale(1.0/EntryCount);
   DivideByBin(HLinearEEC2, LinearBins);
   HLinearEEC2.SetMarkerColor(Colors[0]);
   HLinearEEC2.SetMarkerStyle(20);
   HLinearEEC2.SetLineColor(Colors[0]);
   // to be consistent with opal, need to convert sum on the y-axis to radians
   HLinearEEC2.Scale(180/M_PI);
   HLinearEEC2.SetLineWidth(2); 

   HEEC2.Scale(1.0/EntryCount); 
   DivideByBin(HEEC2, Bins); 
   HEEC2.SetMarkerColor(Colors[0]);
   HEEC2.SetMarkerStyle(20);
   HEEC2.SetLineColor(Colors[0]);

   
   TH1D opalLinear = (TH1D)getOPALPlot(); 
   vector<TH1D*> hists; 
   hists.push_back(&HLinearEEC2); 
   hists.push_back(&opalLinear); 

   vector<TH1D*> histsDoubleLog; 
   histsDoubleLog.push_back(&HEEC2);  
 
  
   // now make the canvas
   MakeCanvas(hists, {"Archived Gen MC", "OPAL"} , "LinearE2C_OpalComp",
    "#theta_{L}", "#frac{1}{N_{event}} #frac{d(Sum E_{i}E_{j}/E^{2})}{d #theta_{L}}", 2e-3, 10, false, false);

   MakeCanvas(hists, {"Archived Gen MC", "OPAL"} , "LinearE2C_OpalCompRatio",
    "#theta_{L}", "#frac{1}{N_{event}} #frac{d(Sum E_{i}E_{j}/E^{2})}{d #theta_{L}}", 2e-3, 10, true, false);

    MakeCanvas(histsDoubleLog, {"Archived Gen MC"} , "DoubleLogE2C_OpalComp",
    "#theta_{L}", "#frac{1}{N_{event}} #frac{d(Sum E_{i}E_{j}/E^{2})}{d #theta_{L}}", 2e-3, 10, false, true);
   
   // cleanup
   InputFile->Close(); 
   delete InputFile; 
   delete MGen; 

    

   return 0;
}


TH1D getOPALPlot(){
   // points for the central values of the theory curves
   double centralvals[100][2] =  {{0.9, 3.24}, {2.7, 1.3}, {4.5, 1.38}, {6.3, 1.197}, {8.1, 0.974}, {9.9, 0.796}, {11.7, 0.661}, {13.5, 0.553}, {15.3, 0.461}, {17.1, 0.390}, {18.9, 0.334}, {20.7, 0.289}, {22.5, 0.255}, {24.3, 0.227}, {26.1, 0.206}, {27.9, 0.188}, {29.7, 0.174}, {31.5, 0.161}, {33.3, 0.150}, {35.1, 0.141}, {36.9, 0.133}, {38.7, 0.127}, {40.5, 0.121}, {42.3, 0.116}, {44.1, 0.110}, {45.9, 0.106}, {47.7, 0.102}, {49.5, 0.099}, {51.3, 0.096}, {53.1, 0.093}, {54.9, 0.091}, {56.7, 0.089}, {58.5, 0.087}, {60.3, 0.085}, {62.1, 0.083}, {63.9, 0.081}, {65.7, 0.081}, {67.5, 0.079}, {69.3, 0.078}, {71.1, 0.077}, {72.9, 0.076}, {74.7, 0.076}, {76.5, 0.075}, {78.3, 0.075}, {80.1, 0.075}, {81.9, 0.075}, {83.7, 0.074}, {85.5, 0.074}, {87.3, 0.074}, {89.1, 0.075}, {90.9, 0.076}, {92.7, 0.076}, {94.5, 0.076}, {96.3, 0.078}, {98.1, 0.078}, {99.9, 0.079}, {101.7, 0.080}, {103.5, 0.082}, {105.3, 0.083}, {107.1, 0.085}, {108.9, 0.087}, {110.7, 0.089}, {112.5, 0.091}, {114.3, 0.094}, {116.1, 0.096}, {117.9, 0.099}, {119.7, 0.102}, {121.5, 0.107}, {123.3, 0.110}, {125.1, 0.116}, {126.9, 0.121}, {128.7, 0.125}, {130.5, 0.131}, {132.3, 0.138}, {134.1, 0.146}, {135.9, 0.155}, {137.7, 0.164}, {139.5, 0.174}, {141.3, 0.186}, {143.1, 0.200}, {144.9, 0.213}, {146.7, 0.230}, {148.5, 0.250}, {150.3, 0.272}, {152.1, 0.299}, {153.9, 0.329}, {155.7, 0.365}, {157.5, 0.410}, {159.3, 0.457}, {161.1, 0.521}, {162.9, 0.595}, {164.7, 0.682}, {166.5, 0.783}, {168.3, 0.906}, {170.1, 1.049}, {171.9, 1.19}, {173.7, 1.31}, {175.5, 1.34}, {177.3, 1.12}, {179.1, 0.46}};

   // set up the binning so that we know where to put the points
   const int BinCount = 50;
   double LinearBins[2*BinCount+1];
   for(int i = 0; i <= 2 * BinCount; i++) LinearBins[i] = 180.0 / (2 * BinCount) * i;
   TH1D HTemp("HTemp", "HTemp",  2 * BinCount, LinearBins); 

   // loop over the points and fill the
   for(int i = 1; i < 101; i++)
   {
      int iGraph = i-1;
      
      // get the central value for the data-point in degrees
      double thetaDegrees = centralvals[iGraph][0]; 
      HTemp.Fill(thetaDegrees, centralvals[iGraph][1]); 
      HTemp.SetBinError(i, 1e-5); 

   }


   static vector<int> Colors = GetCVDColors6();
   HTemp.SetMarkerStyle(20);
   HTemp.SetMarkerColor(Colors[4]);
   HTemp.SetLineColor(Colors[4]);
   HTemp.SetLineWidth(2);

   return HTemp; 

}

TGraphAsymmErrors getOPALGraph(){
   // points for the central values of the theory curves
   double centralvals[100][2] =  {{0.9, 3.24}, {2.7, 1.3}, {4.5, 1.38}, {6.3, 1.197}, {8.1, 0.974}, {9.9, 0.796}, {11.7, 0.661}, {13.5, 0.553}, {15.3, 0.461}, {17.1, 0.390}, {18.9, 0.334}, {20.7, 0.289}, {22.5, 0.255}, {24.3, 0.227}, {26.1, 0.206}, {27.9, 0.188}, {29.7, 0.174}, {31.5, 0.161}, {33.3, 0.150}, {35.1, 0.141}, {36.9, 0.133}, {38.7, 0.127}, {40.5, 0.121}, {42.3, 0.116}, {44.1, 0.110}, {45.9, 0.106}, {47.7, 0.102}, {49.5, 0.099}, {51.3, 0.096}, {53.1, 0.093}, {54.9, 0.091}, {56.7, 0.089}, {58.5, 0.087}, {60.3, 0.085}, {62.1, 0.083}, {63.9, 0.081}, {65.7, 0.081}, {67.5, 0.079}, {69.3, 0.078}, {71.1, 0.077}, {72.9, 0.076}, {74.7, 0.076}, {76.5, 0.075}, {78.3, 0.075}, {80.1, 0.075}, {81.9, 0.075}, {83.7, 0.074}, {85.5, 0.074}, {87.3, 0.074}, {89.1, 0.075}, {90.9, 0.076}, {92.7, 0.076}, {94.5, 0.076}, {96.3, 0.078}, {98.1, 0.078}, {99.9, 0.079}, {101.7, 0.080}, {103.5, 0.082}, {105.3, 0.083}, {107.1, 0.085}, {108.9, 0.087}, {110.7, 0.089}, {112.5, 0.091}, {114.3, 0.094}, {116.1, 0.096}, {117.9, 0.099}, {119.7, 0.102}, {121.5, 0.107}, {123.3, 0.110}, {125.1, 0.116}, {126.9, 0.121}, {128.7, 0.125}, {130.5, 0.131}, {132.3, 0.138}, {134.1, 0.146}, {135.9, 0.155}, {137.7, 0.164}, {139.5, 0.174}, {141.3, 0.186}, {143.1, 0.200}, {144.9, 0.213}, {146.7, 0.230}, {148.5, 0.250}, {150.3, 0.272}, {152.1, 0.299}, {153.9, 0.329}, {155.7, 0.365}, {157.5, 0.410}, {159.3, 0.457}, {161.1, 0.521}, {162.9, 0.595}, {164.7, 0.682}, {166.5, 0.783}, {168.3, 0.906}, {170.1, 1.049}, {171.9, 1.19}, {173.7, 1.31}, {175.5, 1.34}, {177.3, 1.12}, {179.1, 0.46}};
   double errs[100][2] =  {{0.9, 0.18}, {2.7, 0.03}, {4.5, 0.03}, {6.3, 0.016}, {8.1, 0.017}, {9.9, 0.019}, {11.7, 0.019}, {13.5, 0.016}, {15.3, 0.011}, {17.1, 0.011}, {18.9, 0.011}, {20.7, 0.010}, {22.5, 0.011}, {24.3, 0.009}, {26.1, 0.007}, {27.9, 0.007}, {29.7, 0.007}, {31.5, 0.006}, {33.3, 0.006}, {35.1, 0.004}, {36.9, 0.005}, {38.7, 0.004}, {40.5, 0.004}, {42.3, 0.004}, {44.1, 0.003}, {45.9, 0.003}, {47.7, 0.003}, {49.5, 0.002}, {51.3, 0.003}, {53.1, 0.002}, {54.9, 0.003}, {56.7, 0.002}, {58.5, 0.002}, {60.3, 0.003}, {62.1, 0.003}, {63.9, 0.003}, {65.7, 0.003}, {67.5, 0.003}, {69.3, 0.002}, {71.1, 0.002}, {72.9, 0.002}, {74.7, 0.002}, {76.5, 0.002}, {78.3, 0.002}, {80.1, 0.002}, {81.9, 0.002}, {83.7, 0.002}, {85.5, 0.002}, {87.3, 0.001}, {89.1, 0.002}, {90.9, 0.002}, {92.7, 0.002}, {94.5, 0.002}, {96.3, 0.002}, {98.1, 0.002}, {99.9, 0.003}, {101.7, 0.002}, {103.5, 0.002}, {105.3, 0.002}, {107.1, 0.002}, {108.9, 0.002}, {110.7, 0.002}, {112.5, 0.003}, {114.3, 0.002}, {116.1, 0.003}, {117.9, 0.003}, {119.7, 0.002}, {121.5, 0.003}, {123.3, 0.004}, {125.1, 0.004}, {126.9, 0.004}, {128.7, 0.003}, {130.5, 0.003}, {132.3, 0.003}, {134.1, 0.003}, {135.9, 0.004}, {137.7, 0.004}, {139.5, 0.003}, {141.3, 0.003}, {143.1, 0.002}, {144.9, 0.002}, {146.7, 0.002}, {148.5, 0.003}, {150.3, 0.004}, {152.1, 0.003}, {153.9, 0.005}, {155.7, 0.005}, {157.5, 0.005}, {159.3, 0.009}, {161.1, 0.006}, {162.9, 0.010}, {164.7, 0.012}, {166.5, 0.007}, {168.3, 0.008}, {170.1, 0.014}, {171.9, 0.02}, {173.7, 0.02}, {175.5, 0.04}, {177.3, 0.04}, {179.1, 0.02}};

   // set up the binning so that we know where to put the points
   int BinCountDoubleLog = 100; 
   double Bins[2*BinCountDoubleLog+1];
   double BinMin = 0.002;
   double BinMax = M_PI / 2;
   for(int i = 0; i <= BinCountDoubleLog; i++)
   {
      Bins[i] = exp(log(BinMin) + (log(BinMax) - log(BinMin)) / BinCountDoubleLog * i);
      Bins[2*BinCountDoubleLog-i] = BinMax * 2 - exp(log(BinMin) + (log(BinMax) - log(BinMin)) / BinCountDoubleLog * i);
   }
   TGraphAsymmErrors graph(100); 

   TH1D HTempDoubleLog("HTempDoubleLog", "HTempDoubleLog", 2 * BinCountDoubleLog, 0, 2 * BinCountDoubleLog);


   // loop over the points and fill the
   for(int i = 1; i < 101; i++)
   {
      int iGraph = i-1;
      // get the central value for the data-point in degrees
      double thetaDegrees = centralvals[iGraph][0]; 
      // convert theta to radians before looking up the binning
      double thetaRadians = thetaDegrees*(TMath::Pi()/180.0); 
      double thetaRadiansUp; 
      double thetaRadiansDown;

      // handle lower edge case 
      if(i == 1) thetaRadiansDown = 0.000002;
      else thetaRadiansDown = centralvals[iGraph-1][0]*(TMath::Pi()/180.0); 

      // handle upper edge case
      if(i==100) thetaRadiansUp = TMath::Pi() - 0.000002; 
      else thetaRadiansUp = centralvals[iGraph+1][0]*(TMath::Pi()/180.0);  

      int bin     = FindBin(thetaRadians, 2*BinCountDoubleLog, Bins);
      int binUp   = FindBin(thetaRadiansUp, 2*BinCountDoubleLog, Bins);
      int binDown = FindBin(thetaRadiansDown, 2*BinCountDoubleLog, Bins); 

      double binFraction = FindBinFraction(thetaRadians, 2*BinCountDoubleLog, Bins); 
      double binFractionUp = FindBinFraction(thetaRadiansUp, 2*BinCountDoubleLog, Bins); 
      double binFractionDown = FindBinFraction(thetaRadiansDown, 2*BinCountDoubleLog, Bins); 

      // set the point, not worrying about the errors. 
      graph.SetPoint(iGraph, HTempDoubleLog.GetBinLowEdge(bin)+binFraction, centralvals[iGraph][1]); 
      
      // calculate what the error should be
      double binWidthRadiansUp = (HTempDoubleLog.GetBinLowEdge(binUp)+binFractionUp) - (HTempDoubleLog.GetBinLowEdge(bin)+binFraction); 
      double binWidthRadiansDown =  (HTempDoubleLog.GetBinLowEdge(bin)+binFraction) - (HTempDoubleLog.GetBinLowEdge(binDown)+binFractionDown); 
      std::cout << "iGraph: " << iGraph << " with Lower Edge: " << HTempDoubleLog.GetBinLowEdge(bin)+binFraction - binWidthRadiansDown/2.0 << " and upper edge: " <<   HTempDoubleLog.GetBinLowEdge(bin)+binFraction + binWidthRadiansUp/2.0 << std::endl; 
      graph.SetPointError(iGraph, binWidthRadiansDown/2.0 , binWidthRadiansUp/2.0, errs[iGraph][1], errs[iGraph][1]); 

   }


   static vector<int> Colors = GetCVDColors6();
   graph.SetMarkerStyle(20);
   graph.SetMarkerColor(Colors[4]);
   graph.SetLineColor(Colors[4]);
   gStyle->SetEndErrorSize(4);
   graph.SetLineWidth(1.5);

   return graph; 
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

int FindBin(double Value, int NBins, double Bins[])
{
   for(int i = 0; i < NBins; i++)
      if(Value < Bins[i])
         return i - 1;
   return NBins;
}


void MakeCanvas(vector<TH1D *> Histograms, vector<string> Labels, string Output,
   string X, string Y, double WorldMin, double WorldMax, bool DoRatio, bool LogX)
{
   int NLine = Histograms.size();
   int N = Histograms[0]->GetNbinsX();

   double MarginL = 180;
   double MarginR = 90;
   double MarginB = 120;
   double MarginT = 90;

   double WorldXMin = LogX ? 0 : 0;
   double WorldXMax = LogX ? N : 180;
   
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
   // Canvas.SetLogy();
   // Canvas.SetRightMargin(MarginR);
   // Canvas.SetLeftMargin(MarginL);
   // Canvas.SetTopMargin(MarginT);
   // Canvas.SetBottomMargin(MarginB);

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
   for(TH1D *H : Histograms){
      H->Draw("same");
   }

   TGraphAsymmErrors opal; 
   if(LogX){
      opal = (TGraphAsymmErrors)getOPALGraph(); 
      opal.Draw("p same"); 
      opal.Draw("E same"); 
   }



   TGraph G;
   G.SetPoint(0, LogX ? N / 2 : 180.0/2, 0);
   G.SetPoint(1, LogX ? N / 2 : 180.0/2, 1000);
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
         TH1D *H = (TH1D *)Histograms[i]->Clone();
         H->Divide(Histograms[0]);
         H->Draw("same");
      }

      G.Draw("l");

      G2.SetPoint(0, 0, 1);
      G2.SetPoint(1, 99999, 1);
      G2.Draw("l");
   }
   
   double BinMin    = 0.002;
   double BinMiddle = 180 / 2;
   double BinMax    = 180 - 0.002;

   Canvas.cd();
   TGaxis X1(MarginL, MarginB, MarginL + PadWidth / 2, MarginB, BinMin, BinMiddle, 510, "GS");
   TGaxis X2(MarginL + PadWidth, MarginB, MarginL + PadWidth / 2, MarginB, BinMin, BinMiddle, 510, "-GS");
   TGaxis X3(MarginL, MarginB + PadRHeight, MarginL + PadWidth / 2, MarginB + PadRHeight, BinMin, BinMiddle, 510, "+-GS");
   TGaxis X4(MarginL + PadWidth, MarginB + PadRHeight, MarginL + PadWidth / 2, MarginB + PadRHeight, BinMin, BinMiddle, 510, "+-GS");
   TGaxis Y1(MarginL, MarginB, MarginL, MarginB + PadRHeight, WorldRMin, WorldRMax, 505, "");
   TGaxis Y2(MarginL, MarginB + PadRHeight, MarginL, MarginB + PadRHeight + PadHeight, WorldMin, WorldMax, 510, "G");
   
   TGaxis XL1(MarginL, MarginB, MarginL + PadWidth, MarginB, 0, 180, 510, "S");
   TGaxis XL2(MarginL, MarginB + PadRHeight, MarginL + PadWidth, MarginB + PadRHeight, 0, 180, 510, "+-S");

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
   if(LogX) Latex.DrawLatex(MarginL + (1 - MarginR - MarginL) * 0.115, MarginB - 0.01, "0.01");
   if(LogX) Latex.DrawLatex(MarginL + (1 - MarginR - MarginL) * 0.290, MarginB - 0.01, "0.1");
   if(LogX) Latex.DrawLatex(MarginL + (1 - MarginR - MarginL) * 0.465, MarginB - 0.01, "1");
   if(LogX) Latex.DrawLatex(MarginL + (1 - MarginR - MarginL) * 0.535, MarginB - 0.01, "#pi - 1");
   if(LogX) Latex.DrawLatex(MarginL + (1 - MarginR - MarginL) * 0.710, MarginB - 0.01, "#pi - 0.1");
   if(LogX) Latex.DrawLatex(MarginL + (1 - MarginR - MarginL) * 0.885, MarginB - 0.01, "#pi - 0.01");

   Latex.SetTextAlign(12);
   Latex.SetTextAngle(270);
   Latex.SetTextColor(kGray);
   if(!LogX)Latex.DrawLatex(MarginL + (1 - MarginR - MarginL) * 0.5 + 0.0175, 1 - MarginT - 0.015, "#theta_{L} = 180");
   else Latex.DrawLatex(MarginL + (1 - MarginR - MarginL) * 0.5 + 0.0175, 1 - MarginT - 0.015, "#theta_{L} = #pi/2");

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
   Latex.DrawLatex(0.01, 0.01, "Work-in-progress OPAL Check 12/12/2024");

   TLegend Legend(0.15, 0.90, 0.35, 0.90 - 0.06 * min(NLine, 4));
   Legend.SetTextFont(42);
   Legend.SetTextSize(0.035);
   Legend.SetFillStyle(0);
   Legend.SetBorderSize(0);
   for(int i = 0; i < NLine && i < 4; i++)
      Legend.AddEntry(Histograms[i], Labels[i].c_str(), "pl");

   if(LogX)Legend.AddEntry(&opal, "OPAL", "pl"); 
   Legend.Draw();

   TLegend Legend2(0.7, 0.90, 0.9, 0.90 - 0.06 * (NLine - 4));
   Legend2.SetTextFont(42);
   Legend2.SetTextSize(0.035);
   Legend2.SetFillStyle(0);
   Legend2.SetBorderSize(0);
   if(NLine >= 4)
   {
      for(int i = 4; i < NLine; i++)
         Legend2.AddEntry(Histograms[i], Labels[i].c_str(), "pl");
      Legend2.Draw();
   }

   Canvas.SaveAs((Output + ".pdf").c_str());
}

void SetPad(TPad &P)
{
   P.SetLeftMargin(0);
   P.SetTopMargin(0);
   P.SetRightMargin(0);
   P.SetBottomMargin(0);
   P.SetTickx();
   P.SetTicky();
   P.Draw();
}


double FindBinFraction(double Value, int NBins, double Bins[])
{
   int binIndex = FindBin(Value, NBins, Bins);

   // Handle out-of-range values
   if (binIndex < 0 || binIndex >= NBins)
      return -1.0; // Indicates out-of-range

   // Calculate the fraction within the bin
   double BinStart = Bins[binIndex];
   double BinEnd = Bins[binIndex + 1];
   return (Value - BinStart) / (BinEnd - BinStart);
}