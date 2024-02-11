#include <iostream>
using namespace std;

#include "TCanvas.h"
#include "TFile.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TGaxis.h"
#include "TPad.h"

#include "SetStyle.h"
#include "PlotHelper4.h"

int main(int argc, char *argv[]);
void MakePlot(vector<string> FileNames, vector<string> Labels, string Histogram, string Output,
   string X = "", string Y = "", double WorldMin = 1e-3, double WorldMax = 1);
void MakePlot2(vector<string> FileNames, vector<string> Labels, string Histogram1, string Histogram2,
   string Output, string X = "", string Y = "", double WorldMin = 1e-3, double WorldMax = 1);
void MakeCanvasOld(vector<TH1D *> Histograms, vector<string> Labels, string Output,
   string X = "", string Y = "", double WorldMin = 1e-3, double WorldMax = 1);
void MakeCanvas(vector<TH1D *> Histograms, vector<string> Labels, string Output,
   string X = "", string Y = "", double WorldMin = 1e-3, double WorldMax = 1);
double GetMin(TH1D *H);
void SetPad(TPad &P);

int main(int argc, char *argv[])
{
   SetThesisStyle();

   vector<string> FileNames{"PlotGen.root", "PlotReco.root", "PlotPythia8.root"};
   vector<string> Labels{"GEN Archived MC", "RECO Archived MC", "Pythia8"};
   
   MakePlot(FileNames, Labels, "HEEC2", "EEC2.pdf",
      "#theta_{L}", "#frac{1}{N_{event}} #frac{d(Sum E_{i}E_{j}/E^{2})}{d #theta_{L}}", 2e-3, 10);
   MakePlot(FileNames, Labels, "HEEC3", "EEC3.pdf",
      "#theta_{L}", "#frac{1}{N_{event}} #frac{d(Sum E_{i}E_{j}E_{k}/E^{3})}{d #theta_{L}}", 2e-6, 10);
   MakePlot(FileNames, Labels, "HEEC4", "EEC4.pdf",
      "#theta_{L}", "#frac{1}{N_{event}} #frac{d(Sum E_{i}E_{j}E_{k}E_{l}/E^{4})}{d #theta_{L}}", 2e-8, 10);
   
   MakePlot2(FileNames, Labels, "HEEC2", "HEEC3", "EEC3EEC2.pdf", "#theta_{L}", "Ratio E3C/E2C", 2e-4, 10);
   MakePlot2(FileNames, Labels, "HEEC2", "HEEC4", "EEC4EEC2.pdf", "#theta_{L}", "Ratio E4C/E2C", 2e-6, 10);

   return 0;
}

void MakePlot(vector<string> FileNames, vector<string> Labels, string Histogram, string Output,
   string X, string Y, double WorldMin, double WorldMax)
{
   static vector<int> Colors = GetCVDColors6();

   int N = FileNames.size();

   vector<TFile *> Files(N);
   for(int i = 0; i < N; i++)
      Files[i] = new TFile(FileNames[i].c_str());

   vector<TH1D *> Histograms(N);
   for(int i = 0; i < N; i++)
   {
      TH1D *HN = (TH1D *)Files[i]->Get("HN");
      TH1D *H = (TH1D *)Files[i]->Get(Histogram.c_str());
      H->Scale(1 / HN->GetBinContent(1));
      H->SetMarkerColor(Colors[i]);
      H->SetMarkerStyle(20);
      H->SetLineColor(Colors[i]);
      H->SetLineWidth(2);
      Histograms[i] = H;
   }

   MakeCanvas(Histograms, Labels, Output, X, Y, WorldMin, WorldMax);

   for(int i = 0; i < N; i++)
   {
      Files[i]->Close();
      delete Files[i];
   }
}

void MakePlot2(vector<string> FileNames, vector<string> Labels, string Histogram1, string Histogram2,
   string Output, string X, string Y, double WorldMin, double WorldMax)
{
   static vector<int> Colors = GetCVDColors6();
   
   int N = FileNames.size();

   vector<TFile *> Files(N);
   for(int i = 0; i < N; i++)
      Files[i] = new TFile(FileNames[i].c_str());
  
   vector<TH1D *> Histograms(N);
   for(int i = 0; i < N; i++)
   {
      TH1D *HN = (TH1D *)Files[i]->Get("HN");
      TH1D *H1 = (TH1D *)Files[i]->Get(Histogram1.c_str());
      TH1D *H2 = (TH1D *)Files[i]->Get(Histogram2.c_str());

      H1->Scale(1 / HN->GetBinContent(1));
      H2->Scale(1 / HN->GetBinContent(1));

      H1->SetMarkerColor(Colors[i]);
      H1->SetMarkerStyle(20);
      H1->SetLineColor(Colors[i]);
      H1->SetLineWidth(2);
      H2->SetMarkerColor(Colors[i]);
      H2->SetMarkerStyle(20);
      H2->SetLineColor(Colors[i]);
      H2->SetLineWidth(2);

      H2->Divide(H1);

      Histograms[i] = H2;
   }

   MakeCanvas(Histograms, Labels, Output, X, Y, WorldMin, WorldMax);

   for(int i = 0; i < N; i++)
   {
      Files[i]->Close();
      delete Files[i];
   }
}

void MakeCanvasOld(vector<TH1D *> Histograms, vector<string> Labels, string Output,
   string X, string Y, double WorldMin, double WorldMax)
{
   int NLine = Histograms.size();
   int N = Histograms[0]->GetNbinsX();

   double L = 0.12;
   double R = 0.06;
   double B = 0.12;
   double T = 0.09;

   TCanvas Canvas("Canvas", "", 1536, 1024);
   Canvas.SetLogy();
   Canvas.SetRightMargin(R);
   Canvas.SetLeftMargin(L);
   Canvas.SetTopMargin(T);
   Canvas.SetBottomMargin(B);

   TH2D HWorld("HWorld", "", N, 0, N, 100, WorldMin, WorldMax);
   HWorld.GetXaxis()->SetTitle(X.c_str());
   HWorld.GetYaxis()->SetTitle(Y.c_str());
   HWorld.SetStats(0);

   HWorld.GetXaxis()->SetTickLength(0);
   HWorld.GetXaxis()->SetLabelSize(0);
   HWorld.GetXaxis()->CenterTitle();
   HWorld.GetXaxis()->SetTitleOffset(1.4);
   HWorld.GetYaxis()->CenterTitle();
   HWorld.GetYaxis()->SetTitleOffset(1.5);

   HWorld.Draw();
   for(TH1D *H : Histograms)
      H->Draw("same");

   double BinMin = 0.002;
   double BinMiddle = M_PI / 2;
   double BinMax = M_PI - 0.002;

   TGaxis X1(0, WorldMin, N / 2, WorldMin, BinMin, BinMiddle, 510, "G");
   TGaxis X2(N, WorldMin, N / 2, WorldMin, BinMin, BinMiddle, 510, "-G");

   X1.SetLabelSize(0);
   X2.SetLabelSize(0);

   X1.Draw();
   X2.Draw();

   TLatex Latex;
   Latex.SetNDC();
   Latex.SetTextFont(42);
   Latex.SetTextSize(0.035);
   Latex.SetTextAlign(23);
   Latex.DrawLatex(L + (1 - R - L) * 0.115, B - 0.01, "0.01");
   Latex.DrawLatex(L + (1 - R - L) * 0.290, B - 0.01, "0.1");
   Latex.DrawLatex(L + (1 - R - L) * 0.465, B - 0.01, "1");
   Latex.DrawLatex(L + (1 - R - L) * 0.535, B - 0.01, "#pi - 1");
   Latex.DrawLatex(L + (1 - R - L) * 0.710, B - 0.01, "#pi - 0.1");
   Latex.DrawLatex(L + (1 - R - L) * 0.885, B - 0.01, "#pi - 0.01");

   Latex.SetTextAlign(12);
   Latex.SetTextAngle(270);
   Latex.SetTextColor(kGray);
   Latex.DrawLatex(L + (1 - R - L) * 0.5 + 0.015, 1 - T - 0.015, "#theta_{L} = #pi/2");

   TGraph G;
   G.SetPoint(0, N / 2, 0);
   G.SetPoint(1, N / 2, 1000);
   G.SetLineStyle(kDashed);
   G.SetLineColor(kGray);
   G.SetLineWidth(1);
   G.Draw("l");

   TLegend Legend(0.2, 0.85, 0.4, 0.85 - 0.06 * NLine);
   Legend.SetTextFont(42);
   Legend.SetTextSize(0.035);
   Legend.SetFillStyle(0);
   Legend.SetBorderSize(0);
   for(int i = 0; i < NLine; i++)
      Legend.AddEntry(Histograms[i], Labels[i].c_str(), "pl");
   Legend.Draw();

   Canvas.SaveAs(Output.c_str());
}

void MakeCanvas(vector<TH1D *> Histograms, vector<string> Labels, string Output,
   string X, string Y, double WorldMin, double WorldMax)
{
   int NLine = Histograms.size();
   int N = Histograms[0]->GetNbinsX();

   double MarginL = 180;
   double MarginR = 90;
   double MarginB = 120;
   double MarginT = 90;

   double PadWidth = 1200;
   double PadHeight = 640;
   double PadRHeight = 240;

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
   SetPad(PadR);

   Pad.cd();

   TH2D HWorld("HWorld", "", N, 0, N, 100, WorldMin, WorldMax);
   HWorld.SetStats(0);
   HWorld.GetXaxis()->SetTickLength(0);
   HWorld.GetXaxis()->SetLabelSize(0);

   HWorld.Draw("axis");
   for(TH1D *H : Histograms)
      H->Draw("same");

   TGraph G;
   G.SetPoint(0, N / 2, 0);
   G.SetPoint(1, N / 2, 1000);
   G.SetLineStyle(kDashed);
   G.SetLineColor(kGray);
   G.SetLineWidth(1);
   G.Draw("l");

   PadR.cd();

   double WorldRMin = 0.5;
   double WorldRMax = 1.5;
   
   TH2D HWorldR("HWorldR", "", N, 0, N, 100, WorldRMin, WorldRMax);
   HWorldR.SetStats(0);
   HWorldR.GetXaxis()->SetTickLength(0);
   HWorldR.GetXaxis()->SetLabelSize(0);
   HWorldR.GetYaxis()->SetNdivisions(503);

   HWorldR.Draw("axis");
   for(int i = 1; i < NLine; i++)
   {
      TH1D *H = (TH1D *)Histograms[i]->Clone();
      H->Divide(Histograms[0]);
      H->Draw("same");
   }

   G.Draw("l");

   TGraph G2;
   G2.SetPoint(0, 0, 1);
   G2.SetPoint(1, 99999, 1);
   G2.Draw("l");

   double BinMin = 0.002;
   double BinMiddle = M_PI / 2;
   double BinMax = M_PI - 0.002;

   Canvas.cd();
   TGaxis X1(MarginL, MarginB, MarginL + PadWidth / 2, MarginB, BinMin, BinMiddle, 510, "GS");
   TGaxis X2(MarginL + PadWidth, MarginB, MarginL + PadWidth / 2, MarginB, BinMin, BinMiddle, 510, "-GS");
   TGaxis X3(MarginL, MarginB + PadRHeight, MarginL + PadWidth / 2, MarginB + PadRHeight, BinMin, BinMiddle, 510, "+-GS");
   TGaxis X4(MarginL + PadWidth, MarginB + PadRHeight, MarginL + PadWidth / 2, MarginB + PadRHeight, BinMin, BinMiddle, 510, "+-GS");
   TGaxis Y1(MarginL, MarginB, MarginL, MarginB + PadRHeight, WorldRMin, WorldRMax, 503, "");
   TGaxis Y2(MarginL, MarginB + PadRHeight, MarginL, MarginB + PadRHeight + PadHeight, WorldMin, WorldMax, 510, "G");

   Y1.SetLabelFont(42);
   Y2.SetLabelFont(42);

   X1.SetLabelSize(0);
   X2.SetLabelSize(0);
   X3.SetLabelSize(0);
   X4.SetLabelSize(0);

   X1.SetTickSize(0.06);
   X2.SetTickSize(0.06);
   X3.SetTickSize(0.06);
   X4.SetTickSize(0.06);

   X1.Draw();
   X2.Draw();
   X3.Draw();
   X4.Draw();
   Y1.Draw();
   Y2.Draw();

   TLatex Latex;
   Latex.SetNDC();
   Latex.SetTextFont(42);
   Latex.SetTextSize(0.035);
   Latex.SetTextAlign(23);
   Latex.DrawLatex(MarginL + (1 - MarginR - MarginL) * 0.115, MarginB - 0.01, "0.01");
   Latex.DrawLatex(MarginL + (1 - MarginR - MarginL) * 0.290, MarginB - 0.01, "0.1");
   Latex.DrawLatex(MarginL + (1 - MarginR - MarginL) * 0.465, MarginB - 0.01, "1");
   Latex.DrawLatex(MarginL + (1 - MarginR - MarginL) * 0.535, MarginB - 0.01, "#pi - 1");
   Latex.DrawLatex(MarginL + (1 - MarginR - MarginL) * 0.710, MarginB - 0.01, "#pi - 0.1");
   Latex.DrawLatex(MarginL + (1 - MarginR - MarginL) * 0.885, MarginB - 0.01, "#pi - 0.01");

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
   Latex.DrawLatex(MarginL * 0.3, MarginB + PadRHeight * 0.5, "Ratio");
   Latex.DrawLatex(MarginL * 0.3, MarginB + PadRHeight + PadHeight * 0.5, Y.c_str());

   TLegend Legend(0.2, 0.85, 0.4, 0.85 - 0.06 * NLine);
   Legend.SetTextFont(42);
   Legend.SetTextSize(0.035);
   Legend.SetFillStyle(0);
   Legend.SetBorderSize(0);
   for(int i = 0; i < NLine; i++)
      Legend.AddEntry(Histograms[i], Labels[i].c_str(), "pl");
   Legend.Draw();

   Canvas.SaveAs(Output.c_str());
}

double GetMin(TH1D *H)
{
   double Min = -1;

   for(int i = 1; i <= H->GetNbinsX(); i++)
   {
      if(H->GetBinContent(i) == 0)
         continue;

      if(Min < 0 || Min > H->GetBinContent(i))
         Min = H->GetBinContent(i);
   }

   cout << "Min = " << Min << endl;

   return Min;
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

