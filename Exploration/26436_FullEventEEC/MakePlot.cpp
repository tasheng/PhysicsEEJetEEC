#include <iostream>
using namespace std;

#include "TCanvas.h"
#include "TFile.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TGaxis.h"

#include "SetStyle.h"
#include "PlotHelper4.h"

int main(int argc, char *argv[]);
void MakePlot(string GenFile, string RecoFile, string DataFile, string Histogram, string Output,
   string X = "", string Y = "", double WorldMin = 1e-3, double WorldMax = 1);
void MakePlot2(string GenFile, string RecoFile, string DataFile, string Histogram1, string Histogram2,
   string Output, double WorldMin = 1e-3, double WorldMax = 1);
void MakeCanvas(TH1D *H1, TH1D *H2, TH1D *H3, TH1D *HMin, TH1D *HMax, string Output, string X = "", string Y = "", double WorldMin = 1e-3, double WorldMax = 1);
double GetMin(TH1D *H);

int main(int argc, char *argv[])
{
   SetThesisStyle();

   string GenFile = "PlotGen_01.root";
   string RecoFile = "PlotReco_01.root";
   string DataFile = "PlotPythia8.root";
   
   MakePlot(GenFile, RecoFile, DataFile, "HEEC2", "EEC2.pdf",
      "#theta_{L}", "#frac{1}{N_{event}} #frac{d(Sum E_{i}E_{j}/E^{2})}{d #theta_{L}}", 1e-3, 10);
   MakePlot(GenFile, RecoFile, DataFile, "HEEC3", "EEC3.pdf",
      "#theta_{L}", "#frac{1}{N_{event}} #frac{d(Sum E_{i}E_{j}E_{k}/E^{3})}{d #theta_{L}}", 1e-6, 10);
   MakePlot(GenFile, RecoFile, DataFile, "HEEC4", "EEC4.pdf",
      "#theta_{L}", "#frac{1}{N_{event}} #frac{d(Sum E_{i}E_{j}E_{k}E_{l}/E^{4})}{d #theta_{L}}", 1e-8, 10);
   
   MakePlot2(GenFile, RecoFile, DataFile, "HEEC2", "HEEC3", "EEC3EEC2.pdf", 1e-4, 10);
   MakePlot2(GenFile, RecoFile, DataFile, "HEEC2", "HEEC4", "EEC4EEC2.pdf", 1e-6, 10);

   return 0;
}

void MakePlot(string GenFile, string RecoFile, string DataFile, string Histogram, string Output, string X, string Y, double WorldMin, double WorldMax)
{
   TFile F1(GenFile.c_str());
   TFile F2(RecoFile.c_str());
   TFile F3(DataFile.c_str());

   TH1D *N1 = (TH1D *)F1.Get("HN");
   TH1D *N2 = (TH1D *)F2.Get("HN");
   TH1D *N3 = (TH1D *)F3.Get("HN");

   TH1D *H1 = (TH1D *)F1.Get(Histogram.c_str());
   TH1D *H2 = (TH1D *)F2.Get(Histogram.c_str());
   TH1D *H3 = (TH1D *)F3.Get(Histogram.c_str());
   
   TH1D *HMin = (TH1D *)F1.Get("HBinMin");
   TH1D *HMax = (TH1D *)F1.Get("HBinMax");

   H1->Scale(1 / N1->GetBinContent(1));
   H2->Scale(1 / N2->GetBinContent(1));
   H3->Scale(1 / N3->GetBinContent(1));

   H1->SetMarkerColor(kGreen + 3);
   H1->SetMarkerStyle(20);
   H1->SetLineColor(kGreen + 3);
   H1->SetLineWidth(2);
   H2->SetMarkerColor(kRed);
   H2->SetMarkerStyle(20);
   H2->SetLineColor(kRed);
   H2->SetLineWidth(2);
   H3->SetMarkerColor(kBlue);
   H3->SetMarkerStyle(20);
   H3->SetLineColor(kBlue);
   H3->SetLineWidth(2);

   MakeCanvas(H1, H2, H3, HMin, HMax, Output, X, Y, WorldMin, WorldMax);

   F3.Close();
   F2.Close();
   F1.Close();
}

void MakePlot2(string GenFile, string RecoFile, string DataFile, string Histogram1, string Histogram2, string Output, double WorldMin, double WorldMax)
{
   TFile F1(GenFile.c_str());
   TFile F2(RecoFile.c_str());
   TFile F3(DataFile.c_str());
   
   TH1D *N1 = (TH1D *)F1.Get("HN");
   TH1D *N2 = (TH1D *)F2.Get("HN");
   TH1D *N3 = (TH1D *)F3.Get("HN");

   TH1D *H11 = (TH1D *)F1.Get(Histogram1.c_str());
   TH1D *H21 = (TH1D *)F2.Get(Histogram1.c_str());
   TH1D *H31 = (TH1D *)F3.Get(Histogram1.c_str());
   TH1D *H12 = (TH1D *)F1.Get(Histogram2.c_str());
   TH1D *H22 = (TH1D *)F2.Get(Histogram2.c_str());
   TH1D *H32 = (TH1D *)F3.Get(Histogram2.c_str());
   
   TH1D *HMin = (TH1D *)F1.Get("HBinMin");
   TH1D *HMax = (TH1D *)F1.Get("HBinMax");

   H11->Scale(1 / N1->GetBinContent(1));
   H21->Scale(1 / N2->GetBinContent(1));
   H31->Scale(1 / N3->GetBinContent(1));
   H12->Scale(1 / N1->GetBinContent(1));
   H22->Scale(1 / N2->GetBinContent(1));
   H32->Scale(1 / N3->GetBinContent(1));

   H11->SetMarkerColor(kGreen + 3);
   H11->SetMarkerStyle(20);
   H11->SetLineColor(kGreen + 3);
   H11->SetLineWidth(2);
   H21->SetMarkerColor(kRed);
   H21->SetMarkerStyle(20);
   H21->SetLineColor(kRed);
   H21->SetLineWidth(2);
   H31->SetMarkerColor(kBlue);
   H31->SetMarkerStyle(20);
   H31->SetLineColor(kBlue);
   H31->SetLineWidth(2);
   H12->SetMarkerColor(kGreen + 3);
   H12->SetMarkerStyle(20);
   H12->SetLineColor(kGreen + 3);
   H12->SetLineWidth(2);
   H22->SetMarkerColor(kRed);
   H22->SetMarkerStyle(20);
   H22->SetLineColor(kRed);
   H22->SetLineWidth(2);
   H32->SetMarkerColor(kBlue);
   H32->SetMarkerStyle(20);
   H32->SetLineColor(kBlue);
   H32->SetLineWidth(2);

   H12->Divide(H11);
   H22->Divide(H21);
   H32->Divide(H31);

   MakeCanvas(H12, H22, H32, HMin, HMax, Output, "#theta_{L}", Histogram2 + "/" + Histogram1, WorldMin, WorldMax);

   F3.Close();
   F2.Close();
   F1.Close();
}

void MakeCanvas(TH1D *H1, TH1D *H2, TH1D *H3, TH1D *HMin, TH1D *HMax, string Output, string X, string Y, double WorldMin, double WorldMax)
{
   int N = H1->GetNbinsX();

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
   H1->Draw("same");
   H2->Draw("same");
   H3->Draw("same");

   // double BinMin = HMin->GetBinContent(1) / 1;
   // double BinMiddle = M_PI / 2;
   // double BinMax = HMax->GetBinContent(N) / 1;
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

   TLegend Legend(0.2, 0.6, 0.4, 0.8);
   Legend.SetTextFont(42);
   Legend.SetTextSize(0.035);
   Legend.SetFillStyle(0);
   Legend.SetBorderSize(0);
   Legend.AddEntry(H1, "GEN", "pl");
   Legend.AddEntry(H2, "RECO", "pl");
   Legend.AddEntry(H3, "Pythia8", "pl");
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


