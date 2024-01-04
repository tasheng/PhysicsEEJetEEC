#include <iostream>
using namespace std;

#include "TCanvas.h"
#include "TFile.h"
#include "TH1D.h"
#include "TLegend.h"

#include "SetStyle.h"
#include "PlotHelper4.h"

int main(int argc, char *argv[]);
void MakePlot(string GenFile, string RecoFile, string DataFile, string Histogram, string Output,
   string X = "", string Y = "");
void MakePlot2(string GenFile, string RecoFile, string DataFile, string Histogram1, string Histogram2,
   string Output);

int main(int argc, char *argv[])
{
   SetThesisStyle();

   string GenFile = "PlotGen.root";
   string RecoFile = "PlotReco.root";
   string DataFile = "PlotData.root";
   
   MakePlot(GenFile, RecoFile, DataFile, "HEEC2", "EEC2.pdf",
      "R_{L}", "#frac{1}{N_{event}} #frac{d(Sum E_{i}E_{j}/E^{2})}{d R_{L}}");
   MakePlot(GenFile, RecoFile, DataFile, "HEEC3", "EEC3.pdf",
      "R_{L}", "#frac{1}{N_{event}} #frac{d(Sum E_{i}E_{j}E_{k}/E^{3})}{d R_{L}}");
   MakePlot(GenFile, RecoFile, DataFile, "HEEC4", "EEC4.pdf",
      "R_{L}", "#frac{1}{N_{event}} #frac{d(Sum E_{i}E_{j}E_{k}E_{l}/E^{4})}{d R_{L}}");
   MakePlot(GenFile, RecoFile, DataFile, "HEEC5", "EEC5.pdf",
      "R_{L}", "#frac{1}{N_{event}} #frac{d(Sum E_{i}E_{j}E_{k}E_{l}E_{m}/E^{5})}{d R_{L}}");
   MakePlot(GenFile, RecoFile, DataFile, "HEEC2O", "EEC2O.pdf",
      "#pi - R_{L}", "#frac{1}{N_{event}} #frac{d(Sum E_{i}E_{j}/E^{2})}{d R_{L}}");
   MakePlot(GenFile, RecoFile, DataFile, "HEEC3O", "EEC3O.pdf",
      "#pi - R_{L}", "#frac{1}{N_{event}} #frac{d(Sum E_{i}E_{j}E_{k}/E^{3})}{d R_{L}}");
   MakePlot(GenFile, RecoFile, DataFile, "HEEC4O", "EEC4O.pdf",
      "#pi - R_{L}", "#frac{1}{N_{event}} #frac{d(Sum E_{i}E_{j}E_{k}E_{l}/E^{4})}{d R_{L}}");
   MakePlot(GenFile, RecoFile, DataFile, "HEEC5O", "EEC5O.pdf",
      "#pi - R_{L}", "#frac{1}{N_{event}} #frac{d(Sum E_{i}E_{j}E_{k}E_{l}E_{m}/E^{5})}{d R_{L}}");
   
   MakePlot2(GenFile, RecoFile, DataFile, "HEEC2", "HEEC3", "EEC3EEC2.pdf");
   MakePlot2(GenFile, RecoFile, DataFile, "HEEC2", "HEEC4", "EEC4EEC2.pdf");
   MakePlot2(GenFile, RecoFile, DataFile, "HEEC2", "HEEC5", "EEC5EEC2.pdf");
   MakePlot2(GenFile, RecoFile, DataFile, "HEEC2O", "HEEC3O", "EEC3OEEC2O.pdf");
   MakePlot2(GenFile, RecoFile, DataFile, "HEEC2O", "HEEC4O", "EEC4OEEC2O.pdf");
   MakePlot2(GenFile, RecoFile, DataFile, "HEEC2O", "HEEC5O", "EEC5OEEC2O.pdf");

   return 0;
}

void MakePlot(string GenFile, string RecoFile, string DataFile, string Histogram, string Output, string X, string Y)
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

   TCanvas Canvas;
   Canvas.SetLogx();
   Canvas.SetLogy();

   if(X != "")
      H1->GetXaxis()->SetTitle(X.c_str());
   if(Y != "")
      H1->GetYaxis()->SetTitle(Y.c_str());

   H1->Draw();
   H2->Draw("same");
   H3->Draw("same");

   TLegend Legend(0.2, 0.6, 0.4, 0.8);
   Legend.SetTextFont(42);
   Legend.SetTextSize(0.035);
   Legend.SetFillStyle(0);
   Legend.SetBorderSize(0);
   Legend.AddEntry(H1, "GEN", "pl");
   Legend.AddEntry(H2, "RECO", "pl");
   Legend.AddEntry(H3, "Data", "pl");
   Legend.Draw();

   Canvas.SaveAs(Output.c_str());

   F3.Close();
   F2.Close();
   F1.Close();
}

void MakePlot2(string GenFile, string RecoFile, string DataFile, string Histogram1, string Histogram2, string Output)
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

   TCanvas Canvas;
   Canvas.SetLogx();
   Canvas.SetLogy();

   H12->GetXaxis()->SetTitle("R_{L}");

   H12->Draw();
   H22->Draw("same");
   H32->Draw("same");

   TLegend Legend(0.6, 0.2, 0.8, 0.4);
   Legend.SetTextFont(42);
   Legend.SetTextSize(0.035);
   Legend.SetFillStyle(0);
   Legend.SetBorderSize(0);
   Legend.AddEntry(H12, "GEN", "pl");
   Legend.AddEntry(H22, "RECO", "pl");
   Legend.AddEntry(H32, "Data", "pl");
   Legend.Draw();

   Canvas.SaveAs(Output.c_str());

   F3.Close();
   F2.Close();
   F1.Close();
}
