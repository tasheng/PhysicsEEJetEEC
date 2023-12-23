#include <iostream>
using namespace std;

#include "TCanvas.h"
#include "TFile.h"
#include "TH1D.h"
#include "TLegend.h"

#include "SetStyle.h"
#include "PlotHelper4.h"

int main(int argc, char *argv[]);
void MakePlot(string GenFile, string RecoFile, string DataFile, string Histogram, string Output);

int main(int argc, char *argv[])
{
   SetThesisStyle();

   MakePlot("PlotGen.root", "PlotReco.root", "PlotData.root", "HEEC2", "EEC2.pdf");
   MakePlot("PlotGen.root", "PlotReco.root", "PlotData.root", "HEEC3", "EEC3.pdf");
   MakePlot("PlotGen.root", "PlotReco.root", "PlotData.root", "HEEC4", "EEC4.pdf");
   MakePlot("PlotGen.root", "PlotReco.root", "PlotData.root", "HEEC5", "EEC5.pdf");

   return 0;
}

void MakePlot(string GenFile, string RecoFile, string DataFile, string Histogram, string Output)
{
   TFile F1(GenFile.c_str());
   TFile F2(RecoFile.c_str());
   TFile F3(DataFile.c_str());

   TH1D *H1 = (TH1D *)F1.Get(Histogram.c_str());
   TH1D *H2 = (TH1D *)F2.Get(Histogram.c_str());
   TH1D *H3 = (TH1D *)F3.Get(Histogram.c_str());

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

   H1->Draw();
   H2->Draw("same");
   H3->Draw("same");

   TLegend Legend(0.2, 0.2, 0.4, 0.4);
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
