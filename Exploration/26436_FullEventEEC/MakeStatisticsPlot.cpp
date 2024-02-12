#include <iostream>
using namespace std;

#include "TCanvas.h"
#include "TLatex.h"
#include "TGaxis.h"
#include "TFile.h"
#include "TH1D.h"

int main(int argc, char *argv[]);
void MakePlot(string FileName, string Histogram, string Output);

int main(int argc, char *argv[])
{
   MakePlot("PlotData.root", "HEEC2", "DataStatEEC2.pdf");
   MakePlot("PlotData.root", "HEEC3", "DataStatEEC3.pdf");
   MakePlot("PlotData.root", "HEEC4", "DataStatEEC4.pdf");

   return 0;
}

void MakePlot(string FileName, string Histogram, string Output)
{
   double MarginL = 180;
   double MarginR = 90;
   double MarginB = 120;
   double MarginT = 90;

   double PadWidth = 1200;
   double PadHeight = 640 + 240;

   double CanvasWidth = MarginL + PadWidth + MarginR;
   double CanvasHeight = MarginT + PadHeight + MarginB;

   MarginL = MarginL / CanvasWidth;
   MarginR = MarginR / CanvasWidth;
   MarginT = MarginT / CanvasHeight;
   MarginB = MarginB / CanvasHeight;

   PadWidth   = PadWidth / CanvasWidth;
   PadHeight  = PadHeight / CanvasHeight;

   TFile File(FileName.c_str());

   TH1D *HN = (TH1D *)File.Get("HN");
   TH1D *H = (TH1D *)File.Get(Histogram.c_str());

   H->Scale(1 / HN->GetBinContent(1));

   for(int i = 1; i <= H->GetNbinsX(); i++)
      H->SetBinContent(i, 0);

   int N = H->GetNbinsX();

   TCanvas Canvas("Canvas", "", CanvasWidth, CanvasHeight);
   Canvas.SetRightMargin(MarginR);
   Canvas.SetLeftMargin(MarginL);
   Canvas.SetTopMargin(MarginT);
   Canvas.SetBottomMargin(MarginB);

   H->SetStats(0);
   H->Draw();
   H->GetYaxis()->SetRangeUser(-0.002, 0.002);

   H->GetXaxis()->SetTickLength(0);
   H->GetXaxis()->SetLabelSize(0);
   H->GetXaxis()->CenterTitle();
   H->GetXaxis()->SetTitleOffset(1.4);
   H->GetYaxis()->CenterTitle();
   H->GetYaxis()->SetTitleOffset(1.5);

   H->GetXaxis()->SetTitle("#theta_{L}");
   H->GetYaxis()->SetTitle("Raw stat uncertainty");

   TGaxis X1(0, -0.002, N / 2, -0.002, 0.002, M_PI / 2, 510, "G");
   TGaxis X2(N, -0.002, N / 2, -0.002, 0.002, M_PI / 2, 510, "-G");

   X1.SetLabelSize(0);
   X2.SetLabelSize(0);

   X1.Draw();
   X2.Draw();

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

   Canvas.SaveAs(Output.c_str());

   File.Close();
}




