#include <iostream>
using namespace std;

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"

#include "SetStyle.h"
#include "ProgressBar.h"
#include "CommandLine.h"
#include "Messenger.h"
#include "JetCorrector.h"
#include "alephTrkEfficiency.h"

int main(int argc, char *argv[]);
void DivideByBin(TH1D &H, double Bins[]);
int FindBin(double Value, int NBins, double Bins[]);
double GetMax(vector<double> X);

int main(int argc, char *argv[])
{
   SetThesisStyle();

   CommandLine CL(argc, argv);

   string InputFileName    = CL.Get("Input");
   string OutputFileName   = CL.Get("Output", "Plots.root");
   string ParticleTreeName = CL.Get("Particle", "t");
   double MinParticleE     = CL.GetDouble("MinParticleE", 0);
   double MinParticlePT    = CL.GetDouble("MinParticlePT", 0.2);
   double MinTheta         = CL.GetDouble("MinTheta", 0.35);   // cos(0.35) = 0.94
   bool IsReco             = CL.GetBool("IsReco", true);
   bool ChargedOnly        = CL.GetBool("ChargedOnly", true);
   bool DoEENormalize      = CL.GetBool("EENormalize", true);
   bool UseFullEnergy      = CL.GetBool("UseFullEnergy", true);
   bool DoWeight           = CL.GetBool("DoWeight", false);

   TFile OutputFile(OutputFileName.c_str(), "RECREATE");

   const int BinCount = 100;
   double Bins[2*BinCount+1];
   double BinMin = 0.002;
   double BinMax = M_PI / 2;
   for(int i = 0; i <= BinCount; i++)
   {
      Bins[i] = exp(log(BinMin) + (log(BinMax) - log(BinMin)) / BinCount * i);
      Bins[2*BinCount-i] = BinMax * 2 - exp(log(BinMin) + (log(BinMax) - log(BinMin)) / BinCount * i);
   }

   TH1D HN("HN", ";;", 1, 0, 1);
   TH1D HEEC2("HEEC2", ";EEC_{2};", 2 * BinCount, 0, 2 * BinCount);
   TH1D HEEC3("HEEC3", ";EEC_{3};", 2 * BinCount, 0, 2 * BinCount);
   TH1D HEEC4("HEEC4", ";EEC_{4};", 2 * BinCount, 0, 2 * BinCount);
   TH1D HEEC5("HEEC5", ";EEC_{5};", 2 * BinCount, 0, 2 * BinCount);
   TH1D HBinMin("HBinMin", ";EEC;BinMin", 2 * BinCount, 0, 2 * BinCount);
   TH1D HBinMax("HBinMax", ";EEC;BinMax", 2 * BinCount, 0, 2 * BinCount);

   for(int i = 0; i < 2 * BinCount; i++)
   {
      HBinMin.SetBinContent(i + 1, Bins[i]);
      HBinMax.SetBinContent(i + 1, Bins[i+1]);
   }

   HEEC2.SetStats(0);
   HEEC3.SetStats(0);
   HEEC4.SetStats(0);
   HEEC5.SetStats(0);

   TFile File(InputFileName.c_str());

   float NEvent = 0;

   ParticleTreeMessenger MParticle(File, ParticleTreeName.c_str());

   alephTrkEfficiency efficiencyCorrector;

   int EntryCount = MParticle.GetEntries();
   ProgressBar Bar(cout, EntryCount);
   for(int iE = 0; iE < EntryCount; iE++)
   {
      if(EntryCount < 300 || (iE % (EntryCount / 250) == 0))
      {
         Bar.Update(iE);
         Bar.Print();
      }

      MParticle.GetEntry(iE);

      if(IsReco == true && MParticle.PassBaselineCut() == false)
         continue;

      NEvent = NEvent + 1;

      // now we gather the particles
      vector<FourVector> P;
      vector<double> W;
      double TotalE = 0;
      for(int iP = 0; iP < MParticle.nParticle; iP++)
      {
         if(MParticle.P[iP][0] < MinParticleE)
            continue;
         if(MParticle.P[iP].GetPT() < MinParticlePT)
            continue;
         if(ChargedOnly == true && (MParticle.charge[iP] == 0 && MParticle.isCharged[iP] == 0))
            continue;
         if(MParticle.P[iP].GetTheta() < MinTheta || MParticle.P[iP].GetTheta() > M_PI - MinTheta)
            continue;

         FourVector &Momentum = MParticle.P[iP];
         P.push_back(Momentum);
         
         if(DoWeight == true)
         {
            double Efficiency = efficiencyCorrector.efficiency(Momentum.GetTheta(), Momentum.GetPhi(), Momentum.GetPT(), MParticle.nChargedHadronsHP);
            W.push_back((Efficiency > 0) ? (1 / Efficiency) : 0);
         }
         else
            W.push_back(1);
         
         // W.push_back(MParticle.weight[iP]);
         
         TotalE = TotalE + MParticle.P[iP][0];
      }

      if(UseFullEnergy == true)
         TotalE = 91.1876;

      double TotalE2 = DoEENormalize ? (TotalE  * TotalE) : 1;
      double TotalE3 = DoEENormalize ? (TotalE2 * TotalE) : 1;
      double TotalE4 = DoEENormalize ? (TotalE3 * TotalE) : 1;
      double TotalE5 = DoEENormalize ? (TotalE4 * TotalE) : 1;
   
      // Fill EECs
      int N = P.size();
      vector<vector<double>> D(N);
      for(int i = 0; i < N; i++)
      {
         D[i].resize(N);
         for(int j = 0; j < N; j++)
            D[i][j] = GetAngle(P[i], P[j]);
      }

      for(int i1 = 0; i1 < N; i1++)
      {
         for(int i2 = i1 + 1; i2 < N; i2++)
         {
            double Max2 = D[i1][i2];
            int Bin2 = FindBin(Max2, BinCount * 2, Bins);
            HEEC2.Fill(Bin2, P[i1][0] * P[i2][0] / TotalE2 * W[i1] * W[i2]);

            for(int i3 = i2 + 1; i3 < N; i3++)
            {
               double Max3 = GetMax({Max2, D[i1][i3], D[i2][i3]});
               int Bin3 = FindBin(Max3, BinCount * 2, Bins);
               HEEC3.Fill(Bin3, P[i1][0] * P[i2][0] * P[i3][0] / TotalE3 * W[i1] * W[i2] * W[i3]);
           
               for(int i4 = i3 + 1; i4 < N; i4++)
               {
                  double Max4 = GetMax({Max3, D[i1][i4], D[i2][i4], D[i3][i4]});
                  int Bin4 = FindBin(Max4, BinCount * 2, Bins);
                  HEEC4.Fill(Bin4, P[i1][0] * P[i2][0] * P[i3][0] * P[i4][0] / TotalE4 * W[i1] * W[i2] * W[i3] * W[i4]);

                  /*
                  for(int i5 = i4 + 1; i5 < N; i5++)
                  {
                     double Max5 = GetMax({Max4, D[i1][i5], D[i2][i5], D[i3][i5], D[i4][i5]});
                     int Bin5 = FindBin(Max5, BinCount * 2, Bins);
                     HEEC5.Fill(Bin5, P[i1][0] * P[i2][0] * P[i3][0] * P[i4][0] * P[i5][0] / TotalE5 * W[i1] * W[i2] * W[i3] * W[i4] * W[i5]);
                  }
                  */
               }
            }
         }
      }
   }
   Bar.Update(EntryCount);
   Bar.Print();
   Bar.PrintLine();

   File.Close();

   HN.SetBinContent(1, NEvent);
   DivideByBin(HEEC2, Bins);
   DivideByBin(HEEC3, Bins);
   DivideByBin(HEEC4, Bins);
   DivideByBin(HEEC5, Bins);
  
   OutputFile.cd();

   HN.Write();
   HEEC2.Write();
   HEEC3.Write();
   HEEC4.Write();
   HEEC5.Write();
   HBinMin.Write();
   HBinMax.Write();

   OutputFile.Close();

   return 0;
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

double GetMax(vector<double> X)
{
   if(X.size() == 0)
      return -1;

   double Result = X[0];
   for(int i = 1; i < X.size(); i++)
      if(Result < X[i])
         Result = X[i];
   return Result;
}
