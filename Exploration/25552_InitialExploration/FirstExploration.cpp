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

int main(int argc, char *argv[]);
void DivideByBin(TH1D &H);
double GetMax(vector<double> X);

int main(int argc, char *argv[])
{
   SetThesisStyle();

   CommandLine CL(argc, argv);

   string InputFileName    = CL.Get("Input");
   string OutputFileName   = CL.Get("Output", "Plots.root");
   string ParticleTreeName = CL.Get("Particle", "t");
   string JetTreeName      = CL.Get("Jet", "akR4ESchemeJetTree");
   double JetR             = CL.GetDouble("JetR", 0.4);
   double MinJetE          = CL.GetDouble("MinJetE", 20);
   double MinParticleE     = CL.GetDouble("MinParticleE", 0.5);
   bool IsReco             = CL.GetBool("IsReco", true);
   vector<string> JECFiles = CL.GetStringVector("JEC");

   TFile OutputFile(OutputFileName.c_str(), "RECREATE");

   const int BinCount = 100;
   double Bins[BinCount+1];
   double BinMin = 0.0005;
   double BinMax = 1;
   for(int i = 0; i <= BinCount; i++)
      Bins[i] = exp(log(BinMin) + (log(BinMax) - log(BinMin)) / BinCount * i);

   double NJet = 0;
   TH1D HEEC2("HEEC2", ";EEC_{2};", 100, Bins);
   TH1D HEEC3("HEEC3", ";EEC_{3};", 100, Bins);
   TH1D HEEC4("HEEC4", ";EEC_{4};", 100, Bins);
   TH1D HEEC5("HEEC5", ";EEC_{5};", 100, Bins);

   HEEC2.SetStats(0);
   HEEC3.SetStats(0);
   HEEC4.SetStats(0);
   HEEC5.SetStats(0);

   JetCorrector JEC(JECFiles);

   TFile File(InputFileName.c_str());

   ParticleTreeMessenger MParticle(File, ParticleTreeName.c_str());
   JetTreeMessenger      MJet(File, JetTreeName.c_str());

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
      MJet.GetEntry(iE);

      if(IsReco == true && MParticle.PassBaselineCut() == false)
         continue;

      for(int iJ = 0; iJ < MJet.nref; iJ++)
      {
         JEC.SetJetP(MJet.Jet[iJ].GetP());
         JEC.SetJetE(MJet.Jet[iJ][0]);
         JEC.SetJetTheta(MJet.Jet[iJ].GetTheta());
         JEC.SetJetPhi(MJet.Jet[iJ].GetPhi());
         double Correction = JEC.GetCorrection();
         if(Correction > 0)
            MJet.Jet[iJ] = MJet.Jet[iJ] * Correction;
      }

      for(int iJ = 0; iJ < MJet.nref; iJ++)
      {
         // too close to beam pipe
         if(MJet.Jet[iJ].GetTheta() < 0.2 * M_PI || MJet.Jet[iJ].GetTheta() > 0.8 * M_PI)
            continue;

         // too low jet energy
         if(MJet.Jet[iJ][0] < MinJetE)
            continue;

         // now we gather the particles around the jets
         vector<FourVector> P;
         vector<double> W;
         for(int iP = 0; iP < MParticle.nParticle; iP++)
         {
            if(GetAngle(MJet.Jet[iJ], MParticle.P[iP]) > JetR)
               continue;
            if(MParticle.P[iP][0] < MinParticleE)
               continue;
            P.push_back(MParticle.P[iP]);
            W.push_back(MParticle.weight[iP]);
         }

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
               HEEC2.Fill(Max2, P[i1][0] * P[i2][0] / MJet.Jet[iJ][0] / MJet.Jet[iJ][0]);

               for(int i3 = i2 + 1; i3 < N; i3++)
               {
                  double Max3 = GetMax({Max2, D[i1][i3], D[i2][i3]});
                  HEEC3.Fill(Max3, P[i1][0] * P[i2][0] * P[i3][0] / MJet.Jet[iJ][0] / MJet.Jet[iJ][0] / MJet.Jet[iJ][0]);
               
                  for(int i4 = i3 + 1; i4 < N; i4++)
                  {
                     double Max4 = GetMax({Max3, D[i1][i4], D[i2][i4], D[i3][i4]});
                     HEEC4.Fill(Max4, P[i1][0] * P[i2][0] * P[i3][0] * P[i4][0] / MJet.Jet[iJ][0] / MJet.Jet[iJ][0] / MJet.Jet[iJ][0] / MJet.Jet[iJ][0]);
                  
                     for(int i5 = i4 + 1; i5 < N; i5++)
                     {
                        double Max5 = GetMax({Max4, D[i1][i5], D[i2][i5], D[i3][i5], D[i4][i5]});
                        HEEC5.Fill(Max5, P[i1][0] * P[i2][0] * P[i3][0] * P[i4][0] * P[i5][0] / MJet.Jet[iJ][0] / MJet.Jet[iJ][0] / MJet.Jet[iJ][0] / MJet.Jet[iJ][0] / MJet.Jet[iJ][0]);
                     }
                  }
               }
            }
         }

         NJet = NJet + 1;
      }
   }
   Bar.Update(EntryCount);
   Bar.Print();
   Bar.PrintLine();

   File.Close();

   HEEC2.Scale(1 / NJet);   DivideByBin(HEEC2);
   HEEC3.Scale(1 / NJet);   DivideByBin(HEEC3);
   HEEC4.Scale(1 / NJet);   DivideByBin(HEEC4);
   HEEC5.Scale(1 / NJet);   DivideByBin(HEEC5);
  
   OutputFile.cd();

   HEEC2.Write();
   HEEC3.Write();
   HEEC4.Write();
   HEEC5.Write();

   OutputFile.Close();

   return 0;
}

void DivideByBin(TH1D &H)
{
   int N = H.GetNbinsX();
   for(int i = 1; i <= N; i++)
   {
      double L = H.GetXaxis()->GetBinLowEdge(i);
      double R = H.GetXaxis()->GetBinUpEdge(i);
      H.SetBinContent(i, H.GetBinContent(i) / (R - L));
      H.SetBinError(i, H.GetBinError(i) / (R - L));
   }
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
