#include <iostream>
#include <vector>
using namespace std;

#include "TFile.h"

#include "CommandLine.h"
#include "Messenger.h"
#include "alephTrkEfficiency.h"

#define MAX 10000

int main(int argc, char *argv[])
{
   CommandLine CL(argc, argv);

   string InputFileName  = CL.Get("Input");
   string OutputFileName = CL.Get("Output");
   string TreeName       = CL.Get("Tree", "t");
   bool GenLevel         = CL.GetBool("GenLevel", false);
   double MinTheta       = CL.GetDouble("MinTheta", 0.35);
   double Fraction       = CL.GetDouble("Fraction", 1.00);

   TFile InputFile(InputFileName.c_str());
   ParticleTreeMessenger M(InputFile, TreeName);

   cout << M.GetEntries() << endl;

   TFile OutputFile(OutputFileName.c_str(), "RECREATE");
   TTree OutputTree("Tree", "Reduced tree");

   int N;
   float Momentum[MAX], Mass[MAX], Theta[MAX], Phi[MAX], Weight[MAX];
   short Charge[MAX];
   bool PassCut;
   OutputTree.Branch("N", &N, "N/I");
   OutputTree.Branch("Momentum", &Momentum, "Momentum[N]/F");
   OutputTree.Branch("Mass", &Mass, "Mass[N]/F");
   OutputTree.Branch("Theta", &Theta, "Theta[N]/F");
   OutputTree.Branch("Phi", &Phi, "Phi[N]/F");
   OutputTree.Branch("Weight", &Weight, "Weight[N]/F");
   OutputTree.Branch("Charge", &Charge, "Charge[N]/S");
   OutputTree.Branch("PassCut", &PassCut, "PassCut/O");

   alephTrkEfficiency efficiencyCorrector;

   int EntryCount = M.GetEntries() * Fraction;
   for(int iE = 0; iE < EntryCount; iE++)
   {
      M.GetEntry(iE);

      // Event selection here
      PassCut = M.PassBaselineCut();

      N = 0;
      for(int iP = 0; iP < M.nParticle; iP++)
      {
         // Particle selection here
         if(GenLevel == false)
         {
            if(M.highPurity[iP] == false)
               continue;
            if(M.charge[iP] == 0)
               continue;
            if(M.P[iP].GetTheta() < MinTheta || M.P[iP].GetTheta() > M_PI - MinTheta)
               continue;
         }
         else
         {
            if(M.charge[iP] == 0)
               continue;
            if(M.P[iP].GetTheta() < MinTheta || M.P[iP].GetTheta() > M_PI - MinTheta)
               continue;
         }

         double Efficiency = 1;
         if(GenLevel == false)
            efficiencyCorrector.efficiency(M.P[iP].GetTheta(), M.P[iP].GetPhi(), M.P[iP].GetPT(), M.nChargedHadronsHP);

         Momentum[N] = M.P[iP].GetP();
         Mass[N]     = M.P[iP].GetMass();
         Theta[N]    = M.P[iP].GetTheta();
         Phi[N]      = M.P[iP].GetPhi();
         Weight[N]   = (Efficiency > 0) ? (1 / Efficiency) : 0;
         Charge[N]   = M.charge[iP];
         N = N + 1;
      }

      OutputTree.Fill();
   }

   OutputFile.cd();
   OutputTree.Write();
   OutputFile.Close();

   InputFile.Close();

   return 0;
}






