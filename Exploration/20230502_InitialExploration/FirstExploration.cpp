#include <iostream>
using namespace std;

#include "TFile.h"
#include "TTree.h"

#include "ProgressBar.h"
#include "CommandLine.h"
#include "Messenger.h"
#include "JetCorrector.h"

int main(int argc, char *argv[])
{
   CommandLine CL(argc, argv);

   string InputFileName    = CL.Get("Input");
   string JetTreeName      = CL.Get("Jet", "akR4ESchemeJetTree");
   double JetR             = CL.GetDouble("JetR", 0.4);
   vector<string> JECFiles = CL.GetStringVector("JEC");

   JetCorrector JEC(JECFiles);

   TFile File(InputFileName.c_str());

   ParticleTreeMessenger MParticle(File, "t");
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

      if(MParticle.PassBaselineCut() == false)
         continue;

      for(int iJ = 0; iJ < MJet.nref; iJ++)
      {
         JEC.SetJetP(MJet.Jet[iJ].GetP());
         JEC.SetJetE(MJet.Jet[iJ][0]);
         JEC.SetJetTheta(MJet.Jet[iJ].GetTheta());
         JEC.SetJetPhi(MJet.Jet[iJ].GetPhi());
         double Correction = JEC.GetCorrection();
         MJet.Jet[iJ] = MJet.Jet[iJ] * Correction;
      }

      for(int iJ = 0; iJ < MJet.nref; iJ++)
      {
         // too close to beam pipe
         if(MJet.Jet[iJ].GetTheta() < 0.2 * M_PI || MJet.Jet[iJ].GetTheta() > 0.8 * M_PI)
            continue;

         // now we gather the particles around the jets
         for(int iP = 0; iP < MParticle.nParticle; iP++)
         {
            if(GetAngle(MJet.Jet[iJ], MParticle.P[iP]) > JetR)
               continue;
         }
      }
   }
   Bar.Update(EntryCount);
   Bar.Print();
   Bar.PrintLine();

   File.Close();

   return 0;
}



