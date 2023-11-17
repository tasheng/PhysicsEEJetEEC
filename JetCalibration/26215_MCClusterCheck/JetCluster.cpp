#include <vector>
#include <iostream>
using namespace std;

#include "TFile.h"
#include "TTree.h"

#include "fastjet/ClusterSequence.hh"
using namespace fastjet;

#include "TauHelperFunctions3.h"
#include "CommandLine.h"

#include "Messenger.h"

int main(int argc, char *argv[]);

int main(int argc, char *argv[])
{
   CommandLine CL(argc, argv);

   string InputFileName    = CL.Get("Input");
   string OutputFileName   = CL.Get("Output");
   double JetR             = CL.GetDouble("JetR", 0.4);
   string JetTreeName      = CL.Get("Jet", "akR4ESchemeJetTree");
   string ParticleTreeName = CL.Get("P", "t");

   TFile File(InputFileName.c_str());

   JetTreeMessenger MJet(File, JetTreeName);
   ParticleTreeMessenger MP(File, ParticleTreeName);

   int JetCount = 0;
   int MatchedJetCount = 0;

   int EventCount = MJet.GetEntries();
   for(int iE = 0; iE < EventCount; iE++)
   {
      MP.GetEntry(iE);
      MJet.GetEntry(iE);

      // Collect particles and cluster jets
      vector<PseudoJet> FastJetParticles;
      for(int iP = 0; iP < MP.nParticle; iP++)
      {
         // FourVector P(MP.pmag[iP], MP.px[iP], MP.py[iP], MP.pz[iP]);
         // P[0] = sqrt(MP.pmag[iP] * MP.pmag[iP] + MP.mass[iP] * MP.mass[iP]);

         // if(MP.pwflag[iP] < 0 || MP.pwflag[iP] > 5)
         //    continue;

         FourVector P(0, MP.px[iP], MP.py[iP], MP.pz[iP]);
         P[0] = sqrt(P.GetP2() + MP.mass[iP] * MP.mass[iP]);

         // if() continue; // status?
         // if(MP.pid[iP] == 12 || MP.pid[iP] == -12)   continue;
         // if(MP.pid[iP] == 14 || MP.pid[iP] == -14)   continue;
         // if(MP.pid[iP] == 16 || MP.pid[iP] == -16)   continue;

         FastJetParticles.push_back(PseudoJet(P[1], P[2], P[3], P[0]));
      }
      JetDefinition Definition(ee_genkt_algorithm, JetR, -1, RecombinationScheme(E_scheme));
      ClusterSequence Sequence(FastJetParticles, Definition);
      vector<PseudoJet> FastJets = sorted_by_pt(Sequence.inclusive_jets(0.5));

      // Now we compare with existing jets
      // cout << "=== Event ===" << endl;
      // cout << " Stored" << endl;
      // for(FourVector J : MJet.Jet)
      //    if(J[0] > 0.5)
      //       cout << J[0] << " " << J.GetTheta() << " " << J.GetPhi() << endl;
      // cout << " Recluster" << endl;
      // for(PseudoJet J : FastJets)
      //    if(J.e() > 0.5)
      //       cout << J.e() << " " << J.theta() << " " << J.phi_std() << endl;

      int NJet = MJet.nref;
      for(int iJ = 0; iJ < NJet; iJ++)
      {
         FourVector J = MJet.Jet[iJ];

         if(J[0] < 5)
            continue;

         FourVector BestJet;
         PseudoJet BestFastJet;
         double BestAngle = -1;
         for(PseudoJet FJ : FastJets)
         {
            FourVector K(FJ.e(), FJ.px(), FJ.py(), FJ.pz());
            double Angle = GetAngle(J, K);

            if(Angle < BestAngle || BestAngle < 0)
            {
               BestAngle = Angle;
               BestJet = K;
               BestFastJet = FJ;
            }
         }

         JetCount = JetCount + 1;
         if(fabs(J[0] - BestJet[0]) > 0.1 || BestAngle > 0.01)
         {
            cout << "Hmm?" << endl;
            cout << J[0] << " " << BestJet[0] << " " << J[0] - BestJet[0] << " " << BestAngle << endl;
            cout << "Stored      " << J[0] << " " << J.GetTheta() << " " << J.GetPhi() << ", N = " << MJet.jtN[iJ] << endl;
            cout << "Reclustered " << BestJet[0] << " " << BestJet.GetTheta() << " " << BestJet.GetPhi() << ", N = " << BestFastJet.constituents().size() << endl;
         }
         else
         {
            MatchedJetCount = MatchedJetCount + 1;
            // cout << "Close" << endl;
            // cout << J[0] << " " << BestJet[0] << " " << J[0] - BestJet[0] << " " << BestAngle << endl;
            // cout << "Stored      " << J[0] << " " << J.GetTheta() << " " << J.GetPhi() << ", N = " << MJet.jtN[iJ] << endl;
            // cout << "Reclustered " << BestJet[0] << " " << BestJet.GetTheta() << " " << BestJet.GetPhi() << ", N = " << BestFastJet.constituents().size() << endl;
         }
      }
   }

   cout << "Final summary: " << MatchedJetCount << "/" << JetCount << " (" << MatchedJetCount * 100.0 / JetCount << "%) matched" << endl;

   File.Close();

   return 0;
}




