#include <vector>
#include <iostream>
using namespace std;

#include "TFile.h"
#include "TTree.h"

#include "fastjet/ClusterSequence.hh"
using namespace fastjet;

#include "TauHelperFunctions3.h"
#include "CommandLine.h"
#include "ProgressBar.h"

#include "Messenger.h"

#define MAXR 20
#define MAX 1000
#define TOTALGHOST 1e-8

int main(int argc, char *argv[]);
void AddGhosts(vector<PseudoJet> &P, int GhostSpacing);
double GetArea(PseudoJet &J);

int main(int argc, char *argv[])
{
   CommandLine CL(argc, argv);

   string InputFileName             = CL.Get("Input");
   string OutputFileName            = CL.Get("Output");
   vector<double> JetR              = CL.GetDoubleVector("JetR", vector<double>{0.2, 0.4, 0.6, 0.8, 1.0});
   string ParticleTreeName          = CL.Get("Reco", "t");
   string GenParticleTreeName       = CL.Get("Gen", "tgen");
   string GenBeforeParticleTreeName = CL.Get("GenBefore", "tgenBefore");
   bool SkipGen                     = CL.GetBool("SkipGen", false);

   int GhostSpacing                 = CL.GetInt("GhostSpacing", 50);

   TFile InputFile(InputFileName.c_str());

   ParticleTreeMessenger MReco(InputFile, ParticleTreeName);
   ParticleTreeMessenger MGen(InputFile, GenParticleTreeName);
   ParticleTreeMessenger MGenBefore(InputFile, GenBeforeParticleTreeName);

   TFile OutputFile(OutputFileName.c_str(), "RECREATE");

   vector<TTree *> RecoTree, GenTree, GenBeforeTree;

   int NRecoJet[MAXR];
   float RecoJetPT[MAXR][MAX], RecoJetEta[MAXR][MAX], RecoJetPhi[MAXR][MAX], RecoJetM[MAXR][MAX], RecoJetA[MAXR][MAX];
   int RecoJetN[MAXR][MAX];
   int NGenJet[MAXR];
   float GenJetPT[MAXR][MAX], GenJetEta[MAXR][MAX], GenJetPhi[MAXR][MAX], GenJetM[MAXR][MAX], GenJetA[MAXR][MAX];
   int GenJetN[MAXR][MAX];
   int NGenBeforeJet[MAXR];
   float GenBeforeJetPT[MAXR][MAX], GenBeforeJetEta[MAXR][MAX], GenBeforeJetPhi[MAXR][MAX], GenBeforeJetM[MAXR][MAX], GenBeforeJetA[MAXR][MAX];
   int GenBeforeJetN[MAXR][MAX];

   for(int iR = 0; iR < (int)JetR.size(); iR++)
   {
      double R = JetR[iR];
      RecoTree.push_back(     new TTree(Form("RecoR%d",      (int)(R * 10)), "Reclustered trees"));
      if(SkipGen == false)
      {
         GenTree.push_back(      new TTree(Form("GenR%d",       (int)(R * 10)), "Reclustered trees"));
         GenBeforeTree.push_back(new TTree(Form("GenBeforeR%d", (int)(R * 10)), "Reclustered trees"));
      }
      
      RecoTree[iR]->Branch("nref",       &NRecoJet[iR],        "nref/I");
      RecoTree[iR]->Branch("jtpt",       &RecoJetPT[iR],       "jtpt[nref]/F");
      RecoTree[iR]->Branch("jteta",      &RecoJetEta[iR],      "jteta[nref]/F");
      RecoTree[iR]->Branch("jtphi",      &RecoJetPhi[iR],      "jtphi[nref]/F");
      RecoTree[iR]->Branch("jtm",        &RecoJetM[iR],        "jtm[nref]/F");
      RecoTree[iR]->Branch("jtN",        &RecoJetN[iR],        "jtN[nref]/I");
      RecoTree[iR]->Branch("jtA",        &RecoJetA[iR],        "jtA[nref]/F");
      if(SkipGen == false)
      {
         GenTree[iR]->Branch("nref",        &NGenJet[iR],         "nref/I");
         GenTree[iR]->Branch("jtpt",        &GenJetPT[iR],        "jtpt[nref]/F");
         GenTree[iR]->Branch("jteta",       &GenJetEta[iR],       "jteta[nref]/F");
         GenTree[iR]->Branch("jtphi",       &GenJetPhi[iR],       "jtphi[nref]/F");
         GenTree[iR]->Branch("jtm",         &GenJetM[iR],         "jtm[nref]/F");
         GenTree[iR]->Branch("jtN",         &GenJetN[iR],         "jtN[nref]/I");
         GenTree[iR]->Branch("jtA",         &GenJetA[iR],         "jtA[nref]/F");
         GenBeforeTree[iR]->Branch("nref",  &NGenBeforeJet[iR],   "nref/I");
         GenBeforeTree[iR]->Branch("jtpt",  &GenBeforeJetPT[iR],  "jtpt[nref]/F");
         GenBeforeTree[iR]->Branch("jteta", &GenBeforeJetEta[iR], "jteta[nref]/F");
         GenBeforeTree[iR]->Branch("jtphi", &GenBeforeJetPhi[iR], "jtphi[nref]/F");
         GenBeforeTree[iR]->Branch("jtm",   &GenBeforeJetM[iR],   "jtm[nref]/F");
         GenBeforeTree[iR]->Branch("jtN",   &GenBeforeJetN[iR],   "jtN[nref]/I");
         GenBeforeTree[iR]->Branch("jtA",   &GenBeforeJetA[iR],   "jtA[nref]/F");
      }
   }

   int EventCount = MReco.GetEntries();
   ProgressBar Bar(cout, EventCount);
   Bar.SetStyle(-1);
   for(int iE = 0; iE < EventCount; iE++)
   {
      if(EventCount < 500 || iE % (EventCount / 250) == 0)
      {
         Bar.Update(iE);
         Bar.Print();
      }

      MReco.GetEntry(iE);
      if(SkipGen == false)
      {
         MGen.GetEntry(iE);
         MGenBefore.GetEntry(iE);
      }

      // Collect particles and cluster jets
      vector<PseudoJet> RecoFastJetParticles, GenFastJetParticles, GenBeforeFastJetParticles;
      for(int iP = 0; iP < MReco.nParticle; iP++)
      {
         if(MReco.pwflag[iP] < 0 || MReco.pwflag[iP] > 5)
            continue;

         FourVector P(0, MReco.px[iP], MReco.py[iP], MReco.pz[iP]);
         P[0] = sqrt(P.GetP2() + MReco.mass[iP] * MReco.mass[iP]);

         RecoFastJetParticles.push_back(PseudoJet(P[1], P[2], P[3], P[0]));
      }
      AddGhosts(RecoFastJetParticles, GhostSpacing);
      if(SkipGen == false)
      {
         for(int iP = 0; iP < MGen.nParticle; iP++)
         {
            FourVector P(0, MGen.px[iP], MGen.py[iP], MGen.pz[iP]);
            P[0] = sqrt(P.GetP2() + MGen.mass[iP] * MGen.mass[iP]);

            GenFastJetParticles.push_back(PseudoJet(P[1], P[2], P[3], P[0]));
         }
         AddGhosts(GenFastJetParticles, GhostSpacing);
         for(int iP = 0; iP < MGenBefore.nParticle; iP++)
         {
            FourVector P(0, MGenBefore.px[iP], MGenBefore.py[iP], MGenBefore.pz[iP]);
            P[0] = sqrt(P.GetP2() + MGenBefore.mass[iP] * MGenBefore.mass[iP]);

            GenBeforeFastJetParticles.push_back(PseudoJet(P[1], P[2], P[3], P[0]));
         }
         AddGhosts(GenBeforeFastJetParticles, GhostSpacing);
      }

      // Now do all the clustering
      for(int iR = 0; iR < (int)JetR.size(); iR++)
      {
         JetDefinition RecoDefinition(ee_genkt_algorithm, JetR[iR], -1, RecombinationScheme(E_scheme));
         ClusterSequence RecoSequence(RecoFastJetParticles, RecoDefinition);
         vector<PseudoJet> RecoFastJets = sorted_by_pt(RecoSequence.inclusive_jets(0));

         JetDefinition GenDefinition(ee_genkt_algorithm, JetR[iR], -1, RecombinationScheme(E_scheme));
         ClusterSequence GenSequence(GenFastJetParticles, GenDefinition);
         vector<PseudoJet> GenFastJets = sorted_by_pt(GenSequence.inclusive_jets(0));

         JetDefinition GenBeforeDefinition(ee_genkt_algorithm, JetR[iR], -1, RecombinationScheme(E_scheme));
         ClusterSequence GenBeforeSequence(GenBeforeFastJetParticles, GenBeforeDefinition);
         vector<PseudoJet> GenBeforeFastJets = sorted_by_pt(GenBeforeSequence.inclusive_jets(0));

         NRecoJet[iR] = RecoFastJets.size();
         for(int iJ = 0; iJ < (int)RecoFastJets.size(); iJ++)
         {
            RecoJetPT[iR][iJ] = RecoFastJets[iJ].perp();
            RecoJetEta[iR][iJ] = RecoFastJets[iJ].eta();
            RecoJetPhi[iR][iJ] = RecoFastJets[iJ].phi();
            RecoJetM[iR][iJ] = RecoFastJets[iJ].m();
            RecoJetN[iR][iJ] = RecoFastJets[iJ].constituents().size();
            RecoJetA[iR][iJ] = GetArea(RecoFastJets[iJ]);
         }

         if(SkipGen == false)
         {
            NGenJet[iR] = GenFastJets.size();
            for(int iJ = 0; iJ < (int)GenFastJets.size(); iJ++)
            {
               GenJetPT[iR][iJ] = GenFastJets[iJ].perp();
               GenJetEta[iR][iJ] = GenFastJets[iJ].eta();
               GenJetPhi[iR][iJ] = GenFastJets[iJ].phi();
               GenJetM[iR][iJ] = GenFastJets[iJ].m();
               GenJetN[iR][iJ] = GenFastJets[iJ].constituents().size();
               GenJetA[iR][iJ] = GetArea(GenFastJets[iJ]);
            }

            NGenBeforeJet[iR] = GenBeforeFastJets.size();
            for(int iJ = 0; iJ < (int)GenBeforeFastJets.size(); iJ++)
            {
               GenBeforeJetPT[iR][iJ] = GenBeforeFastJets[iJ].perp();
               GenBeforeJetEta[iR][iJ] = GenBeforeFastJets[iJ].eta();
               GenBeforeJetPhi[iR][iJ] = GenBeforeFastJets[iJ].phi();
               GenBeforeJetM[iR][iJ] = GenBeforeFastJets[iJ].m();
               GenBeforeJetN[iR][iJ] = GenBeforeFastJets[iJ].constituents().size();
               GenBeforeJetA[iR][iJ] = GetArea(GenBeforeFastJets[iJ]);
            }
         }

         RecoTree[iR]->Fill();
         if(SkipGen == false)
         {
            GenTree[iR]->Fill();
            GenBeforeTree[iR]->Fill();
         }
      }
   }

   Bar.Update(EventCount);
   Bar.Print();
   Bar.PrintLine();

   for(int iR = 0; iR < (int)JetR.size(); iR++)
   {
      RecoTree[iR]->Write();
      if(SkipGen == false)
      {
         GenTree[iR]->Write();
         GenBeforeTree[iR]->Write();
      }
   }

   OutputFile.Close();

   InputFile.Close();

   return 0;
}

void AddGhosts(vector<PseudoJet> &P, int GhostSpacing)
{
   double TotalArea = 0;

   for(int i = 0; i < GhostSpacing; i++)
   {
      double Theta = M_PI * (i + 0.5) / GhostSpacing;

      int NPhi = (int)(GhostSpacing * sin(Theta)) + 1;
      for(int j = 0; j < NPhi; j++)
      {
         double Phi = 2 * M_PI * (j + 0.5) / NPhi;

         double A = (M_PI / GhostSpacing) * (2 * M_PI * sin(Theta) / NPhi);
         TotalArea = TotalArea + A;

         FourVector G;
         G.SetSizeThetaPhi(TOTALGHOST * A, Theta, Phi);
         P.push_back(PseudoJet(G[1], G[2], G[3], G[0]));
      }
   }

   for(int i = 0; i < (int)P.size(); i++)
   {
      if(P[i].e() < TOTALGHOST)
         P[i] = P[i] / TotalArea;
   }
}

double GetArea(PseudoJet &J)
{
   double TotalGhostE = 0;

   for(PseudoJet &P : J.constituents())
   {
      if(P.e() < TOTALGHOST)
         TotalGhostE = TotalGhostE + P.e();
   }

   return TotalGhostE / TOTALGHOST * 4 * M_PI;
}





