#include <iostream>
#include <vector>
#include <map>
using namespace std;

#include "TTree.h"
#include "TFile.h"

#include "Messenger.h"
#include "CommandLine.h"
#include "Matching.h"
#include "ProgressBar.h"
#include "TauHelperFunctions3.h"

#define MAX 1000
#define MAXPAIR 10000

int main(int argc, char *argv[]);
double MetricAngle(FourVector A, FourVector B);

int main(int argc, char *argv[])
{
   CommandLine CL(argc, argv);

   string InputFileName  = CL.Get("Input");
   string GenTreeName    = CL.Get("Gen", "tgen");
   string RecoTreeName   = CL.Get("Reco", "t");
   string OutputFileName = CL.Get("Output");
   double Fraction       = CL.GetDouble("Fraction", 1.00);

   TFile InputFile(InputFileName.c_str());
   TFile OutputFile(OutputFileName.c_str(), "RECREATE");

   TTree OutputTree("MatchedTree", "");
   TTree OutputPairTree("PairTree", "");

   int NParticle;
   double GenE[MAX], GenX[MAX], GenY[MAX], GenZ[MAX];
   double RecoE[MAX], RecoX[MAX], RecoY[MAX], RecoZ[MAX];
   double Distance[MAX];
   OutputTree.Branch("NParticle", &NParticle, "NParticle/I");
   OutputTree.Branch("GenE", &GenE, "GenE[NParticle]/D");
   OutputTree.Branch("GenX", &GenX, "GenX[NParticle]/D");
   OutputTree.Branch("GenY", &GenY, "GenY[NParticle]/D");
   OutputTree.Branch("GenZ", &GenZ, "GenZ[NParticle]/D");
   OutputTree.Branch("RecoE", &RecoE, "RecoE[NParticle]/D");
   OutputTree.Branch("RecoX", &RecoX, "RecoX[NParticle]/D");
   OutputTree.Branch("RecoY", &RecoY, "RecoY[NParticle]/D");
   OutputTree.Branch("RecoZ", &RecoZ, "RecoZ[NParticle]/D");
   OutputTree.Branch("Distance", &Distance, "Distance[NParticle]/D");

   int NPair;
   double GenE1[MAXPAIR], GenX1[MAXPAIR], GenY1[MAXPAIR], GenZ1[MAXPAIR];
   double GenE2[MAXPAIR], GenX2[MAXPAIR], GenY2[MAXPAIR], GenZ2[MAXPAIR];
   double RecoE1[MAXPAIR], RecoX1[MAXPAIR], RecoY1[MAXPAIR], RecoZ1[MAXPAIR];
   double RecoE2[MAXPAIR], RecoX2[MAXPAIR], RecoY2[MAXPAIR], RecoZ2[MAXPAIR];
   double DistanceGen[MAXPAIR], DistanceReco[MAXPAIR], Distance1[MAXPAIR], Distance2[MAXPAIR];
   OutputPairTree.Branch("NPair", &NPair);
   OutputPairTree.Branch("GenE1", &GenE1, "GenE1[NPair]/D");
   OutputPairTree.Branch("GenX1", &GenX1, "GenX1[NPair]/D");
   OutputPairTree.Branch("GenY1", &GenY1, "GenY1[NPair]/D");
   OutputPairTree.Branch("GenZ1", &GenZ1, "GenZ1[NPair]/D");
   OutputPairTree.Branch("GenE2", &GenE2, "GenE2[NPair]/D");
   OutputPairTree.Branch("GenX2", &GenX2, "GenX2[NPair]/D");
   OutputPairTree.Branch("GenY2", &GenY2, "GenY2[NPair]/D");
   OutputPairTree.Branch("GenZ2", &GenZ2, "GenZ2[NPair]/D");
   OutputPairTree.Branch("RecoE1", &RecoE1, "RecoE1[NPair]/D");
   OutputPairTree.Branch("RecoX1", &RecoX1, "RecoX1[NPair]/D");
   OutputPairTree.Branch("RecoY1", &RecoY1, "RecoY1[NPair]/D");
   OutputPairTree.Branch("RecoZ1", &RecoZ1, "RecoZ1[NPair]/D");
   OutputPairTree.Branch("RecoE2", &RecoE2, "RecoE2[NPair]/D");
   OutputPairTree.Branch("RecoX2", &RecoX2, "RecoX2[NPair]/D");
   OutputPairTree.Branch("RecoY2", &RecoY2, "RecoY2[NPair]/D");
   OutputPairTree.Branch("RecoZ2", &RecoZ2, "RecoZ2[NPair]/D");
   OutputPairTree.Branch("DistanceGen", &DistanceGen, "DistanceGen[NPair]/D");
   OutputPairTree.Branch("DistanceReco", &DistanceReco, "DistanceReco[NPair]/D");
   OutputPairTree.Branch("Distance1", &Distance1, "Distance1[NPair]/D");
   OutputPairTree.Branch("Distance2", &Distance2, "Distance2[NPair]/D");

   ParticleTreeMessenger MGen(InputFile, GenTreeName);
   ParticleTreeMessenger MReco(InputFile, RecoTreeName);

   int EntryCount = MGen.GetEntries() * Fraction;
   ProgressBar Bar(cout, EntryCount);
   Bar.SetStyle(-1);
   for(int iE = 0; iE < EntryCount; iE++)
   {
      if(EntryCount < 500 || iE % (EntryCount / 300) == 0)
      {
         Bar.Update(iE);
         Bar.Print();
      }

      MGen.GetEntry(iE);
      MReco.GetEntry(iE);

      if(MReco.PassBaselineCut() == false)
         continue;

      vector<FourVector> PGen, PReco;
      for(int i = 0; i < MGen.nParticle; i++)
      {
         if(MGen.charge[i] == 0)
            continue;
         if(MGen.highPurity[i] == false)
            continue;
         if(MGen.P[i][0] < 0)   // Here we can implement minimum energy cut
            continue;

         PGen.push_back(MGen.P[i]);
      }
      for(int i = 0; i < MReco.nParticle; i++)
      {
         if(MReco.charge[i] == 0)
            continue;
         if(MReco.highPurity[i] == false)
            continue;
         if(MReco.P[i][0] < 0)   // Here we can implement minimum energy cut
            continue;

         PReco.push_back(MReco.P[i]);
      }

      // cout << PGen.size() << " " << PReco.size() << endl;
      map<int, int> Matching = MatchJetsHungarian(MetricAngle, PGen, PReco);

      int Count = 0;
      NParticle = Matching.size();
      for(auto iter : Matching)
      {
         FourVector Gen = PGen[iter.first];
         FourVector Reco = iter.second >= 0 ? PReco[iter.second] : FourVector(-1, 0, 0, 0);

         GenE[Count] = Gen[0];
         GenX[Count] = Gen[1];
         GenY[Count] = Gen[2];
         GenZ[Count] = Gen[3];
         RecoE[Count] = Reco[0];
         RecoX[Count] = Reco[1];
         RecoY[Count] = Reco[2];
         RecoZ[Count] = Reco[3];
         Distance[Count] = MetricAngle(Gen, Reco);

         Count = Count + 1;
      }

      OutputTree.Fill();

      NPair = 0;
      for(int i = 0; i < NParticle; i++)
      {
         for(int j = i + 1; j < NParticle; j++)
         {
            GenE1[NPair] = GenE[i];
            GenX1[NPair] = GenX[i];
            GenY1[NPair] = GenY[i];
            GenZ1[NPair] = GenZ[i];
            GenE2[NPair] = GenE[j];
            GenX2[NPair] = GenX[j];
            GenY2[NPair] = GenY[j];
            GenZ2[NPair] = GenZ[j];
            RecoE1[NPair] = RecoE[i];
            RecoX1[NPair] = RecoX[i];
            RecoY1[NPair] = RecoY[i];
            RecoZ1[NPair] = RecoZ[i];
            RecoE2[NPair] = RecoE[j];
            RecoX2[NPair] = RecoX[j];
            RecoY2[NPair] = RecoY[j];
            RecoZ2[NPair] = RecoZ[j];

            Distance1[NPair] = Distance[i];
            Distance2[NPair] = Distance[j];

            FourVector Gen1(GenE[i], GenX[i], GenY[i], GenZ[i]);
            FourVector Gen2(GenE[j], GenX[j], GenY[j], GenZ[j]);
            FourVector Reco1(RecoE[i], RecoX[i], RecoY[i], RecoZ[i]);
            FourVector Reco2(RecoE[j], RecoX[j], RecoY[j], RecoZ[j]);

            DistanceGen[NPair] = GetAngle(Gen1, Gen2);
            DistanceReco[NPair] = GetAngle(Reco1, Reco2);

            NPair = NPair + 1;
         }
      }

      OutputPairTree.Fill();
   }
   Bar.Update(EntryCount);
   Bar.Print();
   Bar.PrintLine();

   OutputFile.cd();
   OutputTree.Write();
   OutputPairTree.Write();
   OutputFile.Close();

   OutputFile.Close();
   InputFile.Close();

   return 0;
}

double MetricAngle(FourVector A, FourVector B)
{
   double Angle = GetAngle(A, B);
   double EDiff = A[0] - B[0];

   return Angle;
}


