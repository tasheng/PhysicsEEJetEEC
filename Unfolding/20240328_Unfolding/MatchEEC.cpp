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
#include "alephTrkEfficiency.h"

#include "TCanvas.h"
#include "TH1D.h"

#define MAX 1000
#define MAXPAIR 10000

int main(int argc, char *argv[]);
double MetricAngle(FourVector A, FourVector B);
double MatchingMetric(FourVector A, FourVector B);
int FindBin(double Value, int NBins, double Bins[]); 

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

   //------------------------------------
   // define the binning
   //------------------------------------


   // theta binning
   const int BinCount = 100;
   double Bins[2*BinCount+1];
   double BinMin = 0.002;
   double BinMax = M_PI / 2;

   // z binning
   double zBins[2*BinCount+1]; 
   double zBinMin = (1- cos(0.002))/2; 
   double zBinMax = 0.5;

   // energy binning
   double EnergyBins[BinCount+1];
   double EnergyBinMin = 4e-6; 
   double EnergyBinMax = 0.2; 
   double logMin = std::log10(EnergyBinMin);
   double logMax = std::log10(EnergyBinMax);
   double logStep = (logMax - logMin) / (BinCount);

   for(int i = 0; i <= BinCount; i++){
      // theta double log binning
      Bins[i] = exp(log(BinMin) + (log(BinMax) - log(BinMin)) / BinCount * i);
      Bins[2*BinCount-i] = BinMax * 2 - exp(log(BinMin) + (log(BinMax) - log(BinMin)) / BinCount * i);
     
      // z double log binning
      zBins[i] = exp(log(zBinMin) + (log(zBinMax) - log(zBinMin)) / BinCount * i);
      zBins[2*BinCount-i] = zBinMax * 2 - exp(log(zBinMin) + (log(zBinMax) - log(zBinMin)) / BinCount * i);

      double logValue = logMin + i * logStep;
      EnergyBins[i] =  std::pow(10, logValue);//exp(log(EnergyBinMin) + (log(EnergyBinMax) - log(EnergyBinMin)) / BinCount * i);
      //std::cout << "Adding energy bin " << EnergyBins[i] << std::endl;
   }

   //------------------------------------
   // define the trees
   //------------------------------------

   // tree for single track matching
   TTree OutputTree("MatchedTree", "");
   // tree for matched pairs
   TTree OutputPairTree("PairTree", "");
   // tree that just takes pairs, not matching taken into account
   TTree OutputUnmatchedTree("UnmatchedPairTree", ""); 

   //------------------------------------
   // define the histograms
   //------------------------------------


   // as a function of theta
   TH1D recoUnmatched("recoUnmatched", "recoUnmatched", 2 * BinCount, 0, 2 * BinCount);
   TH1D genUnmatched("genUnmatched", "genUnmatched", 2 * BinCount, 0, 2 * BinCount);
   TH1D recoMatched("recoMatched", "recoMatched", 2 * BinCount, 0, 2 * BinCount);
   TH1D genMatched("genMatched", "genMatched", 2 * BinCount, 0, 2 * BinCount);

   // as a function of z
   TH1D recoUnmatched_z("recoUnmatched_z", "recoUnmatched_z", 2 * BinCount, 0, 2 * BinCount);
   TH1D genUnmatched_z("genUnmatched_z", "genUnmatched_z", 2 * BinCount, 0, 2 * BinCount);
   TH1D recoMatched_z("recoMatched_z", "recoMatched_z", 2 * BinCount, 0, 2 * BinCount);
   TH1D genMatched_z("genMatched_z", "genMatched_z", 2 * BinCount, 0, 2 * BinCount);

   // as a function of e1e2   
   TH1D e1e2RecoUnmatched("e1e2RecoUnmatched", "e1e2RecoUnmatched", BinCount, EnergyBins); 
   TH1D e1e2RecoMatched("e1e2RecoMatched", "e1e2RecoMatched", BinCount, EnergyBins); 
   TH1D e1e2GenUnmatched("e1e2GenUnmatched", "e1e2GenUnmatched", BinCount, EnergyBins); 
   TH1D e1e2GenMatched("e1e2GenMatched", "e1e2GenMatched", BinCount, EnergyBins); 


   alephTrkEfficiency efficiencyCorrector;
   // variables for the matched tree
   int NParticle, eventID; 
   double GenE[MAX], GenX[MAX], GenY[MAX], GenZ[MAX];
   double RecoE[MAX], RecoX[MAX], RecoY[MAX], RecoZ[MAX];
   double Distance[MAX], DeltaPhi[MAX], Metric[MAX], DeltaE[MAX], DeltaTheta[MAX];
   double RecoEta[MAX], RecoPhi[MAX], RecoTheta[MAX], RecoRapidity[MAX];
   int RecoPWFlag[MAX];
   double GenEta[MAX], GenPhi[MAX], GenTheta[MAX], GenRapidity[MAX];
   double RecoEfficiency[MAX];
   OutputTree.Branch("EventID", &eventID, "EventID/I");
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
   OutputTree.Branch("DeltaPhi", &DeltaPhi, "DeltaPhi[NParticle]/D");
   OutputTree.Branch("DeltaTheta", &DeltaTheta, "DeltaTheta[NParticle]/D"); 
   OutputTree.Branch("DeltaE", &DeltaE, "DeltaE[NParticle]/D");
   OutputTree.Branch("Metric", &Metric, "Metric[NParticle]/D");
   OutputTree.Branch("RecoEta", &RecoEta, "RecoEta[NParticle]/D");
   OutputTree.Branch("RecoPhi", &RecoPhi, "RecoPhi[NParticle]/D");
   OutputTree.Branch("RecoTheta", &RecoTheta, "RecoTheta[NParticle]/D");
   OutputTree.Branch("GenEta", &GenEta, "GenEta[NParticle]/D");
   OutputTree.Branch("GenPhi", &GenPhi, "GenPhi[NParticle]/D");
   OutputTree.Branch("GenTheta", &GenTheta, "GenTheta[NParticle]/D");
   OutputTree.Branch("RecoRapidity", &RecoRapidity, "RecoRapidity[NParticle]/D");
   OutputTree.Branch("GenRapidity", &GenRapidity, "GenRapidity[NParticle]/D");
   OutputTree.Branch("RecoPwFlag", &RecoPWFlag, "RecoPWFlag[NParticle]/I");
   OutputTree.Branch("RecoEfficiency", &RecoEfficiency, "RecoEfficiency[NParticle]/D");

   // valiables for the pair tree
   int NPair;
   double GenE1[MAXPAIR], GenX1[MAXPAIR], GenY1[MAXPAIR], GenZ1[MAXPAIR];
   double GenE2[MAXPAIR], GenX2[MAXPAIR], GenY2[MAXPAIR], GenZ2[MAXPAIR];
   double RecoE1[MAXPAIR], RecoX1[MAXPAIR], RecoY1[MAXPAIR], RecoZ1[MAXPAIR];
   double RecoE2[MAXPAIR], RecoX2[MAXPAIR], RecoY2[MAXPAIR], RecoZ2[MAXPAIR];
   double DistanceGen[MAXPAIR], DistanceReco[MAXPAIR], Distance1[MAXPAIR], Distance2[MAXPAIR];
   double E1E2Gen[MAXPAIR], E1E2Reco[MAXPAIR];
   double RecoEfficiency1[MAXPAIR], RecoEfficiency2[MAXPAIR];
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
   OutputPairTree.Branch("E1E2Gen", &E1E2Gen, "E1E2Gen[NPair]/D");
   OutputPairTree.Branch("E1E2Reco", &E1E2Reco, "E1E2Reco[NPair]/D");
   OutputPairTree.Branch("RecoEfficiency1", &RecoEfficiency1, "RecoEfficiency1[NPair]/D");
   OutputPairTree.Branch("RecoEfficiency2", &RecoEfficiency2, "RecoEfficiency2[NPair]/D");

   // variables for the unmatched pair tree
   int NUnmatchedPair;
   double GenE1Unmatched[MAXPAIR], GenX1Unmatched[MAXPAIR], GenY1Unmatched[MAXPAIR], GenZ1Unmatched[MAXPAIR];
   double GenE2Unmatched[MAXPAIR], GenX2Unmatched[MAXPAIR], GenY2Unmatched[MAXPAIR], GenZ2Unmatched[MAXPAIR];
   double RecoE1Unmatched[MAXPAIR], RecoX1Unmatched[MAXPAIR], RecoY1Unmatched[MAXPAIR], RecoZ1Unmatched[MAXPAIR];
   double RecoE2Unmatched[MAXPAIR], RecoX2Unmatched[MAXPAIR], RecoY2Unmatched[MAXPAIR], RecoZ2Unmatched[MAXPAIR];
   double DistanceUnmatchedGen[MAXPAIR], DistanceUnmatchedReco[MAXPAIR], DeltaPhiUnmatchedGen[MAXPAIR], DeltaPhiUnmatchedReco[MAXPAIR], DeltaEUnmatchedReco[MAXPAIR], DeltaEUnmatchedGen[MAXPAIR], DeltaThetaUnmatchedGen[MAXPAIR], DeltaThetaUnmatchedReco[MAXPAIR]; 
   double E1E2GenUnmatched[MAXPAIR], E1E2RecoUnmatched[MAXPAIR];
   OutputUnmatchedTree.Branch("NUnmatchedPair", &NUnmatchedPair);
   OutputUnmatchedTree.Branch("GenE1Unmatched", &GenE1Unmatched, "GenE1Unmatched[NUnmatchedPair]/D");
   OutputUnmatchedTree.Branch("GenE2Unmatched", &GenE2Unmatched, "GenE2Unmatched[NUnmatchedPair]/D");
   OutputUnmatchedTree.Branch("RecoE1Unmatched", &RecoE1Unmatched, "RecoE1Unmatched[NUnmatchedPair]/D");
   OutputUnmatchedTree.Branch("RecoE2Unmatched", &RecoE2Unmatched, "RecoE2Unmatched[NUnmatchedPair]/D");
   OutputUnmatchedTree.Branch("DistanceUnmatchedGen", &DistanceUnmatchedGen, "DistanceUnmatchedGen[NUnmatchedPair]/D");
   OutputUnmatchedTree.Branch("DistanceUnmatchedReco", &DistanceUnmatchedReco, "DistanceUnmatchedReco[NUnmatchedPair]/D");
   OutputUnmatchedTree.Branch("DeltaPhiUnmatchedGen", &DeltaPhiUnmatchedGen, "DeltaPhiUnmatchedGen[NUnmatchedPair]/D");
   OutputUnmatchedTree.Branch("DeltaPhiUnmatchedReco", &DeltaPhiUnmatchedReco, "DeltaPhiUnmatchedReco[NUnmatchedPair]/D");
   OutputUnmatchedTree.Branch("DeltaEUnmatchedGen", &DeltaEUnmatchedGen, "DeltaEUnmatchedGen[NUnmatchedPair]/D");
   OutputUnmatchedTree.Branch("DeltaEUnmatchedReco", &DeltaEUnmatchedReco, "DeltaEUnmatchedReco[NUnmatchedPair]/D");
   OutputUnmatchedTree.Branch("DeltaThetaUnmatchedGen", &DeltaThetaUnmatchedGen, "DeltaThetaUnmatchedGen[NUnmatchedPair]/D");
   OutputUnmatchedTree.Branch("DeltaThetaUnmatchedReco", &DeltaThetaUnmatchedReco, "DeltaThetaUnmatchedReco[NUnmatchedPair]/D");
   OutputUnmatchedTree.Branch("E1E2GenUnmatched", &E1E2GenUnmatched, "E1E2GenUnmatched[NUnmatchedPair]/D");
   OutputUnmatchedTree.Branch("E1E2RecoUnmatched", &E1E2RecoUnmatched, "E1E2RecoUnmatched[NUnmatchedPair]/D");

   ParticleTreeMessenger MGen(InputFile, GenTreeName);
   ParticleTreeMessenger MReco(InputFile, RecoTreeName);

   int EntryCount = MReco.GetEntries() * Fraction;
   int nAcceptedEvents = 0; 
   ProgressBar Bar(cout, EntryCount);
   Bar.SetStyle(-1); 
   for(int iE = 0; iE < EntryCount; iE++) 
   {

      MGen.GetEntry(iE);
      MReco.GetEntry(iE);

      // event selection (not applied for the moment as this is assumed to be applied to the tree)
      /*
      if(MReco.STheta < (7*M_PI)/36 || MReco.STheta >(29*M_PI)/36)continue;
      if(MGen.STheta < (7*M_PI)/36 || MGen.STheta > (29*M_PI)/36)continue;
      if(MGen.passesAll == false) continue;
      if(MReco.passesAll == false) continue;
      if(MReco.PassBaselineCut() == false) continue;
      if(MReco.passesLEP1TwoPC == false) continue;
      */
    

      // place a cut on the number of charged hadrons if you want to make the mult plot
      // if( MGen.nChargedHadronsHP < 40)continue; 
     

      vector<FourVector> PGen, PReco;
      for(int i = 0; i < MGen.nParticle; i++){
         if(MGen.charge[i] == 0) continue;
         if(MGen.highPurity[i] == false) continue;
         PGen.push_back(MGen.P[i]);
      }

      for(int i = 0; i < MReco.nParticle; i++){
         if(MReco.charge[i] == 0) continue;
         if(MReco.highPurity[i] == false) continue;
         // place cut on the reco energy, not included at gen level
         if(MReco.P[i][0] < 0.2) continue;        
         PReco.push_back(MReco.P[i]);
      }


      // now fill the unmatched tree
      int index_counter = 0;
      double TotalE = 91.1876;
      int NunmatchedRecoPairs = 0; 
      for(int i = 0; i < PReco.size(); i++){
         for(int j = i+1; j < PReco.size();j++){
            if(i == j) continue; // don't match particles with themselves
            // fill the tree
            FourVector Reco1 = PReco.at(i);
            FourVector Reco2 = PReco.at(j);
            RecoE1Unmatched[index_counter] = Reco1[0];
            RecoE2Unmatched[index_counter] = Reco2[0];
            DistanceUnmatchedReco[index_counter] = GetAngle(Reco1,Reco2);
            DeltaPhiUnmatchedReco[index_counter] = GetDPhi(Reco1, Reco2); 
            DeltaEUnmatchedReco[index_counter] = Reco1[0]-Reco2[0]; 
            DeltaThetaUnmatchedReco[index_counter] = Reco1.GetTheta() - Reco2.GetTheta(); 
            E1E2RecoUnmatched[index_counter] =  Reco1[0]*Reco2[0]; 

            // fill the theta distributions
            int BinThetaMeasuredMC = FindBin(GetAngle(Reco1,Reco2), 2 * BinCount, Bins);
            recoUnmatched.Fill(BinThetaMeasuredMC,Reco1[0]*Reco2[0]/(TotalE*TotalE));
     
            // fille in the z distributions
            double zUnmatched = (1-cos(GetAngle(Reco1, Reco2)))/2; 
            int BinZMeasured = FindBin(zUnmatched, 2*BinCount, zBins); 
            recoUnmatched_z.Fill(BinZMeasured, Reco1[0]*Reco2[0]/(TotalE*TotalE));

            // fill the energy distributions
            int BinEnergyMeasured = FindBin(Reco1[0]*Reco2[0]/(TotalE*TotalE), BinCount, EnergyBins);
            e1e2RecoUnmatched.Fill(Reco1[0]*Reco2[0]/(TotalE*TotalE));
         
            NunmatchedRecoPairs++; 
            index_counter++;
         }
      }

      int index_gen = 0; 
      int NunmatchedGenPairs = 0; 
      for(int i = 0; i < PGen.size(); i++){
      for(int j = i+1; j < PGen.size();j++){
            FourVector Gen1 = PGen.at(i);
            FourVector Gen2 = PGen.at(j);
            GenE1Unmatched[index_gen] = Gen1[0];
            GenE2Unmatched[index_gen] = Gen2[0];
            DistanceUnmatchedGen[index_gen] = GetAngle(Gen1,Gen2);
            DeltaPhiUnmatchedGen[index_gen] = GetDPhi(Gen1, Gen2); 
            DeltaEUnmatchedGen[index_gen] = Gen1[0]-Gen2[0]; 
            DeltaThetaUnmatchedGen[index_gen] = Gen1.GetTheta() - Gen2.GetTheta(); 
            E1E2GenUnmatched[index_gen] =  Gen1[0]*Gen2[0]; 

            // theta histograms
            int BinThetaGenMC = FindBin(GetAngle(Gen1,Gen2), 2 * BinCount, Bins);
            genUnmatched.Fill(BinThetaGenMC,Gen1[0]*Gen2[0]/(TotalE*TotalE));

            // energy histograms
            int BinEnergyGen = FindBin(Gen1[0]*Gen2[0]/(TotalE*TotalE), BinCount, EnergyBins);
            e1e2GenUnmatched.Fill(Gen1[0]*Gen2[0]/(TotalE*TotalE));

            // z histograms
            double zGenUnmatched = (1-cos(GetAngle(Gen1, Gen2)))/2; 
            int BinZGen = FindBin(zGenUnmatched, 2*BinCount, zBins); 
            genUnmatched_z.Fill(BinZGen, Gen1[0]*Gen2[0]/(TotalE*TotalE));


            index_gen++;
            NunmatchedGenPairs++; 
         }
      }

      if(index_gen > index_counter)NUnmatchedPair = index_gen;
      else NUnmatchedPair = index_counter;


      nAcceptedEvents++; 
      OutputUnmatchedTree.Fill();

      // perform the matching
      map<int, int> Matching = MatchJetsHungarian(MatchingMetric, PGen, PReco);
      int Count = 0;
      NParticle = Matching.size();
      eventID = MGen.EventNo;
      for(auto iter : Matching)
      {
         FourVector Gen = iter.first >= 0 ? PGen[iter.first] : FourVector(-1, 0, 0, 0);
         FourVector Reco = iter.second >= 0 ? PReco[iter.second] : FourVector(-1, 0, 0, 0);
         // if(GetAngle(Gen, Reco) > 1.0) {
         //    Gen = FourVector(-1, 0, 0, 0); 
         //    Reco = FourVector(-1, 0, 0, 0); 
         // }
         GenE[Count] = Gen[0];
         GenX[Count] = Gen[1];
         GenY[Count] = Gen[2];
         GenZ[Count] = Gen[3];
         GenEta[Count] = Gen.GetEta();
         GenRapidity[Count] = Gen.GetRapidity();
         GenPhi[Count] = Gen.GetPhi();
         GenTheta[Count] = Gen.GetTheta();
         RecoE[Count] = Reco[0];
         RecoX[Count] = Reco[1];
         RecoY[Count] = Reco[2];
         RecoZ[Count] = Reco[3];
         RecoRapidity[Count] = Reco.GetRapidity();
         RecoEta[Count] = Reco.GetEta();
         RecoPhi[Count] = Reco.GetPhi();
         RecoTheta[Count] = Reco.GetTheta();
         Distance[Count] = GetAngle(Gen, Reco);
         DeltaPhi[Count] = GetDPhi(Gen,Reco); 
         DeltaE[Count] = Gen[0] - Reco[0]; 
         DeltaTheta[Count] = Gen.GetTheta() - Reco.GetTheta(); 
         double Efficiency; 
         if (Reco.GetPT() <  0.2) Efficiency = 1;
         else Efficiency = efficiencyCorrector.efficiency(Reco.GetTheta(), Reco.GetPhi(), Reco.GetPT(), MReco.nChargedHadronsHP);
         RecoEfficiency[Count] = 1/Efficiency;
         Count = Count + 1;
      }
 
      OutputTree.Fill(); // fill the tree

      // now fill the tree for the matched pairs
      NPair = 0; 
      int NMatchedPairs = 0; 
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
            RecoEfficiency1[NPair] = RecoEfficiency[i];
            RecoEfficiency2[NPair] = RecoEfficiency[j];

            Distance1[NPair] = Distance[i];
            Distance2[NPair] = Distance[j];

            FourVector Gen1(GenE[i], GenX[i], GenY[i], GenZ[i]);
            FourVector Gen2(GenE[j], GenX[j], GenY[j], GenZ[j]);
            FourVector Reco1(RecoE[i], RecoX[i], RecoY[i], RecoZ[i]);
            FourVector Reco2(RecoE[j], RecoX[j], RecoY[j], RecoZ[j]);

            DistanceGen[NPair] = GetAngle(Gen1, Gen2);
            DistanceReco[NPair] = GetAngle(Reco1, Reco2);

            E1E2Gen[NPair] = (GenE[i]*GenE[j])/(TotalE*TotalE);
            E1E2Reco[NPair] = (RecoE[i]*RecoE[j])/(TotalE*TotalE);

            if(RecoE[i] > 0 && RecoE[j] > 0 && GenE[i] > 0 && GenE[j] > 0){
               // theta histograms
               int BinThetaMeasuredMC = FindBin(GetAngle(Reco1,Reco2), 2 * BinCount, Bins);
               recoMatched.Fill(BinThetaMeasuredMC,Reco1[0]*Reco2[0]/(TotalE*TotalE));
               int BinThetaGenMC = FindBin(GetAngle(Gen1,Gen2), 2 * BinCount, Bins);
               genMatched.Fill(BinThetaGenMC,Gen1[0]*Gen2[0]/(TotalE*TotalE));
            
               // energy histograms
               int BinEnergyMeasured = FindBin(Reco1[0]*Reco2[0]/(TotalE*TotalE), BinCount, EnergyBins);
               e1e2RecoMatched.Fill(Reco1[0]*Reco2[0]/(TotalE*TotalE));
               int BinEnergyGen = FindBin(Gen1[0]*Gen2[0]/(TotalE*TotalE), BinCount, EnergyBins);
               e1e2GenMatched.Fill(Gen1[0]*Gen2[0]/(TotalE*TotalE));

               // z histograms
               double zRecoMatched = (1-cos(GetAngle(Reco1, Reco2)))/2; 
               int BinZMeasured = FindBin(zRecoMatched, 2*BinCount, zBins);
               recoMatched_z.Fill(BinZMeasured, Reco1[0]*Reco2[0]/(TotalE*TotalE));

               double zGenMatched = (1-cos(GetAngle(Gen1, Gen2)))/2; 
               int BinZMC = FindBin(zGenMatched, 2*BinCount, zBins); 
               genMatched_z.Fill(BinZMC, Gen1[0]*Gen2[0]/(TotalE*TotalE));

               // increment the number of matched pairs
               NMatchedPairs++; 
            }
   
            NPair = NPair + 1;
         }
      }

      OutputPairTree.Fill();
   }
   std::cout << "Number of accepted events is " << nAcceptedEvents << std::endl;
   Bar.Update(EntryCount);
   Bar.Print();
   Bar.PrintLine();

   OutputFile.cd();

   // write the trees to the output file
   OutputTree.Write();
   OutputPairTree.Write();
   OutputUnmatchedTree.Write(); 

   // -------------------------------------------------------------------
   // histograms for the matching efficiency and fake fraction
   // -------------------------------------------------------------------

   //function of theta
   recoUnmatched.Write(); 
   recoMatched.Write(); 
   genUnmatched.Write(); 
   genMatched.Write(); 

   // function of z
   recoUnmatched_z.Write(); 
   recoMatched_z.Write(); 
   genUnmatched_z.Write(); 
   genMatched_z.Write(); 

   // function of energy
   e1e2GenMatched.Write(); 
   e1e2RecoMatched.Write(); 
   e1e2GenUnmatched.Write(); 
   e1e2RecoUnmatched.Write(); 
   // -------------------------------------------------------------------

   // write the output files
   OutputFile.Close();
   InputFile.Close();

   return 0;
}

// metric angle used as the first guess of the matching metric
double MetricAngle(FourVector A, FourVector B)
{
   double Angle = GetAngle(A, B);
   return Angle; 
}

double MatchingMetric(FourVector A, FourVector B){
   double Angle = GetAngle(A,B);
   double dPhi = GetDPhi(A,B); 
   double Ediff = (A[0] - B[0]); 
   double meanE = (A[0] + B[0])/2;  
   double phiRes = 0.002; // somewhere in the ballpark of 0.002-0.005
   double Eres = 0.1*meanE; // use 10% energy resolution //0.85/sqrt(meanE); // take energy resolution based on the mean energy
   double AngleRes = (0.01*0.01); // take angular resolution based on the detector resolution
   double chi2metric = (Angle/AngleRes)*(Angle/AngleRes)  + (dPhi/phiRes)*(dPhi/phiRes) + (Ediff/Eres)*(Ediff/Eres);
   return chi2metric; 
}

// function for filling the histograms
int FindBin(double Value, int NBins, double Bins[])
{
   for(int i = 0; i < NBins; i++)
      if(Value < Bins[i])
         return i - 1;
   return NBins;
}
