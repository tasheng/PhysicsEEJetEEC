#include <iostream>
#include <vector>

#include "TTree.h"
#include "TFile.h"

#include "TauHelperFunctions3.h"

#define MAXJET 1000
#define MAXPW 6
#define MAXPARTICLE 1000

class JetTreeMessenger;
class ParticleTreeMessenger;
class ReducedTreeMessenger;

class JetTreeMessenger
{
public:
   TTree *Tree;
   int    nref;
   float  jtpt[MAXJET];
   float  jteta[MAXJET];
   float  jtphi[MAXJET];
   float  jtm[MAXJET];
   int    jtN[MAXJET];
   int    jtNPW[MAXJET][MAXPW];
   float  jtptFracPW[MAXJET][MAXPW];
   float  zgJtPt_Beta0p00ZCut0p10[MAXJET];
   float  zgJtPhi_Beta0p00ZCut0p10[MAXJET];
   float  zgJtEta_Beta0p00ZCut0p10[MAXJET];
   float  zg_Beta0p00ZCut0p10[MAXJET];
   float  rg_Beta0p00ZCut0p10[MAXJET];
public:
   std::vector<FourVector> Jet;
public:
   JetTreeMessenger();
   JetTreeMessenger(TFile &file, std::string name);
   JetTreeMessenger(TFile *file, std::string name);
   JetTreeMessenger(TTree *tree);
   bool Initialize(TTree *tree);
   bool Initialize();
   bool GetEntry(int iEntry);
   int GetEntries();
};

class ParticleTreeMessenger
{
public:
   TTree *Tree;
   int           EventNo;
   int           RunNo;
   int           year;
   int           subDir;
   int           process;
   bool          isMC;
   unsigned long long uniqueID;
   float         Energy;
   int           bFlag;
   float         particleWeight;
   float         bx;
   float         by;
   float         ebx;
   float         eby;
   int           nParticle;
   float         px[MAXPARTICLE];
   float         py[MAXPARTICLE];
   float         pz[MAXPARTICLE];
   float         pt[MAXPARTICLE];
   float         pmag[MAXPARTICLE];
   float         rap[MAXPARTICLE];
   float         eta[MAXPARTICLE];
   float         theta[MAXPARTICLE];
   float         phi[MAXPARTICLE]; 
   float         mass[MAXPARTICLE];
   short         charge[MAXPARTICLE];
   bool          isCharged[MAXPARTICLE];
   short         pwflag[MAXPARTICLE];
   int           pid[MAXPARTICLE];
   float         d0[MAXPARTICLE];
   float         z0[MAXPARTICLE];
   bool          highPurity[MAXPARTICLE];
   short         ntpc[MAXPARTICLE];
   short         nitc[MAXPARTICLE];
   short         nvdet[MAXPARTICLE];
   float         vx[MAXPARTICLE];
   float         vy[MAXPARTICLE];
   float         vz[MAXPARTICLE];
   float         weight[MAXPARTICLE];
   float         pt_wrtThr[MAXPARTICLE];
   float         eta_wrtThr[MAXPARTICLE];
   float         rap_wrtThr[MAXPARTICLE];
   float         theta_wrtThr[MAXPARTICLE];
   float         phi_wrtThr[MAXPARTICLE];
   float         pt_wrtThrPerp[MAXPARTICLE];
   float         eta_wrtThrPerp[MAXPARTICLE];
   float         rap_wrtThrPerp[MAXPARTICLE];
   float         theta_wrtThrPerp[MAXPARTICLE];
   float         phi_wrtThrPerp[MAXPARTICLE];
   float         pt_wrtChThr[MAXPARTICLE];
   float         eta_wrtChThr[MAXPARTICLE];
   float         rap_wrtChThr[MAXPARTICLE];
   float         theta_wrtChThr[MAXPARTICLE];
   float         phi_wrtChThr[MAXPARTICLE];
   float         pt_wrtChThrPerp[MAXPARTICLE];
   float         eta_wrtChThrPerp[MAXPARTICLE];
   float         rap_wrtChThrPerp[MAXPARTICLE];
   float         theta_wrtChThrPerp[MAXPARTICLE];
   float         phi_wrtChThrPerp[MAXPARTICLE];
   float         pt_wrtNeuThr[MAXPARTICLE];
   float         eta_wrtNeuThr[MAXPARTICLE];
   float         rap_wrtNeuThr[MAXPARTICLE];
   float         theta_wrtNeuThr[MAXPARTICLE];
   float         phi_wrtNeuThr[MAXPARTICLE];
   float         pt_wrtNeuThrPerp[MAXPARTICLE];
   float         eta_wrtNeuThrPerp[MAXPARTICLE];
   float         theta_wrtNeuThrPerp[MAXPARTICLE];
   float         phi_wrtNeuThrPerp[MAXPARTICLE];
   float         rap_wrtNeuThrPerp[MAXPARTICLE];
   float         pt_wrtThrCorr[MAXPARTICLE];
   float         eta_wrtThrCorr[MAXPARTICLE];
   float         rap_wrtThrCorr[MAXPARTICLE];
   float         theta_wrtThrCorr[MAXPARTICLE];
   float         phi_wrtThrCorr[MAXPARTICLE];
   float         pt_wrtThrCorrPerp[MAXPARTICLE];
   float         eta_wrtThrCorrPerp[MAXPARTICLE];
   float         rap_wrtThrCorrPerp[MAXPARTICLE];
   float         theta_wrtThrCorrPerp[MAXPARTICLE];
   float         phi_wrtThrCorrPerp[MAXPARTICLE];
   float         pt_wrtThrCorrInverse[MAXPARTICLE];
   float         eta_wrtThrCorrInverse[MAXPARTICLE];
   float         rap_wrtThrCorrInverse[MAXPARTICLE];
   float         theta_wrtThrCorrInverse[MAXPARTICLE];
   float         phi_wrtThrCorrInverse[MAXPARTICLE];
   float         pt_wrtThrCorrInversePerp[MAXPARTICLE];
   float         eta_wrtThrCorrInversePerp[MAXPARTICLE];
   float         rap_wrtThrCorrInversePerp[MAXPARTICLE];
   float         theta_wrtThrCorrInversePerp[MAXPARTICLE];
   float         phi_wrtThrCorrInversePerp[MAXPARTICLE];
   float         pt_wrtThrMissP[MAXPARTICLE];
   float         eta_wrtThrMissP[MAXPARTICLE];
   float         rap_wrtThrMissP[MAXPARTICLE];
   float         theta_wrtThrMissP[MAXPARTICLE];
   float         phi_wrtThrMissP[MAXPARTICLE];
   float         pt_wrtThrMissPPerp[MAXPARTICLE];
   float         eta_wrtThrMissPPerp[MAXPARTICLE];
   float         rap_wrtThrMissPPerp[MAXPARTICLE];
   float         theta_wrtThrMissPPerp[MAXPARTICLE];
   float         phi_wrtThrMissPPerp[MAXPARTICLE];
   bool          passesArtificAccept[MAXPARTICLE];
   float         artificAcceptEffCorrection[MAXPARTICLE];
   bool          passesNTupleAfterCut;
   bool          passesTotalChgEnergyMin;
   bool          passesNTrkMin;
   bool          passesSTheta;
   bool          passesMissP;
   bool          passesISR;
   bool          passesWW;
   bool          passesNeuNch;
   bool          passesAll;
   float         missP;
   float         missPt;
   float         missTheta;
   float         missPhi;
   float         missChargedP;
   float         missChargedPt;
   float         missChargedTheta;
   float         missChargedPhi;
   int           nChargedHadrons;
   int           nChargedHadronsHP;
   float         nChargedHadronsHP_Corrected;
   int           nChargedHadrons_GT0p4;
   int           nChargedHadrons_GT0p4Thrust;
   int           nChargedParticleHP;
   float         Thrust;
   float         TTheta;
   float         TPhi;
   float         Thrust_charged;
   float         TTheta_charged;
   float         TPhi_charged;
   float         Thrust_neutral;
   float         TTheta_neutral;
   float         TPhi_neutral;
   float         ThrustCorr;
   float         TThetaCorr;
   float         TPhiCorr;
   float         ThrustCorrInverse;
   float         TThetaCorrInverse;
   float         TPhiCorrInverse;
   float         ThrustWithMissP;
   float         TThetaWithMissP;
   float         TPhiWithMissP;
   float         Sphericity;
   float         STheta;
   float         SPhi;
   float         Aplanarity;
   float         Sphericity_linearized;
   float         STheta_linearized;
   float         SPhi_linearized;
   float         Aplanarity_linearized;
   float         C_linearized;
   float         D_linearized;
   bool          passesLEP1TwoPC;
   float         ThrustWithReco;
   float         TThetaWithReco;
   float         TPhiWithReco;
   float         pt_wrtThrWithReco[MAXPARTICLE];
   float         eta_wrtThrWithReco[MAXPARTICLE];
   float         rap_wrtThrWithReco[MAXPARTICLE];
   float         theta_wrtThrWithReco[MAXPARTICLE];
   float         phi_wrtThrWithReco[MAXPARTICLE];
   float         ThrustWithGenIneff;
   float         TThetaWithGenIneff;
   float         TPhiWithGenIneff;
   float         pt_wrtThrWithGenIneff[MAXPARTICLE];
   float         eta_wrtThrWithGenIneff[MAXPARTICLE];
   float         rap_wrtThrWithGenIneff[MAXPARTICLE];
   float         theta_wrtThrWithGenIneff[MAXPARTICLE];
   float         phi_wrtThrWithGenIneff[MAXPARTICLE];
   float         ThrustWithGenIneffFake;
   float         TThetaWithGenIneffFake;
   float         TPhiWithGenIneffFake;
   float         pt_wrtThrWithGenIneffFake[MAXPARTICLE];
   float         eta_wrtThrWithGenIneffFake[MAXPARTICLE];
   float         rap_wrtThrWithGenIneffFake[MAXPARTICLE];
   float         theta_wrtThrWithGenIneffFake[MAXPARTICLE];
   float         phi_wrtThrWithGenIneffFake[MAXPARTICLE];
   float         ThrustWithRecoCorr;
   float         TThetaWithRecoCorr;
   float         TPhiWithRecoCorr;
   float         pt_wrtThrWithRecoCorr[MAXPARTICLE];
   float         eta_wrtThrWithRecoCorr[MAXPARTICLE];
   float         rap_wrtThrWithRecoCorr[MAXPARTICLE];
   float         theta_wrtThrWithRecoCorr[MAXPARTICLE];
   float         phi_wrtThrWithRecoCorr[MAXPARTICLE];
   float         ThrustWithRecoCorrInverse;
   float         TThetaWithRecoCorrInverse;
   float         TPhiWithRecoCorrInverse;
   float         pt_wrtThrWithRecoCorrInverse[MAXPARTICLE];
   float         eta_wrtThrWithRecoCorrInverse[MAXPARTICLE];
   float         rap_wrtThrWithRecoCorrInverse[MAXPARTICLE];
   float         theta_wrtThrWithRecoCorrInverse[MAXPARTICLE];
   float         phi_wrtThrWithRecoCorrInverse[MAXPARTICLE];
   float         ThrustWithRecoAndMissP;
   float         TThetaWithRecoAndMissP;
   float         TPhiWithRecoAndMissP;
   float         pt_wrtThrWithRecoAndMissP[MAXPARTICLE];
   float         eta_wrtThrWithRecoAndMissP[MAXPARTICLE];
   float         rap_wrtThrWithRecoAndMissP[MAXPARTICLE];
   float         theta_wrtThrWithRecoAndMissP[MAXPARTICLE];
   float         phi_wrtThrWithRecoAndMissP[MAXPARTICLE];
public:
   std::vector<FourVector> P;
public:
   ParticleTreeMessenger();
   ParticleTreeMessenger(TFile &file, std::string name);
   ParticleTreeMessenger(TFile *file, std::string name);
   ParticleTreeMessenger(TTree *tree);
   bool Initialize(TTree *tree);
   bool Initialize();
   bool GetEntry(int iEntry);
   int GetEntries();
   bool PassBaselineCut();
};

class ReducedTreeMessenger
{
public:
   TTree *Tree;
   int    N;
   float  Momentum[MAXPARTICLE];
   float  Mass[MAXPARTICLE];
   float  Theta[MAXPARTICLE];
   float  Phi[MAXPARTICLE];
   float  Weight[MAXPARTICLE];
   short  Charge[MAXPARTICLE];
   bool   PassCut;
public:
   std::vector<FourVector> P;
public:
   ReducedTreeMessenger();
   ReducedTreeMessenger(TFile &file, std::string name);
   ReducedTreeMessenger(TFile *file, std::string name);
   ReducedTreeMessenger(TTree *tree);
   bool Initialize(TTree *tree);
   bool Initialize();
   bool GetEntry(int iEntry);
   int GetEntries();
};


