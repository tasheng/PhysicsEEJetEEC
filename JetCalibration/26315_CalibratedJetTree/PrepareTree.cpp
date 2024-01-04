#include <iostream>
#include <vector>
using namespace std;

#include "TFile.h"
#include "TTree.h"

#include "ProgressBar.h"
#include "CommandLine.h"

int main(int argc, char *argv[])
{
   CommandLine CL(argc, argv);

   string InputFileName     = CL.Get("Input");
   string JetInputFileName  = CL.Get("JetInput");
   string OutputFileName    = CL.Get("Output");

   bool IsGen               = CL.GetBool("IsGen", false);
   string ParticleTreeName  = CL.Get("Particle", "t");

   string InputTreeName     = CL.Get("InputTree", "RecoR4");
   string OutputTreeName    = CL.Get("OutputTree", "akR4ESchemeReclusteredJetTree");

   TFile InputFile(InputFileName.c_str());
   TFile JetInputFile(JetInputFileName.c_str());
   TFile OutputFile(OutputFileName.c_str(), "RECREATE");

   ParticleTreeMessenger MP(InputFile, ParticleTreeName);
   JetTreeMessenger MJ(InputFile, InputTreeName);

   int EntryCount = MPReco.GetEntries();
   for(int iE = 0; iE < EntryCount; iE++)
   {
      MP.GetEntry(iE);
      MJ.GetEntry(iE);
   }

   OutputFile.Close();
   JetInputFile.Close();
   InputFile.Close();

   return 0;
}





