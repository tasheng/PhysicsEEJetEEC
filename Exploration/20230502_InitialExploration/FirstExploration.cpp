#include <iostream>
using namespace std;

#include "TFile.h"
#include "TTree.h"

#include "ProgressBar.h"
#include "CommandLine.h"
#include "Messenger.h"

int main(int argc, char *argv[])
{
   CommandLine CL(argc, argv);

   string InputFileName = CL.Get("Input");
   string JetTreeName   = CL.Get("Jet", "akR4ESchemeJetTree");

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
   }
   Bar.Update(EntryCount);
   Bar.Print();
   Bar.PrintLine();

   File.Close();

   return 0;
}



