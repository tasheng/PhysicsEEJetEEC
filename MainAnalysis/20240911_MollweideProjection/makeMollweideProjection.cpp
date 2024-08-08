#include <iostream>
#include <vector>
#include <map>
using namespace std;

#include "TTree.h"
#include "TFile.h"
#include "TStyle.h"
#include "Messenger.h"
#include "CommandLine.h"
#include "Matching.h"
#include "ProgressBar.h"
#include "TauHelperFunctions3.h"
#include "alephTrkEfficiency.h"


#include "TCanvas.h"
#include "TH2D.h"
#include "TEllipse.h"
#include "TMath.h"
#include "TVector3.h"


int main(int argc, char *argv[]);
double mollweideX(double lambda, double theta); 
double mollweideY(double theta); 

int main(int argc, char *argv[])
{
   CommandLine CL(argc, argv);

   string InputFileName  = CL.Get("Input");
   string GenTreeName    = CL.Get("Gen", "tgen");
   bool isLEP1 = CL.GetBool("DoLEP1", true);


   TFile* InputFile  = new TFile(InputFileName.c_str());


   ParticleTreeMessenger* MGen = new ParticleTreeMessenger(InputFile, GenTreeName);
   const int nBinsX = 200;
   const int nBinsY = 100;
   TH2D *hCMBplot = new TH2D("hCMBplot", "hCMBplotY", nBinsX, -TMath::Pi(), TMath::Pi(), nBinsY, -TMath::Pi()/2, TMath::Pi()/2);


   int EntryCount = MGen->GetEntries();
   ProgressBar Bar(cout, EntryCount);
   Bar.SetStyle(-1); 
   bool foundHighMultEvent = false; 
   for(int iE = 0; iE < EntryCount; iE++) 
   {

      MGen->GetEntry(iE);

      vector<FourVector> PGen;
      for(int i = 0; i < MGen->nParticle; i++){
         
        // charged particle selection 
       if(MGen->charge[i] == 0) continue;
         if(MGen->highPurity[i] == false) continue;
        TVector3* thrustVec = new TVector3(); 
        thrustVec->SetPtEtaPhi(MGen->pt_wrtThrMissP[i], MGen->eta_wrtThrMissP[i], MGen->phi_wrtThrMissP[i]); 
        double px = thrustVec->Px(); 
        double py = thrustVec->Py(); 
        double pz = thrustVec->Pz(); 
        double lambda = TMath::ATan2(py, px); // Azimuthal angle (longitude)
        double r = TMath::Sqrt(px * px + py * py + pz*pz); 
        double theta = TMath::ASin(pz/r); // latitude
        double x = mollweideX(lambda, theta);
        double y = mollweideY(theta);
        hCMBplot->Fill(x, y,MGen->P[i][0]); // Using energy as the weight

      }

   } // end loop over the number of events
 
    // now draw the canvas
    gStyle->SetOptStat(0);
    gStyle->SetPalette(kTemperatureMap);

    // Create a pad for the Mollweide projection
    TCanvas *canvas = new TCanvas("canvas", "CMB Mollweide Plot", 1000, 500);
    canvas->SetFillColor(1);
    canvas->SetLogz(); 

    // Draw the histogram with the Mollweide projection
    hCMBplot->GetZaxis()->SetLabelColor(kWhite); 
    hCMBplot->Draw("COLZ");

    // Create an ellipse to outline the Mollweide projection
    TEllipse *ellipse = new TEllipse(0, 0, TMath::Pi(), TMath::Pi()/2);
    ellipse->SetFillStyle(0);
    ellipse->SetLineColor(kWhite);
    ellipse->SetLineWidth(2);
    ellipse->Draw();   
    // Save the plot to a file
    if(isLEP1)canvas->SaveAs("EEC_MollweidePlot_LEP1.pdf");
    else canvas->SaveAs("EEC_MollweidePlot_LEP2.pdf");

   // cleanup
   InputFile->Close(); 
   delete InputFile; 
   delete MGen; 

    

   return 0;
}


double mollweideX(double lambda, double theta) {
    return 2 * sqrt(2) * lambda * TMath::Cos(theta) / (TMath::Pi());// * sqrt(1 + TMath::Cos(theta)));
}

double mollweideY(double theta) {
    return sqrt(2) * TMath::Sin(theta);// / sqrt(1 + TMath::Cos(theta));
}