#define lifetimeANA_cxx
#include "lifetimeANA.h"

#include <TCanvas.h>
#include <TH1.h>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <string>
#include <TRandom3.h>  

using namespace std;

void lifetimeANA::Loop()
{
//   In a ROOT session, you can do:
//      root> .L lifetimeANA.C
//      root> lifetimeANA t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//



  // Create histograms
  TH1D *histo_data_MKpi    = new TH1D("histo_data_MKpi","",100,1.8,1.95);
  TH1D *histo_mc_MKpi      = new TH1D("histo_mc_MKpi","",100,1.8,1.95);

  histo_data_MKpi        ->Sumw2();
  histo_mc_MKpi          ->Sumw2();
  
  

   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

      // #### My code #########
      if ( id==1){ //DATA
	
	histo_data_MKpi ->Fill(M0_MKpi);
	// Print out the Mpipi values
	//std::cout <<M0_MKpi << std::endl;
	
      }; //end DATA
      
      if (id==13){ //MC
	histo_mc_MKpi ->Fill(M0_MKpi);
	//std::cout <<M0_MKpi << std::endl;
	
      }; //end MC
      
   };// #### end loop over jentry

    //Create Canvas
    TCanvas *c_histo_data_MKpi = new TCanvas("c_histo_data_MKpi","canvas histo",500,500);
    c_histo_data_MKpi->cd();
    histo_data_MKpi->Draw();
    c_histo_data_MKpi->Print("./_fig/c_histo_data_MKpi.pdf");
    c_histo_data_MKpi->Print("./_fig/c_histo_data_MKpi.eps");


    //Create Canvas
    TCanvas *c_histo_mc_MKpi = new TCanvas("c_histo_mc_MKpi","canvas histo",500,500);
    c_histo_mc_MKpi->cd();
    histo_mc_MKpi->Draw();
    c_histo_mc_MKpi->Print("./_fig/c_histo_mc_MKpi.pdf");
    c_histo_mc_MKpi->Print("./_fig/c_histo_mc_MKpi.eps");


    
  //Create a new file to store histograms
   TFile *histo_file = new TFile("./_root/histo_file.root","RECREATE","put a title");
   histo_file->cd();
   histo_data_MKpi->Write();
   histo_mc_MKpi->Write();
   histo_file->Close();

}
