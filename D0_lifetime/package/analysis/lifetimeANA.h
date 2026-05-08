//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Dec  8 17:01:22 2022 by ROOT version 6.26/00
// from TTree t_M0pipi/a Tree with data
// found on file: tree_DKpi_data_S.root
//////////////////////////////////////////////////////////

#ifndef lifetimeANA_h
#define lifetimeANA_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class lifetimeANA {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           id;
   Int_t           event;
   Double_t        h1_px;
   Double_t        h1_py;
   Double_t        h1_pz;
   Double_t        h1_px_true;
   Double_t        h1_py_true;
   Double_t        h1_pz_true;
   Double_t        h1_pt;
   Double_t        h1_p;
   Double_t        h1_eta;
   Double_t        h1_IP;
   Double_t        h1_charge;
   Double_t        h1_thetaC0;
   Double_t        h1_thetaC1;
   Double_t        h1_thetaC2;
   Double_t        h1_DLL;
   Double_t        h2_px;
   Double_t        h2_py;
   Double_t        h2_pz;
   Double_t        h2_px_true;
   Double_t        h2_py_true;
   Double_t        h2_pz_true;
   Double_t        h2_pt;
   Double_t        h2_p;
   Double_t        h2_eta;
   Double_t        h2_IP;
   Double_t        h2_charge;
   Double_t        h2_thetaC0;
   Double_t        h2_thetaC1;
   Double_t        h2_thetaC2;
   Double_t        h2_DLL;
   Double_t        piS_px;
   Double_t        piS_py;
   Double_t        piS_px_true;
   Double_t        piS_py_true;
   Double_t        piS_pt;
   Double_t        piS_pz;
   Double_t        piS_pz_true;
   Double_t        piS_p;
   Double_t        piS_eta;
   Double_t        piS_IP;
   Double_t        piS_charge;
   Double_t        PVx;
   Double_t        PVy;
   Double_t        PVz;
   Double_t        PVx_true;
   Double_t        PVy_true;
   Double_t        PVz_true;
   Double_t        SVx;
   Double_t        SVy;
   Double_t        SVz;
   Double_t        SVx_true;
   Double_t        SVy_true;
   Double_t        SVz_true;
   Double_t        M0_px;
   Double_t        M0_py;
   Double_t        M0_pt;
   Double_t        M0_pz;
   Double_t        M0_p;
   Double_t        M0_eta;
   Double_t        M0_FDx;
   Double_t        M0_FDy;
   Double_t        M0_FDt;
   Double_t        M0_FDz;
   Double_t        M0_FD;
   Double_t        M0_IP;
   Double_t        M0_time;
   Double_t        M0_time_true;
   Double_t        M0_Mpipi;
   Double_t        M0_MKK;
   Double_t        M0_MKpi;
   Double_t        M0_MpiK;

   // List of branches
   TBranch        *b_id;   //!
   TBranch        *b_event;   //!
   TBranch        *b_h1_px;   //!
   TBranch        *b_h1_py;   //!
   TBranch        *b_h1_pz;   //!
   TBranch        *b_h1_px_true;   //!
   TBranch        *b_h1_py_true;   //!
   TBranch        *b_h1_pz_true;   //!
   TBranch        *b_h1_pt;   //!
   TBranch        *b_h1_p;   //!
   TBranch        *b_h1_eta;   //!
   TBranch        *b_h1_IP;   //!
   TBranch        *b_h1_charge;   //!
   TBranch        *b_h1_thetaC0;   //!
   TBranch        *b_h1_thetaC1;   //!
   TBranch        *b_h1_thetaC2;   //!
   TBranch        *b_h1_DLL;   //!
   TBranch        *b_h2_px;   //!
   TBranch        *b_h2_py;   //!
   TBranch        *b_h2_pz;   //!
   TBranch        *b_h2_px_true;   //!
   TBranch        *b_h2_py_true;   //!
   TBranch        *b_h2_pz_true;   //!
   TBranch        *b_h2_pt;   //!
   TBranch        *b_h2_p;   //!
   TBranch        *b_h2_eta;   //!
   TBranch        *b_h2_IP;   //!
   TBranch        *b_h2_charge;   //!
   TBranch        *b_h2_thetaC0;   //!
   TBranch        *b_h2_thetaC1;   //!
   TBranch        *b_h2_thetaC2;   //!
   TBranch        *b_h2_DLL;   //!
   TBranch        *b_piS_px;   //!
   TBranch        *b_piS_py;   //!
   TBranch        *b_piS_px_true;   //!
   TBranch        *b_piS_py_true;   //!
   TBranch        *b_piS_pt;   //!
   TBranch        *b_piS_pz;   //!
   TBranch        *b_piS_pz_true;   //!
   TBranch        *b_piS_p;   //!
   TBranch        *b_piS_eta;   //!
   TBranch        *b_piS_IP;   //!
   TBranch        *b_piS_charge ;   //!
   TBranch        *b_PVx;   //!
   TBranch        *b_PVy;   //!
   TBranch        *b_PVz;   //!
   TBranch        *b_PVx_true;   //!
   TBranch        *b_PVy_true;   //!
   TBranch        *b_PVz_true;   //!
   TBranch        *b_SVx;   //!
   TBranch        *b_SVy;   //!
   TBranch        *b_SVz;   //!
   TBranch        *b_SVx_true;   //!
   TBranch        *b_SVy_true;   //!
   TBranch        *b_SVz_true;   //!
   TBranch        *b_M0_px;   //!
   TBranch        *b_M0_py;   //!
   TBranch        *b_M0_pt;   //!
   TBranch        *b_M0_pz;   //!
   TBranch        *b_M0_p;   //!
   TBranch        *b_M0_eta;   //!
   TBranch        *b_M0_FDx;   //!
   TBranch        *b_M0_FDy;   //!
   TBranch        *b_M0_FDt;   //!
   TBranch        *b_M0_FDz;   //!
   TBranch        *b_M0_FD;   //!
   TBranch        *b_M0_IP;   //!
   TBranch        *b_M0_time;   //!
   TBranch        *b_M0_time_true;   //!
   TBranch        *b_M0_Mpipi;   //!
   TBranch        *b_M0_MKK;   //!
   TBranch        *b_M0_MKpi;   //!
   TBranch        *b_M0_MpiK;   //!

   lifetimeANA(TTree *tree=0);
   virtual ~lifetimeANA();
   //   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef lifetimeANA_cxx
lifetimeANA::lifetimeANA(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   /* if (tree == 0) { */
   /*    TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("tree_DKpi_data_S.root"); */
   /*    if (!f || !f->IsOpen()) { */
   /*       f = new TFile("tree_DKpi_data_S.root"); */
   /*    } */
   /*    f->GetObject("t_M0pipi",tree); */

   /* } */
   /* Init(tree); */

  TChain * chain = new TChain("t_M0pipi","");
  chain->Add("../data/tree_DKPi_data_SB.root");
  chain->Add("../mc/tree_DKPi_mc_410.root");
  
  tree = chain; 
  Init(tree);
  
}

lifetimeANA::~lifetimeANA()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t lifetimeANA::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t lifetimeANA::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void lifetimeANA::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("id", &id, &b_id);
   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("h1_px", &h1_px, &b_h1_px);
   fChain->SetBranchAddress("h1_py", &h1_py, &b_h1_py);
   fChain->SetBranchAddress("h1_pz", &h1_pz, &b_h1_pz);
   fChain->SetBranchAddress("h1_px_true", &h1_px_true, &b_h1_px_true);
   fChain->SetBranchAddress("h1_py_true", &h1_py_true, &b_h1_py_true);
   fChain->SetBranchAddress("h1_pz_true", &h1_pz_true, &b_h1_pz_true);
   fChain->SetBranchAddress("h1_pt", &h1_pt, &b_h1_pt);
   fChain->SetBranchAddress("h1_p", &h1_p, &b_h1_p);
   fChain->SetBranchAddress("h1_eta", &h1_eta, &b_h1_eta);
   fChain->SetBranchAddress("h1_IP", &h1_IP, &b_h1_IP);
   fChain->SetBranchAddress("h1_charge", &h1_charge, &b_h1_charge);
   fChain->SetBranchAddress("h1_thetaC0", &h1_thetaC0, &b_h1_thetaC0);
   fChain->SetBranchAddress("h1_thetaC1", &h1_thetaC1, &b_h1_thetaC1);
   fChain->SetBranchAddress("h1_thetaC2", &h1_thetaC2, &b_h1_thetaC2);
   fChain->SetBranchAddress("h1_DLL", &h1_DLL, &b_h1_DLL);
   fChain->SetBranchAddress("h2_px", &h2_px, &b_h2_px);
   fChain->SetBranchAddress("h2_py", &h2_py, &b_h2_py);
   fChain->SetBranchAddress("h2_pz", &h2_pz, &b_h2_pz);
   fChain->SetBranchAddress("h2_px_true", &h2_px_true, &b_h2_px_true);
   fChain->SetBranchAddress("h2_py_true", &h2_py_true, &b_h2_py_true);
   fChain->SetBranchAddress("h2_pz_true", &h2_pz_true, &b_h2_pz_true);
   fChain->SetBranchAddress("h2_pt", &h2_pt, &b_h2_pt);
   fChain->SetBranchAddress("h2_p", &h2_p, &b_h2_p);
   fChain->SetBranchAddress("h2_eta", &h2_eta, &b_h2_eta);
   fChain->SetBranchAddress("h2_IP", &h2_IP, &b_h2_IP);
   fChain->SetBranchAddress("h2_charge", &h2_charge, &b_h2_charge);
   fChain->SetBranchAddress("h2_thetaC0", &h2_thetaC0, &b_h2_thetaC0);
   fChain->SetBranchAddress("h2_thetaC1", &h2_thetaC1, &b_h2_thetaC1);
   fChain->SetBranchAddress("h2_thetaC2", &h2_thetaC2, &b_h2_thetaC2);
   fChain->SetBranchAddress("h2_DLL", &h2_DLL, &b_h2_DLL);
   fChain->SetBranchAddress("piS_px", &piS_px, &b_piS_px);
   fChain->SetBranchAddress("piS_py", &piS_py, &b_piS_py);
   fChain->SetBranchAddress("piS_px_true", &piS_px_true, &b_piS_px_true);
   fChain->SetBranchAddress("piS_py_true", &piS_py_true, &b_piS_py_true);
   fChain->SetBranchAddress("piS_pt", &piS_pt, &b_piS_pt);
   fChain->SetBranchAddress("piS_pz", &piS_pz, &b_piS_pz);
   fChain->SetBranchAddress("piS_pz_true", &piS_pz_true, &b_piS_pz_true);
   fChain->SetBranchAddress("piS_p", &piS_p, &b_piS_p);
   fChain->SetBranchAddress("piS_eta", &piS_eta, &b_piS_eta);
   fChain->SetBranchAddress("piS_IP", &piS_IP, &b_piS_IP);
   fChain->SetBranchAddress("piS_charge", &piS_charge, &b_piS_charge );
   fChain->SetBranchAddress("PVx", &PVx, &b_PVx);
   fChain->SetBranchAddress("PVy", &PVy, &b_PVy);
   fChain->SetBranchAddress("PVz", &PVz, &b_PVz);
   fChain->SetBranchAddress("PVx_true", &PVx_true, &b_PVx_true);
   fChain->SetBranchAddress("PVy_true", &PVy_true, &b_PVy_true);
   fChain->SetBranchAddress("PVz_true", &PVz_true, &b_PVz_true);
   fChain->SetBranchAddress("SVx", &SVx, &b_SVx);
   fChain->SetBranchAddress("SVy", &SVy, &b_SVy);
   fChain->SetBranchAddress("SVz", &SVz, &b_SVz);
   fChain->SetBranchAddress("SVx_true", &SVx_true, &b_SVx_true);
   fChain->SetBranchAddress("SVy_true", &SVy_true, &b_SVy_true);
   fChain->SetBranchAddress("SVz_true", &SVz_true, &b_SVz_true);
   fChain->SetBranchAddress("M0_px", &M0_px, &b_M0_px);
   fChain->SetBranchAddress("M0_py", &M0_py, &b_M0_py);
   fChain->SetBranchAddress("M0_pt", &M0_pt, &b_M0_pt);
   fChain->SetBranchAddress("M0_pz", &M0_pz, &b_M0_pz);
   fChain->SetBranchAddress("M0_p", &M0_p, &b_M0_p);
   fChain->SetBranchAddress("M0_eta", &M0_eta, &b_M0_eta);
   fChain->SetBranchAddress("M0_FDx", &M0_FDx, &b_M0_FDx);
   fChain->SetBranchAddress("M0_FDy", &M0_FDy, &b_M0_FDy);
   fChain->SetBranchAddress("M0_FDt", &M0_FDt, &b_M0_FDt);
   fChain->SetBranchAddress("M0_FDz", &M0_FDz, &b_M0_FDz);
   fChain->SetBranchAddress("M0_FD", &M0_FD, &b_M0_FD);
   fChain->SetBranchAddress("M0_IP", &M0_IP, &b_M0_IP);
   fChain->SetBranchAddress("M0_time", &M0_time, &b_M0_time);
   fChain->SetBranchAddress("M0_time_true", &M0_time_true, &b_M0_time_true);
   fChain->SetBranchAddress("M0_Mpipi", &M0_Mpipi, &b_M0_Mpipi);
   fChain->SetBranchAddress("M0_MKK", &M0_MKK, &b_M0_MKK);
   fChain->SetBranchAddress("M0_MKpi", &M0_MKpi, &b_M0_MKpi);
   fChain->SetBranchAddress("M0_MpiK", &M0_MpiK, &b_M0_MpiK);
   Notify();
}

Bool_t lifetimeANA::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void lifetimeANA::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
//Int_t lifetimeANA::Cut(Long64_t entry)
//{
//// This function may be called from Loop.
//// returns  1 if entry is accepted.
//// returns -1 otherwise.
//   return 1;
//}
#endif // #ifdef lifetimeANA_cxx
