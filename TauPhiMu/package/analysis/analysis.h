//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Jan  6 13:48:37 2026 by ROOT version 6.34.08
// from TTree t_M0pipi/a Tree with data
// found on file: tree_M0hh_DPLUSPhiPi.root
//////////////////////////////////////////////////////////

#ifndef analysis_h
#define analysis_h

#include <iostream>
#include <stdio.h>
#include <stdlib.h>

#include <TChain.h>
#include <TFile.h>
#include <TH1D.h>
#include <TROOT.h>

using namespace std;

// User-defined struct for storing fit results
struct AuxFitResult
{
    std::vector<Double_t> params;
    std::vector<Double_t> errors;
    std::vector<std::vector<Double_t>> covMatrix;
    bool isValid = false;
};

// Header file for the classes stored in the TTree if any.

class analysis
{
  public:
    TTree *fChain; //! pointer to the analyzed TTree or TChain
    Int_t fCurrent; //! current Tree number in a TChain

    // Fixed size dimensions of array or collections stored in the TTree if any.

    // Declaration of leaf types
    Int_t id;
    Int_t event;
    Double_t weight;
    Double_t h1_px;
    Double_t h1_py;
    Double_t h1_pz;
    Double_t h1_px_true;
    Double_t h1_py_true;
    Double_t h1_pz_true;
    Double_t h1_pt;
    Double_t h1_p;
    Double_t h1_eta;
    Double_t h1_IP;
    Double_t h1_charge;
    Int_t h1_MuonID;
    Double_t h2_px;
    Double_t h2_py;
    Double_t h2_pz;
    Double_t h2_px_true;
    Double_t h2_py_true;
    Double_t h2_pz_true;
    Double_t h2_pt;
    Double_t h2_p;
    Double_t h2_eta;
    Double_t h2_IP;
    Double_t h2_charge;
    Int_t h2_MuonID;
    Double_t h3_px;
    Double_t h3_py;
    Double_t h3_px_true;
    Double_t h3_py_true;
    Double_t h3_pt;
    Double_t h3_pz;
    Double_t h3_pz_true;
    Double_t h3_p;
    Double_t h3_eta;
    Double_t h3_IP;
    Double_t h3_charge;
    Int_t h3_MuonID;
    Double_t PVx;
    Double_t PVy;
    Double_t PVz;
    Double_t PVx_true;
    Double_t PVy_true;
    Double_t PVz_true;
    Double_t DVx;
    Double_t DVy;
    Double_t DVz;
    Double_t DVx_true;
    Double_t DVy_true;
    Double_t DVz_true;
    Double_t M0_px;
    Double_t M0_py;
    Double_t M0_pt;
    Double_t M0_pz;
    Double_t M0_p;
    Double_t M0_eta;
    Double_t M0_FDx;
    Double_t M0_FDy;
    Double_t M0_FDt;
    Double_t M0_FDz;
    Double_t M0_FD;
    Double_t M0_IP;
    Double_t M0_time;
    Double_t M0_time_true;
    Double_t M0_Mpipi;
    Double_t M0_MKK;
    Double_t M0_MKpi;
    Double_t M0_MpiK;
    Double_t D_Vx;
    Double_t D_Vy;
    Double_t D_Vz;
    Double_t D_Vx_true;
    Double_t D_Vy_true;
    Double_t D_Vz_true;
    Double_t D_px;
    Double_t D_py;
    Double_t D_pt;
    Double_t D_pz;
    Double_t D_p;
    Double_t D_eta;
    Double_t D_FDx;
    Double_t D_FDy;
    Double_t D_FDt;
    Double_t D_FDz;
    Double_t D_FD;
    Double_t D_IP;
    Double_t D_time;
    Double_t D_time_true;
    Double_t D_M;

    // List of branches
    TBranch *b_id; //!
    TBranch *b_event; //!
    TBranch *b_weight; //!
    TBranch *b_h1_px; //!
    TBranch *b_h1_py; //!
    TBranch *b_h1_pz; //!
    TBranch *b_h1_px_true; //!
    TBranch *b_h1_py_true; //!
    TBranch *b_h1_pz_true; //!
    TBranch *b_h1_pt; //!
    TBranch *b_h1_p; //!
    TBranch *b_h1_eta; //!
    TBranch *b_h1_IP; //!
    TBranch *b_h1_charge; //!
    TBranch *b_h1_MuonID; //!
    TBranch *b_h2_px; //!
    TBranch *b_h2_py; //!
    TBranch *b_h2_pz; //!
    TBranch *b_h2_px_true; //!
    TBranch *b_h2_py_true; //!
    TBranch *b_h2_pz_true; //!
    TBranch *b_h2_pt; //!
    TBranch *b_h2_p; //!
    TBranch *b_h2_eta; //!
    TBranch *b_h2_IP; //!
    TBranch *b_h2_charge; //!
    TBranch *b_h2_MuonID; //!
    TBranch *b_h3_px; //!
    TBranch *b_h3_py; //!
    TBranch *b_h3_px_true; //!
    TBranch *b_h3_py_true; //!
    TBranch *b_h3_pt; //!
    TBranch *b_h3_pz; //!
    TBranch *b_h3_pz_true; //!
    TBranch *b_h3_p; //!
    TBranch *b_h3_eta; //!
    TBranch *b_h3_IP; //!
    TBranch *b_h3_charge; //!
    TBranch *b_h3_MuonID; //!
    TBranch *b_PVx; //!
    TBranch *b_PVy; //!
    TBranch *b_PVz; //!
    TBranch *b_PVx_true; //!
    TBranch *b_PVy_true; //!
    TBranch *b_PVz_true; //!
    TBranch *b_DVx; //!
    TBranch *b_DVy; //!
    TBranch *b_DVz; //!
    TBranch *b_DVx_true; //!
    TBranch *b_DVy_true; //!
    TBranch *b_DVz_true; //!
    TBranch *b_M0_px; //!
    TBranch *b_M0_py; //!
    TBranch *b_M0_pt; //!
    TBranch *b_M0_pz; //!
    TBranch *b_M0_p; //!
    TBranch *b_M0_eta; //!
    TBranch *b_M0_FDx; //!
    TBranch *b_M0_FDy; //!
    TBranch *b_M0_FDt; //!
    TBranch *b_M0_FDz; //!
    TBranch *b_M0_FD; //!
    TBranch *b_M0_IP; //!
    TBranch *b_M0_time; //!
    TBranch *b_M0_time_true; //!
    TBranch *b_M0_Mpipi; //!
    TBranch *b_M0_MKK; //!
    TBranch *b_M0_MKpi; //!
    TBranch *b_M0_MpiK; //!
    TBranch *b_D_Vx; //!
    TBranch *b_D_Vy; //!
    TBranch *b_D_Vz; //!
    TBranch *b_D_Vx_true; //!
    TBranch *b_D_Vy_true; //!
    TBranch *b_D_Vz_true; //!
    TBranch *b_D_px; //!
    TBranch *b_D_py; //!
    TBranch *b_D_pt; //!
    TBranch *b_D_pz; //!
    TBranch *b_D_p; //!
    TBranch *b_D_eta; //!
    TBranch *b_D_FDx; //!
    TBranch *b_D_FDy; //!
    TBranch *b_D_FDt; //!
    TBranch *b_D_FDz; //!
    TBranch *b_D_FD; //!
    TBranch *b_D_IP; //!
    TBranch *b_D_time; //!
    TBranch *b_D_time_true; //!
    TBranch *b_D_M; //!

    analysis(bool isMC = true);
    virtual ~analysis();
    void LoadDataset(bool isMC);

    // virtual Int_t    Cut(Long64_t entry);
    virtual Int_t GetEntry(Long64_t entry);
    virtual Long64_t LoadTree(Long64_t entry);
    virtual void Init(TTree *tree);
    virtual void Loop();
    virtual bool Notify();
    virtual void Show(Long64_t entry = -1);

    // User functions
    AuxFitResult FitTemplateMass(Int_t mcID);
    AuxFitResult FitCombinatorialBkg();
    void DoFullBlindedUnbinnedFit();

    // Helpers for blind analysis
    inline TH1D *GetBlindedClone(TH1D *h, Double_t blindMin, Double_t blindMax);
    inline void DrawBlindedFunction(
        TF1 *f, Double_t blindMin, Double_t blindMax, Option_t *option = "SAME");
};

#endif

#ifdef analysis_cxx
// Il costruttore ora delega interamente il lavoro a LoadDataset
analysis::analysis(bool isMC)
    : fChain(nullptr)
{
    LoadDataset(isMC);
}

// Distruttore aggiornato per eliminare in sicurezza la TChain
analysis::~analysis()
{
    if(fChain)
    {
        delete fChain;
        fChain = nullptr;
    }
}

void analysis::LoadDataset(bool isMC)
{
    // Se esiste già una catena attiva, eliminala in sicurezza per liberare memoria
    if(fChain)
    {
        delete fChain;
        fChain = nullptr;
    }

    TChain *chain = new TChain("t_M0pipi", "");

    if(isMC)
    {
        chain->Add("../mc/tree_DSPLUS_PhiPi_mc_50M.root");
        chain->Add("../mc/tree_DPLUS_PhiPi_mc_50M.root");
        chain->Add("../mc/tree_DSPLUS_PhiMuNu_mc_5M.root");
        chain->Add("../mc/tree_DSPLUS_TauNu_mc_1M.root");
        cout << "[analysis] Dataset impostato su: MONTE CARLO" << endl;
    }
    else
    {
        chain->Add("../data/tree_data.root");
        cout << "[analysis] Dataset impostato su: DATI REALI" << endl;
    }

    // Esegue l'inizializzazione del nuovo albero/catena
    Init(chain);
}

Int_t analysis::GetEntry(Long64_t entry)
{
    // Read contents of entry.
    if(!fChain)
        return 0;
    return fChain->GetEntry(entry);
}
Long64_t analysis::LoadTree(Long64_t entry)
{
    // Set the environment to read one entry
    if(!fChain)
        return -5;
    Long64_t centry = fChain->LoadTree(entry);
    if(centry < 0)
        return centry;
    if(fChain->GetTreeNumber() != fCurrent)
    {
        fCurrent = fChain->GetTreeNumber();
        Notify();
    }
    return centry;
}

void analysis::Init(TTree *tree)
{
    // The Init() function is called when the selector needs to initialize
    // a new tree or chain. Typically here the branch addresses and branch
    // pointers of the tree will be set.
    // It is normally not necessary to make changes to the generated
    // code, but the routine can be extended by the user if needed.
    // Init() will be called many times when running on PROOF
    // (once per file to be processed).

    // Set branch addresses and branch pointers
    if(!tree)
        return;
    fChain = tree;
    fCurrent = -1;
    fChain->SetMakeClass(1);

    fChain->SetBranchAddress("id", &id, &b_id);
    fChain->SetBranchAddress("event", &event, &b_event);
    fChain->SetBranchAddress("weight", &weight, &b_weight);
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
    fChain->SetBranchAddress("h1_MuonID", &h1_MuonID, &b_h1_MuonID);
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
    fChain->SetBranchAddress("h2_MuonID", &h2_MuonID, &b_h2_MuonID);
    fChain->SetBranchAddress("h3_px", &h3_px, &b_h3_px);
    fChain->SetBranchAddress("h3_py", &h3_py, &b_h3_py);
    fChain->SetBranchAddress("h3_px_true", &h3_px_true, &b_h3_px_true);
    fChain->SetBranchAddress("h3_py_true", &h3_py_true, &b_h3_py_true);
    fChain->SetBranchAddress("h3_pt", &h3_pt, &b_h3_pt);
    fChain->SetBranchAddress("h3_pz", &h3_pz, &b_h3_pz);
    fChain->SetBranchAddress("h3_pz_true", &h3_pz_true, &b_h3_pz_true);
    fChain->SetBranchAddress("h3_p", &h3_p, &b_h3_p);
    fChain->SetBranchAddress("h3_eta", &h3_eta, &b_h3_eta);
    fChain->SetBranchAddress("h3_IP", &h3_IP, &b_h3_IP);
    fChain->SetBranchAddress("h3_charge", &h3_charge, &b_h3_charge);
    fChain->SetBranchAddress("h3_MuonID", &h3_MuonID, &b_h3_MuonID);
    fChain->SetBranchAddress("PVx", &PVx, &b_PVx);
    fChain->SetBranchAddress("PVy", &PVy, &b_PVy);
    fChain->SetBranchAddress("PVz", &PVz, &b_PVz);
    fChain->SetBranchAddress("PVx_true", &PVx_true, &b_PVx_true);
    fChain->SetBranchAddress("PVy_true", &PVy_true, &b_PVy_true);
    fChain->SetBranchAddress("PVz_true", &PVz_true, &b_PVz_true);
    fChain->SetBranchAddress("DVx", &DVx, &b_DVx);
    fChain->SetBranchAddress("DVy", &DVy, &b_DVy);
    fChain->SetBranchAddress("DVz", &DVz, &b_DVz);
    fChain->SetBranchAddress("DVx_true", &DVx_true, &b_DVx_true);
    fChain->SetBranchAddress("DVy_true", &DVy_true, &b_DVy_true);
    fChain->SetBranchAddress("DVz_true", &DVz_true, &b_DVz_true);
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
    fChain->SetBranchAddress("D_Vx", &D_Vx, &b_D_Vx);
    fChain->SetBranchAddress("D_Vy", &D_Vy, &b_D_Vy);
    fChain->SetBranchAddress("D_Vz", &D_Vz, &b_D_Vz);
    fChain->SetBranchAddress("D_Vx_true", &D_Vx_true, &b_D_Vx_true);
    fChain->SetBranchAddress("D_Vy_true", &D_Vy_true, &b_D_Vy_true);
    fChain->SetBranchAddress("D_Vz_true", &D_Vz_true, &b_D_Vz_true);
    fChain->SetBranchAddress("D_px", &D_px, &b_D_px);
    fChain->SetBranchAddress("D_py", &D_py, &b_D_py);
    fChain->SetBranchAddress("D_pt", &D_pt, &b_D_pt);
    fChain->SetBranchAddress("D_pz", &D_pz, &b_D_pz);
    fChain->SetBranchAddress("D_p", &D_p, &b_D_p);
    fChain->SetBranchAddress("D_eta", &D_eta, &b_D_eta);
    fChain->SetBranchAddress("D_FDx", &D_FDx, &b_D_FDx);
    fChain->SetBranchAddress("D_FDy", &D_FDy, &b_D_FDy);
    fChain->SetBranchAddress("D_FDt", &D_FDt, &b_D_FDt);
    fChain->SetBranchAddress("D_FDz", &D_FDz, &b_D_FDz);
    fChain->SetBranchAddress("D_FD", &D_FD, &b_D_FD);
    fChain->SetBranchAddress("D_IP", &D_IP, &b_D_IP);
    fChain->SetBranchAddress("D_time", &D_time, &b_D_time);
    fChain->SetBranchAddress("D_time_true", &D_time_true, &b_D_time_true);
    fChain->SetBranchAddress("D_M", &D_M, &b_D_M);
    Notify();
}

bool analysis::Notify()
{
    // The Notify() function is called when a new file is opened. This
    // can be either for a new TTree in a TChain or when when a new TTree
    // is started when using PROOF. It is normally not necessary to make changes
    // to the generated code, but the routine can be extended by the
    // user if needed. The return value is currently not used.

    return true;
}

void analysis::Show(Long64_t entry)
{
    // Print contents of entry.
    // If entry is not specified, print current entry
    if(!fChain)
        return;
    fChain->Show(entry);
}
// Int_t analysis::Cut(Long64_t entry)
//{
//// This function may be called from Loop.
//// returns  1 if entry is accepted.
//// returns -1 otherwise.
//   return 1;
//}

TH1D *analysis::GetBlindedClone(TH1D *h, Double_t blindMin, Double_t blindMax)
{
    if(!h)
        return nullptr;

    // Crea un clone unico per evitare di sovrascrivere o alterare l'istogramma originale
    TString cloneName = Form("%s_blinded", h->GetName());
    TH1D *h_blind = (TH1D *)h->Clone(cloneName);

    Int_t nBins = h_blind->GetNbinsX();
    for(Int_t i = 1; i <= nBins; ++i)
    {
        Double_t binCenter = h_blind->GetBinCenter(i);
        if(binCenter >= blindMin && binCenter <= blindMax)
        {
            h_blind->SetBinContent(i, 0.0);
            h_blind->SetBinError(i, 0.0);
        }
    }
    return h_blind;
}

void analysis::DrawBlindedFunction(TF1 *f, Double_t blindMin, Double_t blindMax, Option_t *option)
{
    if(!f)
        return;

    // Recupera l'intervallo di definizione originale della funzione
    Double_t xMinOrig, xMaxOrig;
    f->GetRange(xMinOrig, xMaxOrig);

    // Disegna il primo segmento (sinistro)
    f->SetRange(xMinOrig, blindMin);
    f->DrawCopy(option);

    // Disegna il secondo segmento (destro)
    f->SetRange(blindMax, xMaxOrig);
    f->DrawCopy(option);

    // Ripristina l'intervallo originale sulla TF1
    f->SetRange(xMinOrig, xMaxOrig);
}

/*
void Blinding_test()
{
    // 1. Creazione di un istogramma con distribuzione gaussiana (10.000 eventi)
    TH1D *h_test = new TH1D("h_test", "Test Blinding Visivo;X;Entries", 100, -5.0, 5.0);

    TRandom3 rand(42); // Seed fisso per riproducibilità
    for(Int_t i = 0; i < 10000; ++i)
    {
        h_test->Fill(rand.Gaus(0.0, 1.0)); // Genera gaussiana con media 0 e sigma 1
    }

    // 2. Definizione ed esecuzione del Fit (In modalità silenziosa, senza disegnare)
    TF1 *f_gaus = new TF1("f_gaus", "gaus", -5.0, 5.0);

    std::cout << "[INFO] Esecuzione del fit sull'istogramma completo..." << std::endl;
    h_test->Fit(f_gaus, "Q N R"); // Q = Quiet, N = No draw, R = Use range

    // 3. Definizione della regione di blinding manuale
    Double_t blindMin = -0.8;
    Double_t blindMax = 0.8;

    // 4. Fase di plotting: Creazione del Canvas
    TCanvas *c_test = new TCanvas("c_test", "Test Blinding Canvas", 800, 600);
    c_test->cd();

    // Utilizzo dell'Helper 1 per ottenere l'istogramma modificato e disegnarlo
    TH1D *h_visual = GetBlindedClone(h_test, blindMin, blindMax);
    h_visual->SetMarkerStyle(20);
    h_visual->SetMarkerSize(0.8);
    h_visual->SetLineColor(kBlack);
    h_visual->Draw("E");

    // Utilizzo dell'Helper 2 per disegnare la funzione di fit con il vuoto al centro
    f_gaus->SetLineColor(kRed);
    f_gaus->SetLineWidth(3);
    DrawBlindedFunction(f_gaus, blindMin, blindMax, "SAME");

    std::cout << "[SUCCESS] Test completato. La regione [" << blindMin << ", " << blindMax
              << "] e' stata oscurata sul plot." << std::endl;
}
 */

#endif // #ifdef analysis_cxx
