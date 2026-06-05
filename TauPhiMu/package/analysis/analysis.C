#define analysis_cxx
#include <iostream>
#include <stdio.h>
#include <stdlib.h>

#include <TCanvas.h>
#include <TF1.h>
#include <TH1D.h>
#include <TH2.h>
#include <TStyle.h>

#include "analysis.h"
#include "lbrootstyle.hh"

using namespace std;
using namespace TMath;
using namespace lbStyle;

constexpr Bool_t savePlots = false;

// pdf = f_s p(s) + (1-f_s) * (f_1 * p_1 + f_2 * p_2 + f_3 * p_3 + (1 - f_1 - f_2 - f_3) * p_4)

Double_t f_2G_Frac(Double_t *x, Double_t *par)
{
    // par[0] = Yield (Eventi * bin width)
    // par[1] = Mu
    // par[2] = Sigma_1
    // par[3] = Sigma_2
    // par[4] = Frac_1 (Frazione della prima gaussiana)

    double dx1 = (x[0] - par[1]) / par[2];
    double dx2 = (x[0] - par[1]) / par[3];

    double g1 = TMath::Exp(-0.5 * dx1 * dx1) / (par[2] * TMath::Sqrt(TMath::TwoPi()));
    double g2 = TMath::Exp(-0.5 * dx2 * dx2) / (par[3] * TMath::Sqrt(TMath::TwoPi()));

    return par[0] * (par[4] * g1 + (1.0 - par[4]) * g2);
}

void analysis::Loop()
{

    TH1D *histo_D_M = new TH1D("histo_D_M", "", 100, 1.640, 2.100);
    TH1D *histo_D_M_S = new TH1D("histo_D_M_S", "", 100, 1.640, 2.100);
    TH1D *histo_D_time = new TH1D("histo_D_time", "", 100, 0, 5e-12);
    TH1D *histo_D_time_S = new TH1D("histo_D_time_S", "", 100, 0, 5e-12);
    TH1D *histo_D_pt = new TH1D("histo_D_pt", "", 100, 0, 10);
    TH1D *histo_D_pt_S = new TH1D("histo_D_pt_S", "", 100, 0, 10);
    TH1D *histo_h3_p_back = new TH1D("histo_h3_p_back", "", 1000, 0, 200);

    if(fChain == 0)
        return;

    Long64_t nentries = fChain->GetEntriesFast();

    Long64_t nbytes = 0, nb = 0;
    for(Long64_t jentry = 0; jentry < nentries; jentry++)
    {
        Long64_t ientry = LoadTree(jentry);
        if(ientry < 0)
            break;
        nb = fChain->GetEntry(jentry);
        nbytes += nb;
        // if (Cut(ientry) < 0) continue;
        // std::cout << M0_MKK << std::endl;
        histo_h3_p_back->Fill(h3_p, weight);

        // from fit to signal MC
        // Mean    =  1.77699   +/-   1.98392e-05
        // Sigma   =   0.00581298   +/-   1.49607e-05  	 (limited)
        Double_t m_tau_mc;
        Double_t m_sigma_tau_mc;
        m_tau_mc = 1.77699;
        m_sigma_tau_mc = 0.00581298;

        if(TMath::Abs(D_M - m_tau_mc) > 3 * m_sigma_tau_mc)
        {

            if(h3_MuonID == 1)
            {
                // if (1){
                histo_D_M->Fill(D_M);
                histo_D_time->Fill(D_time);
                histo_D_pt->Fill(D_pt);
                // std::cout << D_M << std::endl;
            };
            if(h3_MuonID == 1 && D_time > 0.25e-12 && D_pt > 2.5 && D_FDt > 1e-3)
            {
                // if (h3_MuonID==1 && D_time>0.25e-12 && D_pt>2.5 ){
                // if ( D_time>2.5e-12){
                histo_D_M_S->Fill(D_M);
                histo_D_time_S->Fill(D_time);
                histo_D_pt_S->Fill(D_pt);
            };

        }; // Mass blinding

    }; // End Loop Events

    // Create Canvas
    TCanvas *c_histo_D_M = new TCanvas("c_histo_D_M", "canvas histo", 500, 500);
    c_histo_D_M->cd();
    histo_D_M->SetMinimum(0);
    histo_D_M->Draw();
    histo_D_M_S->SetLineColor(2);
    histo_D_M_S->Draw("same");

    TCanvas *c_histo_D_time = new TCanvas("c_histo_D_time", "canvas histo", 500, 500);
    c_histo_D_time->cd();
    histo_D_time->SetMinimum(0);
    histo_D_time->Draw();
    histo_D_time_S->SetLineColor(2);
    histo_D_time_S->Draw("same");

    TCanvas *c_histo_D_pt = new TCanvas("c_histo_D_pt", "canvas histo", 500, 500);
    c_histo_D_pt->cd();
    histo_D_pt->SetMinimum(0);
    histo_D_pt->Draw();
    histo_D_pt_S->SetLineColor(2);
    histo_D_pt_S->Draw("same");

    // Save canvas in .pdf and .epr format
    c_histo_D_M->Print("./_fig/c_histo.pdf");
    c_histo_D_M->Print("./_fig/c_histo.eps");

    // Create a new file to store histograms
    TFile *histo_file = new TFile("./_root/histo_file.root", "RECREATE", "put a title");
    // TFile *histo_file = new TFile("./_root/histo_p_back_DPLUS_PhiPi.root","RECREATE","put a
    // title");
    histo_file->cd();
    histo_h3_p_back->Write();
    // histo_file->Write();
    histo_D_M->Write();
    histo_D_M_S->Write();
    histo_D_time->Write();
    histo_D_time_S->Write();
    histo_D_pt->Write();
    histo_D_pt_S->Write();
    histo_file->Close();
}

void analysis::FitTemplateMass(Int_t mcID)
{
    SetLBStyle();

    // Controllo di sicurezza sulla validità della TChain
    if(!fChain)
    {
        cout << "[ERROR] fChain is null!" << endl;
        return;
    }

    Long64_t nentries = fChain->GetEntries(); // Più sicuro rispetto a GetEntriesFast()
    if(nentries <= 0)
    {
        cout << "[WARNING] No entries found in fChain!" << endl;
        return;
    }

    std::vector<Double_t> mass_vals;

    // Alloca memoria solo se nentries ha un valore positivo coerente
    mass_vals.reserve(static_cast<size_t>(nentries));

    for(Long64_t jentry = 0; jentry < nentries; jentry++)
    {
        Long64_t ientry = LoadTree(jentry);
        if(ientry < 0)
            break;
        fChain->GetEntry(jentry);

        if(id != mcID)
            continue;
        mass_vals.push_back(D_M);
    }

    if(mass_vals.empty())
    {
        cout << "[WARNING] No events found for mcID = " << mcID << endl;
        return;
    }

    // Calcolo della media e del RMS dei dati reali
    Double_t sum = 0.0;
    for(double m : mass_vals)
        sum += m;
    Double_t meanInit = sum / mass_vals.size();

    Double_t sum_sq = 0.0;
    for(double m : mass_vals)
        sum_sq += (m - meanInit) * (m - meanInit);
    Double_t rmsInit = std::sqrt(sum_sq / mass_vals.size());

    // 1. Definizione della larghezza fissa del bin
    const Double_t binWidthTarget = 5e-4;

    // Intervallo desiderato (circa +/- 4 volte l'RMS)
    Double_t halfRange = 4.0 * rmsInit;

    // Allineamento degli estremi a multipli esatti di 0.005 per eccesso/difetto
    Double_t xMin = std::floor((meanInit - halfRange) / binWidthTarget) * binWidthTarget;
    Double_t xMax = std::ceil((meanInit + halfRange) / binWidthTarget) * binWidthTarget;

    // Calcolo esatto del numero di bin necessari
    Int_t nBins = std::round((xMax - xMin) / binWidthTarget);

    // Sicurezza: se per qualche motivo nBins fosse zero, impostiamo almeno un bin
    if(nBins <= 0)
        nBins = 1;

    // 2. Creazione dell'istogramma con la larghezza del bin fissa a 0.005
    auto h_mass = new TH1D("h_mass", "mass", nBins, xMin, xMax);
    for(double m : mass_vals)
    {
        h_mass->Fill(m);
    }

    // 3. Definizione della funzione di Fit
    TF1 *fFitMass = new TF1("fFitMass", f_2G_Frac, xMin, xMax, 5);

    Double_t binWidth = h_mass->GetBinWidth(1);
    Double_t yieldInit = h_mass->GetEntries() * binWidth; // Area totale del picco

    // Mapping dei parametri di f_2G_Frac:
    // par[0] = Yield (Area totale)
    // par[1] = Mu (Media condivisa)
    // par[2] = Sigma_1
    // par[3] = Sigma_2
    // par[4] = Frac_1 (Frazione della prima gaussiana)
    fFitMass->SetParameters(yieldInit, meanInit, rmsInit * 0.5, rmsInit * 1.5, 0.7);
    fFitMass->SetParNames("Yield", "Shared_Mean", "#sigma_{1}", "#sigma_{2}", "Frac_1");

    // Limiti di sicurezza per evitare che i parametri divergano durante il fit
    fFitMass->SetParLimits(1, xMin, xMax);
    fFitMass->SetParLimits(2, 0.0, rmsInit * 1.5);
    fFitMass->SetParLimits(3, 0.0, rmsInit * 3.0);
    fFitMass->SetParLimits(4, 0.0, 1.0);

    // 4. Creazione del Canvas e disegno
    TCanvas *cMass = new TCanvas("cMass", "Mass Fit MC", 800, 600);
    cMass->cd();
    AddBinSizeOnYTitle(h_mass, "GeV/#it{c}^{2}");

    // 5. Esecuzione del Fit
    cout << "\n--- Fitting Mass for MC (ID: " << mcID << ") ---" << endl;
    h_mass->Fit(fFitMass, "L I R");

    h_mass->Draw("E");

    fFitMass->SetLineWidth(3);
    fFitMass->SetLineColor(kRed);
    fFitMass->Draw("SAME");

    // 6. Disegno separato delle due componenti gaussiane (area-normalizzate "gausn")
    // Componente 1 (Verde)
    TF1 *g1 = new TF1("g1", "gausn", xMin, xMax);
    g1->SetParameters(
        fFitMass->GetParameter(0) * fFitMass->GetParameter(4), // Area = Yield * Frac_1
        fFitMass->GetParameter(1), // Media condivisa
        fFitMass->GetParameter(2) // Sigma 1
    );
    g1->SetLineColor(kGreen + 2);
    g1->SetLineStyle(2);
    g1->SetLineWidth(2);
    g1->Draw("SAME");

    // Componente 2 (Blu)
    TF1 *g2 = new TF1("g2", "gausn", xMin, xMax);
    g2->SetParameters(fFitMass->GetParameter(0)
            * (1.0 - fFitMass->GetParameter(4)), // Area = Yield * (1 - Frac_1)
        fFitMass->GetParameter(1), // Media condivisa
        fFitMass->GetParameter(3) // Sigma 2
    );
    g2->SetLineColor(kBlue);
    g2->SetLineStyle(2);
    g2->SetLineWidth(2);
    g2->Draw("SAME");
}

void analysis::FitCombinatorialBkg() { }
