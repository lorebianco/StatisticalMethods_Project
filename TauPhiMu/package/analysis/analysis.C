#define analysis_cxx
#include <cmath>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>

#include <Math/Factory.h>
#include <Math/Functor.h>
#include <Math/Minimizer.h>
#include <Math/SpecFuncMathCore.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>
#include <TH2.h>
#include <TMath.h>
#include <TMatrixDSym.h>
#include <TRandom3.h>
#include <TStyle.h>

#include "analysis.h"
#include "lbrootstyle.hh"

using namespace std;
using namespace ROOT::Math;
using namespace TMath;
using namespace lbStyle;

constexpr Bool_t savePlots = false;

// pdf = f_s p(s) + (1-f_s) * (f_1 * p_1 + f_2 * p_2 + f_3 * p_3 + (1 - f_1 - f_2 - f_3) * p_4)

inline Double_t f_2G_Frac(Double_t *x, Double_t *par)
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

class ArgusPDF
{
  private:
    double m_L_range;
    double m_H_range;
    double m_binWidth;

  public:
    ArgusPDF(double L_range, double H_range, double binWidth)
        : m_L_range(L_range)
        , m_H_range(H_range)
        , m_binWidth(binWidth)
    {
    }

    double operator()(const double *x, const double *par) const
    {
        const double m = x[0];
        const double Yield = par[0];
        const double m0 = par[1];
        const double c = par[2];
        const double p = par[3];

        if(m >= m0)
        {
            return 0.0;
        }

        const double xL = 1.0 - (m_L_range / m0) * (m_L_range / m0);
        const double xH = 1.0 - (m_H_range / m0) * (m_H_range / m0);

        const double gammaA = ROOT::Math::tgamma(1.0 + p);
        const double dL = gammaA * ROOT::Math::inc_gamma_c(1.0 + p, -c * xL);
        const double dH = gammaA * ROOT::Math::inc_gamma_c(1.0 + p, -c * xH);
        const double norm = (m0 * m0) / (2.0 * c * pow(-c, p)) * (dL - dH);

        const double u = 1.0 - (m / m0) * (m / m0);
        const double pdf = m * pow(u, p) * exp(c * u) / norm;

        // Scalamento per l'area dell'istogramma
        return Yield * m_binWidth * pdf;
    }
};

class Pol1PDF
{
  private:
    double m_L_range;
    double m_H_range;
    double m_binWidth;

  public:
    Pol1PDF(double L_range, double H_range, double binWidth)
        : m_L_range(L_range)
        , m_H_range(H_range)
        , m_binWidth(binWidth)
    {
    }

    double operator()(const double *x, const double *par) const
    {
        const double xx = x[0];
        const double yield = par[0];
        const double slope = par[1];

        const double range = m_H_range - m_L_range;
        if(range <= 0.0)
            return 0.0;

        const double xMid = 0.5 * (m_L_range + m_H_range);

        // PDF lineare normalizzata a 1 nell'intervallo [m_L_range, m_H_range]
        const double pdf = (1.0 / range) * (1.0 + slope * (xx - xMid));

        // Riscalamento per l'area totale dell'istogramma
        return yield * m_binWidth * pdf;
    }
};

class FullUnbinnedPDF
{
  private:
    double m_xMin, m_xMax;
    // Parametri fissati dai fit ausiliari
    double m_sig_mean, m_sig_sigma1, m_sig_sigma2, m_sig_frac1;
    double m_p1_mean, m_p1_sigma1, m_p1_sigma2, m_p1_frac1;
    double m_p2_mean, m_p2_sigma1, m_p2_sigma2, m_p2_frac1;
    double m_argus_m0, m_argus_c, m_argus_p;
    double m_pol1_slope;

  public:
    FullUnbinnedPDF(double xMin, double xMax, const AuxFitResult &sigRes, const AuxFitResult &p1Res,
        const AuxFitResult &p2Res, const AuxFitResult &argusRes, const AuxFitResult &pol1Res)
        : m_xMin(xMin)
        , m_xMax(xMax)
    {
        m_sig_mean = sigRes.params[1];
        m_sig_sigma1 = sigRes.params[2];
        m_sig_sigma2 = sigRes.params[3];
        m_sig_frac1 = sigRes.params[4];

        m_p1_mean = p1Res.params[1];
        m_p1_sigma1 = p1Res.params[2];
        m_p1_sigma2 = p1Res.params[3];
        m_p1_frac1 = p1Res.params[4];

        m_p2_mean = p2Res.params[1];
        m_p2_sigma1 = p2Res.params[2];
        m_p2_sigma2 = p2Res.params[3];
        m_p2_frac1 = p2Res.params[4];

        m_argus_m0 = argusRes.params[1];
        m_argus_c = argusRes.params[2];
        m_argus_p = argusRes.params[3];

        m_pol1_slope = pol1Res.params[1];
    }

    double Eval2G(double x, double mean, double s1, double s2, double f1) const
    {
        double dx1 = (x - mean) / s1;
        double dx2 = (x - mean) / s2;
        double g1 = std::exp(-0.5 * dx1 * dx1) / (s1 * std::sqrt(2.0 * M_PI));
        double g2 = std::exp(-0.5 * dx2 * dx2) / (s2 * std::sqrt(2.0 * M_PI));
        return f1 * g1 + (1.0 - f1) * g2;
    }

    double EvalArgus(double x, double m0, double c, double p) const
    {
        if(x >= m0)
            return 0.0;
        double xL = 1.0 - (m_xMin / m0) * (m_xMin / m0);
        double xH = 1.0 - (m_xMax / m0) * (m_xMax / m0);
        double gammaA = ROOT::Math::tgamma(1.0 + p);
        double dL = gammaA * ROOT::Math::inc_gamma_c(1.0 + p, -c * xL);
        double dH = gammaA * ROOT::Math::inc_gamma_c(1.0 + p, -c * xH);
        double norm = (m0 * m0) / (2.0 * c * std::pow(-c, p)) * (dL - dH);
        double u = 1.0 - (x / m0) * (x / m0);
        return x * std::pow(u, p) * std::exp(c * u) / norm;
    }

    double EvalPol1(double x, double slope) const
    {
        double range = m_xMax - m_xMin;
        if(range <= 0.0)
            return 0.0;
        double xMid = 0.5 * (m_xMin + m_xMax);
        return (1.0 / range) * (1.0 + slope * (x - xMid));
    }

    double operator()(const double *x, const double *par) const
    {
        const double xx = x[0];
        const double f_s = par[0];
        const double f_1 = par[1];
        const double f_2 = par[2];
        const double f_3 = par[3];

        double p_sig = Eval2G(xx, m_sig_mean, m_sig_sigma1, m_sig_sigma2, m_sig_frac1);
        double p_1 = Eval2G(xx, m_p1_mean, m_p1_sigma1, m_p1_sigma2, m_p1_frac1);
        double p_2 = Eval2G(xx, m_p2_mean, m_p2_sigma1, m_p2_sigma2, m_p2_frac1);
        double p_3 = EvalArgus(xx, m_argus_m0, m_argus_c, m_argus_p);
        double p_4 = EvalPol1(xx, m_pol1_slope);

        double p_bkg = f_1 * p_1 + f_2 * p_2 + f_3 * p_3 + (1.0 - f_1 - f_2 - f_3) * p_4;
        return f_s * p_sig + (1.0 - f_s) * p_bkg;
    }
};

// --- NLL for unbinned fit ---
std::vector<double> g_data_events;
FullUnbinnedPDF *g_pdf_unbinned = nullptr;

double UnbinnedNLL(const double *par)
{
    double nll = 0.0;
    for(int i = 0; i < (int)g_data_events.size(); i++)
    {
        double x = g_data_events[i];
        double val = (*g_pdf_unbinned)(&x, par);
        if(!std::isfinite(val) || val <= 0.0)
            val = 1e-10;
        nll -= std::log(val);
    }
    return nll;
}

// ========================================================================

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
        // cout << M0_MKK << endl;
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
                // cout << D_M << endl;
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

AuxFitResult analysis::FitTemplateMass(Int_t mcID)
{
    AuxFitResult res;

    if(mcID != 34 && mcID != 41 && mcID != 42 && mcID != 44)
    {
        cout << "[ERROR] Invalid mcID = " << mcID << ". Allowed IDs are: 34, 41, 42, 44." << endl;
        return res;
    }

    if(!fChain)
    {
        cout << "[ERROR] fChain is null!" << endl;
        return res;
    }

    Long64_t nentries = fChain->GetEntries();
    if(nentries <= 0)
    {
        cout << "[WARNING] No entries found in fChain!" << endl;
        return res;
    }

    // 1. ALLOCHIAMO SUBITO IL CANVAS CORRETTO (Evita i crash grafici interni di ROOT durante GetEntry)
    SetLBStyle();
    
    // --- MODIFICA 1: Sovrascriviamo le impostazioni dello stile globale per abilitare la Box ---
    gStyle->SetOptFit(1111); // Forza la comparsa di: Chi2/ndf, Parametri, Errori, Status

    TCanvas *cMass = new TCanvas(Form("cMass_%d", mcID), Form("Mass Fit MC %d", mcID), 800, 600);
    cMass->cd();

    std::vector<Double_t> mass_vals;
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
        return res;
    }

    // ... [Il codice intermedio di calcolo di meanInit, rmsInit, xMin, xMax rimane identico] ...
    Double_t sum = 0.0;
    for(double m : mass_vals) sum += m;
    Double_t meanInit = sum / mass_vals.size();

    Double_t sum_sq = 0.0;
    for(double m : mass_vals) sum_sq += (m - meanInit) * (m - meanInit);
    Double_t rmsInit = std::sqrt(sum_sq / mass_vals.size());

    const Double_t binWidthTarget = 5e-4;
    Double_t halfRange = 4.0 * rmsInit;

    Double_t xMin = std::floor((meanInit - halfRange) / binWidthTarget) * binWidthTarget;
    Double_t xMax = std::ceil((meanInit + halfRange) / binWidthTarget) * binWidthTarget;

    Int_t nBins = std::round((xMax - xMin) / binWidthTarget);
    if(nBins <= 0) nBins = 1;

    auto h_mass = new TH1D(Form("h_mass_MC_%d", mcID), "mass", nBins, xMin, xMax);
    h_mass->SetDirectory(nullptr);

    for(double m : mass_vals) h_mass->Fill(m);

    Double_t binWidth = h_mass->GetBinWidth(1);
    TFitResultPtr r;

    cMass->cd();
    AddBinSizeOnYTitle(h_mass, "GeV/#it{c}^{2}");
    
    TString particleName = "";
    if (mcID == 34)      particleName = "D^{+} #rightarrow #phi#pi^{+}";
    else if (mcID == 41) particleName = "D_{s}^{+} #rightarrow #phi#pi^{+}";
    else if (mcID == 42) particleName = "D_{s}^{+} #rightarrow #phi#mu^{+}#nu_{#mu}";
    else if (mcID == 44) particleName = "D_{s}^{+} #rightarrow #tau^{+}#nu_{#tau}";
    
    h_mass->GetXaxis()->SetTitle(Form("M(%s) [GeV/#it{c}^{2}]", particleName.Data()));
    
    // --- MODIFICA 2: Diciamo esplicitamente all'istogramma che VOGLIAMO le statistiche ---
    h_mass->SetStats(kTRUE);

    if(mcID != 42)
    {
        TF1 *fFitMass = new TF1(Form("fFitMass_%d", mcID), f_2G_Frac, xMin, xMax, 5);
        Double_t yieldInit = h_mass->GetEntries() * binWidth;

        fFitMass->SetParameters(yieldInit, meanInit, rmsInit * 0.5, rmsInit * 1.5, 0.7);
        fFitMass->SetParNames("Yield", "Shared_Mean", "#sigma_{1}", "#sigma_{2}", "Frac_1");

        fFitMass->SetParLimits(1, xMin, xMax);
        fFitMass->SetParLimits(2, 0.0, rmsInit * 1.5);
        fFitMass->SetParLimits(3, 0.0, rmsInit * 3.0);
        fFitMass->SetParLimits(4, 0.0, 1.0);

        cout << "\n--- Fitting Mass for MC (ID: " << mcID << ") with Double Gaussian ---" << endl;
        
        // 1. Disegna PRIMA per allocare la grafica
        h_mass->Draw("E");
        
        // 2. Fai il Fit
        r = h_mass->Fit(fFitMass, "L I R S Q");
        
        // 3. Genera la stat box in memoria
        cMass->Modified();
        cMass->Update();
        
        // 4. Trova la box e spostala
        TPaveStats *st = (TPaveStats*)h_mass->FindObject("stats");
        if(st) {
            st->SetOptFit(1111);
            
            st->SetX1NDC(0.653); // Bordo sinistro
            st->SetY1NDC(0.652); // Bordo inferiore
            st->SetX2NDC(0.962); // Bordo destro
            st->SetY2NDC(0.972); // Bordo superiore
        } else {
            cout << "[WARNING] Impossibile trovare la stat box!" << endl;
        }
        
        // 5. Disegna solo le funzioni. L'istogramma è GIÀ disegnato.
        fFitMass->SetLineWidth(3);
        fFitMass->SetLineColor(kRed);
        fFitMass->Draw("SAME");

        TF1 *g1 = new TF1(Form("g1_%d", mcID), "gausn", xMin, xMax);
        g1->SetParameters(fFitMass->GetParameter(0) * fFitMass->GetParameter(4),
            fFitMass->GetParameter(1), fFitMass->GetParameter(2));
        g1->SetLineColor(kGreen + 2);
        g1->SetLineStyle(2);
        g1->SetLineWidth(2);
        g1->Draw("SAME");

        TF1 *g2 = new TF1(Form("g2_%d", mcID), "gausn", xMin, xMax);
        g2->SetParameters(fFitMass->GetParameter(0) * (1.0 - fFitMass->GetParameter(4)),
            fFitMass->GetParameter(1), fFitMass->GetParameter(3));
        g2->SetLineColor(kBlue);
        g2->SetLineStyle(2);
        g2->SetLineWidth(2);
        g2->Draw("SAME");
    }
    else
    {
        Double_t L_taumass_range = 1.6;
        Double_t H_taumass_range = 2.0;

        ArgusPDF argus_func(L_taumass_range, H_taumass_range, binWidth);
        TF1 *fFitMass = new TF1(Form("fFitMass_%d", mcID), argus_func, L_taumass_range, H_taumass_range, 4);

        fFitMass->SetParameters(95200, 1.978, -3.84, 1.417);
        fFitMass->SetParNames("Yield", "m_{0}", "c", "p");

        cout << "\n--- Fitting Mass for MC (ID: " << mcID << ") with ARGUS ---" << endl;
        
        // 1. Disegna PRIMA per allocare la struttura grafica sul Canvas
        h_mass->Draw("E");
        
        // 2. Fai il Fit
        r = h_mass->Fit(fFitMass, "L I R S Q");

        // 3. AGGIORNA IL CANVAS (Questo passaggio mancava e generava un puntatore nullo su 'stats')
        cMass->Modified();
        cMass->Update();

        // 4. Ora che è in memoria, estrai e sposta il pannello
        TPaveStats *st = (TPaveStats*)h_mass->FindObject("stats");
        if(st) {
            st->SetOptFit(1111);
            
            st->SetX1NDC(0.199); // Bordo sinistro
            st->SetY1NDC(0.195); // Bordo inferiore
            st->SetX2NDC(0.560); // Bordo destro
            st->SetY2NDC(0.475); // Bordo superiore
        } else {
            cout << "[WARNING] Impossibile trovare la stat box per ARGUS!" << endl;
        }
        
        // 5. Disegna la linea di Fit sopra l'istogramma esistente
        fFitMass->SetLineWidth(3);
        fFitMass->SetLineColor(kRed);
        fFitMass->Draw("SAME");
    }
    
    // Un ultimo refresh globale prima del salvataggio definitivo
    cMass->Modified();
    cMass->Update();

    if (savePlots) {
        cMass->Print(Form("./_fig/FitMass_MC_%d.pdf", mcID));
        cMass->Print(Form("./_fig/FitMass_MC_%d.root", mcID));
    }
    
    if(r.Get())
    {
        res.isValid = r->IsValid();
        Int_t nPar = r->NPar();
        res.params.resize(nPar);
        res.errors.resize(nPar);
        for(Int_t i = 0; i < nPar; ++i)
        {
            res.params[i] = r->Parameter(i);
            res.errors[i] = r->ParError(i);
        }
        res.covMatrix.resize(nPar, std::vector<Double_t>(nPar, 0.0));
        TMatrixDSym cov = r->GetCovarianceMatrix();
        for(Int_t i = 0; i < nPar; ++i)
        {
            for(Int_t j = 0; j < nPar; ++j)
            {
                res.covMatrix[i][j] = cov(i, j);
            }
        }
    }
    return res;
}

AuxFitResult analysis::FitCombinatorialBkg()
{
    AuxFitResult res;
    SetLBStyle();
    
    // --- MODIFICA 1: Abilitiamo lo stile globale per il pannello delle statistiche ---
    gStyle->SetOptFit(1111);

    LoadDataset(0);

    Double_t xMin = 2.0;
    Double_t xMax = 2.09;

    auto h_mass = new TH1D("h_mass_bkg", ";Invariant Mass [GeV/#it{c}^{2}];Entries", 100, xMin, xMax);
    h_mass->SetDirectory(nullptr);
    AddBinSizeOnYTitle(h_mass, "GeV/#it{c}^{2}");

    for(Long64_t jentry = 0; jentry < fChain->GetEntriesFast(); jentry++)
    {
        if(LoadTree(jentry) < 0)
            break;
        fChain->GetEntry(jentry);

        if(id != 0)
        {
            cerr << "[ERROR] Loaded data is MC, not real data!" << endl;
            return res;
        }

        h_mass->Fill(D_M);
    }

    Double_t binWidth = h_mass->GetBinWidth(1);

    Pol1PDF pol1_func(xMin, xMax, binWidth);
    TF1 *f_pol1 = new TF1("f_pol1", pol1_func, xMin, xMax, 2);

    Double_t yieldInit = h_mass->GetEntries();
    f_pol1->SetParameters(yieldInit, 0.0);
    f_pol1->SetParNames("Yield", "Slope");

    cout << "\n--- Fitting Combinatorial Background with Pol1PDF Functor ---" << endl;
    
    // --- MODIFICA 2: Rimosse opzione "N" per permettere a ROOT di salvare i dati del fit nell'istogramma ---
    TFitResultPtr r = h_mass->Fit(f_pol1, "L I R S Q");

    Double_t blindMin = 1.777 - 3 * 0.0058;
    Double_t blindMax = 1.777 + 3 * 0.0058;

    TCanvas *c_bkg = new TCanvas("c_bkg", "Combinatorial Background Fit", 800, 600);
    c_bkg->cd();

    // --- MODIFICA 3: Riorganizzazione della logica di disegno e gestione unificata della Stat Box ---
    TH1D *h_to_draw = h_mass;
    bool isBlinded = (blindMin > xMin && blindMax < xMax);

    if(isBlinded)
    {
        // Se siamo nel range di blinding, creiamo il clone oscurato
        h_to_draw = analysis::GetBlindedClone(h_mass, blindMin, blindMax);
    }

    // Diciamo all'istogramma (normale o blindato) di mostrare le statistiche e lo disegnamo
    h_to_draw->SetStats(kTRUE);
    h_to_draw->SetMinimum(0.0);
    h_to_draw->SetMaximum(80.0);
    h_to_draw->Draw("E");

    // Forza ROOT a generare la stat box in memoria prima di cercarla
    c_bkg->Modified();
    c_bkg->Update();

    // Cerchiamo la box direttamente dall'istogramma appena disegnato
    TPaveStats *st = (TPaveStats*)h_to_draw->FindObject("stats");
    if(st) {
        st->SetOptFit(1111);
        st->SetX1NDC(0.653); // Bordo sinistro
        st->SetY1NDC(0.734); // Bordo inferiore
        st->SetX2NDC(0.961); // Bordo destro
        st->SetY2NDC(0.972); // Bordo superiore
        // Sotto inteso: NESSUN st->Draw() per evitare il bug del PDF vuoto!
    } else {
        cout << "[WARNING] Impossibile trovare la stat box del background!" << endl;
    }

    // Ora disegnamo la funzione sopra l'istogramma
    f_pol1->SetLineColor(kRed);
    f_pol1->SetLineWidth(3);
    if(isBlinded)
    {
        analysis::DrawBlindedFunction(f_pol1, blindMin, blindMax, "SAME");
    }
    else
    {
        f_pol1->Draw("SAME");
    }
    
    // Refresh finale della grafica per incorporare il posizionamento della box
    c_bkg->Modified();
    c_bkg->Update();

    if (savePlots) {
        c_bkg->Print("./_fig/FitCombinatorialBkg.pdf");
        c_bkg->Print("./_fig/FitCombinatorialBkg.root");
    }

    if(r.Get())
    {
        res.isValid = r->IsValid();
        Int_t nPar = r->NPar();
        res.params.resize(nPar);
        res.errors.resize(nPar);
        for(Int_t i = 0; i < nPar; ++i)
        {
            res.params[i] = r->Parameter(i);
            res.errors[i] = r->ParError(i);
        }
        res.covMatrix.resize(nPar, std::vector<Double_t>(nPar, 0.0));
        TMatrixDSym cov = r->GetCovarianceMatrix();
        for(Int_t i = 0; i < nPar; ++i)
        {
            for(Int_t j = 0; j < nPar; ++j)
            {
                res.covMatrix[i][j] = cov(i, j);
            }
        }
    }
    return res;
}

void analysis::DoFullBlindedUnbinnedFit()
{
    SetLBStyle();

    // 1. Esecuzione dei fit ausiliari (MC e fondo combinatorio)
    cout << "\n=== [STEP 1/3] Running MC Auxiliary Fits ===" << endl;
    LoadDataset(1); // Carica Monte Carlo
    AuxFitResult res_sig = FitTemplateMass(44);
    AuxFitResult res_p1 = FitTemplateMass(34);
    AuxFitResult res_p2 = FitTemplateMass(41);
    AuxFitResult res_arg = FitTemplateMass(42);

    cout << "\n=== [STEP 2/3] Running Combinatorial Background Fit ===" << endl;
    AuxFitResult res_pol = FitCombinatorialBkg();

    if(!res_sig.isValid || !res_p1.isValid || !res_p2.isValid || !res_arg.isValid
        || !res_pol.isValid)
    {
        cerr << "[ERROR] Auxiliary fits failed! Aborting unbinned fit." << endl;
        return;
    }

    // 2. Caricamento dei dati reali per il Fit Unbinned
    cout << "\n=== [STEP 3/3] Preparing Real Data for Unbinned Fit ===" << endl;
    LoadDataset(0); // Carica Dati Reali

    Double_t xMin = 1.65;
    Double_t xMax = 2.09;

    g_data_events.clear();
    Long64_t nentries = fChain->GetEntries();
    g_data_events.reserve(static_cast<size_t>(nentries));

    // Istogramma binned per la sola visualizzazione finale
    auto h_mass = new TH1D("h_mass_unbinned", ";Invariant Mass [GeV/#it{c}^{2}];Entries", 100, xMin, xMax);
    h_mass->SetDirectory(nullptr);
    
    for(Long64_t jentry = 0; jentry < fChain->GetEntriesFast(); jentry++)
    {
        if(LoadTree(jentry) < 0)
            break;
        fChain->GetEntry(jentry);
        if(id != 0)
            continue;

        if(D_M >= xMin && D_M <= xMax)
        {
            g_data_events.push_back(D_M);
            h_mass->Fill(D_M); // Riempiamo anche l'istogramma di plot
        }
    }

    Double_t binWidth = h_mass->GetBinWidth(1);
    Double_t nEntries = h_mass->GetEntries();

    // 3. Setup del Functor della PDF totale e di MINUIT
    g_pdf_unbinned = new FullUnbinnedPDF(xMin, xMax, res_sig, res_p1, res_p2, res_arg, res_pol);

    ROOT::Math::Minimizer *minimizer = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
    minimizer->SetMaxFunctionCalls(100000);
    minimizer->SetTolerance(0.01);

    // IMPORTANTE: Silenzia completamente l'output interno di MINUIT durante le iterazioni
    minimizer->SetPrintLevel(0);

    ROOT::Math::Functor fNLL(&UnbinnedNLL, 4);
    minimizer->SetFunction(fNLL);

    // Definizione delle frazioni e limiti [0, 1]
    minimizer->SetLimitedVariable(0, "f_s", 0.02, 0.005, 0.0, 1.0);
    minimizer->SetLimitedVariable(1, "f_1", 0.15, 0.010, 0.0, 1.0);
    minimizer->SetLimitedVariable(2, "f_2", 0.15, 0.010, 0.0, 1.0);
    minimizer->SetLimitedVariable(3, "f_3", 0.30, 0.020, 0.0, 1.0);

    cout << "\n--- Minimizing Unbinned Likelihood with MINUIT (Silenced) ---" << endl;
    minimizer->Minimize();

    // 4. STAMPA USER PERSONALIZZATA (BLINDED)
    cout << "\n=======================================================" << endl;
    cout << "   UNBINNED MAXIMUM LIKELIHOOD FIT RESULTS (BLINDED)  " << endl;
    cout << "=======================================================" << endl;
    cout << "  Fit Status: " << (minimizer->Status() == 0 ? "CONVERGED" : "FAILED") << endl;
    cout << "  Edm:        " << minimizer->Edm() << endl;
    cout << "-------------------------------------------------------" << endl;

    const double *xs = minimizer->X();
    const double *errs = minimizer->Errors();

    for(unsigned int i = 0; i < minimizer->NDim(); ++i)
    {
        std::string varName = minimizer->VariableName(i);

        // Se il parametro è la frazione di segnale "f_s", nascondiamo i valori
        if(varName == "f_s" || i == 0)
        {
            cout << "  " << varName << " \t=  [BLINDED]  +/-  [BLINDED]" << endl;
        }
        else
        {
            cout << "  " << varName << " \t=  " << xs[i] << "  +/-  " << errs[i] << endl;
        }
    }
    cout << "=======================================================" << endl;

    // 5. Fase di Plotting Visivo con Blinding
    Double_t blindMin = 1.777 - 3 * 0.0058;
    Double_t blindMax = 1.777 + 3 * 0.0058;

    TCanvas *c_unbinned = new TCanvas("c_unbinned", "Unbinned Maximum Likelihood Fit", 800, 600);
    c_unbinned->cd();
    AddBinSizeOnYTitle(h_mass, "GeV/#it{c}^{2}");

    // Istogramma accecato
    TH1D *h_mass_blind = analysis::GetBlindedClone(h_mass, blindMin, blindMax);
    h_mass_blind->Draw("E");

    // Creiamo una funzione di disegno per riscalare la PDF unbinned all'altezza dell'istogramma
    FullUnbinnedPDF *pdf_for_draw = g_pdf_unbinned;
    Double_t nEntries_copy = nEntries;
    Double_t binWidth_copy = binWidth;

    auto scale_pdf_lambda = [pdf_for_draw, nEntries_copy, binWidth_copy](double *x, double *par)
    {
        double xx = x[0];
        double pdf_val = (*pdf_for_draw)(&xx, par);
        return nEntries_copy * binWidth_copy * pdf_val;
    };

    TF1 *f_draw = new TF1("f_draw", scale_pdf_lambda, xMin, xMax, 4);
    f_draw->SetParameters(
        minimizer->X()[0], minimizer->X()[1], minimizer->X()[2], minimizer->X()[3]);
    f_draw->SetNpx(10000);
    f_draw->SetLineColor(kRed);
    f_draw->SetLineWidth(3);

    // Disegniamo la curva con il vuoto al centro
    analysis::DrawBlindedFunction(f_draw, blindMin, blindMax, "SAME");

    if (savePlots) {
        c_unbinned->Print("./_fig/FitFullUnbinned.pdf");
        c_unbinned->Print("./_fig/FitFullUnbinned.root");
    }
    
    // Pulizia della memoria - NOTA: g_pdf_unbinned NON deve essere deallocata qui
    // perché la lambda in f_draw la riferisce ancora!
    delete minimizer;
    // delete g_pdf_unbinned; // COMENTATO: evita accessi a memoria non valida
    g_pdf_unbinned = nullptr;
}
