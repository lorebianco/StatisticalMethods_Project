#include "RtypesCore.h"
#define analysis_cxx
#include <algorithm>
#include <cmath>
#include <fstream>
#include <future>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <thread>
#include <vector>

#include <Math/Factory.h>
#include <Math/Functor.h>
#include <Math/Minimizer.h>
#include <Math/SpecFuncMathCore.h>
#include <TAxis.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TH2.h>
#include <TLegend.h>
#include <TLine.h>
#include <TMath.h>
#include <TMatrixDSym.h>
#include <TPaveStats.h>
#include <TRandom3.h>
#include <TStyle.h>

#include "analysis.h"
#include "lbrootstyle.hh"

using namespace std;
using namespace ROOT::Math;
using namespace TMath;
using namespace lbStyle;

constexpr Bool_t savePlots = false;
constexpr UInt_t THE_SEED = 0;

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

class ExpoPDF
{
  private:
    double m_L_range;
    double m_H_range;
    double m_binWidth;

  public:
    ExpoPDF(double L_range, double H_range, double binWidth)
        : m_L_range(L_range)
        , m_H_range(H_range)
        , m_binWidth(binWidth)
    {
    }

    double operator()(const double *x, const double *par) const
    {
        const double xx = x[0];
        const double yield = par[0];
        const double lambda = par[1]; // Esponente della funzione e^(slope * x)

        // Calcolo analitico dell'integrale per la normalizzazione
        double norm = 0.0;
        if(std::abs(lambda) < 1e-6)
        {
            // Se lo slope è quasi zero, la funzione è piatta (evita divisione per zero)
            norm = m_H_range - m_L_range;
        }
        else
        {
            norm = (std::exp(lambda * m_H_range) - std::exp(lambda * m_L_range)) / lambda;
        }

        if(norm <= 0.0)
            return 0.0;

        const double pdf = std::exp(lambda * xx) / norm;

        // Scalamento per l'area dell'istogramma
        return yield * m_binWidth * pdf;
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
  public:
    double m_xMin, m_xMax;
    double m_sig_mean, m_sig_sigma1, m_sig_sigma2, m_sig_frac1;
    double m_p1_mean, m_p1_sigma1, m_p1_sigma2, m_p1_frac1;
    double m_p2_mean, m_p2_sigma1, m_p2_sigma2, m_p2_frac1;
    double m_argus_m0, m_argus_c, m_argus_p;
    bool m_usePol1Bkg;

    FullUnbinnedPDF(double xMin, double xMax, const AuxFitResult &sigRes, const AuxFitResult &p1Res,
        const AuxFitResult &p2Res, const AuxFitResult &argusRes, bool usePol1Bkg = false)
        : m_xMin(xMin)
        , m_xMax(xMax)
        , m_usePol1Bkg(usePol1Bkg)
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
    }

    // Calcolo esatto dell'integrale della gaussiana tra xMin e xMax
    double NormGaussian(double mean, double sigma) const
    {
        double zMin = (m_xMin - mean) / (sigma * M_SQRT2);
        double zMax = (m_xMax - mean) / (sigma * M_SQRT2);
        return 0.5 * (std::erf(zMax) - std::erf(zMin));
    }

    double Eval2G(double x, double mean, double s1, double s2, double f1) const
    {
        // Ottieni la normalizzazione esatta
        double norm1 = NormGaussian(mean, s1);
        double norm2 = NormGaussian(mean, s2);

        double dx1 = (x - mean) / s1;
        double dx2 = (x - mean) / s2;

        // Dividi per la normalizzazione di range oltre a quella standard
        double g1 = std::exp(-0.5 * dx1 * dx1) / (s1 * std::sqrt(2.0 * M_PI) * norm1);
        double g2 = std::exp(-0.5 * dx2 * dx2) / (s2 * std::sqrt(2.0 * M_PI) * norm2);

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

    double EvalExpo(double x, double slope) const
    {
        double norm = (std::exp(slope * m_xMax) - std::exp(slope * m_xMin)) / slope;
        return std::exp(slope * x) / norm;
    }

    double EvalPol1(double x, double slope) const
    {
        double range = m_xMax - m_xMin;
        double xMid = 0.5 * (m_xMin + m_xMax);
        return (1.0 / range) * (1.0 + slope * (x - xMid));
    }

    // --- OPERATORE CHIAMATO DA MINUIT ---
    double operator()(const double *x, const double *par) const
    {
        const double xx = x[0];
        const double f_s = par[0];
        const double f_1 = par[1];
        const double f_2 = par[2];
        const double f_3 = par[3];
        const double slope = par[4];

        double p_sig = Eval2G(xx, par[5], par[6], par[7], par[8]);
        double p_1 = Eval2G(xx, par[9], par[10], par[11], par[12]);
        double p_2 = Eval2G(xx, par[13], par[14], par[15], par[16]);
        double p_3 = EvalArgus(xx, par[17], par[18], par[19]);
        double p_4 = m_usePol1Bkg ? EvalPol1(xx, slope) : EvalExpo(xx, slope);

        double p_bkg = f_1 * p_1 + f_2 * p_2 + f_3 * p_3 + (1.0 - f_1 - f_2 - f_3) * p_4;

        double pdf_val = f_s * p_sig + (1.0 - f_s) * p_bkg;
        return pdf_val;
    }
};

// --- NLL for unbinned fit ---
thread_local std::vector<double> g_data_events;
thread_local FullUnbinnedPDF *g_pdf_unbinned = nullptr;

inline double Unbinned2NLL(const double *par)
{
    double nll = 0.0;
    for(int i = 0; i < (int)g_data_events.size(); i++)
    {
        double x = g_data_events[i];
        double val = (*g_pdf_unbinned)(&x, par);

        nll -= std::log(val);
    }
    return 2.0 * nll;
}

/*
    // --- Likelihood ratio ordering
    inline Double_t LROrdering(Double_t *x, Double_t *par)
    {
    double xx = x[0];
    double mu = par[0];
    double sigma = par[1];

    if(xx >= 0)
        return exp(-0.5 * (xx - mu) * (xx - mu) / (sigma * sigma));
    else
        return exp(
            (-0.5 * (xx - mu) * (xx - mu) / (sigma * sigma)) + 0.5 * (xx * xx) / (sigma * sigma));
}
*/

// --- Likelihood ratio ordering (Feldman-Cousins con sigma dinamica)
inline Double_t LROrdering(Double_t *x, Double_t *par)
{
    double xx = x[0];
    double mu = par[0];
    double sigma0 = par[1];
    double alpha = par[2];

    // Calcolo la sigma per l'ipotesi mu (numeratore)
    double sigma_mu = sigma0 + alpha * mu;

    // Calcolo il best-fit e la sua rispettiva sigma (denominatore)
    double mu_best = std::max(0.0, xx); // Limite fisico: mu_best >= 0
    double sigma_best = sigma0 + alpha * mu_best;

    // Prefattore matematico derivante dal rapporto 1/sigma_mu diviso 1/sigma_best
    double prefactor = sigma_best / sigma_mu;

    if(xx >= 0)
    {
        // Se xx >= 0, allora mu_best = xx. Il termine esponenziale al denominatore è exp(0) = 1.
        return prefactor * exp(-0.5 * (xx - mu) * (xx - mu) / (sigma_mu * sigma_mu));
    }
    else
    {
        // Se xx < 0, allora mu_best = 0. Il denominatore ha un esponenziale valutato in 0.
        return prefactor
            * exp((-0.5 * (xx - mu) * (xx - mu) / (sigma_mu * sigma_mu))
                + (0.5 * (xx * xx) / (sigma_best * sigma_best)));
    }
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
    c_histo_D_M->SaveAs("./_fig/c_histo.pdf");
    c_histo_D_M->SaveAs("./_fig/c_histo.eps");

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

    // 1. ALLOCHIAMO SUBITO IL CANVAS CORRETTO (Evita i crash grafici interni di ROOT durante
    // GetEntry)
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
    for(double m : mass_vals)
        sum += m;
    Double_t meanInit = sum / mass_vals.size();

    Double_t sum_sq = 0.0;
    for(double m : mass_vals)
        sum_sq += (m - meanInit) * (m - meanInit);
    Double_t rmsInit = std::sqrt(sum_sq / mass_vals.size());

    Double_t binWidthTarget = 5e-4;
    Double_t halfRange = 4.0 * rmsInit;

    Double_t xMin = std::floor((meanInit - halfRange) / binWidthTarget) * binWidthTarget;
    Double_t xMax = std::ceil((meanInit + halfRange) / binWidthTarget) * binWidthTarget;

    // Forza il range fisso tra 1.6 e 2.0 solo per il canale 42 (ARGUS)
    if(mcID == 42)
    {
        xMin = 1.60;
        xMax = 1.98;

        binWidthTarget = 5e-3;
    }

    Int_t nBins = std::round((xMax - xMin) / binWidthTarget);
    if(nBins <= 0)
        nBins = 1;

    auto h_mass = new TH1D(Form("h_mass_MC_%d", mcID), "mass", nBins, xMin, xMax);
    h_mass->SetDirectory(nullptr);

    for(double m : mass_vals)
        h_mass->Fill(m);

    Double_t binWidth = h_mass->GetBinWidth(1);
    TFitResultPtr r;

    cMass->cd();
    AddBinSizeOnYTitle(h_mass, "GeV/#it{c}^{2}");

    TString particleName = "";
    if(mcID == 34)
        particleName = "D^{+}#rightarrow#phi#pi^{+}";
    else if(mcID == 41)
        particleName = "D_{s}^{+}#rightarrow#phi#pi^{+}";
    else if(mcID == 42)
        particleName = "D_{s}^{+}#rightarrow#phi#mu^{+}#nu_{#mu}";
    else if(mcID == 44)
        particleName = "D_{s}^{+}#rightarrow#tau^{+}#nu_{#tau}";

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
        TPaveStats *st = (TPaveStats *)h_mass->FindObject("stats");
        if(st)
        {
            st->SetOptFit(1111);

            st->SetX1NDC(0.653); // Bordo sinistro
            st->SetY1NDC(0.652); // Bordo inferiore
            st->SetX2NDC(0.962); // Bordo destro
            st->SetY2NDC(0.972); // Bordo superiore
        }
        else
        {
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
        TF1 *fFitMass
            = new TF1(Form("fFitMass_%d", mcID), argus_func, L_taumass_range, H_taumass_range, 4);

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
        TPaveStats *st = (TPaveStats *)h_mass->FindObject("stats");
        if(st)
        {
            st->SetOptFit(1111);

            st->SetX1NDC(0.199); // Bordo sinistro
            st->SetY1NDC(0.195); // Bordo inferiore
            st->SetX2NDC(0.560); // Bordo destro
            st->SetY2NDC(0.475); // Bordo superiore
        }
        else
        {
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

    if(savePlots)
    {
        cMass->SaveAs(Form("./_fig/FitMass_MC_%d.pdf", mcID));
        cMass->SaveAs(Form("./_root/FitMass_MC_%d.root", mcID));
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

AuxFitResult analysis::FitCombinatorialBkg(bool usePol1)
{
    AuxFitResult res;
    SetLBStyle();
    gStyle->SetOptFit(1111);

    LoadDataset(0);

    Double_t xMin = 2.0;
    Double_t xMax = 2.09;

    auto h_mass
        = new TH1D("h_mass_bkg", ";Invariant Mass [GeV/#it{c}^{2}];Entries", 100, xMin, xMax);
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

    TF1 *f_bkg = nullptr;
    Double_t yieldInit = h_mass->GetEntries();

    if(usePol1)
    {
        Pol1PDF pol1_func(xMin, xMax, binWidth);
        f_bkg = new TF1("f_bkg", pol1_func, xMin, xMax, 2);
        f_bkg->SetParameters(yieldInit, 0.0);
        f_bkg->SetParNames("Yield", "Slope");
        cout << "\n--- Fitting Combinatorial Background with Pol1PDF Functor ---" << endl;
    }
    else
    {
        ExpoPDF expo_func(xMin, xMax, binWidth);
        f_bkg = new TF1("f_bkg", expo_func, xMin, xMax, 2);
        f_bkg->SetParameters(yieldInit, -2.0);
        f_bkg->SetParNames("Yield", "Lambda");
        cout << "\n--- Fitting Combinatorial Background with ExpoPDF Functor ---" << endl;
    }

    TFitResultPtr r = h_mass->Fit(f_bkg, "L I R S Q");

    Double_t blindMin = 1.777 - 3 * 0.0058;
    Double_t blindMax = 1.777 + 3 * 0.0058;

    TCanvas *c_bkg = new TCanvas("c_bkg", "Combinatorial Background Fit", 800, 600);
    c_bkg->cd();

    TH1D *h_to_draw = h_mass;
    bool isBlinded = (blindMin > xMin && blindMax < xMax);

    if(isBlinded)
    {
        h_to_draw = analysis::GetBlindedClone(h_mass, blindMin, blindMax);
    }

    h_to_draw->SetStats(kTRUE);
    h_to_draw->SetMinimum(0.0);
    h_to_draw->SetMaximum(80.0);
    h_to_draw->Draw("E");

    c_bkg->Modified();
    c_bkg->Update();

    TPaveStats *st = (TPaveStats *)h_to_draw->FindObject("stats");
    if(st)
    {
        st->SetOptFit(1111);
        st->SetX1NDC(0.653);
        st->SetY1NDC(0.734);
        st->SetX2NDC(0.961);
        st->SetY2NDC(0.972);
    }

    f_bkg->SetLineColor(kRed);
    f_bkg->SetLineWidth(3);
    if(isBlinded)
    {
        analysis::DrawBlindedFunction(f_bkg, blindMin, blindMax, "SAME");
    }
    else
    {
        f_bkg->Draw("SAME");
    }

    c_bkg->Modified();
    c_bkg->Update();

    if(savePlots)
    {
        c_bkg->SaveAs("./_fig/FitCombinatorialBkg.pdf");
        c_bkg->SaveAs("./_root/FitCombinatorialBkg.root");
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

    constexpr bool usePol1Bkg = false;

    // --- SANITY CHECKS ---
    // Impostare a 'true' per rilasciare i rispettivi gruppi di parametri nel fit
    constexpr bool floatSignalShape = false; // Leave it false
    constexpr bool floatP1Shape = false;
    constexpr bool floatP2Shape = false;
    constexpr bool floatArgusShape = false; // Leave it false

    // 1. Esecuzione dei fit ausiliari (MC e fondo combinatorio)
    cout << "\n=== [STEP 1/2] Running MC Auxiliary Fits ===" << endl;
    LoadDataset(1);
    AuxFitResult res_sig = FitTemplateMass(44);
    AuxFitResult res_p1 = FitTemplateMass(34);
    AuxFitResult res_p2 = FitTemplateMass(41);
    AuxFitResult res_arg = FitTemplateMass(42);

    if(!res_sig.isValid || !res_p1.isValid || !res_p2.isValid || !res_arg.isValid)
    {
        cerr << "[ERROR] Auxiliary fits failed! Aborting unbinned fit." << endl;
        return;
    }

    // 2. Caricamento dei dati reali per il Fit Unbinned
    cout << "\n=== [STEP 2/2] Preparing Real Data for Unbinned Fit ===" << endl;
    LoadDataset(0);

    Double_t xMin = 1.65;
    Double_t xMax = 2.09;

    g_data_events.clear();
    Long64_t nentries = fChain->GetEntries();
    g_data_events.reserve(static_cast<size_t>(nentries));

    auto h_mass
        = new TH1D("h_mass_unbinned", ";Invariant Mass [GeV/#it{c}^{2}];Entries", 100, xMin, xMax);
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
            h_mass->Fill(D_M);
        }
    }

    Double_t binWidth = h_mass->GetBinWidth(1);
    Double_t nEntries = h_mass->GetEntries();

    // 3. Setup della PDF e di MINUIT a 20 parametri
    g_pdf_unbinned = new FullUnbinnedPDF(xMin, xMax, res_sig, res_p1, res_p2, res_arg, usePol1Bkg);

    ROOT::Math::Minimizer *minimizer = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
    minimizer->SetMaxFunctionCalls(100000);
    minimizer->SetTolerance(0.01);
    minimizer->SetPrintLevel(0);

    ROOT::Math::Functor fNLL(&Unbinned2NLL, 20); // Impostato a 20 parametri
    minimizer->SetFunction(fNLL);

    // Parametri primari (Frazioni e pendenza)
    minimizer->SetVariable(0, "f_s", 0., 0.005);
    minimizer->SetVariable(1, "f_1", 0.025, 0.010);
    minimizer->SetVariable(2, "f_2", 0.051, 0.010);
    minimizer->SetVariable(3, "f_3", 0.625, 0.020);

    double init_slope = usePol1Bkg ? -1.14501 : -1.14501;
    if(usePol1Bkg)
        minimizer->SetVariable(4, "pol1_slope", init_slope, 0.05);
    else
        minimizer->SetVariable(4, "expo_slope", init_slope, 0.1);

    // Parametri 5-8: Segnale (Double Gauss)
    minimizer->SetVariable(5, "sig_mean", g_pdf_unbinned->m_sig_mean, 0.001);
    minimizer->SetVariable(6, "sig_sigma1", g_pdf_unbinned->m_sig_sigma1, 0.0005);
    minimizer->SetVariable(7, "sig_sigma2", g_pdf_unbinned->m_sig_sigma2, 0.001);
    minimizer->SetVariable(8, "sig_frac1", g_pdf_unbinned->m_sig_frac1, 0.05);

    // Parametri 9-12: Fondo P1 (Double Gauss)
    minimizer->SetVariable(9, "p1_mean", g_pdf_unbinned->m_p1_mean, 0.001);
    minimizer->SetVariable(10, "p1_sigma1", g_pdf_unbinned->m_p1_sigma1, 0.0005);
    minimizer->SetVariable(11, "p1_sigma2", g_pdf_unbinned->m_p1_sigma2, 0.001);
    minimizer->FixVariable(11);
    minimizer->SetVariable(12, "p1_frac1", g_pdf_unbinned->m_p1_frac1, 0.05);
    minimizer->FixVariable(12);

    // Parametri 13-16: Fondo P2 (Double Gauss)
    minimizer->SetVariable(13, "p2_mean", g_pdf_unbinned->m_p2_mean, 0.001);
    minimizer->SetVariable(14, "p2_sigma1", g_pdf_unbinned->m_p2_sigma1, 0.0005);
    minimizer->SetVariable(15, "p2_sigma2", g_pdf_unbinned->m_p2_sigma2, 0.001);
    minimizer->FixVariable(15);
    minimizer->SetVariable(16, "p2_frac1", g_pdf_unbinned->m_p2_frac1, 0.05);
    minimizer->FixVariable(16);

    // Parametri 17-19: Argus (Fondo 3)
    minimizer->SetVariable(17, "argus_m0", g_pdf_unbinned->m_argus_m0, 0.005);
    minimizer->SetVariable(18, "argus_c", g_pdf_unbinned->m_argus_c, 0.1);
    minimizer->SetVariable(19, "argus_p", g_pdf_unbinned->m_argus_p, 0.05);

    // Congela i parametri se i rispettivi flag di sanity check sono false
    if(!floatSignalShape)
    {
        minimizer->FixVariable(5);
        minimizer->FixVariable(6);
        minimizer->FixVariable(7);
        minimizer->FixVariable(8);
    }
    if(!floatP1Shape)
    {
        minimizer->FixVariable(9);
        minimizer->FixVariable(10);
    }
    if(!floatP2Shape)
    {
        minimizer->FixVariable(13);
        minimizer->FixVariable(14);
    }
    if(!floatArgusShape)
    {
        minimizer->FixVariable(17);
        minimizer->FixVariable(18);
        minimizer->FixVariable(19);
    }

    cout << "\n--- Minimizing Unbinned Likelihood with MINUIT (Silenced) ---" << endl;
    minimizer->Minimize();
    minimizer->Hesse();

    // 4. Stampa dei risultati
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
        bool isFixed = minimizer->IsFixedVariable(i);
        if(varName == "f_s" || i == 0)
        {
            cout << "  " << varName << " \t=  [BLINDED]  +/-  [BLINDED]" << endl;
        }
        else
        {
            cout << "  " << varName << " \t=  " << xs[i] << "  +/-  " << errs[i];
            if(isFixed)
                cout << "  (FIXED)";
            cout << endl;
        }
    }
    cout << "=======================================================" << endl;

    // 5. Plotting e decomposizione dei canali
    Double_t blindMin = 1.777 - 3 * 0.0058;
    Double_t blindMax = 1.777 + 3 * 0.0058;

    TCanvas *c_unbinned
        = new TCanvas("c_unbinned", "Unbinned Maximum Likelihood Fit with Pulls", 900, 900);
    double splitPoint = 0.30;

    TPad *pad1 = new TPad("pad1", "Main Fit Pad", 0.0, splitPoint, 1.0, 1.0);
    pad1->SetBottomMargin(0.02);
    pad1->Draw();
    pad1->cd();

    h_mass->GetYaxis()->SetTitle("Entries");
    AddBinSizeOnYTitle(h_mass, "GeV/#it{c}^{2}");
    h_mass->GetXaxis()->SetLabelSize(0);
    h_mass->GetXaxis()->SetTitleSize(0);

    TH1D *h_mass_blind = analysis::GetBlindedClone(h_mass, blindMin, blindMax);
    h_mass_blind->SetMinimum(0.0);
    h_mass_blind->Draw("E");

    double best_f_s = xs[0];
    double best_f_1 = xs[1];
    double best_f_2 = xs[2];
    double best_f_3 = xs[3];
    double best_slope = xs[4];

    double best_sig_mean = xs[5];
    double best_sig_sigma1 = xs[6];
    double best_sig_sigma2 = xs[7];
    double best_sig_frac1 = xs[8];

    double best_p1_mean = xs[9];
    double best_p1_sigma1 = xs[10];
    double best_p1_sigma2 = xs[11];
    double best_p1_frac1 = xs[12];

    double best_p2_mean = xs[13];
    double best_p2_sigma1 = xs[14];
    double best_p2_sigma2 = xs[15];
    double best_p2_frac1 = xs[16];

    double best_argus_m0 = xs[17];
    double best_argus_c = xs[18];
    double best_argus_p = xs[19];

    FullUnbinnedPDF *pdf_for_draw = g_pdf_unbinned;

    // Aggiorniamo le lambda expression affinché usino le variabili post-fit locali 'best_...'
    auto scale_pdf_lambda = [pdf_for_draw, nEntries, binWidth, xs](double *x, double *par)
    {
        double xx = x[0];
        return nEntries * binWidth * (*pdf_for_draw)(&xx, xs);
    };

    auto scale_sig_lambda
        = [pdf_for_draw, nEntries, binWidth, best_f_s, best_sig_mean, best_sig_sigma1,
              best_sig_sigma2, best_sig_frac1](double *x, double *par)
    {
        double xx = x[0];
        double p_sig = pdf_for_draw->Eval2G(
            xx, best_sig_mean, best_sig_sigma1, best_sig_sigma2, best_sig_frac1);
        return nEntries * binWidth * (best_f_s * p_sig);
    };

    auto scale_comb_lambda = [pdf_for_draw, nEntries, binWidth, best_f_s, best_f_1, best_f_2,
                                 best_f_3, best_slope, usePol1Bkg](double *x, double *par)
    {
        double xx = x[0];
        double p_4 = usePol1Bkg ? pdf_for_draw->EvalPol1(xx, best_slope)
                                : pdf_for_draw->EvalExpo(xx, best_slope);
        double frac_comb = (1.0 - best_f_s) * (1.0 - best_f_1 - best_f_2 - best_f_3);
        return nEntries * binWidth * (frac_comb * p_4);
    };

    auto scale_p1_lambda
        = [pdf_for_draw, nEntries, binWidth, best_f_s, best_f_1, best_p1_mean, best_p1_sigma1,
              best_p1_sigma2, best_p1_frac1](double *x, double *par)
    {
        double xx = x[0];
        double p_1
            = pdf_for_draw->Eval2G(xx, best_p1_mean, best_p1_sigma1, best_p1_sigma2, best_p1_frac1);
        return nEntries * binWidth * ((1.0 - best_f_s) * best_f_1 * p_1);
    };

    auto scale_p2_lambda
        = [pdf_for_draw, nEntries, binWidth, best_f_s, best_f_2, best_p2_mean, best_p2_sigma1,
              best_p2_sigma2, best_p2_frac1](double *x, double *par)
    {
        double xx = x[0];
        double p_2
            = pdf_for_draw->Eval2G(xx, best_p2_mean, best_p2_sigma1, best_p2_sigma2, best_p2_frac1);
        return nEntries * binWidth * ((1.0 - best_f_s) * best_f_2 * p_2);
    };

    auto scale_argus_lambda = [pdf_for_draw, nEntries, binWidth, best_f_s, best_f_3, best_argus_m0,
                                  best_argus_c, best_argus_p](double *x, double *par)
    {
        double xx = x[0];
        double p_3 = pdf_for_draw->EvalArgus(xx, best_argus_m0, best_argus_c, best_argus_p);
        return nEntries * binWidth * ((1.0 - best_f_s) * best_f_3 * p_3);
    };

    TF1 *f_draw = new TF1("f_draw", scale_pdf_lambda, xMin, xMax, 0);
    f_draw->SetNpx(10000);
    f_draw->SetLineColor(kBlue);
    f_draw->SetLineWidth(3);

    TF1 *f_sig = new TF1("f_sig", scale_sig_lambda, xMin, xMax, 0);
    f_sig->SetNpx(10000);
    f_sig->SetLineColor(kRed);
    f_sig->SetLineStyle(2);
    f_sig->SetLineWidth(2);

    TF1 *f_comb = new TF1("f_comb", scale_comb_lambda, xMin, xMax, 0);
    f_comb->SetNpx(10000);
    f_comb->SetLineColor(kGray + 2);
    f_comb->SetLineStyle(3);
    f_comb->SetLineWidth(2);

    TF1 *f_p1 = new TF1("f_p1", scale_p1_lambda, xMin, xMax, 0);
    f_p1->SetNpx(10000);
    f_p1->SetLineColor(kOrange + 1);
    f_p1->SetLineStyle(7);
    f_p1->SetLineWidth(2);

    TF1 *f_p2 = new TF1("f_p2", scale_p2_lambda, xMin, xMax, 0);
    f_p2->SetNpx(10000);
    f_p2->SetLineColor(kMagenta);
    f_p2->SetLineStyle(7);
    f_p2->SetLineWidth(2);

    TF1 *f_argus = new TF1("f_argus", scale_argus_lambda, xMin, xMax, 0);
    f_argus->SetNpx(10000);
    f_argus->SetLineColor(kCyan + 1);
    f_argus->SetLineStyle(5);
    f_argus->SetLineWidth(2);

    analysis::DrawBlindedFunction(f_draw, blindMin, blindMax, "SAME");
    analysis::DrawBlindedFunction(f_sig, blindMin, blindMax, "SAME");
    analysis::DrawBlindedFunction(f_comb, blindMin, blindMax, "SAME");
    analysis::DrawBlindedFunction(f_p1, blindMin, blindMax, "SAME");
    analysis::DrawBlindedFunction(f_p2, blindMin, blindMax, "SAME");
    analysis::DrawBlindedFunction(f_argus, blindMin, blindMax, "SAME");

    TH1D *hPull = (TH1D *)h_mass->Clone("hPull_unbinned");
    hPull->Reset();
    hPull->SetStats(0);
    for(int i = 1; i <= hPull->GetNbinsX(); i++)
    {
        hPull->SetBinContent(i, -999.0);
        hPull->SetBinError(i, 0.0);
    }

    double chi2 = 0.0;
    int nBinsUsed = 0;
    for(int i = h_mass->FindBin(xMin); i <= h_mass->FindBin(xMax); i++)
    {
        double x = h_mass->GetBinCenter(i);

        if(x >= blindMin && x <= blindMax)
            continue;

        double obs = h_mass->GetBinContent(i);
        double err = h_mass->GetBinError(i);
        if(err > 0)
        {
            double val = f_draw->Eval(x);
            double pull = (obs - val) / err;
            hPull->SetBinContent(i, pull);
            hPull->SetBinError(i, 0.0);
            chi2 += pull * pull;
            nBinsUsed++;
        }
    }

    int ndf = nBinsUsed - minimizer->NFree();

    TPaveText *pave = new TPaveText(0.18, 0.62, 0.52, 0.88, "NDC");
    pave->SetBorderSize(0);
    pave->SetFillStyle(0);
    pave->SetTextFont(42);
    pave->SetTextSize(0.033);

    pave->AddText(Form("#chi^{2} / ndf = %.1f / %d", chi2, ndf));
    pave->AddText(Form("Prob = %.1f%%", TMath::Prob(chi2, ndf) * 100.0));
    pave->Draw();

    TLegend *leg = new TLegend(0.55, 0.50, 0.95, 0.88);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextFont(42);
    leg->SetTextSize(0.028);
    leg->AddEntry(h_mass_blind, "Data (Blinded)", "ep");
    leg->AddEntry(f_draw, "Total Fit", "l");
    leg->AddEntry(f_sig, "Signal (D_{s}^{+}#rightarrow#tau^{+}#nu_{#tau})", "l");
    leg->AddEntry(f_p1, "D^{+}#rightarrow#phi#pi^{+} Bkg", "l");
    leg->AddEntry(f_p2, "D_{s}^{+}#rightarrow#phi#pi^{+} Bkg", "l");
    leg->AddEntry(f_argus, "D_{s}^{+}#rightarrow#phi#mu^{+}#nu_{#mu} (Argus) Bkg", "l");
    if(usePol1Bkg)
        leg->AddEntry(f_comb, "Combinatorial Bkg (Pol1)", "l");
    else
        leg->AddEntry(f_comb, "Combinatorial Bkg (Expo)", "l");
    leg->Draw("SAME");

    // Pad Inferiore (Residui)
    c_unbinned->cd();
    TPad *pad2 = new TPad("pad2", "Pull Pad", 0.0, 0.0, 1.0, splitPoint);
    pad2->SetTopMargin(0.02);
    pad2->SetBottomMargin(0.35);
    pad2->SetGridy();
    pad2->Draw();
    pad2->cd();

    hPull->GetYaxis()->SetTitle("Pull");
    hPull->GetYaxis()->SetTitleSize(gStyle->GetTitleSize("Y"));
    hPull->GetYaxis()->SetLabelSize(gStyle->GetLabelSize("Y"));
    hPull->GetYaxis()->SetTitleOffset(gStyle->GetTitleOffset("Y"));
    hPull->GetYaxis()->SetRangeUser(-5, 5);
    hPull->GetXaxis()->SetTitle("M(D_{s}^{+}) [GeV/#it{c}^{2}]");
    hPull->GetXaxis()->SetTitleSize(gStyle->GetTitleSize("X"));
    hPull->GetXaxis()->SetLabelSize(gStyle->GetLabelSize("X"));
    hPull->GetXaxis()->SetTitleOffset(gStyle->GetTitleOffset("X"));
    hPull->SetMarkerStyle(20);
    hPull->SetMarkerSize(0.8);
    hPull->SetMarkerColor(kBlack);
    hPull->SetLineColor(kBlack);
    hPull->Draw("P");

    TLine *line0 = new TLine(xMin, 0.0, xMax, 0.0);
    line0->SetLineColor(kRed);
    line0->SetLineWidth(2);
    line0->Draw();

    c_unbinned->Update();

    if(savePlots)
    {
        c_unbinned->SaveAs("./_fig/FitFullUnbinned_Pulls.pdf");
        c_unbinned->SaveAs("./_root/FitFullUnbinned_Pulls.root");
    }

    // ========================================================================
    // POINT 7: VALIDATION OF THE PROCEDURE (MEASUREMENT OF R_phi_pi)
    // ========================================================================
    cout << "\n=======================================================" << endl;
    cout << "   VALIDATION CHECK: MEASUREMENT OF R_phi_pi" << endl;
    cout << "=======================================================" << endl;

    // 1. Calcolo esatto delle efficienze binomiali dai MC
    // ID 34 = D+ -> phi pi, ID 41 = Ds+ -> phi pi (Tabella 1: 50M generati ciascuno)
    const double n_gen_Dplus = 50.0e6;
    const double n_gen_Dsplus = 50.0e6;

    LoadDataset(1); // Carichiamo il file MC
    long long n_pass_Dplus = 0;
    long long n_pass_Dsplus = 0;

    for(Long64_t jentry = 0; jentry < fChain->GetEntriesFast(); jentry++)
    {
        if(LoadTree(jentry) < 0)
            break;
        fChain->GetEntry(jentry);
        if(id == 34)
            n_pass_Dplus++;
        if(id == 41)
            n_pass_Dsplus++;
    }

    // Efficienze fisiche
    double eff_Dplus = (double)n_pass_Dplus / n_gen_Dplus;
    double eff_Dsplus = (double)n_pass_Dsplus / n_gen_Dsplus;

    // Rapporto delle efficienze
    double ratio_eff = eff_Dsplus / eff_Dplus;

    // --- PROPAGAZIONE BINOMIALE ESATTA ---
    double rel_var_Dsplus = (1.0 - eff_Dsplus) / n_pass_Dsplus;
    double rel_var_Dplus = (1.0 - eff_Dplus) / n_pass_Dplus;

    double err_ratio_eff = ratio_eff * std::sqrt(rel_var_Dsplus + rel_var_Dplus);

    // 2. Calcolo del rapporto delle frazioni dal fit ai dati (f1 / f2)
    double f1_val = xs[1];
    double f2_val = xs[2];
    double f1_err = errs[1];
    double f2_err = errs[2];
    double cov_f1_f2 = minimizer->CovMatrix(1, 2);

    double ratio_yield = f1_val / f2_val;

    // Propagazione dell'errore esatta per una divisione usando la covarianza
    double err_ratio_yield = ratio_yield
        * std::sqrt((f1_err * f1_err) / (f1_val * f1_val) + (f2_err * f2_err) / (f2_val * f2_val)
            - (2.0 * cov_f1_f2) / (f1_val * f2_val));

    // 3. Valore Finale Misurato
    double R_phipi_meas = ratio_yield * ratio_eff;

    // Errore totale (errore fit e errore MC si sommano in quadratura come incertezze relative)
    double err_R_phipi_meas = R_phipi_meas
        * std::sqrt((err_ratio_yield / ratio_yield) * (err_ratio_yield / ratio_yield)
            + (err_ratio_eff / ratio_eff) * (err_ratio_eff / ratio_eff));

    // 4. Valore Teorico dal PDF (Pagina 4)
    double sigma_Dplus = 834.0;
    double sigma_Dsplus = 353.0;
    double BR_Dplus = 2.69e-3;
    double BR_Dsplus = 2.25e-2;

    // Nota: Il testo dice che l'errore grande sulle sezioni d'urto si cancella ampiamente nel
    // rapporto. Propaghiamo solo gli errori statistici decorrelati dei Branching Ratios.
    double err_BR_Dplus = 0.08e-3;
    double err_BR_Dsplus = 0.05e-2;

    double R_phipi_ref = (sigma_Dplus * BR_Dplus) / (sigma_Dsplus * BR_Dsplus);
    double err_R_phipi_ref = R_phipi_ref
        * std::sqrt((err_BR_Dplus / BR_Dplus) * (err_BR_Dplus / BR_Dplus)
            + (err_BR_Dsplus / BR_Dsplus) * (err_BR_Dsplus / BR_Dsplus));

    // 5. Confronto (Pull)
    double diff = std::abs(R_phipi_meas - R_phipi_ref);
    double combined_err
        = std::sqrt(err_R_phipi_meas * err_R_phipi_meas + err_R_phipi_ref * err_R_phipi_ref);
    double pull = diff / combined_err;

    cout << Form("  Eff Ratio (Ds / D+) : %.4f +/- %.4f", ratio_eff, err_ratio_eff) << endl;
    cout << Form("  Yield Ratio (f1/f2) : %.4f +/- %.4f", ratio_yield, err_ratio_yield) << endl;
    cout << "-------------------------------------------------------" << endl;
    cout << Form("  Measured R_phi_pi   : %.4f +/- %.4f", R_phipi_meas, err_R_phipi_meas) << endl;
    cout << Form("  Expected R_phi_pi   : %.4f +/- %.4f", R_phipi_ref, err_R_phipi_ref) << endl;
    cout << "-------------------------------------------------------" << endl;
    cout << Form("  Compatibility       : %.2f sigma", pull) << endl;
    if(pull < 3.0)
        cout << "  => SUCCESS: The analysis procedure is validated!" << endl;
    else
        cout << "  => WARNING: Discrepancy observed. Check fit model!" << endl;
    cout << "=======================================================\n" << endl;

    // Riporta LoadDataset a 0 per sicurezza se devi fare altre operazioni dopo
    LoadDataset(0);

    delete minimizer;
    g_pdf_unbinned = nullptr;
}

// ================================================================================
// ToyMC Simulation
// ================================================================================

struct ToyResult
{
    bool converged = false;
    double fs_val = 0, fs_err = 0;
    double f1_val = 0, f1_err = 0;
    double f2_val = 0, f2_err = 0;
    double f3_val = 0, f3_err = 0;
    double slope_val = 0, slope_err = 0;
    double cov_fs_f3 = 0; // Covarianza tra f_s (parametro 0) e f_3 (parametro 3)
};

// Funzione worker isolata per il singolo thread
ToyResult RunSingleToy(int toyId, int nEvents, const FullUnbinnedPDF &templatePdf,
    const std::vector<double> &gen_pars, bool usePol1Bkg, double xMin, double xMax)
{
    ToyResult res;

    // 1. Generatore di numeri casuali locale al thread (evita conflitti su gRandom)
    TRandom3 threadRandom(THE_SEED * (toyId + 1));
    // cout << Form("Toy %d: seed = %d", toyId, threadRandom.GetSeed()) << endl;

    // 2. Copia locale della PDF per questo thread
    FullUnbinnedPDF localPdf = templatePdf;
    g_pdf_unbinned = &localPdf; // Assegna il puntatore thread-local

    // 3. Wrapper locale TF1 per la generazione unbinned
    auto gen_lambda = [&localPdf, &gen_pars](double *x, double *par)
    {
        double xx = x[0];
        return localPdf(&xx, gen_pars.data());
    };
    TF1 fGen(Form("fGen_%d", toyId), gen_lambda, xMin, xMax, 0);
    fGen.SetNpx(10000);

    // 4. Generazione del dataset locale al thread
    g_data_events.clear();
    // nEvents = threadRandom->Poisson(nEvents);
    g_data_events.reserve(nEvents);
    for(int ev = 0; ev < nEvents; ++ev)
    {
        g_data_events.push_back(fGen.GetRandom(&threadRandom));
    }

    // 5. Setup locale del Minimizzatore
    ROOT::Math::Minimizer *minimizer = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
    minimizer->SetMaxFunctionCalls(50000);
    minimizer->SetTolerance(0.01);
    minimizer->SetPrintLevel(-1); // Silenzioso

    ROOT::Math::Functor fNLL(&Unbinned2NLL, 20);
    minimizer->SetFunction(fNLL);

    // =========================================================================
    //   STRESS TEST: RANDOM MULTISTART INITIALIZATION
    // =========================================================================

    // 1. f_s (segnale) fluttua casualmente tra -0.05 e +0.05 attorno al valore vero
    double start_fs = gen_pars[0] + threadRandom.Uniform(-0.02, 0.02);

    // 2. Le frazioni dei fondi f_1, f_2, f_3 variano casualmente del +/- 10%
    //    rispetto al loro valore vero (es. moltiplicate per un fattore tra 0.8 e 1.2)
    double start_f1 = gen_pars[1] * threadRandom.Uniform(0.9, 1.1);
    double start_f2 = gen_pars[2] * threadRandom.Uniform(0.9, 1.1);
    double start_f3 = gen_pars[3] * threadRandom.Uniform(0.9, 1.1);

    // 3. La pendenza del fondo combinatorio varia casualmente del +/- 10%
    double start_slope = gen_pars[4] * threadRandom.Uniform(0.9, 1.1);

    // Impostiamo le variabili di partenza di Minuit con questi valori casuali
    minimizer->SetVariable(0, "f_s", start_fs, 0.005);
    minimizer->SetVariable(1, "f_1", start_f1, 0.005);
    minimizer->SetVariable(2, "f_2", start_f2, 0.005);
    minimizer->SetVariable(3, "f_3", start_f3, 0.005);

    // Definizione delle variabili locali
    // minimizer->SetVariable(0, "f_s", gen_pars[0], 0.005);
    // minimizer->FixVariable(0);
    // minimizer->SetVariable(1, "f_1", gen_pars[1], 0.005);
    // minimizer->SetVariable(2, "f_2", gen_pars[2], 0.005);
    // minimizer->FixVariable(2);
    // minimizer->SetVariable(3, "f_3", gen_pars[3], 0.005);
    // minimizer->FixVariable(3);

    if(usePol1Bkg)
    {
        minimizer->SetVariable(4, "pol1_slope", start_slope, 0.05);
    }
    else
    {
        minimizer->SetVariable(4, "expo_slope", start_slope, 0.005);
    }

    // Congelamento dei parametri di forma (fissati ai parametri nominali di generazione)
    minimizer->SetFixedVariable(5, "sig_mean", gen_pars[5]);
    minimizer->SetFixedVariable(6, "sig_sigma1", gen_pars[6]);
    minimizer->SetFixedVariable(7, "sig_sigma2", gen_pars[7]);
    minimizer->SetFixedVariable(8, "sig_frac1", gen_pars[8]);

    minimizer->SetFixedVariable(9, "p1_mean", gen_pars[9]);
    minimizer->SetFixedVariable(10, "p1_sigma1", gen_pars[10]);
    minimizer->SetFixedVariable(11, "p1_sigma2", gen_pars[11]);
    minimizer->SetFixedVariable(12, "p1_frac1", gen_pars[12]);

    minimizer->SetFixedVariable(13, "p2_mean", gen_pars[13]);
    minimizer->SetFixedVariable(14, "p2_sigma1", gen_pars[14]);
    minimizer->SetFixedVariable(15, "p2_sigma2", gen_pars[15]);
    minimizer->SetFixedVariable(16, "p2_frac1", gen_pars[16]);

    minimizer->SetFixedVariable(17, "argus_m0", gen_pars[17]);
    minimizer->SetFixedVariable(18, "argus_c", gen_pars[18]);
    minimizer->SetFixedVariable(19, "argus_p", gen_pars[19]);

    minimizer->Minimize();
    // minimizer->Hesse(); // Abilitato per ottenere la matrice di covarianza accurata

    int status = minimizer->Status();

    // Se il fit ha fallito (status != 0)
    if(status != 0)
    {
        // Creiamo un file di log specifico per questo singolo toy fallito
        std::ofstream logFile(Form("./_fig/failed_toy_%d_debug.txt", toyId));

        logFile << "=========================================\n";
        logFile << "   DIAGNOSTIC LOG FOR FAILED TOY #" << toyId << "\n";
        logFile << "=========================================\n";
        logFile << "Minimizer Status: " << status << "\n";
        logFile << "EDM:              " << minimizer->Edm() << "\n\n";

        // 1. Salviamo lo stato dei parametri al momento del fallimento
        logFile << "Parameter values at failure:\n";
        for(unsigned int p = 0; p < minimizer->NDim(); ++p)
        {
            logFile << "  " << minimizer->VariableName(p) << " = " << minimizer->X()[p];
            if(minimizer->IsFixedVariable(p))
                logFile << " (FIXED)";
            logFile << "\n";
        }
        logFile << "\n";

        // 2. Controlliamo se ci sono eventi che mandano la PDF a zero o negativa
        logFile << "Checking event PDF values at failure point:\n";
        int negative_pdf_count = 0;
        for(size_t i = 0; i < g_data_events.size(); ++i)
        {
            double x = g_data_events[i];
            double pdf_val = (*g_pdf_unbinned)(&x, minimizer->X());
            if(pdf_val <= 0.0)
            {
                negative_pdf_count++;
                logFile << "  [WARNING] Event #" << i << " (x=" << x << ") has PDF = " << pdf_val
                        << " (<= 0!)\n";
            }
        }
        logFile << "\nTotal events with PDF <= 0: " << negative_pdf_count << "\n";

        // 3. Salva anche i dati di questo toy specifico così puoi rifittarlo da solo
        logFile << "\nDataset events:\n";
        for(double x : g_data_events)
        {
            logFile << x << "\n";
        }

        logFile.close();
    }

    if(minimizer->Status() == 0)
    {
        res.converged = true;
        res.fs_val = minimizer->X()[0];
        res.fs_err = minimizer->Errors()[0];
        res.f1_val = minimizer->X()[1];
        res.f1_err = minimizer->Errors()[1];
        res.f2_val = minimizer->X()[2];
        res.f2_err = minimizer->Errors()[2];
        res.f3_val = minimizer->X()[3];
        res.f3_err = minimizer->Errors()[3];
        res.slope_val = minimizer->X()[4];
        res.slope_err = minimizer->Errors()[4];

        // Estrazione dell'elemento cov(0, 3) della matrice di errore
        res.cov_fs_f3 = minimizer->CovMatrix(0, 3);
    }

    delete minimizer;
    g_pdf_unbinned = nullptr; // Reset
    return res;
}

std::pair<double, double> analysis::RunToyMC(int nToys, double true_fs, bool useBR)
{
    auto start = std::chrono::high_resolution_clock::now();

    SetLBStyle();
    constexpr bool usePol1Bkg = false;
    // Generatore di numeri casuali locale al thread principale
    TRandom3 mainRandom(THE_SEED * (static_cast<int>(std::thread::hardware_concurrency())) + 1);

    TF1::DefaultAddToGlobalList(kFALSE);

    // --- [STEP 1] Fit ausiliari nominali ---
    cout << "\n=== [TOY MC] Running Auxiliary Fits ===" << endl;
    LoadDataset(1);
    AuxFitResult res_sig = FitTemplateMass(44);
    AuxFitResult res_p1 = FitTemplateMass(34);
    AuxFitResult res_p2 = FitTemplateMass(41);
    AuxFitResult res_arg = FitTemplateMass(42);

    if(!res_sig.isValid || !res_p1.isValid || !res_p2.isValid || !res_arg.isValid)
    {
        cerr << "[ERROR] Auxiliary fits failed! Aborting Toy MC." << endl;
        return std::make_pair(0.0, 0.0);
    }

    // --- [STEP 2] Stima delle efficienze dal Monte Carlo ---
    cout << "\n--> Estimating Selection Efficiencies from MC..." << endl;
    const double n_gen_sig = 1.0e6;
    const double n_gen_norm = 5.0e6;

    long long n_pass_sig = 0;
    long long n_pass_norm = 0;

    Long64_t nentries_mc = fChain->GetEntries();
    for(Long64_t jentry = 0; jentry < nentries_mc; jentry++)
    {
        if(LoadTree(jentry) < 0)
            break;
        fChain->GetEntry(jentry);
        if(id == 44)
            n_pass_sig++; // Ds+ -> tau+ nu_tau (Segnale)
        if(id == 42)
            n_pass_norm++; // Ds+ -> phi mu+ nu_mu (Normalizzazione)
    }

    double eff_sig = (double)n_pass_sig / n_gen_sig;
    double eff_norm = (double)n_pass_norm / n_gen_norm;
    double r_eff = eff_norm / eff_sig;

    // --- Formula binomiale esatta per l'incertezza del rapporto ---
    double r_eff_err
        = r_eff * std::sqrt((1.0 - eff_sig) / n_pass_sig + (1.0 - eff_norm) / n_pass_norm);

    cout << Form("  N_pass (sig)  = %lld, Eff (sig)  = %.6f", n_pass_sig, eff_sig) << endl;
    cout << Form("  N_pass (norm) = %lld, Eff (norm) = %.6f", n_pass_norm, eff_norm) << endl;
    cout << Form("  Ratio eff_norm / eff_sig (Binomial) = %.4f +/- %.4f", r_eff, r_eff_err) << endl;

    // --- [STEP 3] Definizione costanti esterne (PDG) ---
    const double br_taunu_nom = 5.39e-2;
    const double br_taunu_err = 0.09e-2;
    const double br_phimunu_nom = 2.24e-2;
    const double br_phimunu_err = 0.11e-2;

    double k_factor_nom = r_eff * (br_phimunu_nom / br_taunu_nom);

    const double true_fs_val = true_fs;
    const double true_f1 = 0.026;
    const double true_f2 = 0.051;
    const double true_f3 = 0.625;
    const double true_slope = -1.145;

    double true_est = true_fs_val / ((1.0 - true_fs_val) * true_f3);
    double true_br = true_est * k_factor_nom;

    std::vector<double> gen_pars = { true_fs_val, true_f1, true_f2, true_f3, true_slope,
        res_sig.params[1], res_sig.params[2], res_sig.params[3], res_sig.params[4],
        res_p1.params[1], res_p1.params[2], res_p1.params[3], res_p1.params[4], res_p2.params[1],
        res_p2.params[2], res_p2.params[3], res_p2.params[4], res_arg.params[1], res_arg.params[2],
        res_arg.params[3] };

    Double_t xMin = 1.65;
    Double_t xMax = 2.09;

    FullUnbinnedPDF templatePdf(xMin, xMax, res_sig, res_p1, res_p2, res_arg, usePol1Bkg);

    LoadDataset(0);
    int nEvents = 0;
    for(Long64_t jentry = 0; jentry < fChain->GetEntriesFast(); jentry++)
    {
        if(LoadTree(jentry) < 0)
            break;
        fChain->GetEntry(jentry);
        if(id == 0 && D_M >= xMin && D_M <= xMax)
        {
            nEvents++;
        }
    }
    if(nEvents == 0)
        nEvents = 10000;

    int nCores = static_cast<int>(std::thread::hardware_concurrency());
    if(nCores == 0)
        nCores = 4;
    cout << "[INFO] Launching Toy MC using " << nCores << " parallel threads." << endl;

    std::vector<ToyResult> toyResults;
    toyResults.reserve(nToys);

    for(int i = 0; i < nToys; i += nCores)
    {
        std::vector<std::future<ToyResult>> futures;
        unsigned int currentBatchSize = 0;

        for(int t = 0; (t < nCores) && (i + t) < nToys; ++t)
        {
            int toyId = i + t;
            futures.push_back(std::async(std::launch::async, RunSingleToy, toyId, nEvents,
                std::ref(templatePdf), std::cref(gen_pars), usePol1Bkg, xMin, xMax));
            currentBatchSize++;
        }

        for(auto &f : futures)
        {
            toyResults.push_back(f.get());
        }
        cout << "  Completed toys: " << toyResults.size() << " / " << nToys << "..." << endl;
    }

    // --- [STEP 4] Booking Istogrammi ---
    // Istogrammi f_s
    auto h_fit_fs = new TH1D("h_fit_fs", "Fitted f_{s};#Delta f_{s};Toys", 50, -0, 0);
    auto h_pull_fs = new TH1D("h_pull_fs", "Pull f_{s};Pull;Toys", 50, -5.0, 5.0);

    // Istogrammi per l'estimatore di resa grezzo Y = f_s / ((1 - f_s) * f_3)
    auto h_fit_est = new TH1D("h_fit_est", "Fitted Yield Ratio Y;#Delta Y;Toys", 50, -0, 0);
    auto h_pull_est = new TH1D("h_pull_est", "Pull Yield Ratio Y;Pull;Toys", 50, -5.0, 5.0);

    // Istogrammi per il Branching Ratio fisico
    auto h_fit_br = new TH1D(
        "h_fit_br", "Fitted #f{B}(#tau^{+}#rightarrow#phi#mu^{+});#Delta#f{B};Toys", 50, -0, 0);
    auto h_pull_br = new TH1D(
        "h_pull_br", "Pull #f{B}(#tau^{+}#rightarrow#phi#mu^{+});Pull;Toys", 50, -5.0, 5.0);

    // Altri parametri di background
    auto h_fit_f1 = new TH1D("h_fit_f1", "Fitted f_{1};#Delta f_{1};Toys", 50, -0, 0);
    auto h_fit_f2 = new TH1D("h_fit_f2", "Fitted f_{2};#Delta f_{2};Toys", 50, -0, 0);
    auto h_fit_f3 = new TH1D("h_fit_f3", "Fitted f_{3};#Delta f_{3};Toys", 50, -0, 0);
    auto h_fit_slope = new TH1D("h_fit_slope", "Fitted Slope;#Delta Slope;Toys", 50, -0, 0);

    auto h_pull_f1 = new TH1D("h_pull_f1", "Pull f_{1};Pull;Toys", 50, -5.0, 5.0);
    auto h_pull_f2 = new TH1D("h_pull_f2", "Pull f_{2};Pull;Toys", 50, -5.0, 5.0);
    auto h_pull_f3 = new TH1D("h_pull_f3", "Pull f_{3};Pull;Toys", 50, -5.0, 5.0);
    auto h_pull_slope = new TH1D("h_pull_slope", "Pull Slope;Pull;Toys", 50, -5.0, 5.0);

    // --- [STEP 5] Analisi risultati dei Toy con Smearing dei parametri esterni ---
    int convergedToys = 0;
    for(const auto &res : toyResults)
    {
        if(res.converged)
        {
            convergedToys++;

            h_fit_fs->Fill(res.fs_val - true_fs_val);
            h_fit_f1->Fill(res.f1_val - true_f1);
            h_fit_f2->Fill(res.f2_val - true_f2);
            h_fit_f3->Fill(res.f3_val - true_f3);
            h_fit_slope->Fill(res.slope_val - true_slope);

            if(res.fs_err > 0)
                h_pull_fs->Fill((res.fs_val - true_fs_val) / res.fs_err);
            if(res.f1_err > 0)
                h_pull_f1->Fill((res.f1_val - true_f1) / res.f1_err);
            if(res.f2_err > 0)
                h_pull_f2->Fill((res.f2_val - true_f2) / res.f2_err);
            if(res.f3_err > 0)
                h_pull_f3->Fill((res.f3_val - true_f3) / res.f3_err);
            if(res.slope_err > 0)
                h_pull_slope->Fill((res.slope_val - true_slope) / res.slope_err);

            double u = res.fs_val;
            double v = res.f3_val;
            double du = res.fs_err;
            double dv = res.f3_err;
            double cov_uv = res.cov_fs_f3;

            if(std::abs((1.0 - u) * v) > 1e-9)
            {
                double est_val = u / ((1.0 - u) * v);
                h_fit_est->Fill(est_val - true_est);

                // Calcolo errore propagato dell'estimatore Y
                double dF_du = 1.0 / ((1.0 - u) * (1.0 - u) * v);
                double dF_dv = -u / ((1.0 - u) * v * v);
                double variance_est = (dF_du * dF_du * du * du) + (dF_dv * dF_dv * dv * dv)
                    + (2.0 * dF_du * dF_dv * cov_uv);

                if(variance_est > 0.0)
                {
                    h_pull_est->Fill((est_val - true_est) / std::sqrt(variance_est));
                }

                // Smearing dei parametri esterni per questo toy (Incertezze statistiche e PDG)
                double r_eff_toy = mainRandom.Gaus(r_eff, r_eff_err);
                double br_phimunu_toy = mainRandom.Gaus(br_phimunu_nom, br_phimunu_err);
                double br_taunu_toy = mainRandom.Gaus(br_taunu_nom, br_taunu_err);

                double k_factor_toy = r_eff_toy * (br_phimunu_toy / br_taunu_toy);

                // Calcolo del Branching Ratio per il toy corrente
                double br_val = est_val * k_factor_toy;
                h_fit_br->Fill(br_val - true_br);

                // Incertezza relativa al quadrato del fattore K (usa l'errore binomiale esatto
                // r_eff_err)
                double rel_err_k2 = (r_eff_err / r_eff) * (r_eff_err / r_eff)
                    + (br_phimunu_err / br_phimunu_nom) * (br_phimunu_err / br_phimunu_nom)
                    + (br_taunu_err / br_taunu_nom) * (br_taunu_err / br_taunu_nom);
                double err_k = k_factor_nom * std::sqrt(rel_err_k2);

                // Incertezza totale sul BR (composizione dell'errore di Y e di K, indipendenti)
                double variance_br = (k_factor_nom * k_factor_nom * variance_est)
                    + (est_val * est_val * err_k * err_k);

                if(variance_br > 0.0)
                {
                    double br_err = std::sqrt(variance_br);
                    h_pull_br->Fill((br_val - true_br) / br_err);
                }
            }
        }
    }

    cout << "\n[TOY MC RESULTS] Converged: " << convergedToys << " / " << nToys << endl;

    TF1::DefaultAddToGlobalList(kTRUE);

    // --- [STEP 6] Disegno delle due Canvas richieste ---
    auto drawResult = [](TVirtualPad *pad, TH1D *h, double trueVal)
    {
        pad->cd();
        h->SetStats(kTRUE);
        gStyle->SetOptStat("emr");
        h->Draw();
        h->Fit("gaus", "Q L I");
        gStyle->SetOptFit(111);
    };

    // Canvas 1 aggiornata: fs, Yield Ratio (Y) e Branching Ratio (BR) disposti su 3 colonne
    TCanvas *cMainEst
        = new TCanvas("cMainEst", "Main Estimators (fs, Yield Ratio, BR)", 1500, 1000);
    cMainEst->Divide(3, 2);

    drawResult(cMainEst->GetPad(1), h_fit_fs, 0.0);
    drawResult(cMainEst->GetPad(2), h_fit_est, 0.0);
    drawResult(cMainEst->GetPad(3), h_fit_br, 0.0);

    drawResult(cMainEst->GetPad(4), h_pull_fs, 0.0);
    drawResult(cMainEst->GetPad(5), h_pull_est, 0.0);
    drawResult(cMainEst->GetPad(6), h_pull_br, 0.0);

    cMainEst->Update();
    if(savePlots)
    {
        cMainEst->SaveAs("./_fig/ToyMC_MainEstimators.pdf");
        cMainEst->SaveAs("./_root/ToyMC_MainEstimators.root");
    }

    // Canvas 2: Altri parametri di fit (f1, f2, f3, slope)
    TCanvas *cOtherPars = new TCanvas("cOtherPars", "Other Fit Parameters", 1600, 800);
    cOtherPars->Divide(4, 2);

    drawResult(cOtherPars->GetPad(1), h_fit_f1, 0.0);
    drawResult(cOtherPars->GetPad(2), h_fit_f2, 0.0);
    drawResult(cOtherPars->GetPad(3), h_fit_f3, 0.0);
    drawResult(cOtherPars->GetPad(4), h_fit_slope, 0.0);

    drawResult(cOtherPars->GetPad(5), h_pull_f1, 0.0);
    drawResult(cOtherPars->GetPad(6), h_pull_f2, 0.0);
    drawResult(cOtherPars->GetPad(7), h_pull_f3, 0.0);
    drawResult(cOtherPars->GetPad(8), h_pull_slope, 0.0);

    cOtherPars->Update();
    if(savePlots)
    {
        cOtherPars->SaveAs("./_fig/ToyMC_OtherParameters.pdf");
        cOtherPars->SaveAs("./_root/ToyMC_OtherParameters.root");
    }

    // ================================================================================
    // Sezione diagnostica finale (singolo toy)
    // ================================================================================
    cout << "\n=== [TOY MC] Generating and fitting a single Toy example for diagnostics ==="
         << endl;

    // Assegniamo la PDF di riferimento
    g_pdf_unbinned = &templatePdf;

    auto gen_lambda_diag = [&templatePdf, &gen_pars](double *x, double *par)
    {
        double xx = x[0];
        return templatePdf(&xx, gen_pars.data());
    };
    TF1 fGenDiag("fGenDiag", gen_lambda_diag, xMin, xMax, 0);
    fGenDiag.SetNpx(2000);

    // Vettore locale e istogramma per il Toy di esempio
    std::vector<double> diag_events;
    diag_events.reserve(nEvents);
    auto h_single_toy = new TH1D("h_single_toy",
        "Single Toy Pseudo-Data Fit;Mass [GeV/#it{c}^{2}];Entries", 100, xMin, xMax);
    h_single_toy->SetDirectory(nullptr);
    AddBinSizeOnYTitle(h_single_toy, "GeV/#it{c}^{2}");

    for(int ev = 0; ev < nEvents; ++ev)
    {
        double val = fGenDiag.GetRandom(&mainRandom);
        diag_events.push_back(val);
        h_single_toy->Fill(val);
    }

    // Eseguiamo il fit di questo dataset di esempio nel thread principale
    g_data_events = diag_events;

    ROOT::Math::Minimizer *minDiag = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
    minDiag->SetMaxFunctionCalls(50000);
    minDiag->SetTolerance(0.01);
    minDiag->SetPrintLevel(0); // Mostra l'output di questo singolo fit di controllo

    ROOT::Math::Functor fNLLDiag(&Unbinned2NLL, 20);
    minDiag->SetFunction(fNLLDiag);

    // Inizializzazione parametri (vicini ai valori di generazione)
    minDiag->SetVariable(0, "f_s", 0.01, 0.005);
    minDiag->SetVariable(1, "f_1", 0.026, 0.005);
    minDiag->SetVariable(2, "f_2", 0.051, 0.005);
    minDiag->SetVariable(3, "f_3", 0.625, 0.005);

    if(usePol1Bkg)
        minDiag->SetVariable(4, "pol1_slope", gen_pars[4], 0.005);
    else
        minDiag->SetVariable(4, "expo_slope", gen_pars[4], 0.005);

    // Fissiamo tutte le forme ai parametri veri
    for(int p = 5; p < 20; ++p)
    {
        minDiag->SetVariable(p, Form("p_%d", p), gen_pars[p], 0.005);
        minDiag->FixVariable(p);
    }

    minDiag->Minimize();
    minDiag->Hesse();

    const double *xsDiag = minDiag->X();

    // Disegniamo l'istogramma e vi sovrapponiamo le componenti del fit
    TCanvas *cSingleToy = new TCanvas("cSingleToy", "Single Toy Fit Example", 800, 600);
    cSingleToy->cd();
    h_single_toy->SetMinimum(0.0);
    h_single_toy->Draw("E");

    double binWidthDiag = h_single_toy->GetBinWidth(1);
    double nEntriesDiag = h_single_toy->GetEntries();

    // --- DEFINIZIONE DELLE LAMBDA (PRIMA dei TF1) ---
    auto scale_pdf_diag = [&templatePdf, nEntriesDiag, binWidthDiag, xsDiag](double *x, double *par)
    {
        double xx = x[0];
        return nEntriesDiag * binWidthDiag * templatePdf(&xx, xsDiag);
    };

    auto scale_sig_diag = [&templatePdf, nEntriesDiag, binWidthDiag, xsDiag](double *x, double *par)
    {
        double xx = x[0];
        double p_sig = templatePdf.Eval2G(xx, xsDiag[5], xsDiag[6], xsDiag[7], xsDiag[8]);
        return nEntriesDiag * binWidthDiag * (xsDiag[0] * p_sig);
    };

    auto scale_comb_diag
        = [&templatePdf, nEntriesDiag, binWidthDiag, xsDiag, usePol1Bkg](double *x, double *par)
    {
        double xx = x[0];
        double p_4 = usePol1Bkg ? templatePdf.EvalPol1(xx, xsDiag[4])
                                : templatePdf.EvalExpo(xx, xsDiag[4]);
        double frac_comb = (1.0 - xsDiag[0]) * (1.0 - xsDiag[1] - xsDiag[2] - xsDiag[3]);
        return nEntriesDiag * binWidthDiag * (frac_comb * p_4);
    };

    auto scale_p1_diag = [&templatePdf, nEntriesDiag, binWidthDiag, xsDiag](double *x, double *par)
    {
        double xx = x[0];
        double p_1 = templatePdf.Eval2G(xx, xsDiag[9], xsDiag[10], xsDiag[11], xsDiag[12]);
        return nEntriesDiag * binWidthDiag * ((1.0 - xsDiag[0]) * xsDiag[1] * p_1);
    };

    auto scale_p2_diag = [&templatePdf, nEntriesDiag, binWidthDiag, xsDiag](double *x, double *par)
    {
        double xx = x[0];
        double p_2 = templatePdf.Eval2G(xx, xsDiag[13], xsDiag[14], xsDiag[15], xsDiag[16]);
        return nEntriesDiag * binWidthDiag * ((1.0 - xsDiag[0]) * xsDiag[2] * p_2);
    };

    auto scale_argus_diag
        = [&templatePdf, nEntriesDiag, binWidthDiag, xsDiag](double *x, double *par)
    {
        double xx = x[0];
        double p_3 = templatePdf.EvalArgus(xx, xsDiag[17], xsDiag[18], xsDiag[19]);
        return nEntriesDiag * binWidthDiag * ((1.0 - xsDiag[0]) * xsDiag[3] * p_3);
    };

    // --- COSTRUZIONE DEI TF1 ---
    TF1 *f_draw_diag = new TF1("f_draw_diag", scale_pdf_diag, xMin, xMax, 0);
    f_draw_diag->SetNpx(5000);
    f_draw_diag->SetLineColor(kBlue);
    f_draw_diag->SetLineWidth(3);
    f_draw_diag->Draw("SAME");

    TF1 *f_sig_diag = new TF1("f_sig_diag", scale_sig_diag, xMin, xMax, 0);
    f_sig_diag->SetNpx(5000);
    f_sig_diag->SetLineColor(kRed);
    f_sig_diag->SetLineStyle(2);
    f_sig_diag->SetLineWidth(2);
    f_sig_diag->Draw("SAME");

    TF1 *f_comb_diag = new TF1("f_comb_diag", scale_comb_diag, xMin, xMax, 0);
    f_comb_diag->SetNpx(5000);
    f_comb_diag->SetLineColor(kGray + 2);
    f_comb_diag->SetLineStyle(3);
    f_comb_diag->SetLineWidth(2);
    f_comb_diag->Draw("SAME");

    TF1 *f_p1_diag = new TF1("f_p1_diag", scale_p1_diag, xMin, xMax, 0);
    f_p1_diag->SetNpx(5000);
    f_p1_diag->SetLineColor(kOrange + 1);
    f_p1_diag->SetLineStyle(7);
    f_p1_diag->SetLineWidth(2);
    f_p1_diag->Draw("SAME");

    TF1 *f_p2_diag = new TF1("f_p2_diag", scale_p2_diag, xMin, xMax, 0);
    f_p2_diag->SetNpx(5000);
    f_p2_diag->SetLineColor(kMagenta);
    f_p2_diag->SetLineStyle(7);
    f_p2_diag->SetLineWidth(2);
    f_p2_diag->Draw("SAME");

    TF1 *f_arg_diag = new TF1("f_arg_diag", scale_argus_diag, xMin, xMax, 0);
    f_arg_diag->SetNpx(5000);
    f_arg_diag->SetLineColor(kCyan + 1);
    f_arg_diag->SetLineStyle(5);
    f_arg_diag->SetLineWidth(2);
    f_arg_diag->Draw("SAME");

    // Legenda descrittiva
    TLegend *legDiag = new TLegend(0.55, 0.50, 0.95, 0.88);
    legDiag->SetBorderSize(0);
    legDiag->SetFillStyle(0);
    legDiag->SetTextFont(42);
    legDiag->SetTextSize(0.028);
    legDiag->AddEntry(h_single_toy, "Toy Pseudo-Data", "ep");
    legDiag->AddEntry(f_draw_diag, "Total Unbinned Fit", "l");
    legDiag->AddEntry(f_sig_diag, "Signal Component", "l");
    legDiag->AddEntry(f_p1_diag, "D^{+} Bkg Component", "l");
    legDiag->AddEntry(f_p2_diag, "D_{s}^{+} Bkg Component", "l");
    legDiag->AddEntry(f_arg_diag, "Argus Bkg Component", "l");
    legDiag->AddEntry(f_comb_diag, "Combinatorial Bkg", "l");
    legDiag->Draw("SAME");

    cSingleToy->Update();
    cSingleToy->SaveAs("./_fig/ToyMC_SingleFit_Diagnostic.pdf");
    cSingleToy->SaveAs("./_root/ToyMC_SingleFit_Diagnostic.root");

    // Pulizia
    delete minDiag;
    delete cSingleToy;
    delete f_draw_diag;
    delete f_sig_diag;
    delete f_comb_diag;
    delete f_p1_diag;
    delete f_p2_diag;
    delete f_arg_diag;
    delete h_single_toy;
    g_pdf_unbinned = nullptr;

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    std::cout << "Toy MC completed in " << elapsed.count() << " seconds." << std::endl;

    // --- INTERRUTTORE DI RITORNO ---
    if(useBR)
    {
        double mean_br_fit = 0.0;
        double sigma_br_fit = 0.0;
        TFitResultPtr fitResBr = h_fit_br->Fit("gaus", "S L Q I N");
        if(fitResBr.Get() && fitResBr->IsValid())
        {
            mean_br_fit = fitResBr->Parameter(1);
            sigma_br_fit = fitResBr->Parameter(2);
        }
        else
        {
            mean_br_fit = h_fit_br->GetMean();
            sigma_br_fit = h_fit_br->GetRMS();
        }
        cout << Form("[TOY MC RESULT - BR mode] Mean Shift = %.5e, Sigma (Res) = %.5e", mean_br_fit,
            sigma_br_fit)
             << endl;
        return std::make_pair(mean_br_fit, sigma_br_fit);
    }
    else
    {
        double mean_fs_fit = 0.0;
        double sigma_fs_fit = 0.0;
        TFitResultPtr fitResFs = h_fit_fs->Fit("gaus", "S L Q I N");
        if(fitResFs.Get() && fitResFs->IsValid())
        {
            mean_fs_fit = fitResFs->Parameter(1);
            sigma_fs_fit = fitResFs->Parameter(2);
        }
        else
        {
            mean_fs_fit = h_fit_fs->GetMean();
            sigma_fs_fit = h_fit_fs->GetRMS();
        }
        cout << Form("[TOY MC RESULT - fs mode] Mean Shift = %.5e, Sigma (Res) = %.5e", mean_fs_fit,
            sigma_fs_fit)
             << endl;
        return std::make_pair(mean_fs_fit, sigma_fs_fit);
    }
}

// ================================================================================
// Belt Construction
// ================================================================================
// Struttura d'appoggio per ordinare l'array in base al Likelihood Ratio
struct FCPoint
{
    double x; // Valore misurato
    double prob; // Probabilità P(x|mu) * dx
    double R; // Likelihood Ratio

    // Operatore per ordinare in senso DECRESCENTE rispetto a R
    bool operator<(const FCPoint &other) const
    {
        return R > other.R;
    }
};

void analysis::ConstructBelt(double sigma0, double alpha, double max_val, bool useBR)
{
    SetLBStyle();

    // Fissiamo univocamente il Confidence Level al 90% come richiesto dalle istruzioni
    constexpr double target_CL = 0.90;

    // Range dinamico basato sul massimo valore esplorato nei Toy
    double mu_start = 0.0;
    double mu_end = max_val;
    int n_steps = 200; // Numero di punti per rendere la curva liscia
    double mu_step = (mu_end - mu_start) / n_steps;

    std::vector<double> vec_mu;
    std::vector<double> vec_x_lower;
    std::vector<double> vec_x_upper;

    cout << "\n--> Constructing Feldman-Cousins Belt Boundaries (" << (target_CL * 100)
         << "% CL)..." << endl;

    for(double mu = mu_start; mu <= mu_end; mu += mu_step)
    {
        double sigma_mu = sigma0 + alpha * mu;

        std::vector<FCPoint> points;
        double dx = sigma_mu / 200.0;
        double x_start = mu - 5.0 * sigma_mu;
        double x_end = mu + 5.0 * sigma_mu;

        for(double x = x_start; x <= x_end; x += dx)
        {
            double par[3] = { mu, sigma0, alpha }; // Ora passi anche alpha e sigma0!
            double prob_density = TMath::Gaus(x, mu, sigma_mu, kTRUE);
            double R = LROrdering(&x, par);

            points.push_back({ x, prob_density * dx, R });
        }
        /*
            double sigma = sigma0 + alpha * mu;

            std::vector<FCPoint> points;
            double dx = sigma / 200.0;
            double x_start = mu - 5.0 * sigma;
            double x_end = mu + 5.0 * sigma;

            for(double x = x_start; x <= x_end; x += dx)
            {
            double par[2] = { mu, sigma };
            double prob_density = TMath::Gaus(x, mu, sigma, true);
            double R = LROrdering(&x, par);

            points.push_back({ x, prob_density * dx, R });
            }
        */
        std::sort(points.begin(), points.end());

        double sum_prob = 0.0;
        double x_min = 999.0;
        double x_max = -999.0;

        for(const auto &pt : points)
        {
            sum_prob += pt.prob;
            if(pt.x < x_min)
                x_min = pt.x;
            if(pt.x > x_max)
                x_max = pt.x;
            if(sum_prob >= target_CL)
                break;
        }

        vec_mu.push_back(mu);
        vec_x_lower.push_back(x_min);
        vec_x_upper.push_back(x_max);
    }

    // --- DISEGNO DELLA BANDA ---
    int nPoints = vec_mu.size();
    std::vector<double> x_closed;
    std::vector<double> mu_closed;
    x_closed.reserve(2 * nPoints);
    mu_closed.reserve(2 * nPoints);

    for(int i = 0; i < nPoints; ++i)
    {
        x_closed.push_back(vec_x_lower[i]);
        mu_closed.push_back(vec_mu[i]);
    }
    for(int i = nPoints - 1; i >= 0; --i)
    {
        x_closed.push_back(vec_x_upper[i]);
        mu_closed.push_back(vec_mu[i]);
    }

    TGraph *g_belt_filled = new TGraph(x_closed.size(), &x_closed[0], &mu_closed[0]);
    if(useBR)
        g_belt_filled->SetFillColorAlpha(kGreen - 9, 0.35);
    else
        g_belt_filled->SetFillColorAlpha(kBlue - 9, 0.35);
    g_belt_filled->SetLineWidth(0);

    TGraph *g_upperEdge = new TGraph(nPoints, &vec_x_lower[0], &vec_mu[0]);
    TGraph *g_lowerEdge = new TGraph(nPoints, &vec_x_upper[0], &vec_mu[0]);
    if(useBR)
    {
        g_upperEdge->SetLineColor(kGreen + 1);
        g_lowerEdge->SetLineColor(kGreen + 1);
    }
    else
    {
        g_upperEdge->SetLineColor(kBlue + 1);
        g_lowerEdge->SetLineColor(kBlue + 1);
    }
    g_upperEdge->SetLineWidth(3);
    g_lowerEdge->SetLineWidth(3);

    TString canvasName = useBR ? "cBeltBR" : "cBelt_fs";
    TCanvas *cBelt = new TCanvas(canvasName, "Feldman-Cousins Confidence Belt", 800, 800);
    cBelt->cd();
    cBelt->SetGrid();

    double plot_x_min = -3.0 * sigma0;
    double plot_x_max = mu_end + 3.0 * sigma0;

    TString titleX = useBR ? "Measured #hat{#font[12]{B}}(#tau^{+}#rightarrow#phi#mu^{+})"
                           : "Measured #hat{f}_{s}";
    TString titleY = useBR ? "True #font[12]{B}(#tau^{+}#rightarrow#phi#mu^{+})" : "True f_{s}";
    TString frameTitle
        = Form("Feldman-Cousins Confidence Belt (90%% CL);%s;%s", titleX.Data(), titleY.Data());

    TH2F *hFrame = new TH2F("hFrame", frameTitle, 100, plot_x_min, plot_x_max, 100, 0.0, mu_end);
    hFrame->SetStats(0);
    hFrame->Draw();

    g_belt_filled->Draw("F SAME");
    g_upperEdge->Draw("L SAME");
    g_lowerEdge->Draw("L SAME");

    TLine *diag = new TLine(0.0, 0.0, mu_end, mu_end);
    diag->SetLineStyle(2);
    diag->SetLineColor(kGray + 2);
    diag->Draw("SAME");

    TLine *vert_zero = new TLine(0.0, 0.0, 0.0, mu_end);
    vert_zero->SetLineStyle(3);
    vert_zero->SetLineColor(kBlack);
    vert_zero->Draw("SAME");

    cBelt->Update();

    TString saveName_fig
        = useBR ? "./_fig/FeldmanCousinsBelt_BR_Filled" : "./_fig/FeldmanCousinsBelt_fs_Filled";
    cBelt->SaveAs(saveName_fig + ".pdf");
    TString saveName_root
        = useBR ? "./_root/FeldmanCousinsBelt_BR_Filled" : "./_root/FeldmanCousinsBelt_fs_Filled";
    cBelt->SaveAs(saveName_root + ".root");

    cout << "Filled Belt (90% CL) successfully generated and saved!" << endl;
}

void analysis::RunFeldmanCousinsPipeline(int nToysPerPoint, bool useBR)
{
    cout << "\n=======================================================" << endl;
    cout << Form("   STARTING AUTOMATED FELDMAN-COUSINS PIPELINE (%s mode)", useBR ? "BR" : "fs")
         << endl;
    cout << "=======================================================" << endl;

    // Fattori nominali costanti
    double r_eff = 0.2212;
    const double br_taunu_nom = 5.39e-2;
    const double br_phimunu_nom = 2.24e-2;
    double k_factor_nom = r_eff * (br_phimunu_nom / br_taunu_nom);
    double true_f3 = 0.625;

    // Definiamo i punti nominali su f_s (generiamo SEMPRE a partire da qui)
    std::vector<double> fs_points = { 0.0, 0.0025, 0.0050, 0.0075, 0.0100 };

    std::vector<double> x_true;
    std::vector<double> y_sigma;
    std::vector<double> y_bias; // Salveremo il bias (mean shift) di ciascun punto

    for(double fs : fs_points)
    {
        // 1. Calcolo esatto del BR equivalente
        double true_br = (fs / ((1.0 - fs) * true_f3)) * k_factor_nom;

        if(useBR)
            cout << Form("\n--> Running Toy: true_fs = %.4f => true_BR = %.2e", fs, true_br)
                 << endl;
        else
            cout << Form("\n--> Running Toy: true_fs = %.4f", fs) << endl;

        // 2. Lancio del Toy (restituisce la coppia <mean_shift, sigma>)
        auto results = RunToyMC(nToysPerPoint, fs, useBR);
        double bias_misurato = results.first;
        double sigma_misurata = results.second;

        // 3. Salviamo le coordinate per i fit e i controlli
        if(useBR)
        {
            x_true.push_back(true_br);
        }
        else
        {
            x_true.push_back(fs);
        }
        y_sigma.push_back(sigma_misurata);
        y_bias.push_back(bias_misurato);
    }

    // Valore massimo dell'asse X per disegnare correttamente la banda
    double max_x_val = x_true.back();

    // --- FIT DELLA RISOLUZIONE (USANDO LE VARIABILI RISCALATE PER MINUIT) ---
    double scale = useBR ? 1e7 : 1e3;
    std::vector<double> x_scaled(x_true.size());
    std::vector<double> y_scaled(y_sigma.size());

    for(size_t i = 0; i < x_true.size(); ++i)
    {
        x_scaled[i] = x_true[i] * scale;
        y_scaled[i] = y_sigma[i] * scale;
    }

    TGraph *g_res = new TGraph(x_scaled.size(), &x_scaled[0], &y_scaled[0]);
    TF1 *f_linear = new TF1("f_linear", "[0] + [1]*x", 0.0, max_x_val * scale * 1.2);

    f_linear->SetParameters(y_scaled[0], 0.0);
    cout << "\n--> Fitting resolution dependency..." << endl;
    g_res->Fit(f_linear, "Q");

    double sigma0 = f_linear->GetParameter(0) / scale;
    double alpha = f_linear->GetParameter(1);

    cout << "=======================================================" << endl;
    cout << Form("   FIT RESULTS FOR RESOLUTION FUNCTION sigma(%s):", useBR ? "BR" : "fs") << endl;
    cout << Form("   sigma(X) = %.5e + %.5f * X", sigma0, alpha) << endl;
    cout << "=======================================================" << endl;

    // =========================================================================
    // SANITY CHECK 1: COMPORTAMENTO DELLA RISOLUZIONE (PLOT FISICO)
    // =========================================================================
    TCanvas *cResCheck = new TCanvas("cResCheck", "Resolution Sanity Check", 800, 600);
    cResCheck->SetGrid();

    TGraphErrors *g_res_physical = new TGraphErrors(x_true.size());
    for(size_t i = 0; i < x_true.size(); ++i)
    {
        g_res_physical->SetPoint(i, x_true[i], y_sigma[i]);
        // Incertezza statistica della sigma stimata dai Toy: sigma / sqrt(2*N)
        double err_y = y_sigma[i] / std::sqrt(2.0 * nToysPerPoint);
        g_res_physical->SetPointError(i, 0.0, err_y);
    }

    g_res_physical->SetMarkerStyle(20);
    g_res_physical->SetMarkerSize(1.2);
    g_res_physical->SetMarkerColor(kBlue + 1);
    g_res_physical->SetLineColor(kBlue + 1);

    TString titleX = useBR ? "True B(#tau^{+}#rightarrow#phi#mu^{+})" : "True f_{s}";
    TString titleY = useBR ? "#sigma(B)" : "#sigma(f_{s})";
    g_res_physical->SetTitle(Form("Resolution Fit;%s;%s", titleX.Data(), titleY.Data()));
    g_res_physical->GetXaxis()->SetMaxDigits(3);
    g_res_physical->GetYaxis()->SetMaxDigits(3);
    g_res_physical->Draw("AP");

    // Sovrapponiamo la funzione lineare con i parametri de-scalati corretti
    TF1 *f_phys = new TF1("f_phys", "[0] + [1]*x", 0.0, max_x_val * 1.2);
    f_phys->SetParameters(sigma0, alpha);
    f_phys->SetLineColor(kRed);
    f_phys->SetLineWidth(3);
    f_phys->Draw("SAME");

    TPaveText *paveRes = new TPaveText(0.15, 0.72, 0.45, 0.88, "NDC");
    paveRes->SetBorderSize(1);
    paveRes->SetFillColor(kWhite);
    paveRes->SetTextFont(42);
    paveRes->SetTextSize(0.03);
    paveRes->AddText(Form("#sigma_{0} = %.3e", sigma0));
    paveRes->AddText(Form("#alpha = %.4f", alpha));
    paveRes->Draw();

    cResCheck->Update();
    cResCheck->SaveAs(
        useBR ? "./_fig/SanityCheck_Resolution_BR.pdf" : "./_fig/SanityCheck_Resolution_fs.pdf");

    // =========================================================================
    // SANITY CHECK 2: COMPORTAMENTO DEL BIAS (COERENTE CON ZERO)
    // =========================================================================
    TCanvas *cBiasCheck = new TCanvas("cBiasCheck", "Bias Sanity Check", 800, 600);
    cBiasCheck->SetGrid();

    TGraphErrors *g_bias_physical = new TGraphErrors(x_true.size());
    for(size_t i = 0; i < x_true.size(); ++i)
    {
        g_bias_physical->SetPoint(i, x_true[i], y_bias[i]);
        // Errore statistico del valor medio (errore sulla media dei Toy): sigma_toy / sqrt(N_toy)
        double err_y = y_sigma[i] / std::sqrt(nToysPerPoint);
        g_bias_physical->SetPointError(i, 0.0, err_y);
    }

    g_bias_physical->SetMarkerStyle(21);
    g_bias_physical->SetMarkerSize(1.2);
    g_bias_physical->SetMarkerColor(kBlack);
    g_bias_physical->SetLineColor(kBlack);

    TString titleY_bias = useBR ? "Bias: Fitted B - True B" : "Bias: Fitted f_{s} - True f_{s}";
    g_bias_physical->SetTitle(Form("Bias Check;%s;%s", titleX.Data(), titleY_bias.Data()));
    g_bias_physical->GetXaxis()->SetMaxDigits(3);
    g_bias_physical->GetYaxis()->SetMaxDigits(3);

    // Impostiamo l'asse Y simmetrico intorno a 0 per una visualizzazione ottimale del bias
    double max_err = 0.0;
    for(size_t i = 0; i < x_true.size(); ++i)
    {
        double val = std::abs(y_bias[i]) + 3.0 * (y_sigma[i] / std::sqrt(nToysPerPoint));
        if(val > max_err)
            max_err = val;
    }
    g_bias_physical->SetMinimum(-max_err);
    g_bias_physical->SetMaximum(max_err);
    g_bias_physical->Draw("AP");

    // Disegnamo la linea dello zero teorico (nessun bias)
    TLine *lZero = new TLine(0.0, 0.0, max_x_val * 1.2, 0.0);
    lZero->SetLineStyle(2);
    lZero->SetLineColor(kRed);
    lZero->SetLineWidth(2);
    lZero->Draw("SAME");

    cBiasCheck->Update();
    cBiasCheck->SaveAs(useBR ? "./_fig/SanityCheck_Bias_BR.pdf" : "./_fig/SanityCheck_Bias_fs.pdf");

    // 4. Costruzione geometrica finale della Belt
    ConstructBelt(sigma0, alpha, max_x_val, useBR);

    delete g_res;
    delete f_linear;
    delete g_res_physical;
    delete f_phys;
    delete cResCheck;
    delete g_bias_physical;
    delete lZero;
    delete cBiasCheck;
}

// ================================================================================
// Wilks Validation
// ================================================================================
struct WilksResult
{
    bool converged = false;
    double delta_F = 0.0; // Questo conterrà il valore di t0 = F_0 - F_min
    double fs_fitted = 0.0;
};

// Funzione worker isolata per calcolare la differenza di Likelihood su singolo Toy
WilksResult RunSingleWilksToy(int toyId, int nEvents, const FullUnbinnedPDF &templatePdf,
    const std::vector<double> &gen_pars, bool usePol1Bkg, double xMin, double xMax)
{
    WilksResult res;

    // 1. Generatore di numeri casuali locale (true_fs = 0.0 per l'ipotesi nulla)
    TRandom3 threadRandom(THE_SEED * (toyId + 1));

    FullUnbinnedPDF localPdf = templatePdf;
    g_pdf_unbinned = &localPdf;

    std::vector<double> null_gen_pars = gen_pars;
    null_gen_pars[0] = 0.0; // Generazione rigorosamente sotto H0 (f_s = 0)

    auto gen_lambda = [&localPdf, &null_gen_pars](double *x, double *par)
    {
        double xx = x[0];
        return localPdf(&xx, null_gen_pars.data());
    };
    TF1 fGen(Form("fGenWilks_%d", toyId), gen_lambda, xMin, xMax, 0);
    fGen.SetNpx(10000);

    g_data_events.clear();
    g_data_events.reserve(nEvents);
    for(int ev = 0; ev < nEvents; ++ev)
    {
        g_data_events.push_back(fGen.GetRandom(&threadRandom));
    }

    // =========================================================================
    // FIT 1: IPOTESI ALTERNATIVA (f_s COMPLETAMENTE LIBERO - UNCONSTRAINED)
    // =========================================================================
    ROOT::Math::Minimizer *minAlt = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
    minAlt->SetMaxFunctionCalls(50000);
    minAlt->SetTolerance(0.01);
    minAlt->SetPrintLevel(-1);

    ROOT::Math::Functor fNLL(&Unbinned2NLL, 20);
    minAlt->SetFunction(fNLL);

    // Definiamo f_s senza limiti: può fluttuare liberamente nel negativo
    minAlt->SetVariable(0, "f_s", 0.0, 0.005);
    // minAlt->SetVariableLimits(0, 0.0, 1.0); // <-- RIMOZIONE DEL LIMITE!

    minAlt->SetVariable(1, "f_1", gen_pars[1], 0.005);
    minAlt->SetVariable(2, "f_2", gen_pars[2], 0.005);
    minAlt->SetVariable(3, "f_3", gen_pars[3], 0.005);

    if(usePol1Bkg)
        minAlt->SetVariable(4, "pol1_slope", gen_pars[4], 0.05);
    else
        minAlt->SetVariable(4, "expo_slope", gen_pars[4], 0.005);

    for(int p = 5; p < 20; ++p)
    {
        minAlt->SetVariable(p, Form("p_%d", p), gen_pars[p], 0.005);
        minAlt->FixVariable(p);
    }

    minAlt->Minimize();

    if(minAlt->Status() != 0)
    {
        delete minAlt;
        g_pdf_unbinned = nullptr;
        return res;
    }

    double F_min = minAlt->MinValue();
    double fs_fitted = minAlt->X()[0];

    // =========================================================================
    // FIT 2: IPOTESI NULLA (f_s FISSO A ZERO)
    // =========================================================================
    ROOT::Math::Minimizer *minNull = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
    minNull->SetMaxFunctionCalls(50000);
    minNull->SetTolerance(0.01);
    minNull->SetPrintLevel(-1);
    minNull->SetFunction(fNLL);

    minNull->SetVariable(0, "f_s", 0.0, 0.005);
    minNull->FixVariable(0);

    minNull->SetVariable(1, "f_1", gen_pars[1], 0.005);
    minNull->SetVariable(2, "f_2", gen_pars[2], 0.005);
    minNull->SetVariable(3, "f_3", gen_pars[3], 0.005);

    if(usePol1Bkg)
        minNull->SetVariable(4, "pol1_slope", gen_pars[4], 0.05);
    else
        minNull->SetVariable(4, "expo_slope", gen_pars[4], 0.005);

    for(int p = 5; p < 20; ++p)
    {
        minNull->SetVariable(p, Form("p_%d", p), gen_pars[p], 0.005);
        minNull->FixVariable(p);
    }

    minNull->Minimize();

    if(minNull->Status() != 0)
    {
        delete minAlt;
        delete minNull;
        g_pdf_unbinned = nullptr;
        return res;
    }

    double F_0 = minNull->MinValue();

    res.converged = true;
    res.fs_fitted = fs_fitted;
    res.delta_F = F_0 - F_min;

    if(res.delta_F < 0.0)
        res.delta_F = 0.0;

    delete minAlt;
    delete minNull;
    g_pdf_unbinned = nullptr;
    return res;
}

void analysis::VerifyWilksTheorem(int nToys)
{
    auto start = std::chrono::high_resolution_clock::now();
    SetLBStyle();
    constexpr bool usePol1Bkg = false;

    // --- SALVAVITA MULTI-THREADING PER ROOT ---
    TF1::DefaultAddToGlobalList(kFALSE);

    // --- [STEP 1] Fit ausiliari per le forme ---
    cout << "\n=== [WILKS CHECK] Running Auxiliary Fits ===" << endl;
    LoadDataset(1);
    AuxFitResult res_sig = FitTemplateMass(44);
    AuxFitResult res_p1 = FitTemplateMass(34);
    AuxFitResult res_p2 = FitTemplateMass(41);
    AuxFitResult res_arg = FitTemplateMass(42);

    if(!res_sig.isValid || !res_p1.isValid || !res_p2.isValid || !res_arg.isValid)
    {
        cerr << "[ERROR] Auxiliary fits failed! Aborting Wilks check." << endl;
        TF1::DefaultAddToGlobalList(kTRUE);
        return;
    }

    std::vector<double> gen_pars = { 0.0, 0.026, 0.051, 0.625, -1.145, res_sig.params[1],
        res_sig.params[2], res_sig.params[3], res_sig.params[4], res_p1.params[1], res_p1.params[2],
        res_p1.params[3], res_p1.params[4], res_p2.params[1], res_p2.params[2], res_p2.params[3],
        res_p2.params[4], res_arg.params[1], res_arg.params[2], res_arg.params[3] };

    Double_t xMin = 1.65;
    Double_t xMax = 2.09;
    FullUnbinnedPDF templatePdf(xMin, xMax, res_sig, res_p1, res_p2, res_arg, usePol1Bkg);

    // Conteggio eventi totali dai dati reali
    LoadDataset(0);
    int nEvents = 0;
    for(Long64_t jentry = 0; jentry < fChain->GetEntriesFast(); jentry++)
    {
        if(LoadTree(jentry) < 0)
            break;
        fChain->GetEntry(jentry);
        if(id == 0 && D_M >= xMin && D_M <= xMax)
            nEvents++;
    }
    if(nEvents == 0)
        nEvents = 10000;

    int nCores = static_cast<int>(std::thread::hardware_concurrency());
    if(nCores == 0)
        nCores = 4;
    cout << "[INFO] Launching " << nToys << " Unconstrained Wilks Toys on " << nCores << " cores."
         << endl;

    std::vector<WilksResult> wilksResults;
    wilksResults.reserve(nToys);

    for(int i = 0; i < nToys; i += nCores)
    {
        std::vector<std::future<WilksResult>> futures;
        for(int t = 0; (t < nCores) && (i + t) < nToys; ++t)
        {
            int toyId = i + t;
            futures.push_back(std::async(std::launch::async, RunSingleWilksToy, toyId, nEvents,
                std::ref(templatePdf), std::cref(gen_pars), usePol1Bkg, xMin, xMax));
        }
        for(auto &f : futures)
        {
            wilksResults.push_back(f.get());
        }
        cout << "\r  Completed Wilks toys: " << wilksResults.size() << " / " << nToys << "..."
             << flush;
    }

    // --- [STEP 2] Booking Istogramma per Delta_F ---
    TH1D *h_delta_F = new TH1D(
        "h_delta_F", "Wilks Theorem Check;#Delta F = F_{0} - F_{min};Toys", 50, 0.0, 10.0);
    h_delta_F->SetDirectory(nullptr);

    int convergedToys = 0;
    for(const auto &res : wilksResults)
    {
        if(res.converged)
        {
            convergedToys++;
            h_delta_F->Fill(res.delta_F);
        }
    }

    cout << "\n[WILKS RESULTS] Converged: " << convergedToys << " / " << nToys << endl;

    // --- [STEP 3] Plotting e Overplot Teorico con Residui (Pulls) ---
    TCanvas *cWilks = new TCanvas("cWilks", "Wilks Theorem Validation with Pulls", 900, 900);
    double splitPoint = 0.30;

    // Pad Superiore (Plot Principale)
    TPad *pad1 = new TPad("pad1", "Main Pad", 0.0, splitPoint, 1.0, 1.0);
    pad1->SetBottomMargin(0.02); // Tocca il pad inferiore senza spazio bianco
    pad1->Draw();
    pad1->cd();
    pad1->SetGrid();

    h_delta_F->SetMinimum(0.0);
    h_delta_F->GetXaxis()->SetLabelSize(0); // Nascondiamo l'asse X superiore
    h_delta_F->GetXaxis()->SetTitleSize(0);
    h_delta_F->SetMarkerStyle(20);
    h_delta_F->SetMarkerSize(1.0);
    h_delta_F->Draw("E P");

    // 1. Definisci la PDF pura (SENZA la larghezza del bin nella normalizzazione)
    // Nota: parto da > 0 per evitare il div by zero nei calcoli interni di ROOT,
    // ma l'integrale ROOT lo gestirà benissimo
    TF1 *f_chi2_theory
        = new TF1("f_chi2_theory", "[0] * exp(-0.5*x) / sqrt(2.0 * TMath::Pi() * x)", 1e-6, 10.0);
    f_chi2_theory->SetNpx(1000);

    // Normalizzazione = solo numero totale di toy
    f_chi2_theory->SetParameter(0, convergedToys);

    // 2. Per disegnare la curva continua riscalata ai bin (solo per la vista)
    TF1 *f_chi2_draw = (TF1 *)f_chi2_theory->Clone("f_chi2_draw");
    f_chi2_draw->SetParameter(0, convergedToys * h_delta_F->GetBinWidth(1));
    f_chi2_draw->SetLineColor(kRed);
    f_chi2_draw->SetLineWidth(3);
    f_chi2_draw->Draw("SAME");

    TLegend *leg = new TLegend(0.45, 0.65, 0.88, 0.85);
    leg->SetBorderSize(1);
    leg->SetFillColor(kWhite);
    leg->SetTextFont(42);
    leg->SetTextSize(0.03);
    leg->AddEntry(h_delta_F, Form("Toy MC under H_{0} (%d toys)", convergedToys), "ep");
    leg->AddEntry(f_chi2_theory, "Pure Wilks Theory: #chi^{2}_{1}(#Delta F)", "l");
    leg->Draw("SAME");

    // Pad Inferiore (Plot dei Residui / Pulls)
    cWilks->cd();
    TPad *pad2 = new TPad("pad2", "Pull Pad", 0.0, 0.0, 1.0, splitPoint);
    pad2->SetTopMargin(0.02);
    pad2->SetBottomMargin(0.35); // Spazio per i titoli dell'asse X
    pad2->SetGridy();
    pad2->Draw();
    pad2->cd();

    TH1D *hPull = (TH1D *)h_delta_F->Clone("hPull_wilks");
    hPull->Reset();
    for(int i = 1; i <= h_delta_F->GetNbinsX(); ++i)
    {
        double obs = h_delta_F->GetBinContent(i);
        double err = h_delta_F->GetBinError(i);
        double binLow = h_delta_F->GetXaxis()->GetBinLowEdge(i);
        double binUp = h_delta_F->GetXaxis()->GetBinUpEdge(i);

        // 1. Calcolo analitico esatto della frazione tramite la CDF (senza integrazione numerica)
        double fraction = TMath::Erf(std::sqrt(binUp / 2.0)) - TMath::Erf(std::sqrt(binLow / 2.0));
        double exp = convergedToys * fraction;

        if(err > 0.0)
        {
            double pull = (obs - exp) / err;
            hPull->SetBinContent(i, pull);
            hPull->SetBinError(i, 0.0); // Nessuna barra d'errore sul punto del pull
        }
        else
        {
            hPull->SetBinContent(i, -999.0);
        }
    }

    hPull->GetYaxis()->SetTitle("Pull");
    hPull->GetYaxis()->SetTitleSize(gStyle->GetTitleSize("Y") * 0.8);
    hPull->GetYaxis()->SetLabelSize(gStyle->GetLabelSize("Y") * 0.8);
    hPull->GetYaxis()->SetTitleOffset(gStyle->GetTitleOffset("Y") * 1.2);
    hPull->GetYaxis()->SetRangeUser(-5.0, 5.0); // I pull oscillano tipicamente tra -3 e 3

    hPull->GetXaxis()->SetTitle("#Delta F = F_{0} - F_{min}");
    hPull->GetXaxis()->SetTitleSize(gStyle->GetTitleSize("X") * 0.8);
    hPull->GetXaxis()->SetLabelSize(gStyle->GetLabelSize("X") * 0.8);
    hPull->GetXaxis()->SetTitleOffset(gStyle->GetTitleOffset("X") * 1.1);

    hPull->SetMarkerStyle(20);
    hPull->SetMarkerSize(0.8);
    hPull->SetMarkerColor(kBlack);
    hPull->SetLineColor(kBlack);
    hPull->Draw("P");

    // Linea di riferimento a Zero
    TLine *line0 = new TLine(0.0, 0.0, 10.0, 0.0);
    line0->SetLineColor(kRed);
    line0->SetLineWidth(2);
    line0->Draw("SAME");

    cWilks->Update();
    cWilks->SaveAs("./_fig/SanityCheck_Wilks_Pure.pdf");
    cWilks->SaveAs("./_root/SanityCheck_Wilks_Pure.root");

    while(gROOT->GetListOfCanvases()->FindObject("cWilks"))
    {
        gSystem->ProcessEvents(); // Gestisce i movimenti del mouse, zoom, click, ecc.
        gSystem->Sleep(50); // Dorme 50 millisecondi per non sovraccaricare la CPU
    }

    // =========================================================================
    // CALCOLO DELLA CDF ED ESECUZIONE DEL TEST DI KOLMOGOROV-SMIRNOV (UNBINNED)
    // =========================================================================

    // 1. Raccogliamo i Delta F dei toy convertiti e ordiniamoli
    std::vector<double> toy_deltas;
    toy_deltas.reserve(convergedToys);
    for(const auto &res : wilksResults)
    {
        if(res.converged)
        {
            toy_deltas.push_back(res.delta_F);
        }
    }
    std::sort(toy_deltas.begin(), toy_deltas.end());

    // 2. Calcolo manuale e rigoroso della statistica KS (D_max)
    double max_D = 0.0;
    int N = toy_deltas.size();

    TGraph *g_cdf_emp = new TGraph(N); // Grafico per la CDF sperimentale (dei toy)

    for(int i = 0; i < N; ++i)
    {
        double x = toy_deltas[i];
        double f_emp = (double)(i + 1) / N; // CDF empirica: frazione di toy <= x
        double f_theo = TMath::Erf(std::sqrt(x / 2.0)); // CDF teorica del Chi2(1)

        g_cdf_emp->SetPoint(i, x, f_emp);

        double diff = std::abs(f_emp - f_theo);
        if(diff > max_D)
        {
            max_D = diff;
        }
    }

    // Calcoliamo il p-value di Kolmogorov-Smirnov usando la libreria di ROOT
    double ks_p_value = TMath::KolmogorovProb(max_D * std::sqrt(N));

    std::cout << "\n=======================================================" << endl;
    std::cout << "   KOLMOGOROV-SMIRNOV TEST RESULTS (UNBINNED)" << endl;
    std::cout << "=======================================================" << endl;
    std::cout << Form("  Number of Toys (N):  %d", N) << endl;
    std::cout << Form("  KS Distance (d_max): %.4f", max_D) << endl;
    std::cout << Form("  KS p-value:          %.4f (%.1f%%)", ks_p_value, ks_p_value * 100.0)
              << endl;
    std::cout << "=======================================================\n" << endl;

    // 3. Creazione del Canvas per il Plot della CDF
    TCanvas *cCDF = new TCanvas("cCDF", "Cumulative Distribution Function & KS Test", 800, 600);
    cCDF->cd();
    cCDF->SetGrid();

    // Creiamo un frame per gli assi (asse Y va rigorosamente da 0 a 1 per una probabilità)
    TH2F *hFrameCDF = new TH2F("hFrameCDF",
        "Wilks Theorem Validation (CDF);#Delta F = F_{0} - F_{min};Cumulative Probability", 100,
        0.0, 10.0, 100, 0.0, 1.05);
    hFrameCDF->SetStats(0);
    hFrameCDF->Draw();

    // Disegniamo la CDF teorica (Chi2 a 1 grado di libertà)
    TF1 *f_cdf_theo = new TF1("f_cdf_theo", "TMath::Erf(std::sqrt(x/2.0))", 0.0, 10.0);
    f_cdf_theo->SetNpx(1000);
    f_cdf_theo->SetLineColor(kRed);
    f_cdf_theo->SetLineWidth(3);
    f_cdf_theo->Draw("SAME");

    // Disegniamo la CDF empirica dei tuoi Toy
    g_cdf_emp->SetMarkerStyle(20);
    g_cdf_emp->SetMarkerSize(0.6);
    g_cdf_emp->SetMarkerColor(kBlack);
    g_cdf_emp->SetLineColor(kBlack);
    g_cdf_emp->SetLineWidth(1);
    g_cdf_emp->Draw("P SAME");

    // Legenda con font e dimensione forzati a valori relativi sicuri
    TLegend *legCDF = new TLegend(0.15, 0.70, 0.55, 0.85);
    legCDF->SetBorderSize(1);
    legCDF->SetFillColor(kWhite);
    legCDF->SetTextFont(42); // Font 42 (dimensione relativa, non in pixel)
    legCDF->SetTextSize(0.035); // 3.5% dell'altezza della finestra (molto sicuro)
    legCDF->AddEntry(g_cdf_emp, Form("Toy MC under H_{0} (%d toys)", N), "ep");
    legCDF->AddEntry(f_cdf_theo, "Pure Wilks Theory: F_{#chi^{2}_{1}}(#Delta F)", "l");
    legCDF->Draw("SAME");

    // Box di testo KS con font e dimensione forzati
    TPaveText *paveKS = new TPaveText(0.55, 0.15, 0.88, 0.32, "NDC");
    paveKS->SetBorderSize(1);
    paveKS->SetFillColor(kWhite);
    paveKS->SetTextFont(42); // Font 42
    paveKS->SetTextSize(0.035); // Dimensione del testo 3.5%
    paveKS->SetTextAlign(12); // Allineamento a sinistra
    paveKS->AddText(Form("KS d_{max} = %.4f", max_D));
    paveKS->AddText(Form("KS p-value = %.4f", ks_p_value));
    paveKS->Draw();

    cCDF->Update();
    cCDF->SaveAs("./_fig/SanityCheck_Wilks_CDF.pdf");
    cCDF->SaveAs("./_root/SanityCheck_Wilks_CDF.root");

    // --- CICLO DI ATTESA PER TENERLO APERTO ---
    std::cout << "--> Grafico CDF pronto! Chiudi la finestra del CDF per proseguire." << std::endl;
    while(gROOT->GetListOfCanvases()->FindObject("cCDF"))
    {
        gSystem->ProcessEvents();
        gSystem->Sleep(50);
    }

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    std::cout << "Wilks validation completed in " << elapsed.count() << " seconds." << std::endl;

    // --- RIPRISTINO DELLO STATO GLOBALE DI ROOT ---
    TF1::DefaultAddToGlobalList(kTRUE);

    delete h_delta_F;
    delete f_chi2_theory;
    delete leg;
    delete hPull;
    delete line0;
    delete cWilks;
}
