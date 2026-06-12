#define analysis_cxx
#include <cmath>
#include <future>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <thread>

#include <Math/Factory.h>
#include <Math/Functor.h>
#include <Math/Minimizer.h>
#include <Math/SpecFuncMathCore.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>
#include <TH2.h>
#include <TLegend.h>
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
  public: // Membri pubblici per compatibilità con il plotting esistente
    double m_xMin, m_xMax;
    double m_sig_mean, m_sig_sigma1, m_sig_sigma2, m_sig_frac1;
    double m_p1_mean, m_p1_sigma1, m_p1_sigma2, m_p1_frac1;
    double m_p2_mean, m_p2_sigma1, m_p2_sigma2, m_p2_frac1;
    double m_argus_m0, m_argus_c, m_argus_p;
    bool m_usePol1Bkg;

  private:
    // --- VARIABILI DI CACHE ---
    mutable std::vector<double> m_last_pars; // Memorizza gli ultimi parametri valutati

    // Costanti di normalizzazione pre-calcolate
    mutable double m_norm_argus;
    mutable double m_norm_expo;
    mutable double m_norm_pol1;
    mutable double m_norm_sig_g1, m_norm_sig_g2;
    mutable double m_norm_p1_g1, m_norm_p1_g2;
    mutable double m_norm_p2_g1, m_norm_p2_g2;

    // Metodo privato per aggiornare le normalizzazioni (chiamato solo al cambio dei parametri)
    void UpdateNormalizations(const double *par) const
    {
        const double slope = par[4];

        // Normalizzazione Esponenziale
        if(std::abs(slope) < 1e-6)
        {
            m_norm_expo = m_xMax - m_xMin;
        }
        else
        {
            m_norm_expo = (std::exp(slope * m_xMax) - std::exp(slope * m_xMin)) / slope;
        }

        // Normalizzazione Pol1
        m_norm_pol1 = m_xMax - m_xMin;

        // Prefattori Gaussiane: 1 / (sigma * sqrt(2*pi))
        m_norm_sig_g1 = 1.0 / (par[6] * std::sqrt(2.0 * M_PI));
        m_norm_sig_g2 = 1.0 / (par[7] * std::sqrt(2.0 * M_PI));

        m_norm_p1_g1 = 1.0 / (par[10] * std::sqrt(2.0 * M_PI));
        m_norm_p1_g2 = 1.0 / (par[11] * std::sqrt(2.0 * M_PI));

        m_norm_p2_g1 = 1.0 / (par[14] * std::sqrt(2.0 * M_PI));
        m_norm_p2_g2 = 1.0 / (par[15] * std::sqrt(2.0 * M_PI));

        // Normalizzazione dell'Argus (Molto costosa)
        const double m0 = par[17];
        const double c = par[18];
        const double p = par[19];

        double xL = 1.0 - (m_xMin / m0) * (m_xMin / m0);
        double xH = 1.0 - (m_xMax / m0) * (m_xMax / m0);
        double gammaA = ROOT::Math::tgamma(1.0 + p);
        double dL = gammaA * ROOT::Math::inc_gamma_c(1.0 + p, -c * xL);
        double dH = gammaA * ROOT::Math::inc_gamma_c(1.0 + p, -c * xH);
        m_norm_argus = (m0 * m0) / (2.0 * c * std::pow(-c, p)) * (dL - dH);
    }

    // Controlla se i parametri del fit sono cambiati rispetto all'ultimo evento valutato
    bool ParametersChanged(const double *par) const
    {
        if(m_last_pars.size() < 20)
            return true;
        for(int i = 0; i < 20; ++i)
        {
            if(std::abs(par[i] - m_last_pars[i]) > 1e-12)
                return true;
        }
        return false;
    }

    // Metodi interni veloci che sfruttano la cache
    double Eval2G_Cached(double x, double mean, double s1, double s2, double f1, double norm_g1,
        double norm_g2) const
    {
        double dx1 = (x - mean) / s1;
        double dx2 = (x - mean) / s2;
        double g1 = std::exp(-0.5 * dx1 * dx1) * norm_g1;
        double g2 = std::exp(-0.5 * dx2 * dx2) * norm_g2;
        return f1 * g1 + (1.0 - f1) * g2;
    }

    double EvalArgus_Cached(double x, double m0, double c, double p) const
    {
        if(x >= m0 || m_norm_argus <= 0.0)
            return 0.0;
        double u = 1.0 - (x / m0) * (x / m0);
        return x * std::pow(u, p) * std::exp(c * u) / m_norm_argus;
    }

    double EvalExpo_Cached(double x, double slope) const
    {
        if(m_norm_expo <= 0.0)
            return 0.0;
        return std::exp(slope * x) / m_norm_expo;
    }

    double EvalPol1_Cached(double x, double slope) const
    {
        if(m_norm_pol1 <= 0.0)
            return 0.0;
        double xMid = 0.5 * (m_xMin + m_xMax);
        return (1.0 / m_norm_pol1) * (1.0 + slope * (x - xMid));
    }

  public:
    FullUnbinnedPDF(double xMin, double xMax, const AuxFitResult &sigRes, const AuxFitResult &p1Res,
        const AuxFitResult &p2Res, const AuxFitResult &argusRes, const AuxFitResult &expo1Res,
        bool usePol1Bkg = false)
        : m_xMin(xMin)
        , m_xMax(xMax)
        , m_usePol1Bkg(usePol1Bkg)
        , m_norm_argus(1.0)
        , m_norm_expo(1.0)
        , m_norm_pol1(1.0)
        , m_norm_sig_g1(1.0)
        , m_norm_sig_g2(1.0)
        , m_norm_p1_g1(1.0)
        , m_norm_p1_g2(1.0)
        , m_norm_p2_g1(1.0)
        , m_norm_p2_g2(1.0)
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

    // Manteniamo le funzioni di Eval pubbliche e non-cached originali
    // per non rompere i puntatori o le lambda usate nel plotting
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

    double EvalExpo(double x, double slope) const
    {
        double norm = 0.0;
        if(std::abs(slope) < 1e-6)
        {
            norm = m_xMax - m_xMin;
        }
        else
        {
            norm = (std::exp(slope * m_xMax) - std::exp(slope * m_xMin)) / slope;
        }
        if(norm <= 0.0)
            return 0.0;
        return std::exp(slope * x) / norm;
    }

    double EvalPol1(double x, double slope) const
    {
        double range = m_xMax - m_xMin;
        if(range <= 0.0)
            return 0.0;
        double xMid = 0.5 * (m_xMin + m_xMax);
        return (1.0 / range) * (1.0 + slope * (x - xMid));
    }

    // --- OPERATORE CHIAMATO DA MINUIT ---
    double operator()(const double *x, const double *par) const
    {
        // Se i parametri correnti sono diversi da quelli dell'ultima iterazione,
        // ricalcoliamo le costanti di normalizzazione (una sola volta per step!)
        if(ParametersChanged(par))
        {
            UpdateNormalizations(par);
            m_last_pars.assign(par, par + 20);
        }

        const double xx = x[0];
        const double f_s = par[0];
        const double f_1 = par[1];
        const double f_2 = par[2];
        const double f_3 = par[3];
        const double slope = par[4];

        // Valutazione ad altissima velocità tramite le funzioni cached
        double p_sig
            = Eval2G_Cached(xx, par[5], par[6], par[7], par[8], m_norm_sig_g1, m_norm_sig_g2);
        double p_1
            = Eval2G_Cached(xx, par[9], par[10], par[11], par[12], m_norm_p1_g1, m_norm_p1_g2);
        double p_2
            = Eval2G_Cached(xx, par[13], par[14], par[15], par[16], m_norm_p2_g1, m_norm_p2_g2);
        double p_3 = EvalArgus_Cached(xx, par[17], par[18], par[19]);
        double p_4 = m_usePol1Bkg ? EvalPol1_Cached(xx, slope) : EvalExpo_Cached(xx, slope);

        double p_bkg = f_1 * p_1 + f_2 * p_2 + f_3 * p_3 + (1.0 - f_1 - f_2 - f_3) * p_4;
        return f_s * p_sig + (1.0 - f_s) * p_bkg;
    }
};

// --- NLL for unbinned fit ---
thread_local std::vector<double> g_data_events;
thread_local FullUnbinnedPDF *g_pdf_unbinned = nullptr;

double Unbinned2NLL(const double *par)
{
    double nll = 0.0;
    for(int i = 0; i < (int)g_data_events.size(); i++)
    {
        double x = g_data_events[i];
        double val = (*g_pdf_unbinned)(&x, par);
        // if(!std::isfinite(val) || val <= 0.0)
        //     val = 1e-10;
        nll -= std::log(val);
    }
    return 2.0 * nll;
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
        particleName = "D^{+} #rightarrow #phi#pi^{+}";
    else if(mcID == 41)
        particleName = "D_{s}^{+} #rightarrow #phi#pi^{+}";
    else if(mcID == 42)
        particleName = "D_{s}^{+} #rightarrow #phi#mu^{+}#nu_{#mu}";
    else if(mcID == 44)
        particleName = "D_{s}^{+} #rightarrow #tau^{+}#nu_{#tau}";

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

    constexpr bool usePol1Bkg = false;

    // --- SANITY CHECKS ---
    // Impostare a 'true' per rilasciare i rispettivi gruppi di parametri nel fit
    constexpr bool floatSignalShape = false; // Leave it false
    constexpr bool floatP1Shape = false;
    constexpr bool floatP2Shape = false;
    constexpr bool floatArgusShape = false; // Leave it false

    // 1. Esecuzione dei fit ausiliari (MC e fondo combinatorio)
    cout << "\n=== [STEP 1/3] Running MC Auxiliary Fits ===" << endl;
    LoadDataset(1);
    AuxFitResult res_sig = FitTemplateMass(44);
    AuxFitResult res_p1 = FitTemplateMass(34);
    AuxFitResult res_p2 = FitTemplateMass(41);
    AuxFitResult res_arg = FitTemplateMass(42);

    cout << "\n=== [STEP 2/3] Running Combinatorial Background Fit ===" << endl;
    AuxFitResult res_bkg = FitCombinatorialBkg(usePol1Bkg);

    if(!res_sig.isValid || !res_p1.isValid || !res_p2.isValid || !res_arg.isValid
        || !res_bkg.isValid)
    {
        cerr << "[ERROR] Auxiliary fits failed! Aborting unbinned fit." << endl;
        return;
    }

    // 2. Caricamento dei dati reali per il Fit Unbinned
    cout << "\n=== [STEP 3/3] Preparing Real Data for Unbinned Fit ===" << endl;
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
    g_pdf_unbinned
        = new FullUnbinnedPDF(xMin, xMax, res_sig, res_p1, res_p2, res_arg, res_bkg, usePol1Bkg);

    ROOT::Math::Minimizer *minimizer = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
    minimizer->SetMaxFunctionCalls(100000);
    minimizer->SetTolerance(0.01);
    minimizer->SetPrintLevel(0);

    ROOT::Math::Functor fNLL(&Unbinned2NLL, 20); // Impostato a 20 parametri
    minimizer->SetFunction(fNLL);

    // Parametri primari (Frazioni e pendenza)
    minimizer->SetVariable(0, "f_s", 0., 0.005);
    // minimizer->SetVariableLimits(0, 0.0, 1.0);
    minimizer->SetVariable(1, "f_1", 0.025, 0.010);
    // minimizer->SetVariableLimits(1, 0.0, 1.0);
    minimizer->SetVariable(2, "f_2", 0.051, 0.010);
    // minimizer->SetVariableLimits(2, 0.0, 1.0);
    minimizer->SetVariable(3, "f_3", 0.625, 0.020);
    // minimizer->SetVariableLimits(3, 0.0, 1.0);

    double init_slope = res_bkg.params[1];
    if(usePol1Bkg)
    {
        minimizer->SetVariable(4, "pol1_slope", init_slope, 0.05);
        // minimizer->SetVariableLimits(4, -10.0, 10.0);
    }
    else
    {
        minimizer->SetVariable(4, "expo_slope", init_slope, 0.1);
        // minimizer->SetVariableLimits(4, -20.0, 0.0);
    }

    // Parametri 5-8: Segnale (Double Gauss)
    minimizer->SetVariable(5, "sig_mean", g_pdf_unbinned->m_sig_mean, 0.001);
    minimizer->SetVariable(6, "sig_sigma1", g_pdf_unbinned->m_sig_sigma1, 0.0005);
    // minimizer->SetVariableLimits(6, 1e-4, 0.05);
    minimizer->SetVariable(7, "sig_sigma2", g_pdf_unbinned->m_sig_sigma2, 0.001);
    // minimizer->SetVariableLimits(7, 1e-4, 0.10);
    minimizer->SetVariable(8, "sig_frac1", g_pdf_unbinned->m_sig_frac1, 0.05);
    // minimizer->SetVariableLimits(8, 0.0, 1.0);

    // Parametri 9-12: Fondo P1 (Double Gauss)
    minimizer->SetVariable(9, "p1_mean", g_pdf_unbinned->m_p1_mean, 0.001);
    minimizer->SetVariable(10, "p1_sigma1", g_pdf_unbinned->m_p1_sigma1, 0.0005);
    // minimizer->SetVariableLimits(10, 1e-4, 0.05);
    minimizer->SetVariable(11, "p1_sigma2", g_pdf_unbinned->m_p1_sigma2, 0.001);
    minimizer->FixVariable(11);
    minimizer->SetVariable(12, "p1_frac1", g_pdf_unbinned->m_p1_frac1, 0.05);
    minimizer->FixVariable(12);

    // Parametri 13-16: Fondo P2 (Double Gauss)
    minimizer->SetVariable(13, "p2_mean", g_pdf_unbinned->m_p2_mean, 0.001);
    minimizer->SetVariable(14, "p2_sigma1", g_pdf_unbinned->m_p2_sigma1, 0.0005);
    // minimizer->SetVariableLimits(14, 1e-4, 0.05);
    minimizer->SetVariable(15, "p2_sigma2", g_pdf_unbinned->m_p2_sigma2, 0.001);
    minimizer->FixVariable(15);
    // minimizer->SetVariableLimits(15, 1e-4, 0.10);
    minimizer->SetVariable(16, "p2_frac1", g_pdf_unbinned->m_p2_frac1, 0.05);
    minimizer->FixVariable(16);
    // minimizer->SetVariableLimits(16, 0.0, 1.0);

    // Parametri 17-19: Argus (Fondo 3)
    minimizer->SetVariable(17, "argus_m0", g_pdf_unbinned->m_argus_m0, 0.005);
    // minimizer->SetVariableLimits(17, 1.90, 2.05);
    minimizer->SetVariable(18, "argus_c", g_pdf_unbinned->m_argus_c, 0.1);
    minimizer->SetVariable(19, "argus_p", g_pdf_unbinned->m_argus_p, 0.05);
    // minimizer->SetVariableLimits(19, 0.0, 5.0);

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
    if(usePol1Bkg)
    {
        pave->AddText(Form("Pol1 Slope = %.3f #pm %.3f", best_slope, errs[4]));
    }
    else
    {
        pave->AddText(Form("Lambda Slope = %.3f #pm %.3f", best_slope, errs[4]));
    }
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
    leg->AddEntry(f_sig, "Signal (D_{s}^{+} #rightarrow #tau^{+}#nu_{#tau})", "l");
    leg->AddEntry(f_p1, "D^{+} #rightarrow #phi#pi^{+} Bkg", "l");
    leg->AddEntry(f_p2, "D_{s}^{+} #rightarrow #phi#pi^{+} Bkg", "l");
    leg->AddEntry(f_argus, "D_{s}^{+} #rightarrow #phi#mu^{+}#nu_{#mu} (Argus) Bkg", "l");
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
        c_unbinned->Print("./_fig/FitFullUnbinned_Pulls.pdf");
        c_unbinned->Print("./_fig/FitFullUnbinned_Pulls.root");
    }

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
};

// Funzione worker isolata per il singolo thread
ToyResult RunSingleToy(int toyId, int nEvents, const FullUnbinnedPDF &templatePdf,
    const std::vector<double> &gen_pars, bool usePol1Bkg, double xMin, double xMax)
{
    ToyResult res;

    // 1. Generatore di numeri casuali locale al thread (evita conflitti su gRandom)
    TRandom3 threadRandom(12345 + toyId);
    gRandom = &threadRandom;

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
    fGen.SetNpx(2000);

    // 4. Generazione del dataset locale al thread
    g_data_events.clear();
    // nEvents = gRandom->Poisson(nEvents);
    g_data_events.reserve(nEvents);
    for(int ev = 0; ev < nEvents; ++ev)
    {
        g_data_events.push_back(fGen.GetRandom());
    }

    // 5. Setup locale del Minimizzatore
    ROOT::Math::Minimizer *minimizer = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
    minimizer->SetMaxFunctionCalls(50000);
    minimizer->SetTolerance(0.01);
    minimizer->SetPrintLevel(-1); // Silenzioso

    ROOT::Math::Functor fNLL(&Unbinned2NLL, 20);
    minimizer->SetFunction(fNLL);

    // Definizione delle variabili locali
    minimizer->SetVariable(0, "f_s", 0., 0.005);
    minimizer->SetVariable(1, "f_1", 0.025, 0.010);
    minimizer->SetVariable(2, "f_2", 0.051, 0.010);
    minimizer->SetVariable(3, "f_3", 0.625, 0.020);

    if(usePol1Bkg)
    {
        minimizer->SetVariable(4, "pol1_slope", gen_pars[4], 0.05);
    }
    else
    {
        minimizer->SetVariable(4, "expo_slope", gen_pars[4], 0.01);
        minimizer->FixVariable(4);
    }

    // Congelamento dei parametri di forma (fissati ai parametri nominali di generazione)
    minimizer->SetVariable(5, "sig_mean", gen_pars[5], 0.001);
    minimizer->FixVariable(5);
    minimizer->SetVariable(6, "sig_sigma1", gen_pars[6], 0.0005);
    minimizer->FixVariable(6);
    minimizer->SetVariable(7, "sig_sigma2", gen_pars[7], 0.001);
    minimizer->FixVariable(7);
    minimizer->SetVariable(8, "sig_frac1", gen_pars[8], 0.05);
    minimizer->FixVariable(8);

    minimizer->SetVariable(9, "p1_mean", gen_pars[9], 0.001);
    minimizer->FixVariable(9);
    minimizer->SetVariable(10, "p1_sigma1", gen_pars[10], 0.0005);
    minimizer->FixVariable(10);
    minimizer->SetVariable(11, "p1_sigma2", gen_pars[11], 0.001);
    minimizer->FixVariable(11);
    minimizer->SetVariable(12, "p1_frac1", gen_pars[12], 0.05);
    minimizer->FixVariable(12);

    minimizer->SetVariable(13, "p2_mean", gen_pars[13], 0.001);
    minimizer->FixVariable(13);
    minimizer->SetVariable(14, "p2_sigma1", gen_pars[14], 0.0005);
    minimizer->FixVariable(14);
    minimizer->SetVariable(15, "p2_sigma2", gen_pars[15], 0.001);
    minimizer->FixVariable(15);
    minimizer->SetVariable(16, "p2_frac1", gen_pars[16], 0.05);
    minimizer->FixVariable(16);

    minimizer->SetVariable(17, "argus_m0", gen_pars[17], 0.005);
    minimizer->FixVariable(17);
    minimizer->SetVariable(18, "argus_c", gen_pars[18], 0.1);
    minimizer->FixVariable(18);
    minimizer->SetVariable(19, "argus_p", gen_pars[19], 0.05);
    minimizer->FixVariable(19);

    minimizer->Minimize();
    minimizer->Hesse();

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
    }

    delete minimizer;
    g_pdf_unbinned = nullptr; // Reset
    return res;
}

void analysis::RunToyMC(int nToys)
{
    auto start = std::chrono::high_resolution_clock::now();

    SetLBStyle();
    constexpr bool usePol1Bkg = false;

    TF1::DefaultAddToGlobalList(kFALSE);

    // 1. Esecuzione dei fit ausiliari nominali
    cout << "\n=== [TOY MC] Running Auxiliary Fits ===" << endl;
    LoadDataset(1);
    AuxFitResult res_sig = FitTemplateMass(44);
    AuxFitResult res_p1 = FitTemplateMass(34);
    AuxFitResult res_p2 = FitTemplateMass(41);
    AuxFitResult res_arg = FitTemplateMass(42);
    AuxFitResult res_bkg = FitCombinatorialBkg(usePol1Bkg);

    if(!res_sig.isValid || !res_p1.isValid || !res_p2.isValid || !res_arg.isValid
        || !res_bkg.isValid)
    {
        cerr << "[ERROR] Auxiliary fits failed! Aborting Toy MC." << endl;
        return;
    }

    // Frazioni e parametri fisici di generazione
    const double true_fs = 0.0;
    const double true_f1 = 0.026;
    const double true_f2 = 0.051;
    const double true_f3 = 0.625;
    const double true_slope = res_bkg.params[1];

    std::vector<double> gen_pars = { true_fs, true_f1, true_f2, true_f3, true_slope,
        res_sig.params[1], res_sig.params[2], res_sig.params[3], res_sig.params[4],
        res_p1.params[1], res_p1.params[2], res_p1.params[3], res_p1.params[4], res_p2.params[1],
        res_p2.params[2], res_p2.params[3], res_p2.params[4], res_arg.params[1], res_arg.params[2],
        res_arg.params[3] };

    Double_t xMin = 1.65;
    Double_t xMax = 2.09;

    // PDF modello base per i thread
    FullUnbinnedPDF templatePdf(xMin, xMax, res_sig, res_p1, res_p2, res_arg, res_bkg, usePol1Bkg);

    // Conteggio eventi sui dati reali
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

    // Rilevamento automatico del numero di CPU disponibili
    unsigned int nCores = std::thread::hardware_concurrency();
    if(nCores == 0)
        nCores = 4; // Fallback generico
    cout << "[INFO] Launching Toy MC using " << nCores << " parallel threads." << endl;

    std::vector<ToyResult> toyResults;
    toyResults.reserve(nToys);

    // 2. Lancio dei thread in blocchi di dimensione pari a 'nCores'
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

        // Raccogliamo i risultati di questa tranche (questo sincronizza i thread)
        for(auto &f : futures)
        {
            toyResults.push_back(f.get());
        }
        cout << "  Completed toys: " << toyResults.size() << " / " << nToys << "..." << endl;
    }

    // 3. Riempimento degli istogrammi post-elaborazione parallela
    auto h_fit_fs = new TH1D("h_fit_fs", "Fitted f_{s};f_{s};Toys", 100, 0., 0.);
    auto h_fit_f1 = new TH1D("h_fit_f1", "Fitted f_{1};f_{1};Toys", 100, 0., 0.);
    auto h_fit_f2 = new TH1D("h_fit_f2", "Fitted f_{2};f_{2};Toys", 100, 0., 0.);
    auto h_fit_f3 = new TH1D("h_fit_f3", "Fitted f_{3};f_{3};Toys", 100, 0., 0.);

    auto h_pull_fs = new TH1D("h_pull_fs", "Pull f_{s};Pull;Toys", 50, -5.0, 5.0);
    auto h_pull_f1 = new TH1D("h_pull_f1", "Pull f_{1};Pull;Toys", 50, -5.0, 5.0);
    auto h_pull_f2 = new TH1D("h_pull_f2", "Pull f_{2};Pull;Toys", 50, -5.0, 5.0);
    auto h_pull_f3 = new TH1D("h_pull_f3", "Pull f_{3};Pull;Toys", 50, -5.0, 5.0);

    int convergedToys = 0;
    for(const auto &res : toyResults)
    {
        if(res.converged)
        {
            convergedToys++;
            h_fit_fs->Fill(res.fs_val - true_fs);
            h_fit_f1->Fill(res.f1_val - true_f1);
            h_fit_f2->Fill(res.f2_val - true_f2);
            h_fit_f3->Fill(res.f3_val - true_f3);

            if(res.fs_err > 0)
                h_pull_fs->Fill((res.fs_val - true_fs) / res.fs_err);
            if(res.f1_err > 0)
                h_pull_f1->Fill((res.f1_val - true_f1) / res.f1_err);
            if(res.f2_err > 0)
                h_pull_f2->Fill((res.f2_val - true_f2) / res.f2_err);
            if(res.f3_err > 0)
                h_pull_f3->Fill((res.f3_val - true_f3) / res.f3_err);
        }
    }

    cout << "\n[TOY MC RESULTS] Converged: " << convergedToys << " / " << nToys << endl;

    // 4. Visualizzazione e salvataggio dei Canvas
    TF1::DefaultAddToGlobalList(kTRUE);

    TCanvas *cToys = new TCanvas("cToys", "Toy MC Study Results", 1600, 800);
    cToys->Divide(4, 2);

    auto drawResult = [](TVirtualPad *pad, TH1D *h, double trueVal, bool isPull)
    {
        pad->cd();
        h->SetStats(kTRUE);
        gStyle->SetOptStat("emr");
        h->Draw();

        if(!isPull)
        {
            TLine *line = new TLine(trueVal, 0, trueVal, h->GetMaximum() * 1.05);
            line->SetLineColor(kRed);
            line->SetLineWidth(2);
            line->SetLineStyle(2);
            line->Draw();
        }
        else
        {
            h->Fit("gaus", "Q L");
            gStyle->SetOptFit(111);
        }
    };

    drawResult(cToys->GetPad(1), h_fit_fs, true_fs, false);
    drawResult(cToys->GetPad(2), h_fit_f1, true_f1, false);
    drawResult(cToys->GetPad(3), h_fit_f2, true_f2, false);
    drawResult(cToys->GetPad(4), h_fit_f3, true_f3, false);

    drawResult(cToys->GetPad(5), h_pull_fs, 0.0, true);
    drawResult(cToys->GetPad(6), h_pull_f1, 0.0, true);
    drawResult(cToys->GetPad(7), h_pull_f2, 0.0, true);
    drawResult(cToys->GetPad(8), h_pull_f3, 0.0, true);

    cToys->Update();

    if(savePlots)
    {
        cToys->Print("./_fig/ToyMC_Results.pdf");
        cToys->Print("./_fig/ToyMC_Results.root");
    }

    // ================================================================================
    // 5. Generazione e Plot di un Toy di Esempio (Diagnostica nel thread principale)
    // ================================================================================
    cout << "\n=== [TOY MC] Generating and fitting a single Toy example for diagnostics ==="
         << endl;

    // Generatore di numeri casuali locale al thread principale
    TRandom3 mainRandom(123456);
    gRandom = &mainRandom;

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
        double val = fGenDiag.GetRandom();
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
    minDiag->SetVariable(1, "f_1", 0.026, 0.010);
    minDiag->SetVariable(2, "f_2", 0.051, 0.010);
    minDiag->SetVariable(3, "f_3", 0.625, 0.020);

    if(usePol1Bkg)
    {
        minDiag->SetVariable(4, "pol1_slope", gen_pars[4], 0.05);
    }
    else
    {
        minDiag->SetVariable(4, "expo_slope", gen_pars[4], 0.1);
    }

    // Fissiamo tutte le forme ai parametri veri
    for(int p = 5; p < 20; ++p)
    {
        minDiag->SetVariable(p, Form("p_%d", p), gen_pars[p], 0.01);
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

    // Creazione ed overlay dei TF1 di fit sul Canvas
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
    cSingleToy->Print("./_fig/ToyMC_SingleFit_Diagnostic.pdf");
    cSingleToy->Print("./_fig/ToyMC_SingleFit_Diagnostic.root");

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
}
