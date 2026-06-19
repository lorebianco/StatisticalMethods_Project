#include "RtypesCore.h"
#define analysis_cxx
#include <algorithm>
#include <cmath>
#include <fstream>
#include <future>
#include <iomanip>
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
#include <TGraphErrors.h>
#include <TGraphSmooth.h>
#include <TH2.h>
#include <TLegend.h>
#include <TLine.h>
#include <TMarker.h>
#include <TMath.h>
#include <TMatrixDSym.h>
#include <TPaveStats.h>
#include <TRandom3.h>
#include <TSpline.h>
#include <TStyle.h>
#include <TTreeFormula.h>

#include "analysis.h"
#include "lbrootstyle.hh"

using namespace std;
using namespace ROOT::Math;
using namespace TMath;
using namespace lbStyle;

constexpr Bool_t savePlots = true;
constexpr Bool_t USE_OPTIMIZED_CUTS = true;
constexpr Bool_t USE_UNBINNED = false;
constexpr Bool_t USE_EXTENDED = true;
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

inline Double_t Pol2Blinded(Double_t *x, Double_t *par)
{
    // Definizione dei limiti di blinding (3 sigma attorno alla massa del tau)
    const double blindMin = 1.777 - 3 * 0.0058; // 1.7596 GeV
    const double blindMax = 1.777 + 3 * 0.0058; // 1.7944 GeV

    if(x[0] > blindMin && x[0] < blindMax)
    {
        TF1::RejectPoint(); // Esclude la blinding region dal fit
        return 0.0;
    }

    double m_shift = x[0] - 1.85; // Shift spaziale per stabilità numerica del fit
    return par[0] + par[1] * m_shift + par[2] * m_shift * m_shift;
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
thread_local TH1D *g_data_hist = nullptr;

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

inline double Binned2NLL(const double *par)
{
    double nll = 0.0;
    int nBins = g_data_hist->GetNbinsX();

    for(int i = 1; i <= nBins; i++)
    {
        double n_obs = g_data_hist->GetBinContent(i);

        if(n_obs > 0)
        {
            // Ottieni i bordi e il centro del bin
            double x_low = g_data_hist->GetXaxis()->GetBinLowEdge(i);
            double x_up = g_data_hist->GetXaxis()->GetBinUpEdge(i);
            double x_mid = g_data_hist->GetBinCenter(i);

            // Valuta la PDF in 3 punti
            double f_low = (*g_pdf_unbinned)(&x_low, par);
            double f_mid = (*g_pdf_unbinned)(&x_mid, par);
            double f_up = (*g_pdf_unbinned)(&x_up, par);

            // Regola di Simpson per l'integrale del bin: (dx / 6) * (f(a) + 4f(m) + f(b))
            double bin_width = x_up - x_low;
            double p_i = (bin_width / 6.0) * (f_low + 4.0 * f_mid + f_up);

            nll -= n_obs * std::log(p_i);
        }
    }
    return 2.0 * nll;
}

inline double BinnedExtended2NLL(const double *par)
{
    double nll = 0.0;
    int nBins = g_data_hist->GetNbinsX();

    // Il parametro 20 è il numero totale atteso di eventi (N_tot)
    double N_tot = par[20];

    for(int i = 1; i <= nBins; i++)
    {
        double n_obs = g_data_hist->GetBinContent(i);

        // Integrazione del bin con la regola di Simpson
        double x_low = g_data_hist->GetXaxis()->GetBinLowEdge(i);
        double x_up = g_data_hist->GetXaxis()->GetBinUpEdge(i);
        double x_mid = g_data_hist->GetBinCenter(i);

        double f_low = (*g_pdf_unbinned)(&x_low, par);
        double f_mid = (*g_pdf_unbinned)(&x_mid, par);
        double f_up = (*g_pdf_unbinned)(&x_up, par);

        double bin_width = x_up - x_low;
        double p_i = (bin_width / 6.0) * (f_low + 4.0 * f_mid + f_up); // Probabilità nel bin

        double mu_i = N_tot * p_i; // Eventi attesi in questo bin

        if(n_obs > 0)
            nll += (mu_i - n_obs * std::log(mu_i));
        else
            nll += mu_i; // Se il bin è vuoto, rimane solo il termine lineare mu_i
    }
    return 2.0 * nll;
}

// --- Likelihood ratio ordering (Feldman-Cousins con sigma costante per fetta mu)
inline Double_t LROrdering(Double_t *x, Double_t *par)
{
    double xx = x[0];
    double mu = par[0];
    double sigma = par[1]; // Sigma costante per questa specifica fetta della banda

    if(xx >= 0)
    {
        // Caso x >= 0: mu_best = x (il denominatore della verosimiglianza è pari a 1)
        return exp(-0.5 * (xx - mu) * (xx - mu) / (sigma * sigma));
    }
    else
    {
        // Caso x < 0: mu_best = 0 (il denominatore è valutato in 0)
        return exp(
            (-0.5 * (xx - mu) * (xx - mu) / (sigma * sigma)) + (0.5 * (xx * xx) / (sigma * sigma)));
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

double analysis::GetMCEfficiencyRatio(bool useOptimized)
{
    // Numero di eventi generati alla sorgente (Tabella 1)
    const double n_gen_sig = 1.0e6; // ID 44 (Signal Ds -> tau nu)
    const double n_gen_norm = 5.0e6; // ID 42 (Norm Ds -> phi mu nu)

    LoadDataset(1); // Carichiamo il file MC
    if(!fChain)
        return 0.0;

    TString cutString = GetCutString(useOptimized);
    TTreeFormula formula("cut", cutString.Data(), fChain);
    Int_t currentTreeNumber = -1;

    long long n_pass_sig = 0;
    long long n_pass_norm = 0;

    Long64_t nentries = fChain->GetEntries();
    for(Long64_t jentry = 0; jentry < nentries; jentry++)
    {
        if(LoadTree(jentry) < 0)
            break;
        if(fChain->GetTreeNumber() != currentTreeNumber)
        {
            currentTreeNumber = fChain->GetTreeNumber();
            formula.UpdateFormulaLeaves();
        }
        fChain->GetEntry(jentry);

        // Conta gli eventi sopravvissuti ai tagli attivi
        if(formula.EvalInstance() > 0)
        {
            if(id == 44)
                n_pass_sig++;
            if(id == 42)
                n_pass_norm++;
        }
    }

    double eff_sig = (double)n_pass_sig / n_gen_sig;
    double eff_norm = (double)n_pass_norm / n_gen_norm;

    if(eff_sig <= 0.0)
        return 0.0;
    return eff_norm / eff_sig; // Ritorna il rapporto r_eff dinamico
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

    SetLBStyle();
    gStyle->SetOptFit(1111);

    TCanvas *cMass = new TCanvas(Form("cMass_%d", mcID), Form("Mass Fit MC %d", mcID), 800, 600);
    cMass->cd();

    std::vector<Double_t> mass_vals;

    // Setup della formula di taglio attiva
    TString cutString = GetCutString(USE_OPTIMIZED_CUTS);
    TTreeFormula formula("cut", cutString.Data(), fChain);
    Int_t currentTreeNumber = -1;

    Long64_t nentries = fChain->GetEntries();
    for(Long64_t jentry = 0; jentry < nentries; jentry++)
    {
        if(LoadTree(jentry) < 0)
            break;
        if(fChain->GetTreeNumber() != currentTreeNumber)
        {
            currentTreeNumber = fChain->GetTreeNumber();
            formula.UpdateFormulaLeaves();
        }
        fChain->GetEntry(jentry);

        if(id != mcID)
            continue;

        // Applica il filtro se attivo
        if(formula.EvalInstance() > 0)
        {
            mass_vals.push_back(D_M);
        }
    }

    if(mass_vals.empty())
    {
        cout << "[WARNING] No events found for mcID = " << mcID << endl;
        return res;
    }

    // ... [Il codice intermedio di calcolo di meanInit, rmsInit, xMin, xMax rimane identico]
    // ...
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
        xMax = 2.0;

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

        // 3. AGGIORNA IL CANVAS (Questo passaggio mancava e generava un puntatore nullo su
        // 'stats')
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

void analysis::DoFullBlindedFit()
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

    Double_t xMin = 1.6;
    Double_t xMax = 2.10;

    g_data_events.clear();
    Long64_t nentries = fChain->GetEntries();
    g_data_events.reserve(static_cast<size_t>(nentries));

    auto h_mass
        = new TH1D("h_mass_unbinned", ";Invariant Mass [GeV/#it{c}^{2}];Entries", 200, xMin, xMax);
    h_mass->SetDirectory(nullptr);

    // Inizializza la formula basata sul constexpr globale
    TString cutString = GetCutString(USE_OPTIMIZED_CUTS);
    TTreeFormula formula("cut", cutString.Data(), fChain);
    Int_t currentTreeNumber = -1;

    for(Long64_t jentry = 0; jentry < fChain->GetEntriesFast(); jentry++)
    {
        if(LoadTree(jentry) < 0)
            break;
        if(fChain->GetTreeNumber() != currentTreeNumber)
        {
            currentTreeNumber = fChain->GetTreeNumber();
            formula.UpdateFormulaLeaves();
        }
        fChain->GetEntry(jentry);
        if(id != 0)
            continue;

        // Filtra i dati reali applicando l'interruttore
        if(formula.EvalInstance() > 0)
        {
            if(D_M >= xMin && D_M <= xMax)
            {
                g_data_events.push_back(D_M);
                h_mass->Fill(D_M);
            }
        }
    }
    cout << "NEntries = " << h_mass->GetEntries() << endl;
    cout << endl;

    g_data_hist = h_mass;

    Double_t binWidth = h_mass->GetBinWidth(1);
    Double_t nEntries = h_mass->GetEntries();

    // 3. Setup della PDF e di MINUIT a 20 parametri
    g_pdf_unbinned = new FullUnbinnedPDF(xMin, xMax, res_sig, res_p1, res_p2, res_arg, usePol1Bkg);

    ROOT::Math::Minimizer *minimizer = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
    minimizer->SetMaxFunctionCalls(100000);
    minimizer->SetTolerance(0.01);
    minimizer->SetPrintLevel(0);

    // Sceglie dinamicamente la NLL in base al flag constexpr
    // // Imposta la dimensione del functor a 21 se extended, altrimenti 20
    int nPars = USE_EXTENDED ? 21 : 20;
    ROOT::Math::Functor fNLL(
        USE_EXTENDED ? &BinnedExtended2NLL : (USE_UNBINNED ? &Unbinned2NLL : &Binned2NLL), nPars);
    minimizer->SetFunction(fNLL);

    // Parametri primari (Frazioni e pendenza)
    minimizer->SetVariable(0, "f_s", 0., 0.001);
    minimizer->SetVariable(1, "f_1", 0.021, 0.001);
    minimizer->SetVariable(2, "f_2", 0.043, 0.001);
    minimizer->SetVariable(3, "f_3", 0.656, 0.001);

    double init_slope = -0.932563;
    if(usePol1Bkg)
        minimizer->SetVariable(4, "pol1_slope", init_slope, 0.001);
    else
        minimizer->SetVariable(4, "expo_slope", init_slope, 0.001);

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

    // Se extended, aggiungiamo N_tot come parametro libero
    if(USE_EXTENDED)
    {
        double nEntries = h_mass->GetEntries();
        minimizer->SetVariable(20, "N_tot", nEntries, std::sqrt(nEntries));
    }

    cout << "\n--- Minimizing Unbinned Likelihood with MINUIT (Silenced) ---" << endl;
    minimizer->Minimize();
    minimizer->Hesse();

    // --- NUOVA PARTE: Salvataggio dei parametri per i Toy ---
    if(minimizer->Status() == 0)
    {
        const double *xs = minimizer->X();
        unsigned int nDim = minimizer->NDim();
        m_fitted_pars.assign(xs, xs + nDim);
        m_fit_done = true;
    }
    else
    {
        cerr << "[WARNING] Full unbinned fit did not converge. Using fallback parameters." << endl;
        m_fit_done = false;
    }

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

    // Calcolo del fattore di scala corretto (se il fit è esteso usa N_tot fittato, altrimenti le
    // entrate dell'istogramma)
    double nTotFit = USE_EXTENDED ? xs[20] : nEntries;

    // Aggiorniamo le lambda expression affinché usino nTotFit e le variabili post-fit locali
    // 'best_...'
    auto scale_pdf_lambda = [pdf_for_draw, nTotFit, binWidth, xs](double *x, double *par)
    {
        double xx = x[0];
        return nTotFit * binWidth * (*pdf_for_draw)(&xx, xs);
    };

    auto scale_sig_lambda
        = [pdf_for_draw, nTotFit, binWidth, best_f_s, best_sig_mean, best_sig_sigma1,
              best_sig_sigma2, best_sig_frac1](double *x, double *par)
    {
        double xx = x[0];
        double p_sig = pdf_for_draw->Eval2G(
            xx, best_sig_mean, best_sig_sigma1, best_sig_sigma2, best_sig_frac1);
        return nTotFit * binWidth * (best_f_s * p_sig);
    };

    auto scale_comb_lambda = [pdf_for_draw, nTotFit, binWidth, best_f_s, best_f_1, best_f_2,
                                 best_f_3, best_slope, usePol1Bkg](double *x, double *par)
    {
        double xx = x[0];
        double p_4 = usePol1Bkg ? pdf_for_draw->EvalPol1(xx, best_slope)
                                : pdf_for_draw->EvalExpo(xx, best_slope);
        double frac_comb = (1.0 - best_f_s) * (1.0 - best_f_1 - best_f_2 - best_f_3);
        return nTotFit * binWidth * (frac_comb * p_4);
    };

    auto scale_p1_lambda
        = [pdf_for_draw, nTotFit, binWidth, best_f_s, best_f_1, best_p1_mean, best_p1_sigma1,
              best_p1_sigma2, best_p1_frac1](double *x, double *par)
    {
        double xx = x[0];
        double p_1
            = pdf_for_draw->Eval2G(xx, best_p1_mean, best_p1_sigma1, best_p1_sigma2, best_p1_frac1);
        return nTotFit * binWidth * ((1.0 - best_f_s) * best_f_1 * p_1);
    };

    auto scale_p2_lambda
        = [pdf_for_draw, nTotFit, binWidth, best_f_s, best_f_2, best_p2_mean, best_p2_sigma1,
              best_p2_sigma2, best_p2_frac1](double *x, double *par)
    {
        double xx = x[0];
        double p_2
            = pdf_for_draw->Eval2G(xx, best_p2_mean, best_p2_sigma1, best_p2_sigma2, best_p2_frac1);
        return nTotFit * binWidth * ((1.0 - best_f_s) * best_f_2 * p_2);
    };

    auto scale_argus_lambda = [pdf_for_draw, nTotFit, binWidth, best_f_s, best_f_3, best_argus_m0,
                                  best_argus_c, best_argus_p](double *x, double *par)
    {
        double xx = x[0];
        double p_3 = pdf_for_draw->EvalArgus(xx, best_argus_m0, best_argus_c, best_argus_p);
        return nTotFit * binWidth * ((1.0 - best_f_s) * best_f_3 * p_3);
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
    const double n_gen_Dplus = 50.0e6;
    const double n_gen_Dsplus = 50.0e6;

    LoadDataset(1); // Carichiamo il file MC
    long long n_pass_Dplus = 0;
    long long n_pass_Dsplus = 0;

    // Usiamo nomi univoci (cutStringVal, formulaVal, currentTreeNumberVal) per evitare conflitti
    TString cutStringVal = GetCutString(USE_OPTIMIZED_CUTS);
    TTreeFormula formulaVal("cutVal", cutStringVal.Data(), fChain);
    Int_t currentTreeNumberVal = -1;

    for(Long64_t jentry = 0; jentry < fChain->GetEntriesFast(); jentry++)
    {
        if(LoadTree(jentry) < 0)
            break;
        if(fChain->GetTreeNumber() != currentTreeNumberVal)
        {
            currentTreeNumberVal = fChain->GetTreeNumber();
            formulaVal.UpdateFormulaLeaves();
        }
        fChain->GetEntry(jentry);

        // Conta solo gli eventi che superano i tagli attivi
        if(formulaVal.EvalInstance() > 0)
        {
            if(id == 34)
                n_pass_Dplus++;
            if(id == 41)
                n_pass_Dsplus++;
        }
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
    delete g_pdf_unbinned; // Dealloca l'oggetto prima di azzerare il puntatore
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
    TH1D *local_hist = nullptr; // Va distrutto alla fine del worker

    if(USE_UNBINNED)
    {
        g_data_events.clear();
        g_data_events.reserve(nEvents);
        for(int ev = 0; ev < nEvents; ++ev)
        {
            g_data_events.push_back(fGen.GetRandom(&threadRandom));
        }
    }
    else
    {
        int nBins = std::round((xMax - xMin) / 0.0025);
        local_hist = new TH1D(Form("h_toy_%d", toyId), "", nBins, xMin, xMax);
        g_data_hist = local_hist;

        // FLUTTUAZIONE POISSONIANA DI N_tot (Proprietà fondamentale dell'Extended Fit)
        int nEventsFluct = USE_EXTENDED ? threadRandom.Poisson(nEvents) : nEvents;
        for(int ev = 0; ev < nEventsFluct; ++ev)
        {
            double val = fGen.GetRandom(&threadRandom);
            local_hist->Fill(val);
        }
    }

    // 5. Setup locale del Minimizzatore
    ROOT::Math::Minimizer *minimizer = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
    minimizer->SetMaxFunctionCalls(50000);
    minimizer->SetTolerance(0.01);
    minimizer->SetPrintLevel(-1); // Silenzioso

    int nPars = USE_EXTENDED ? 21 : 20;
    ROOT::Math::Functor fNLL(
        USE_EXTENDED ? &BinnedExtended2NLL : (USE_UNBINNED ? &Unbinned2NLL : &Binned2NLL), nPars);
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

    if(USE_EXTENDED)
    {
        minimizer->SetVariable(20, "N_tot", nEvents, std::sqrt(nEvents));
    }

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
    if(!USE_UNBINNED && local_hist)
        delete local_hist; // Rimuove il TH1 allocato per questo toy
    g_pdf_unbinned = nullptr; // Reset
    return res;
}

ToyMetrics analysis::RunToyMC(int nToys, double true_fs)
{
    auto start = std::chrono::high_resolution_clock::now();

    ROOT::EnableThreadSafety();
    TH1::AddDirectory(kFALSE);

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
        return ToyMetrics {};
    }

    // --- [STEP 2] Stima delle efficienze dal Monte Carlo ---
    cout << "\n--> Estimating Selection Efficiencies from MC..." << endl;
    const double n_gen_sig = 1.0e6;
    const double n_gen_norm = 5.0e6;

    long long n_pass_sig = 0;
    long long n_pass_norm = 0;

    // Prepariamo la formula di selezione per il MC con nomi univoci (cutStringMC, formulaMC)
    TString cutStringMC = GetCutString(USE_OPTIMIZED_CUTS);
    TTreeFormula formulaMC("cutMC", cutStringMC.Data(), fChain);
    Int_t treeNumMC = -1;

    Long64_t nentries_mc = fChain->GetEntries();
    for(Long64_t jentry = 0; jentry < nentries_mc; jentry++)
    {
        if(LoadTree(jentry) < 0)
            break;
        if(fChain->GetTreeNumber() != treeNumMC)
        {
            treeNumMC = fChain->GetTreeNumber();
            formulaMC.UpdateFormulaLeaves();
        }
        fChain->GetEntry(jentry);

        // Conta solo se l'evento MC supera i tagli attivi
        if(formulaMC.EvalInstance() > 0)
        {
            if(id == 44)
                n_pass_sig++; // Ds+ -> tau+ nu_tau (Segnale)
            if(id == 42)
                n_pass_norm++; // Ds+ -> phi mu+ nu_mu (Normalizzazione)
        }
    }

    double eff_sig = (double)n_pass_sig / n_gen_sig;
    double eff_norm = (double)n_pass_norm / n_gen_norm;

    // --- Calcolo errori binomiali individuali ---
    double eff_sig_err = std::sqrt(eff_sig * (1.0 - eff_sig) / n_gen_sig);
    double eff_norm_err = std::sqrt(eff_norm * (1.0 - eff_norm) / n_gen_norm);

    double r_eff = eff_norm / eff_sig;
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

    // --- Controllo dinamico del Fit reale ---
    if(!m_fit_done || m_fitted_pars.size() < 5)
    {
        cout << "[INFO] Parametri del fit reale non trovati. Esecuzione automatica di "
                "DoFullBlindedFit()..."
             << endl;
        DoFullBlindedFit();
    }

    // Verifichiamo nuovamente se adesso il fit è presente e valido
    double true_fs_val = true_fs;
    double true_f1, true_f2, true_f3, true_slope;

    if(m_fit_done && m_fitted_pars.size() >= 5)
    {
        true_f1 = m_fitted_pars[1];
        true_f2 = m_fitted_pars[2];
        true_f3 = m_fitted_pars[3];
        true_slope = m_fitted_pars[4];
        cout << "[INFO] Toy MC inizializzato correttamente con i parametri del fit reale:" << endl;
        cout << Form("  f1 = %.4f | f2 = %.4f | f3 = %.4f | slope = %.4f", true_f1, true_f2,
            true_f3, true_slope)
             << endl;
    }
    else
    {
        cerr << "[ERROR] Impossibile eseguire i Toy MC: il fit automatico sui dati reali ha "
                "fallito o non è valido."
             << endl;
        return ToyMetrics {}; // Interrompe l'esecuzione restituendo una struttura vuota
    }

    double true_est = true_fs_val / ((1.0 - true_fs_val) * true_f3);
    double true_br = true_est * k_factor_nom;
    double true_rtau = true_est * r_eff;

    std::vector<double> gen_pars = { true_fs_val, true_f1, true_f2, true_f3, true_slope,
        res_sig.params[1], res_sig.params[2], res_sig.params[3], res_sig.params[4],
        res_p1.params[1], res_p1.params[2], res_p1.params[3], res_p1.params[4], res_p2.params[1],
        res_p2.params[2], res_p2.params[3], res_p2.params[4], res_arg.params[1], res_arg.params[2],
        res_arg.params[3] };

    Double_t xMin = 1.60;
    Double_t xMax = 2.10;

    FullUnbinnedPDF templatePdf(xMin, xMax, res_sig, res_p1, res_p2, res_arg, usePol1Bkg);

    LoadDataset(0);
    int nEvents = 0;

    // Prepariamo la formula di selezione per i dati reali con nomi univoci (cutStringData,
    // formulaData)
    TString cutStringData = GetCutString(USE_OPTIMIZED_CUTS);
    TTreeFormula formulaData("cutData", cutStringData.Data(), fChain);
    Int_t treeNumData = -1;

    for(Long64_t jentry = 0; jentry < fChain->GetEntriesFast(); jentry++)
    {
        if(LoadTree(jentry) < 0)
            break;
        if(fChain->GetTreeNumber() != treeNumData)
        {
            treeNumData = fChain->GetTreeNumber();
            formulaData.UpdateFormulaLeaves();
        }
        fChain->GetEntry(jentry);
        if(id == 0 && D_M >= xMin && D_M <= xMax)
        {
            // Applica la selezione attiva prima di incrementare il conteggio eventi reali
            if(formulaData.EvalInstance() > 0)
            {
                nEvents++;
            }
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
        cout << "\r  Completed toys: " << toyResults.size() << " / " << nToys << "..." << flush;
    }

    // --- [STEP 4] Booking Istogrammi ---
    // Istogrammi f_s
    auto h_fit_fs = new TH1D("h_fit_fs", "Fitted f_{s};#Deltaf_{s};Toys", 50, -0, 0);
    auto h_pull_fs = new TH1D("h_pull_fs", "Pull f_{s};Pull;Toys", 50, -5.0, 5.0);

    // Istogrammi per l'estimatore di resa grezzo Y = f_s / ((1 - f_s) * f_3)
    auto h_fit_est = new TH1D("h_fit_est", "Fitted Yield Ratio Y;#DeltaY;Toys", 50, -0, 0);
    auto h_pull_est = new TH1D("h_pull_est", "Pull Yield Ratio Y;Pull;Toys", 50, -5.0, 5.0);

    // Istogrammi per il Branching Ratio fisico
    auto h_fit_br = new TH1D(
        "h_fit_br", "Fitted #f{B}(#tau^{+}#rightarrow#phi#mu^{+});#Delta#f{B};Toys", 50, -0, 0);
    auto h_pull_br = new TH1D(
        "h_pull_br", "Pull #f{B}(#tau^{+}#rightarrow#phi#mu^{+});Pull;Toys", 50, -5.0, 5.0);

    // Istogrammi per R_tau
    auto h_fit_rtau = new TH1D("h_fit_rtau", "Fitted R_{#tau};#Delta R_{#tau};Toys", 50, -0, 0);
    auto h_pull_rtau = new TH1D("h_pull_rtau", "Pull R_{#tau};Pull;Toys", 50, -5.0, 5.0);

    // Altri parametri di background
    auto h_fit_f1 = new TH1D("h_fit_f1", "Fitted f_{1};#Deltaf_{1};Toys", 50, -0, 0);
    auto h_fit_f2 = new TH1D("h_fit_f2", "Fitted f_{2};#Deltaf_{2};Toys", 50, -0, 0);
    auto h_fit_f3 = new TH1D("h_fit_f3", "Fitted f_{3};#Deltaf_{3};Toys", 50, -0, 0);
    auto h_fit_slope = new TH1D("h_fit_slope", "Fitted Slope;#DeltaSlope;Toys", 50, -0, 0);

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
                double eff_sig_toy = mainRandom.Gaus(eff_sig, eff_sig_err);
                double eff_norm_toy = mainRandom.Gaus(eff_norm, eff_norm_err);
                // Impedisci fluttuazioni fisicamente impossibili sotto o uguale a zero
                if(eff_sig_toy <= 0.0)
                    eff_sig_toy = eff_sig;
                if(eff_norm_toy <= 0.0)
                    eff_norm_toy = eff_norm;

                double r_eff_toy = eff_norm_toy / eff_sig_toy;
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

                // ==================== CALCOLO R_TAU ====================
                double rtau_val = est_val * r_eff_toy;
                h_fit_rtau->Fill(rtau_val - true_rtau);

                // Varianza R_tau (ignora le B(D) del PDG!)
                double variance_rtau
                    = (r_eff * r_eff * variance_est) + (est_val * est_val * r_eff_err * r_eff_err);
                if(variance_rtau > 0.0)
                    h_pull_rtau->Fill((rtau_val - true_rtau) / std::sqrt(variance_rtau));
            }
        }
    }

    cout << "\n[TOY MC RESULTS] Converged: " << convergedToys << " / " << nToys << endl;

    TF1::DefaultAddToGlobalList(kTRUE);

    // --- [STEP 6] Disegno delle due Canvas richieste ---
    auto drawResult = [](TVirtualPad *pad, TH1D *h)
    {
        pad->cd();
        h->SetStats(kTRUE);
        gStyle->SetOptStat("emr");
        h->Draw();
        h->Fit("gaus", "Q L I");
        gStyle->SetOptFit(111);
    };

    // Canvas 1 aggiornata
    TCanvas *cMainEst = new TCanvas("cMainEst", "Estimators (fs, Y, BR, R_tau)", 1800, 900);
    cMainEst->Divide(4, 2);

    drawResult(cMainEst->GetPad(1), h_fit_fs);
    drawResult(cMainEst->GetPad(2), h_fit_est);
    drawResult(cMainEst->GetPad(3), h_fit_br);
    drawResult(cMainEst->GetPad(4), h_fit_rtau);

    drawResult(cMainEst->GetPad(5), h_pull_fs);
    drawResult(cMainEst->GetPad(6), h_pull_est);
    drawResult(cMainEst->GetPad(7), h_pull_br);
    drawResult(cMainEst->GetPad(8), h_pull_rtau);

    cMainEst->Update();
    if(savePlots)
    {
        cMainEst->SaveAs("./_fig/ToyMC_MainEstimators.pdf");
        cMainEst->SaveAs("./_root/ToyMC_MainEstimators.root");
    }

    // Canvas 2: Altri parametri di fit (f1, f2, f3, slope)
    TCanvas *cOtherPars = new TCanvas("cOtherPars", "Other Fit Parameters", 1600, 800);
    cOtherPars->Divide(4, 2);

    drawResult(cOtherPars->GetPad(1), h_fit_f1);
    drawResult(cOtherPars->GetPad(2), h_fit_f2);
    drawResult(cOtherPars->GetPad(3), h_fit_f3);
    drawResult(cOtherPars->GetPad(4), h_fit_slope);

    drawResult(cOtherPars->GetPad(5), h_pull_f1);
    drawResult(cOtherPars->GetPad(6), h_pull_f2);
    drawResult(cOtherPars->GetPad(7), h_pull_f3);
    drawResult(cOtherPars->GetPad(8), h_pull_slope);

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
        "Single Toy Pseudo-Data Fit;Mass [GeV/#it{c}^{2}];Entries", 200, xMin, xMax);
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

    // Estrazione dei fit per l'algoritmo FC
    ToyMetrics metrics;
    auto getMetrics = [](TH1D *h, double &mean, double &sigma)
    {
        TF1 *fit = h->GetFunction("gaus");
        if(fit)
        {
            mean = fit->GetParameter(1);
            sigma = fit->GetParameter(2);
        }
        else
        {
            cerr << "Warning: No fit found for histogram " << h->GetName() << endl;
            mean = h->GetMean();
            sigma = h->GetRMS();
        }
    };

    getMetrics(h_fit_fs, metrics.bias_fs, metrics.res_fs);
    getMetrics(h_fit_br, metrics.bias_br, metrics.res_br);
    getMetrics(h_fit_rtau, metrics.bias_rtau, metrics.res_rtau);

    cout << Form("[TOY MC] f_s: Bias=%.2e, Res=%.2e | BR: Bias=%.2e, Res=%.2e | R_tau: "
                 "Bias=%.2e, Res=%.2e",
        metrics.bias_fs, metrics.res_fs, metrics.bias_br, metrics.res_br, metrics.bias_rtau,
        metrics.res_rtau)
         << endl;

    return metrics;
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

void analysis::ConstructBelt(double sigma0, double alpha, double max_val, int mode,
    const std::vector<double> &discrete_mu, const std::vector<double> &discrete_sigma)
{
    SetLBStyle();
    constexpr double target_CL = 0.90;

    double mu_start = 0.0;
    double mu_end = max_val;
    int n_steps = 200;
    double mu_step = (mu_end - mu_start) / n_steps;

    std::vector<double> vec_mu, vec_x_lower, vec_x_upper;

    cout << "\n--> Constructing Feldman-Cousins Continuous Belt boundaries (90% CL)..." << endl;

    // --- 1. CALCOLO BANDA CONTINUA (INTERPOLATA) ---
    for(double mu = mu_start; mu <= mu_end; mu += mu_step)
    {
        double sigma_mu = sigma0 + alpha * mu;

        std::vector<FCPoint> points;
        double dx = sigma_mu / 200.0;
        double x_start = mu - 5.0 * sigma_mu;
        double x_end = mu + 5.0 * sigma_mu;

        for(double x = x_start; x <= x_end; x += dx)
        {
            double par[2] = { mu, sigma_mu };
            double prob_density = TMath::Gaus(x, mu, sigma_mu, kTRUE);
            double R = LROrdering(&x, par);

            points.push_back({ x, prob_density * dx, R });
        }
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

    // --- 2. CALCOLO LIMITI SUI PUNTI DISCRETI DEI TOY ---
    std::vector<double> vec_disc_x_lower, vec_disc_x_upper;
    bool has_discrete = (!discrete_mu.empty() && discrete_mu.size() == discrete_sigma.size());

    if(has_discrete)
    {
        cout << "--> Calculating Feldman-Cousins Discrete boundaries directly from Toy points..."
             << endl;
        for(size_t i = 0; i < discrete_mu.size(); ++i)
        {
            double mu = discrete_mu[i];
            double sig = discrete_sigma[i];

            std::vector<FCPoint> points;
            double dx = sig / 200.0;
            double x_start = mu - 5.0 * sig;
            double x_end = mu + 5.0 * sig;

            for(double x = x_start; x <= x_end; x += dx)
            {
                double par[2] = { mu, sig };
                double prob_density = TMath::Gaus(x, mu, sig, kTRUE);
                double R = LROrdering(&x, par);

                points.push_back({ x, prob_density * dx, R });
            }
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
            vec_disc_x_lower.push_back(x_min);
            vec_disc_x_upper.push_back(x_max);
        }
    }

    // --- COSTRUZIONE DEI GRAFICI ---
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

    // Grafici continui
    TGraph *g_belt = new TGraph(x_closed.size(), &x_closed[0], &mu_closed[0]);
    TGraph *g_up = new TGraph(nPoints, &vec_x_lower[0], &vec_mu[0]);
    TGraph *g_low = new TGraph(nPoints, &vec_x_upper[0], &vec_mu[0]);

    // Grafici discreti (vengono istanziati ma non disegnati sul Canvas corrente)
    TGraph *g_up_discrete = nullptr;
    TGraph *g_low_discrete = nullptr;
    if(has_discrete)
    {
        g_up_discrete = new TGraph(discrete_mu.size(), &vec_disc_x_lower[0], &discrete_mu[0]);
        g_low_discrete = new TGraph(discrete_mu.size(), &vec_disc_x_upper[0], &discrete_mu[0]);
    }

    TString titleX, titleY, nameSuffix;
    int colBase;

    if(mode == 0)
    {
        titleX = "Measured #hat{f}_{s}";
        titleY = "True f_{s}";
        nameSuffix = "fs";
        colBase = kBlue;
    }
    else if(mode == 1)
    {
        titleX = "Measured #hat{#font[12]{B}}(#tau^{+}#rightarrow#phi#mu^{+})";
        titleY = "True #font[12]{B}";
        nameSuffix = "BR";
        colBase = kGreen;
    }
    else
    {
        titleX = "Measured #hat{R}_{#tau}";
        titleY = "True R_{#tau}";
        nameSuffix = "Rtau";
        colBase = kOrange + 1;
    }

    // Stile Grafici Continui
    g_belt->SetFillColorAlpha(colBase - 9, 0.35);
    g_belt->SetLineWidth(0);
    g_up->SetLineColor(colBase + 1);
    g_low->SetLineColor(colBase + 1);
    g_up->SetLineWidth(3);
    g_low->SetLineWidth(3);

    // --- DISEGNO CANVAS PULITO (SOLO CONTURNO CONTINUO) ---
    TCanvas *cBelt
        = new TCanvas(Form("cBelt_%s", nameSuffix.Data()), "Feldman-Cousins Belt", 800, 800);
    cBelt->cd();
    cBelt->SetGrid();

    double plot_x_min = -3.0 * sigma0;
    double plot_x_max = mu_end + 3.0 * sigma0;

    TH2F *hFrame = new TH2F(Form("hFrame_%s", nameSuffix.Data()),
        Form("FC Belt 90%% CL (%s);%s;%s", nameSuffix.Data(), titleX.Data(), titleY.Data()), 100,
        plot_x_min, plot_x_max, 100, 0.0, mu_end);
    hFrame->SetStats(0);
    hFrame->Draw();

    g_belt->Draw("F SAME");
    g_up->Draw("L SAME");
    g_low->Draw("L SAME");

    TLine *diag = new TLine(0.0, 0.0, mu_end, mu_end);
    diag->SetLineStyle(2);
    diag->SetLineColor(kGray + 2);
    diag->Draw("SAME");
    TLine *vert = new TLine(0.0, 0.0, 0.0, mu_end);
    vert->SetLineStyle(3);
    vert->SetLineColor(kBlack);
    vert->Draw("SAME");

    cBelt->Update();

    // --- SALVATAGGIO REALE DEI GRAFICI INDIPENDENTI SU FILE ROOT ---
    TFile *fOut
        = TFile::Open(Form("./_root/FC_Belt_Output_%s.root", nameSuffix.Data()), "RECREATE");
    if(fOut && !fOut->IsZombie())
    {
        fOut->cd();

        // Salvataggio dei grafici continui
        g_belt->SetName("g_belt_interpolated");
        g_belt->Write();
        g_up->SetName("g_up_interpolated");
        g_up->Write();
        g_low->SetName("g_low_interpolated");
        g_low->Write();

        // Salvataggio dei grafici discreti (se presenti)
        if(has_discrete)
        {
            g_up_discrete->SetName("g_up_discrete");
            g_up_discrete->Write();
            g_low_discrete->SetName("g_low_discrete");
            g_low_discrete->Write();
        }

        cBelt->Write("cBelt_Canvas");
        fOut->Close();
        cout << "[INFO] File ROOT salvato con successo: "
             << Form("./_root/FC_Belt_Output_%s.root", nameSuffix.Data()) << endl;
    }

    if(savePlots)
    {
        cBelt->SaveAs(Form("./_fig/FeldmanCousinsBelt_%s.pdf", nameSuffix.Data()));
    }
}

void analysis::RunFeldmanCousinsPipeline(int nToysPerPoint)
{
    cout << "\n=======================================================" << endl;
    cout << "      STARTING AUTOMATED FELDMAN-COUSINS PIPELINE      " << endl;
    cout << "=======================================================" << endl;
    // --- Controllo dinamico del Fit reale ---
    if(!m_fit_done || m_fitted_pars.size() < 4)
    {
        cout << "[INFO] Parametri del fit reale non trovati per la pipeline. Esecuzione "
                "automatica "
                "di DoFullBlindedFit()..."
             << endl;
        DoFullBlindedFit();
    }

    double true_f3;
    if(m_fit_done && m_fitted_pars.size() >= 4)
    {
        true_f3 = m_fitted_pars[3];
        cout << "[INFO] Pipeline Feldman-Cousins configurata con f3 reale = " << true_f3 << endl;
    }
    else
    {
        cerr << "[ERROR] Impossibile avviare la pipeline FC: il fit automatico ha fallito." << endl;
        return; // Interrompe la pipeline
    }

    // Calcolo dinamico di r_eff basato sul constexpr globale USE_OPTIMIZED_CUTS
    double r_eff = GetMCEfficiencyRatio(USE_OPTIMIZED_CUTS);

    const double br_taunu_nom = 5.39e-2;
    const double br_phimunu_nom = 2.24e-2;
    double k_factor_nom = r_eff * (br_phimunu_nom / br_taunu_nom);

    cout << Form("[INFO] r_eff dinamico calcolato dal MC: %.4f", r_eff) << endl;

    // Definiamo i punti nominali su f_s dinamicamente con un ciclo
    std::vector<double> fs_points;

    double fs_start = 0.0;
    double fs_end = 0.01;
    double fs_step = 0.0002;
    int n_steps = std::round((fs_end - fs_start) / fs_step) + 1;
    for(int i = 0; i < n_steps; ++i)
        fs_points.push_back(fs_start + i * fs_step);

    std::vector<double> xt_fs, xt_br, xt_rtau;
    std::vector<double> yb_fs, ys_fs, yb_br, ys_br, yb_rtau, ys_rtau;

    for(double fs : fs_points)
    {
        double true_est = fs / ((1.0 - fs) * true_f3);
        double true_br = true_est * k_factor_nom;
        double true_rtau = true_est * r_eff;

        cout << Form("\n--> Initiating Toy Block for true_fs = %.4f (BR = %.2e, R_tau = %.2e)", fs,
            true_br, true_rtau)
             << endl;

        ToyMetrics m = RunToyMC(nToysPerPoint, fs); // <--- Genera una sola volta per punto!

        xt_fs.push_back(fs);
        yb_fs.push_back(m.bias_fs);
        ys_fs.push_back(m.res_fs);
        xt_br.push_back(true_br);
        yb_br.push_back(m.bias_br);
        ys_br.push_back(m.res_br);
        xt_rtau.push_back(true_rtau);
        yb_rtau.push_back(m.bias_rtau);
        ys_rtau.push_back(m.res_rtau);
    }

    // Lambda helper per elaborare le singole liste generando Fit, Control Plot e Belt
    auto processBeltMetrics
        = [&](const std::vector<double> &x, const std::vector<double> &b,
              const std::vector<double> &s, int mode, double scale, const char *name)
    {
        cout << "\n--- Processing Belts for " << name << " ---" << endl;
        double max_x = x.back();

        std::vector<double> xs(x.size()), ys(s.size());
        for(size_t i = 0; i < x.size(); ++i)
        {
            xs[i] = x[i] * scale;
            ys[i] = s[i] * scale;
        }

        TGraph *g_res = new TGraph(xs.size(), &xs[0], &ys[0]);
        TF1 *f_lin = new TF1(Form("f_lin_%d", mode), "[0] + [1]*x", 0.0, max_x * scale * 1.2);
        f_lin->SetParameters(ys[0], 0.0);
        g_res->Fit(f_lin, "Q");

        double sigma0 = f_lin->GetParameter(0) / scale;
        double alpha = f_lin->GetParameter(1);

        cout << Form("    Resolution: sigma(X) = %.5e + %.5f * X", sigma0, alpha) << endl;

        // 1) Plot Sanity Risoluzione
        TCanvas *cRes = new TCanvas(Form("cRes_%d", mode), Form("Res %s", name), 600, 500);
        cRes->SetGrid();
        TGraphErrors *ge_res = new TGraphErrors(x.size());
        for(size_t i = 0; i < x.size(); ++i)
        {
            ge_res->SetPoint(i, x[i], s[i]);
            ge_res->SetPointError(i, 0, s[i] / std::sqrt(2.0 * nToysPerPoint));
        }
        ge_res->SetMarkerStyle(20);
        ge_res->SetTitle(Form("Resolution %s", name));
        ge_res->Draw("AP");
        TF1 *f_phys = new TF1(Form("f_phys_%d", mode), "[0] + [1]*x", 0.0, max_x * 1.2);
        f_phys->SetParameters(sigma0, alpha);
        f_phys->SetLineColor(kRed);
        f_phys->Draw("SAME");
        if(savePlots)
            cRes->SaveAs(Form("./_fig/SanityCheck_Res_%s.pdf", name));

        // 2) Plot Sanity Bias
        TCanvas *cBias = new TCanvas(Form("cBias_%d", mode), Form("Bias %s", name), 600, 500);
        cBias->SetGrid();
        TGraphErrors *ge_bias = new TGraphErrors(x.size());
        for(size_t i = 0; i < x.size(); ++i)
        {
            ge_bias->SetPoint(i, x[i], b[i]);
            ge_bias->SetPointError(i, 0, s[i] / std::sqrt(nToysPerPoint));
        }
        ge_bias->SetMarkerStyle(21);
        ge_bias->SetTitle(Form("Bias %s", name));
        ge_bias->Draw("AP");
        TLine *lz = new TLine(0.0, 0.0, max_x * 1.2, 0.0);
        lz->SetLineColor(kRed);
        lz->Draw("SAME");
        if(savePlots)
            cBias->SaveAs(Form("./_fig/SanityCheck_Bias_%s.pdf", name));

        // 3) Costruzione FC
        ConstructBelt(sigma0, alpha, max_x, mode, x, s);
    };

    // Chiama la Pipeline per tutte e 3 le metriche (impostando le scale di conversione corrette
    // per Minuit)
    processBeltMetrics(xt_fs, yb_fs, ys_fs, 0, 1e3, "fs");
    processBeltMetrics(xt_br, yb_br, ys_br, 1, 1e7, "BR");
    processBeltMetrics(xt_rtau, yb_rtau, ys_rtau, 2, 1e4, "Rtau");

    cout << "\n=======================================================" << endl;
    cout << "  PIPELINE FELDMAN COUSINS (FS, BR, R_TAU) COMPLETED!  " << endl;
    cout << "=======================================================" << endl;
}

void analysis::InterpolateDiscreteBelt(
    const TString &nameSuffix, double smoothingSpan, bool overlayContinuous)
{
    // --- 1. CARICAMENTO DELLO STILE DI LORENZO ---
    lbStyle::SetLBStyle();

    TString fileName = Form("./_root/FC_Belt_Output_%s.root", nameSuffix.Data());
    TFile *fIn = TFile::Open(fileName, "READ");
    if(!fIn || fIn->IsZombie())
    {
        std::cerr << "[ERROR] Impossibile aprire il file ROOT: " << fileName << std::endl;
        return;
    }

    TGraph *g_up_discrete_raw = (TGraph *)fIn->Get("g_up_discrete");
    TGraph *g_low_discrete_raw = (TGraph *)fIn->Get("g_low_discrete");

    if(!g_up_discrete_raw || !g_low_discrete_raw)
    {
        std::cerr << "[ERROR] Grafici discreti non trovati nel file ROOT!" << std::endl;
        fIn->Close();
        return;
    }

    // --- 2. RETRIEVAL LIMITI ASSI DALLA CANVAS SALVATA ---
    double plot_x_min = 1e9;
    double plot_x_max = -1e9;
    bool foundSavedLimits = false;

    TCanvas *cSaved = (TCanvas *)fIn->Get("cBelt_Canvas");
    if(cSaved)
    {
        TObject *obj = cSaved->FindObject(Form("hFrame_%s", nameSuffix.Data()));
        if(!obj)
            obj = cSaved->GetPrimitive(Form("hFrame_%s", nameSuffix.Data()));

        if(obj && obj->InheritsFrom("TH2"))
        {
            TH2 *hFrameSaved = (TH2 *)obj;
            plot_x_min = hFrameSaved->GetXaxis()->GetXmin();
            plot_x_max = hFrameSaved->GetXaxis()->GetXmax();
            foundSavedLimits = true;
            std::cout << "[INFO] Limiti asse X sincronizzati con il grafico continuo: ["
                      << plot_x_min << ", " << plot_x_max << "]" << std::endl;
        }
    }

    // --- 3. CARICAMENTO BANDA CONTINUA, LOG DI DIAGNOSTICA E AUTO-RISCALAMENTO ---
    TGraph *g_up_interp = nullptr;
    TGraph *g_low_interp = nullptr;
    if(overlayContinuous)
    {
        TGraph *g_up_interp_raw = (TGraph *)fIn->Get("g_up_interpolated");
        TGraph *g_low_interp_raw = (TGraph *)fIn->Get("g_low_interpolated");
        if(g_up_interp_raw)
            g_up_interp = (TGraph *)g_up_interp_raw->Clone();
        if(g_low_interp_raw)
            g_low_interp = (TGraph *)g_low_interp_raw->Clone();

        if(g_up_interp && g_low_interp)
        {
            double max_x_val = -1e9;
            double max_y_val = -1e9;
            for(int i = 0; i < g_up_interp->GetN(); ++i)
            {
                double tx, ty;
                g_up_interp->GetPoint(i, tx, ty);
                if(tx > max_x_val)
                    max_x_val = tx;
                if(ty > max_y_val)
                    max_y_val = ty;
            }

            // Stima del fattore di riscalamento dinamico asimmetrico
            double factor_x = 1.0;
            if(max_x_val > 1.0)
            {
                if(max_x_val > 1000.0)
                    factor_x = 1e7;
                else if(max_x_val > 10.0)
                    factor_x = 1e3;
                else
                    factor_x = 1e4;
            }

            double factor_y = 1.0;
            if(max_y_val > 1.0)
            {
                if(max_y_val > 1000.0)
                    factor_y = 1e7;
                else if(max_y_val > 10.0)
                    factor_y = 1e3;
                else
                    factor_y = 1e4;
            }

            if(factor_x > 1.0 || factor_y > 1.0)
            {
                std::cout << "[INFO] Riscalamento dinamico rilevato -> Fattore X: " << factor_x
                          << ", Fattore Y: " << factor_y << std::endl;
                for(int i = 0; i < g_up_interp->GetN(); ++i)
                {
                    double tx_up, ty_up, tx_low, ty_low;
                    g_up_interp->GetPoint(i, tx_up, ty_up);
                    g_low_interp->GetPoint(i, tx_low, ty_low);

                    g_up_interp->SetPoint(i, tx_up / factor_x, ty_up / factor_y);
                    g_low_interp->SetPoint(i, tx_low / factor_x, ty_low / factor_y);
                }
            }

            // SCRITTURA FILE DI LOG DIAGNOSTICO SU DISCO
            std::ofstream debugFile("./_fig/FC_Belt_Debug_Points.txt");
            if(debugFile.is_open())
            {
                debugFile << "=====================================================\n";
                debugFile << "      DIAGNOSTIC LOG FOR CONTINUOUS RED CURVE POINTS \n";
                debugFile << "=====================================================\n";
                debugFile << "Name Suffix:    " << nameSuffix.Data() << "\n";
                debugFile << "Number of points: " << g_up_interp->GetN() << "\n";
                debugFile << "Max X found (before scaling): " << max_x_val << "\n";
                debugFile << "Max Y found (before scaling): " << max_y_val << "\n\n";
                debugFile << "LIST OF POINTS (after scaling applied in memory):\n";
                for(int i = 0; i < g_up_interp->GetN(); ++i)
                {
                    double x_up, y_up, x_low, y_low;
                    g_up_interp->GetPoint(i, x_up, y_up);
                    g_low_interp->GetPoint(i, x_low, y_low);
                    debugFile << Form(
                        "  Point %3d:  UPPER (x = %e, y = %e)  |  LOWER (x = %e, y = %e)\n", i,
                        x_up, y_up, x_low, y_low);
                }
                debugFile.close();
                std::cout
                    << "[DEBUG] Log delle coordinate salvato in: ./_fig/FC_Belt_Debug_Points.txt"
                    << std::endl;
            }
        }
    }

    // Cloniamo i grafici discreti per svincolarli dal file ROOT
    TGraph *g_up_discrete = (TGraph *)g_up_discrete_raw->Clone();
    TGraph *g_low_discrete = (TGraph *)g_low_discrete_raw->Clone();
    fIn->Close(); // Chiudiamo subito il file di lettura per evitare conflitti

    int nPoints = g_up_discrete->GetN();

    // --- 4. COSTRUZIONE DEI GRAFICI PARAMETRICI (Scambio X <-> Y) ---
    TGraph *g_up_param = new TGraph(nPoints);
    TGraph *g_low_param = new TGraph(nPoints);

    for(int i = 0; i < nPoints; ++i)
    {
        double x_meas_up, y_true_up;
        g_up_discrete->GetPoint(i, x_meas_up, y_true_up);
        g_up_param->SetPoint(i, y_true_up, x_meas_up);

        double x_meas_low, y_true_low;
        g_low_discrete->GetPoint(i, x_meas_low, y_true_low);
        g_low_param->SetPoint(i, y_true_low, x_meas_low);
    }

    g_up_param->Sort();
    g_low_param->Sort();

    // --- 5. LEVIGATURA LOWESS ---
    TGraphSmooth *gs = new TGraphSmooth("fc_smoother");
    TGraph *tmp_up = gs->SmoothLowess(g_up_param, "raw", smoothingSpan);
    TGraph *g_up_param_smooth = (TGraph *)tmp_up->Clone("g_up_param_smooth");

    TGraph *tmp_low = gs->SmoothLowess(g_low_param, "raw", smoothingSpan);
    TGraph *g_low_param_smooth = (TGraph *)tmp_low->Clone("g_low_param_smooth");

    // Blending morbido per raccordare i bordi ed eliminare lo scalino
    int kBlend = 5;
    if(nPoints < 15)
        kBlend = nPoints / 3;

    for(int i = 0; i < nPoints; ++i)
    {
        if(i < kBlend)
        {
            double weight = (double)i / kBlend;
            double x_raw_up, y_raw_up, x_sm_up, y_sm_up;
            g_up_param->GetPoint(i, x_raw_up, y_raw_up);
            g_up_param_smooth->GetPoint(i, x_sm_up, y_sm_up);
            g_up_param_smooth->SetPoint(i, x_sm_up, (1.0 - weight) * y_raw_up + weight * y_sm_up);

            double x_raw_low, y_raw_low, x_sm_low, y_sm_low;
            g_low_param->GetPoint(i, x_raw_low, y_raw_low);
            g_low_param_smooth->GetPoint(i, x_sm_low, y_sm_low);
            g_low_param_smooth->SetPoint(
                i, x_sm_low, (1.0 - weight) * y_raw_low + weight * y_sm_low);
        }
        else if(i >= nPoints - kBlend)
        {
            double weight = (double)(nPoints - 1 - i) / kBlend;
            double x_raw_up, y_raw_up, x_sm_up, y_sm_up;
            g_up_param->GetPoint(i, x_raw_up, y_raw_up);
            g_up_param_smooth->GetPoint(i, x_sm_up, y_sm_up);
            g_up_param_smooth->SetPoint(i, x_sm_up, (1.0 - weight) * y_raw_up + weight * y_sm_up);

            double x_raw_low, y_raw_low, x_sm_low, y_sm_low;
            g_low_param->GetPoint(i, x_raw_low, y_raw_low);
            g_low_param_smooth->GetPoint(i, x_sm_low, y_sm_low);
            g_low_param_smooth->SetPoint(
                i, x_sm_low, (1.0 - weight) * y_raw_low + weight * y_sm_low);
        }
    }

    g_up_param_smooth->Sort();
    g_low_param_smooth->Sort();

    // --- 6. SPLINE CUBICA ---
    TSpline3 spline_up("spline_up", g_up_param_smooth);
    TSpline3 spline_low("spline_low", g_low_param_smooth);

    // --- 7. CAMPIONAMENTO BANDA LISCIA ---
    const int nFine = 1000;
    TGraph *g_up_smooth = new TGraph(nFine);
    g_up_smooth->SetName("g_up_discrete_interpolated");

    TGraph *g_low_smooth = new TGraph(nFine);
    g_low_smooth->SetName("g_low_discrete_interpolated");

    double mu_min = g_up_param->GetX()[0];
    double mu_max = g_up_param->GetX()[nPoints - 1];
    double dmu = (mu_max - mu_min) / (nFine - 1);

    for(int i = 0; i < nFine; ++i)
    {
        double mu = mu_min + i * dmu;
        g_up_smooth->SetPoint(i, spline_up.Eval(mu), mu);
        g_low_smooth->SetPoint(i, spline_low.Eval(mu), mu);
    }

    // Poligono della banda colorata
    TGraph *g_belt_smooth = new TGraph(2 * nFine);
    g_belt_smooth->SetName("g_belt_discrete_interpolated");

    for(int i = 0; i < nFine; ++i)
    {
        double x, y;
        g_up_smooth->GetPoint(i, x, y);
        g_belt_smooth->SetPoint(i, x, y);
    }
    for(int i = 0; i < nFine; ++i)
    {
        double x, y;
        g_low_smooth->GetPoint(nFine - 1 - i, x, y);
        g_belt_smooth->SetPoint(nFine + i, x, y);
    }

    int colBase = kBlue;
    if(nameSuffix == "BR")
        colBase = kGreen;
    else if(nameSuffix == "Rtau")
        colBase = kOrange + 1;

    g_belt_smooth->SetFillColorAlpha(colBase - 9, 0.45);
    g_belt_smooth->SetLineWidth(0);
    g_up_smooth->SetLineColor(colBase + 1);
    g_up_smooth->SetLineWidth(3);
    g_low_smooth->SetLineColor(colBase + 1);
    g_low_smooth->SetLineWidth(3);

    // --- 8. SALVATAGGIO DEI RISULTATI SU FILE ---
    TFile *fOut = TFile::Open(fileName, "UPDATE");
    if(fOut && !fOut->IsZombie())
    {
        g_up_smooth->Write("g_up_discrete_interpolated", TObject::kOverwrite);
        g_low_smooth->Write("g_low_discrete_interpolated", TObject::kOverwrite);
        g_belt_smooth->Write("g_belt_discrete_interpolated", TObject::kOverwrite);
        fOut->Close();
    }

    // --- 9. DISEGNO CANVAS DI CONFRONTO ---
    TCanvas *cCheck = new TCanvas(
        Form("cCheck_Smooth_%s", nameSuffix.Data()), "Feldman-Cousins Smooth Belt", 800, 800);
    cCheck->cd();
    cCheck->SetGrid();

    if(!foundSavedLimits)
    {
        double x_plot_min = 1e9;
        double x_plot_max = -1e9;
        for(int i = 0; i < nFine; ++i)
        {
            double x_up_val, x_low_val, y;
            g_up_smooth->GetPoint(i, x_up_val, y);
            g_low_smooth->GetPoint(i, x_low_val, y);
            if(x_up_val < x_plot_min)
                x_plot_min = x_up_val;
            if(x_up_val > x_plot_max)
                x_plot_max = x_up_val;
            if(x_low_val < x_plot_min)
                x_plot_min = x_low_val;
            if(x_low_val > x_plot_max)
                x_plot_max = x_low_val;
        }
        plot_x_min = x_plot_min - 0.15 * (x_plot_max - x_plot_min);
        plot_x_max = x_plot_max + 0.15 * (x_plot_max - x_plot_min);
    }

    TString titleX = "Measured Value";
    if(nameSuffix == "fs")
        titleX = "Measured #hat{f}_{s}";
    else if(nameSuffix == "BR")
        titleX = "Measured #hat{#font[12]{B}}(#tau^{+}#rightarrow#phi#mu^{+})";
    else if(nameSuffix == "Rtau")
        titleX = "Measured #hat{R}_{#tau}";

    TString titleY = "True Value";
    if(nameSuffix == "fs")
        titleY = "True f_{s}";
    else if(nameSuffix == "BR")
        titleY = "True #font[12]{B}";
    else if(nameSuffix == "Rtau")
        titleY = "True R_{#tau}";

    TH2F *hFrame = new TH2F(Form("hFrame_smooth_%s", nameSuffix.Data()),
        Form(";%s;%s", titleX.Data(), titleY.Data()), 100, plot_x_min, plot_x_max, 100, 0.0,
        mu_max);

    hFrame->SetStats(0);

    // Configurazione font vettoriale 42
    hFrame->GetXaxis()->SetLabelFont(42);
    hFrame->GetXaxis()->SetLabelSize(0.04);
    hFrame->GetXaxis()->SetTitleFont(42);
    hFrame->GetXaxis()->SetTitleSize(0.045);
    hFrame->GetXaxis()->SetTitleOffset(1.2);

    hFrame->GetYaxis()->SetLabelFont(42);
    hFrame->GetYaxis()->SetLabelSize(0.04);
    hFrame->GetYaxis()->SetTitleFont(42);
    hFrame->GetYaxis()->SetTitleSize(0.045);
    hFrame->GetYaxis()->SetTitleOffset(1.5);

    hFrame->Draw();

    // 1. Disegna la banda discreta riempita
    g_belt_smooth->Draw("F SAME");
    g_up_smooth->Draw("L SAME");
    g_low_smooth->Draw("L SAME");

    // 2. Disegna i marker o sovrappone la curva continua basandosi sul booleano
    if(overlayContinuous && g_up_interp && g_low_interp)
    {
        // FIX: Linea ROSSA SOLIDA spessa (Stile 1) + Marker Cerchio vuoto rosso per sicurezza
        // assoluta
        g_up_interp->SetLineColor(kRed);
        g_up_interp->SetLineStyle(1);
        g_up_interp->SetLineWidth(4);
        g_up_interp->SetMarkerColor(kRed);
        g_up_interp->SetMarkerStyle(24); // Cerchio vuoto
        g_up_interp->SetMarkerSize(0.7);
        g_up_interp->Draw("LP SAME"); // Disegna Linea + Punti

        g_low_interp->SetLineColor(kRed);
        g_low_interp->SetLineStyle(1);
        g_low_interp->SetLineWidth(4);
        g_low_interp->SetMarkerColor(kRed);
        g_low_interp->SetMarkerStyle(24);
        g_low_interp->SetMarkerSize(0.7);
        g_low_interp->Draw("LP SAME");

        // Aggiorna la legenda mostrando la linea solida rossa
        TLegend *leg = new TLegend(0.18, 0.72, 0.58, 0.88);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->SetTextFont(42);
        leg->SetTextSize(0.035);
        leg->AddEntry(g_belt_smooth, "Discrete Belt (Toy MC + Spline)", "f");
        leg->AddEntry(g_up_interp, "Continuous Belt (Resolution Extrap.)", "l");
        leg->Draw("SAME");
    }
    else
    {
        // Mostra i marker se non facciamo il confronto
        g_up_discrete->SetMarkerStyle(20);
        g_up_discrete->SetMarkerColor(colBase + 1);
        g_up_discrete->SetMarkerSize(0.9);
        g_up_discrete->Draw("P SAME");

        g_low_discrete->SetMarkerStyle(20);
        g_low_discrete->SetMarkerColor(colBase + 1);
        g_low_discrete->SetMarkerSize(0.9);
        g_low_discrete->Draw("P SAME");
    }

    gPad->RedrawAxis();
    cCheck->Update();

    cCheck->SaveAs(Form("./_fig/FeldmanCousinsBelt_Smooth_%s.pdf", nameSuffix.Data()));

    // Cleanup delle strutture temporanee
    delete g_up_param;
    delete g_low_param;
    delete g_up_param_smooth;
    delete g_low_param_smooth;
    delete gs;
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
    null_gen_pars[0] = 0.0; // Generazione sotto H0 (f_s = 0)

    auto gen_lambda = [&localPdf, &null_gen_pars](double *x, double *par)
    {
        double xx = x[0];
        return localPdf(&xx, null_gen_pars.data());
    };
    TF1 fGen(Form("fGenWilks_%d", toyId), gen_lambda, xMin, xMax, 0);
    fGen.SetNpx(10000);

    TH1D *local_hist = nullptr; // Istogramma locale al thread

    if(USE_UNBINNED)
    {
        g_data_events.clear();
        g_data_events.reserve(nEvents);
        for(int ev = 0; ev < nEvents; ++ev)
        {
            g_data_events.push_back(fGen.GetRandom(&threadRandom));
        }
    }
    else
    {
        int nBins = std::round((xMax - xMin) / 0.0025);
        local_hist = new TH1D(Form("h_toy_wilks_%d", toyId), "", nBins, xMin, xMax);
        g_data_hist = local_hist; // Passa il puntatore a Minuit

        // FLUTTUAZIONE POISSONIANA DI N_tot (Proprietà fondamentale dell'Extended Fit)
        int nEventsFluct = USE_EXTENDED ? threadRandom.Poisson(nEvents) : nEvents;

        // Generazione esatta evento per evento
        for(int ev = 0; ev < nEventsFluct; ++ev)
        {
            double val = fGen.GetRandom(&threadRandom);
            local_hist->Fill(val);
        }
    }

    // =========================================================================
    // FIT 1: IPOTESI ALTERNATIVA (f_s COMPLETAMENTE LIBERO - UNCONSTRAINED)
    // =========================================================================
    ROOT::Math::Minimizer *minAlt = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
    minAlt->SetMaxFunctionCalls(50000);
    minAlt->SetTolerance(0.01);
    minAlt->SetPrintLevel(-1);

    int nPars = USE_EXTENDED ? 21 : 20;
    ROOT::Math::Functor fNLL(
        USE_EXTENDED ? &BinnedExtended2NLL : (USE_UNBINNED ? &Unbinned2NLL : &Binned2NLL), nPars);
    minAlt->SetFunction(fNLL);

    minAlt->SetVariable(0, "f_s", 0.0, 0.005); // f_s fluttua liberamente (anche nel negativo)

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

    if(USE_EXTENDED)
        minAlt->SetVariable(20, "N_tot", nEvents, std::sqrt(nEvents));

    minAlt->Minimize();

    if(minAlt->Status() != 0)
    {
        delete minAlt;
        if(!USE_UNBINNED && local_hist)
            delete local_hist; // <--- AGGIUNGI QUESTA
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

    if(USE_EXTENDED)
        minNull->SetVariable(20, "N_tot", nEvents, std::sqrt(nEvents));

    minNull->Minimize();

    if(minNull->Status() != 0)
    {
        delete minAlt;
        delete minNull;
        if(!USE_UNBINNED && local_hist)
            delete local_hist; // <--- AGGIUNGI QUESTA
        g_pdf_unbinned = nullptr;
        return res;
    }

    double F_0 = minNull->MinValue();

    res.converged = true;
    res.fs_fitted = fs_fitted;

    // =========================================================================
    // TRUCCO UNILATERALE (One-Sided Wilks)
    // =========================================================================
    if(fs_fitted >= 0.0)
    {
        res.delta_F = F_0 - F_min;
    }
    else
    {
        res.delta_F = 0.0; // Se f_s fittato < 0, il minimo fisico vincolato a f_s >= 0 è a 0,
                           // quindi F_0 - F_min = 0
    }

    if(res.delta_F < 0.0)
        res.delta_F = 0.0;

    delete minAlt;
    delete minNull;
    if(!USE_UNBINNED && local_hist)
        delete local_hist; // <--- AGGIUNGI QUESTA
    g_pdf_unbinned = nullptr;
    return res;
}

void analysis::VerifyWilksTheorem(int nToys)
{
    auto start = std::chrono::high_resolution_clock::now();
    ROOT::EnableThreadSafety();
    TH1::AddDirectory(kFALSE);

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

    // --- [STEP 2] Caricamento dinamico dei parametri dal fit reale ---
    if(!m_fit_done || m_fitted_pars.size() < 5)
    {
        cout << "[INFO] Parametri del fit reale non trovati per Wilks. Esecuzione automatica "
                "di DoFullBlindedFit()..."
             << endl;
        DoFullBlindedFit();
    }

    double true_f1, true_f2, true_f3, true_slope;

    if(m_fit_done && m_fitted_pars.size() >= 5)
    {
        true_f1 = m_fitted_pars[1];
        true_f2 = m_fitted_pars[2];
        true_f3 = m_fitted_pars[3];
        true_slope = m_fitted_pars[4];
        cout << "[INFO] Validazione di Wilks configurata con i parametri del fit reale:" << endl;
        cout << Form("  f1 = %.4f | f2 = %.4f | f3 = %.4f | slope = %.4f", true_f1, true_f2,
            true_f3, true_slope)
             << endl;
    }
    else
    {
        cerr << "[ERROR] Impossibile verificare il teorema di Wilks: il fit di calibrazione "
                "ha restituito parametri non validi."
             << endl;
        TF1::DefaultAddToGlobalList(kTRUE);
        return;
    }

    // fs_gen = 0.0 per generare sotto l'ipotesi nulla H0
    std::vector<double> gen_pars = { 0.0, true_f1, true_f2, true_f3, true_slope, res_sig.params[1],
        res_sig.params[2], res_sig.params[3], res_sig.params[4], res_p1.params[1], res_p1.params[2],
        res_p1.params[3], res_p1.params[4], res_p2.params[1], res_p2.params[2], res_p2.params[3],
        res_p2.params[4], res_arg.params[1], res_arg.params[2], res_arg.params[3] };

    Double_t xMin = 1.60;
    Double_t xMax = 2.10;
    FullUnbinnedPDF templatePdf(xMin, xMax, res_sig, res_p1, res_p2, res_arg, usePol1Bkg);

    // Conteggio eventi totali dai dati reali sotto i tagli attivi
    LoadDataset(0);
    int nEvents = 0;

    TString cutString = GetCutString(USE_OPTIMIZED_CUTS);
    TTreeFormula formula("cut", cutString.Data(), fChain);
    Int_t currentTreeNumber = -1;

    for(Long64_t jentry = 0; jentry < fChain->GetEntriesFast(); jentry++)
    {
        if(LoadTree(jentry) < 0)
            break;
        if(fChain->GetTreeNumber() != currentTreeNumber)
        {
            currentTreeNumber = fChain->GetTreeNumber();
            formula.UpdateFormulaLeaves();
        }
        fChain->GetEntry(jentry);
        if(id == 0 && D_M >= xMin && D_M <= xMax)
        {
            if(formula.EvalInstance() > 0)
            {
                nEvents++;
            }
        }
    }
    if(nEvents == 0)
        nEvents = 10000;

    // --- [STEP 3] Esecuzione parallela dei Toy pseudo-esperimenti ---
    int nCores = static_cast<int>(std::thread::hardware_concurrency());
    if(nCores == 0)
        nCores = 4;
    cout << "[INFO] Launching " << nToys << " One-Sided Wilks Toys on " << nCores << " cores."
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

    // --- [STEP 4] Booking e riempimento Istogramma ---
    TH1D *h_delta_F = new TH1D("h_delta_F",
        "One-Sided Wilks Theorem Check;#Delta F = F_{0} - F_{min};Toys", 50, 0.0, 10.0);
    h_delta_F->SetDirectory(nullptr);

    int convergedToys = 0;
    int zeroToys = 0;
    for(const auto &res : wilksResults)
    {
        if(res.converged)
        {
            convergedToys++;
            h_delta_F->Fill(res.delta_F);
            if(res.delta_F <= 0.0)
                zeroToys++;
        }
    }

    cout << "\n[WILKS RESULTS] Converged: " << convergedToys << " / " << nToys << endl;
    cout << Form("  Toys with Delta_F = 0 (Physical boundary): %d (%.1f%%, Expected ~50%%)",
        zeroToys, (double)zeroToys / convergedToys * 100.0)
         << endl;

    // --- [STEP 5] Plotting e Overplot Teorico Unilaterale con Residui ---
    TCanvas *cWilks
        = new TCanvas("cWilks", "One-Sided Wilks Theorem Validation with Pulls", 900, 900);
    double splitPoint = 0.30;

    TPad *pad1 = new TPad("pad1", "Main Pad", 0.0, splitPoint, 1.0, 1.0);
    pad1->SetBottomMargin(0.02);
    pad1->Draw();
    pad1->cd();
    pad1->SetGrid();

    h_delta_F->SetMinimum(0.0);
    h_delta_F->GetXaxis()->SetLabelSize(0);
    h_delta_F->GetXaxis()->SetTitleSize(0);
    h_delta_F->SetMarkerStyle(20);
    h_delta_F->SetMarkerSize(1.0);
    h_delta_F->Draw("E P");

    // 1. PDF teorica unilaterale (Chernoff) per il calcolo matematico
    TF1 *f_chi2_theory = new TF1(
        "f_chi2_theory",
        [](double *x, double *p)
        {
            if(x[0] <= 1e-6)
                return 0.0;
            return 0.5 * p[0] * std::exp(-0.5 * x[0]) / std::sqrt(2.0 * TMath::Pi() * x[0]);
        },
        1e-6, 10.0, 1);
    f_chi2_theory->SetParameter(0, convergedToys);

    // 2. Curva teorica indipendente per l'istogramma, riscalata per la larghezza del bin
    TF1 *f_chi2_draw = new TF1(
        "f_chi2_draw",
        [](double *x, double *p)
        {
            if(x[0] <= 1e-6)
                return 0.0;
            return 0.5 * p[0] * std::exp(-0.5 * x[0]) / std::sqrt(2.0 * TMath::Pi() * x[0]);
        },
        1e-6, 10.0, 1);

    // Impostiamo esplicitamente p[0] = N * BinWidth
    f_chi2_draw->SetParameter(0, convergedToys * h_delta_F->GetBinWidth(1));
    f_chi2_draw->SetNpx(5000); // Campionamento fit ad alta risoluzione
    f_chi2_draw->SetLineColor(kRed);
    f_chi2_draw->SetLineWidth(3);
    f_chi2_draw->Draw("SAME");

    TLegend *leg = new TLegend(0.45, 0.65, 0.88, 0.85);
    leg->SetBorderSize(1);
    leg->SetFillColor(kWhite);
    leg->SetTextFont(42);
    leg->SetTextSize(0.03);
    leg->AddEntry(h_delta_F, Form("Toy MC under H_{0} (%d toys)", convergedToys), "ep");
    leg->AddEntry(f_chi2_theory,
        "One-Sided Theory: #frac{1}{2}#delta(0) + #frac{1}{2}#chi^{2}_{1}(#Delta F)", "l");
    leg->Draw("SAME");

    // Pad Inferiore (Residui / Pull)
    cWilks->cd();
    TPad *pad2 = new TPad("pad2", "Pull Pad", 0.0, 0.0, 1.0, splitPoint);
    pad2->SetTopMargin(0.02);
    pad2->SetBottomMargin(0.35);
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

        double exp = 0.0;
        if(binLow <= 0.0)
        {
            // Primo bin: contiene la Delta di Dirac (50%) + l'integrale della coda da 0 a binUp
            double fraction = 0.5 + 0.5 * std::erf(std::sqrt(binUp / 2.0));
            exp = convergedToys * fraction;
        }
        else
        {
            // Bin successivi: contengono solo la coda della chi2_1 riscalata di 0.5
            double fraction
                = 0.5 * (std::erf(std::sqrt(binUp / 2.0)) - std::erf(std::sqrt(binLow / 2.0)));
            exp = convergedToys * fraction;
        }

        if(err > 0.0)
        {
            double pull = (obs - exp) / err;
            hPull->SetBinContent(i, pull);
            hPull->SetBinError(i, 0.0);
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
    hPull->GetYaxis()->SetRangeUser(-5.0, 5.0);

    hPull->GetXaxis()->SetTitle("#Delta F = F_{0} - F_{min}");
    hPull->GetXaxis()->SetTitleSize(gStyle->GetTitleSize("X") * 0.8);
    hPull->GetXaxis()->SetLabelSize(gStyle->GetLabelSize("X") * 0.8);
    hPull->GetXaxis()->SetTitleOffset(gStyle->GetTitleOffset("X") * 1.1);

    hPull->SetMarkerStyle(20);
    hPull->SetMarkerSize(0.8);
    hPull->SetMarkerColor(kBlack);
    hPull->SetLineColor(kBlack);
    hPull->Draw("P");

    TLine *line0 = new TLine(0.0, 0.0, 10.0, 0.0);
    line0->SetLineColor(kRed);
    line0->SetLineWidth(2);
    line0->Draw("SAME");

    cWilks->Update();
    cWilks->SaveAs("./_fig/SanityCheck_Wilks_OneSided.pdf");
    cWilks->SaveAs("./_root/SanityCheck_Wilks_OneSided.root");

    // =========================================================================
    // [STEP 6] CDF CON CORRETTA GESTIONE DEL FATTORE 1/2 E DEI TIES (TEST KS COMPLETO)
    // =========================================================================
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

    double max_D = 0.0;
    int N = toy_deltas.size();

    TGraph *g_cdf_emp = new TGraph(N);

    for(int i = 0; i < N; ++i)
    {
        double x = toy_deltas[i];

        // --- GESTIONE DEI "TIES" (OSSERVAZIONI RIPETUTE A ZERO) ---
        int last_idx = i;
        while(last_idx + 1 < N && std::abs(toy_deltas[last_idx + 1] - x) < 1e-9)
        {
            last_idx++;
        }

        double f_emp = (double)(last_idx + 1) / N;
        g_cdf_emp->SetPoint(i, x, f_emp);

        // CDF Teorica Unilaterale Completa (Chernoff)
        double f_theo = 0.5 + 0.5 * std::erf(std::sqrt(x / 2.0));

        double diff = std::abs(f_emp - f_theo);
        if(diff > max_D)
        {
            max_D = diff;
        }
    }

    // Calcolo corretto del p-value del test KS
    double ks_p_value = TMath::KolmogorovProb(max_D * std::sqrt(N));
    std::cout << "\n=======================================================" << endl;
    std::cout << "   KOLMOGOROV-SMIRNOV TEST RESULTS (ONE-SIDED WITH 1/2)" << endl;
    std::cout << "=======================================================" << endl;
    std::cout << Form("  Number of Toys (N):  %d", N) << endl;
    std::cout << Form("  KS Distance (d_max): %.4f", max_D) << endl;
    std::cout << Form("  KS p-value:          %.4f (%.1f%%)", ks_p_value, ks_p_value * 100.0)
              << endl;
    std::cout << "=======================================================\n" << endl;

    TCanvas *cCDF
        = new TCanvas("cCDF", "Cumulative Distribution Function & KS Test (One-Sided)", 800, 600);
    cCDF->cd();
    cCDF->SetGrid();

    TH2F *hFrameCDF = new TH2F("hFrameCDF",
        "One-Sided Wilks Validation (CDF);#Delta F = F_{0} - F_{min};Cumulative Probability", 100,
        0.0, 10.0, 100, 0.0, 1.05);
    hFrameCDF->SetStats(0);
    hFrameCDF->Draw();

    // Disegniamo la CDF teorica unilaterale tramite Lambda C++ compilata
    TF1 *f_cdf_theo = new TF1(
        "f_cdf_theo",
        [](double *x, double *p) { return 0.5 + 0.5 * std::erf(std::sqrt(x[0] / 2.0)); }, 0.0, 10.0,
        0);
    f_cdf_theo->SetLineColor(kRed);
    f_cdf_theo->SetLineWidth(3);
    f_cdf_theo->Draw("SAME");

    g_cdf_emp->SetMarkerStyle(20);
    g_cdf_emp->SetMarkerSize(0.6);
    g_cdf_emp->SetMarkerColor(kBlack);
    g_cdf_emp->SetLineColor(kBlack);
    g_cdf_emp->SetLineWidth(1);
    g_cdf_emp->Draw("P SAME");

    TLegend *legCDF = new TLegend(0.15, 0.70, 0.55, 0.85);
    legCDF->SetBorderSize(1);
    legCDF->SetFillColor(kWhite);
    legCDF->SetTextFont(42);
    legCDF->SetTextSize(0.035);
    legCDF->AddEntry(g_cdf_emp, Form("Toy MC under H_{0} (%d toys)", N), "ep");
    legCDF->AddEntry(f_cdf_theo, "One-Sided Theory: F_{#chi^{2}_{1s}}(\\Delta F)", "l");
    legCDF->Draw("SAME");

    TPaveText *paveKS = new TPaveText(0.55, 0.15, 0.88, 0.32, "NDC");
    paveKS->SetBorderSize(1);
    paveKS->SetFillColor(kWhite);
    paveKS->SetTextFont(42);
    paveKS->SetTextSize(0.035);
    paveKS->SetTextAlign(12);
    paveKS->AddText(Form("KS d_{max} = %.4f", max_D));
    paveKS->AddText(Form("KS p-value = %.4f", ks_p_value));
    paveKS->Draw();

    cCDF->Update();
    cCDF->SaveAs("./_fig/SanityCheck_Wilks_CDF_OneSided.pdf");
    cCDF->SaveAs("./_root/SanityCheck_Wilks_CDF_OneSided.root");

    std::cout << "--> Grafico CDF pronto! Chiudi la finestra del CDF per proseguire." << std::endl;
    while(gROOT->GetListOfCanvases()->FindObject("cCDF"))
    {
        gSystem->ProcessEvents();
        gSystem->Sleep(50);
    }

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    std::cout << "Wilks validation completed in " << elapsed.count() << " seconds." << std::endl;

    TF1::DefaultAddToGlobalList(kTRUE);

    delete h_delta_F;
    delete f_chi2_theory;
    delete f_chi2_draw;
    delete leg;
    delete hPull;
    delete line0;
    delete cWilks;
}

// ================================================================================
// CheckCut: Utility to test offline selection cuts on BLINDED Real Data
// ================================================================================

const TString cut_kaons = "(h1_pt > 0.450 && h2_pt > 0.450) && "
                          "(h1_p > 3.0 && h2_p > 3.0) && "
                          "(h1_eta > 2.0 && h1_eta < 4.2) && (h2_eta > 2.0 && h2_eta < 4.2) && "
                          "(h1_IP > 60e-6 && h2_IP > 60e-6)";
const TString cut_muon = "h3_pt > 0.350 && h3_p > 3.0 && "
                         "(h3_eta > 2.0 && h3_eta < 4.2) && "
                         "h3_IP > 90e-6 && h3_MuonID == 1";
const TString cut_phi = "M0_pt > 0.9";
const TString cut_parent = "D_pt > 2.5 && D_time > 0.25e-12 && (D_M > 1.6 && D_M < 2.1)";
// --- TAGLIO COMBINATO (BASELINE DEL PDF TAB. 2) ---
const TString cut_baseline = cut_kaons + " && " + cut_muon + " && " + cut_phi + " && " + cut_parent;

void analysis::CheckCut(const TString &cutString)
{
    SetLBStyle(); // Applica il tuo stile

    cout << "\n=======================================================" << endl;
    cout << "   TESTING CUT: " << cutString << endl;
    cout << "=======================================================" << endl;

    // 1. Carica i dati reali
    LoadDataset(0);
    if(!fChain)
    {
        cerr << "[ERROR] fChain is null!" << endl;
        return;
    }

    // 2. Definisci il range e la regione di blinding
    Double_t xMin = 1.60;
    Double_t xMax = 2.10;
    Double_t blindMin = 1.777 - 3 * 0.0058;
    Double_t blindMax = 1.777 + 3 * 0.0058;

    auto h_before = new TH1D(
        "h_before_cut", "Before Cut;M(D_{s}^{+}) [GeV/#it{c}^{2}];Entries", 100, xMin, xMax);
    auto h_after = new TH1D(
        "h_after_cut", "After Cut;M(D_{s}^{+}) [GeV/#it{c}^{2}];Entries", 100, xMin, xMax);

    h_before->SetDirectory(nullptr);
    h_after->SetDirectory(nullptr);
    AddBinSizeOnYTitle(h_before, "GeV/#it{c}^{2}");

    // 3. TTreeFormula per interpretare la stringa di taglio al volo
    // TTreeFormula richiede un const char*, quindi usiamo .Data() sulla nostra TString
    TTreeFormula formula("cutFormula", cutString.Data(), fChain);
    if(formula.GetNdim() == 0)
    {
        cerr << "[ERROR] Invalid cut string provided!" << endl;
        return;
    }

    // 4. Loop sugli eventi
    Int_t currentTreeNumber = -1;
    Long64_t nentries = fChain->GetEntriesFast();

    int evts_before_sideband = 0;
    int evts_after_sideband = 0;

    for(Long64_t jentry = 0; jentry < nentries; jentry++)
    {
        if(LoadTree(jentry) < 0)
            break;

        // Necessario per TTreeFormula quando si usa una TChain (file multipli)
        if(fChain->GetTreeNumber() != currentTreeNumber)
        {
            currentTreeNumber = fChain->GetTreeNumber();
            formula.UpdateFormulaLeaves();
        }

        fChain->GetEntry(jentry);

        // Assicuriamoci di guardare solo i dati reali
        if(id != 0)
            continue;

        // Limita l'analisi al range di massa di interesse
        if(D_M < xMin || D_M > xMax)
            continue;

        bool is_sideband = (D_M < blindMin || D_M > blindMax);

        // Prima del taglio
        h_before->Fill(D_M);
        if(is_sideband)
            evts_before_sideband++;

        // Valuta la condizione di taglio
        if(formula.EvalInstance() > 0)
        {
            h_after->Fill(D_M);
            if(is_sideband)
                evts_after_sideband++;
        }
    }

    // 5. Applica la maschera di Blinding
    TH1D *h_before_blind = GetBlindedClone(h_before, blindMin, blindMax);
    TH1D *h_after_blind = GetBlindedClone(h_after, blindMin, blindMax);

    // 6. Stile grafico
    h_before_blind->SetLineColor(kBlack);
    h_before_blind->SetMarkerColor(kBlack);
    h_before_blind->SetMarkerStyle(20);
    h_before_blind->SetMarkerSize(0.8);

    h_after_blind->SetLineColor(kRed);
    h_after_blind->SetMarkerColor(kRed);
    h_after_blind->SetMarkerStyle(21);
    h_after_blind->SetMarkerSize(0.8);

    // 7. Disegno
    TCanvas *cCut = new TCanvas("cCut", "Cut Effect Validation", 800, 600);
    cCut->cd();

    double max_y = h_before_blind->GetMaximum() * 1.3;
    h_before_blind->SetMaximum(max_y);
    h_before_blind->SetMinimum(0.0);
    h_before_blind->SetStats(0); // Rimuove la stat box per pulizia visiva

    h_before_blind->Draw("E");
    h_after_blind->Draw("HIST SAME"); // Istogramma colorato
    h_after_blind->Draw("E SAME"); // Più i punti con errore

    // Disegna le linee rosse per delimitare l'area blinded
    TLine *l1 = new TLine(blindMin, 0, blindMin, max_y);
    TLine *l2 = new TLine(blindMax, 0, blindMax, max_y);
    l1->SetLineStyle(2);
    l1->SetLineColor(kGray + 2);
    l2->SetLineStyle(2);
    l2->SetLineColor(kGray + 2);
    l1->Draw("SAME");
    l2->Draw("SAME");

    // 8. Legenda e metriche
    // Uso una stringa troncata per la legenda se il taglio è troppo lungo (evita che sbordi)
    TString legendCutName = cutString;
    if(legendCutName.Length() > 40)
    {
        legendCutName.Remove(37);
        legendCutName += "...";
    }

    TLegend *leg = new TLegend(0.15, 0.75, 0.55, 0.88);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextFont(42);
    leg->AddEntry(h_before_blind, "Data (Before Cut)", "lep");
    leg->AddEntry(h_after_blind, Form("Cut: %s", legendCutName.Data()), "flep");
    leg->Draw("SAME"); // Ora de-commentato e sicuro!

    double bkg_retention = (evts_before_sideband > 0)
        ? (double)evts_after_sideband / evts_before_sideband * 100.0
        : 0.0;

    TPaveText *pt = new TPaveText(0.58, 0.70, 0.92, 0.88, "NDC");
    pt->SetBorderSize(1);
    pt->SetFillColor(kWhite);
    pt->SetTextFont(42);
    pt->SetTextAlign(12);
    pt->AddText(Form("Sideband Evts (Before): %d", evts_before_sideband));
    pt->AddText(Form("Sideband Evts (After):  %d", evts_after_sideband));
    pt->AddText(Form("Bkg Retention:          %.1f%%", bkg_retention));
    pt->Draw("SAME");

    cCut->Update();

    // Pulizia stringa per il nome del file usando i metodi nativi TString
    TString safeName = cutString;
    safeName.ReplaceAll(" ", "");
    safeName.ReplaceAll(">", "GT");
    safeName.ReplaceAll("<", "LT");
    safeName.ReplaceAll("=", "EQ");
    safeName.ReplaceAll("&", "AND");
    safeName.ReplaceAll("|", "OR");
    safeName.ReplaceAll(".", "p");
    safeName.ReplaceAll("-", "m");
    safeName.ReplaceAll("(", "");
    safeName.ReplaceAll(")", "");

    // Taglia a 50 caratteri per evitare nomi file illegali per il Sistema Operativo
    if(safeName.Length() > 50)
        safeName.Remove(50);

    if(savePlots)
    {
        cCut->SaveAs(Form("./_fig/CheckCut_%s.pdf", safeName.Data()));
    }

    cout << "  Cut efficiency on Background (Sidebands): " << bkg_retention << "%" << endl;
    cout << "=======================================================\n" << endl;

    // Cleanup mem
    delete h_before;
    delete h_after;
}

double analysis::EvaluateFOM(const TString &cutString, double a_param)
{
    Double_t xMin = 1.777 - 12 * 0.0058;
    Double_t xMax = 1.777 + 12 * 0.0058;
    Double_t blindMin = 1.777 - 3 * 0.0058;
    Double_t blindMax = 1.777 + 3 * 0.0058;

    // =========================================================================
    // 1. STIMA EFFICIENZA SEGNALE (S) DAL MONTE CARLO
    // =========================================================================
    LoadDataset(1);
    TTreeFormula formulaMC("cutMC", cutString.Data(), fChain);

    int S_before = 0;
    int S_after = 0;
    Int_t currentTreeNumber = -1;

    for(Long64_t jentry = 0; jentry < fChain->GetEntriesFast(); jentry++)
    {
        if(LoadTree(jentry) < 0)
            break;
        if(fChain->GetTreeNumber() != currentTreeNumber)
        {
            currentTreeNumber = fChain->GetTreeNumber();
            formulaMC.UpdateFormulaLeaves();
        }
        fChain->GetEntry(jentry);

        if(id != 44)
            continue;

        S_before++;
        if(formulaMC.EvalInstance() > 0)
        {
            S_after++;
        }
    }

    double eff_sig = (S_before > 0) ? (double)S_after / S_before : 0.0;
    if(eff_sig <= 0.0)
        return 0.0;

    // =========================================================================
    // 2. STIMA DEL FONDO INTERPOLATO (B) DAI DATI REALI
    // =========================================================================
    LoadDataset(0);
    TTreeFormula formulaData("cutData", cutString.Data(), fChain);

    TH1D *h_mass_bkg = new TH1D("h_mass_bkg", "Mass Bkg", 50, xMin, xMax);
    h_mass_bkg->SetDirectory(nullptr);

    currentTreeNumber = -1;
    for(Long64_t jentry = 0; jentry < fChain->GetEntriesFast(); jentry++)
    {
        if(LoadTree(jentry) < 0)
            break;
        if(fChain->GetTreeNumber() != currentTreeNumber)
        {
            currentTreeNumber = fChain->GetTreeNumber();
            formulaData.UpdateFormulaLeaves();
        }
        fChain->GetEntry(jentry);

        if(id != 0)
            continue;

        if(formulaData.EvalInstance() > 0)
        {
            h_mass_bkg->Fill(D_M);
        }
    }

    // Eseguiamo il fit con la funzione accecata
    TF1 *f_bkg_fit = new TF1("f_bkg_fit", Pol2Blinded, xMin, xMax, 3);
    f_bkg_fit->SetParameters(h_mass_bkg->GetEntries() / 50.0, 0.0, 0.0);
    h_mass_bkg->Fit(f_bkg_fit, "R Q N");

    // TRUCCO RISOLUTIVO: Copiamo i parametri su una funzione analitica CONTINUA
    // per calcolare l'integrale reale senza il blocco "return 0.0"
    TF1 f_proj("f_proj", "[0] + [1]*(x-1.85) + [2]*(x-1.85)*(x-1.85)", xMin, xMax);
    f_proj.SetParameters(f_bkg_fit->GetParameters());

    double binWidth = h_mass_bkg->GetBinWidth(1);
    double integral = f_proj.Integral(blindMin, blindMax); // Integrazione corretta non nulla!

    double B_est = integral / binWidth;

    if(B_est < 0.0 || std::isnan(B_est))
    {
        B_est = 0.0;
    }

    // Cleanup memoria
    delete h_mass_bkg;
    delete f_bkg_fit;

    // =========================================================================
    // 3. CALCOLO DELLA PUNZI FOM
    // =========================================================================
    double punzi_fom = eff_sig / ((a_param / 2.0) + std::sqrt(B_est));

    return punzi_fom;
}

TGraph *analysis::ScanVariable(const TString &baseline, const TString &varFormula, double start,
    double stop, double step, double baseline_fom, double a_param)
{
    TString varName = varFormula;
    int idx = varName.Index("X");
    if(idx > 0)
        varName.Remove(idx);
    varName.ReplaceAll(">", "");
    varName.ReplaceAll("<", "");
    varName.ReplaceAll(" ", "");
    varName.ReplaceAll("(", "");
    varName.ReplaceAll(")", "");

    std::cout << "\n--- Scanning: " << varFormula.Data() << " ---" << std::endl;

    std::vector<double> x_vals;
    std::vector<double> y_foms;

    double best_fom = -1.0;
    double best_val = 0.0;

    int n_steps = std::abs((stop - start) / step) + 1;

    for(int i = 0; i < n_steps; ++i)
    {
        double current_val = start + i * step;

        TString test_cond = varFormula;
        test_cond.ReplaceAll("X",
            Form("%e",
                current_val)); // Usa notazione scientifica %e per gestire i metri (es. IP)

        TString full_cut = baseline + " && " + test_cond;

        // Calcoliamo la FOM con il nuovo metodo blindato
        double fom = EvaluateFOM(full_cut, a_param);

        x_vals.push_back(current_val);
        y_foms.push_back(fom);

        if(fom > best_fom)
        {
            best_fom = fom;
            best_val = current_val;
        }
    }

    std::cout << ">>> BEST CUT for " << varName.Data() << " is: " << best_val
              << " (FOM = " << best_fom << ")" << std::endl;

    TGraph *g_fom = new TGraph(x_vals.size(), &x_vals[0], &y_foms[0]);
    g_fom->SetName(Form("g_%s", varName.Data()));
    g_fom->SetTitle(Form("%s;Cut Value;Punzi FOM", varName.Data()));
    g_fom->SetMarkerStyle(20);
    g_fom->SetMarkerSize(1.0);
    g_fom->SetMarkerColor(kBlue + 1);
    g_fom->SetLineColor(kBlue + 1);
    g_fom->SetLineWidth(2);

    return g_fom;
}

void analysis::OptimizeParentCuts(double a_param)
{
    std::cout << "\n=======================================================" << std::endl;
    std::cout << "   STARTING OPTIMIZATION FOR SELECTED 4 PARENT VARIABLES" << std::endl;
    std::cout << Form("   PUNZI FOM SIGNIFICANCE TARGET: a = %.1f sigma", a_param) << std::endl;
    std::cout << "=======================================================" << std::endl;

    TString baseline_daughters
        = "(h1_pt > 0.450 && h2_pt > 0.450) && (h1_p > 3.0 && h2_p > 3.0) && "
          "(h1_eta > 2.0 && h1_eta < 4.2) && (h2_eta > 2.0 && h2_eta < 4.2) && "
          "(h1_IP > 60e-6 && h2_IP > 60e-6) && "
          "h3_pt > 0.350 && h3_p > 3.0 && (h3_eta > 2.0 && h3_eta < 4.2) && "
          "h3_IP > 90e-6 && h3_MuonID == 1 && "
          "M0_pt > 0.9 && (D_M > 1.6 && D_M < 2.1)";

    TString trigger_tau = baseline_daughters + " && D_pt > 2.5 && D_time > 0.25e-12";
    double baseline_fom = EvaluateFOM(trigger_tau, a_param);

    std::cout << "  Initial Trigger-Level FOM (Baseline) = " << baseline_fom << std::endl;
    std::cout << "-------------------------------------------------------" << std::endl;

    std::vector<TGraph *> graphs;

    // -- 1. D_p (Momento totale del padre) --
    graphs.push_back(ScanVariable(trigger_tau, "D_p > X", 10.0, 50.0, 5.0, baseline_fom, a_param));

    // -- 2. D_time (Tempo di volo proprio ct) --
    graphs.push_back(ScanVariable(
        trigger_tau, "D_time > X", 0.25e-12, 1.45e-12, 0.15e-12, baseline_fom, a_param));

    // -- 3. D_IP (Impact Parameter del padre - scansione inversa con passo fine) --
    // Partiamo da 150 micrometri (plateau) e scendiamo a 10 micrometri con passo fine di -10
    // micrometri
    graphs.push_back(
        ScanVariable(trigger_tau, "D_IP < X", 150e-6, 10e-6, -10e-6, baseline_fom, a_param));

    // -- 4. D_FD (Distanza di volo 3D) --
    graphs.push_back(
        ScanVariable(trigger_tau, "D_FD > X", 0.0, 0.01, 0.001, baseline_fom, a_param));

    // ==========================================================
    // DISEGNO DEL DASHBOARD RIASSUNTIVO (GRIGLIA 2x2)
    // ==========================================================
    TCanvas *c_parent_opt
        = new TCanvas("c_parent_opt", "Selected Parent Cuts Optimization", 1200, 1000);
    c_parent_opt->Divide(2, 2);

    TString titles[4] = { "p(D) Cut [GeV/#it{c}]", "#it{ct}(D) Cut [seconds]", "IP(D) Cut [meters]",
        "FD(D) Cut [meters]" };

    for(size_t i = 0; i < graphs.size(); ++i)
    {
        c_parent_opt->cd(i + 1);
        gPad->SetGrid();
        gPad->SetBottomMargin(0.15);
        gPad->SetLeftMargin(0.15);

        graphs[i]->GetXaxis()->SetTitle(titles[i].Data());
        graphs[i]->GetYaxis()->SetTitle("Punzi FOM");
        graphs[i]->Draw("APL");

        // Disegna la linea della FOM di riferimento iniziale (Baseline del trigger)
        double start = graphs[i]->GetX()[0];
        double stop = graphs[i]->GetX()[graphs[i]->GetN() - 1];
        TLine *l_base = new TLine(start, baseline_fom, stop, baseline_fom);
        l_base->SetLineStyle(2);
        l_base->SetLineColor(kGray + 2);
        l_base->SetLineWidth(2);
        l_base->Draw("SAME");
    }

    c_parent_opt->Update();

    std::cout << "\n=======================================================" << std::endl;
    std::cout << "   OPTIMIZATION DASHBOARD GENERATED SUCCESSFULLY!" << std::endl;
    std::cout << "=======================================================\n" << std::endl;
}

void analysis::StudyVertexCorrelations()
{
    SetLBStyle(); // Applica il tuo stile di pubblicazione (font 43, 26px, no titoli)
    gStyle->SetOptStat(0);
    gStyle->SetPaintTextFormat("+.2f"); // Forza la scrittura dei numeri con segno e 2 decimali

    const int nVars = 5;
    double blindMin = 1.777 - 3 * 0.0058;
    double blindMax = 1.777 + 3 * 0.0058;

    std::vector<std::vector<double>> data_Signal(nVars);
    std::vector<std::vector<double>> data_Bkg(nVars);

    // =========================================================================
    // 1. LETTURA MC SEGNALE (ID 44)
    // =========================================================================
    LoadDataset(1);
    if(!fChain)
    {
        cerr << "[ERROR] fChain nullo per il MC!" << endl;
        return;
    }

    cout << "--> Lettura in corso: MC Segnale (ID 44)..." << endl;
    Long64_t nEntriesMC = fChain->GetEntries();
    for(Long64_t jentry = 0; jentry < nEntriesMC; jentry++)
    {
        if(LoadTree(jentry) < 0)
            break;
        fChain->GetEntry(jentry);
        if(id != 44)
            continue;

        // Lettura diretta super-veloce e crash-proof delle variabili membro
        data_Signal[0].push_back(D_IP);
        data_Signal[1].push_back(D_FD);
        data_Signal[2].push_back(D_FDt);
        data_Signal[3].push_back(D_FDz);
        data_Signal[4].push_back(D_time);
    }

    // =========================================================================
    // 2. LETTURA DATI REALI (FONDO)
    // =========================================================================
    LoadDataset(0);
    if(!fChain)
    {
        cerr << "[ERROR] fChain nullo per i Dati Reali!" << endl;
        return;
    }

    cout << "--> Lettura in corso: Dati Reali Blindati (ID 0)..." << endl;
    Long64_t nEntriesData = fChain->GetEntries();
    for(Long64_t jentry = 0; jentry < nEntriesData; jentry++)
    {
        if(LoadTree(jentry) < 0)
            break;
        fChain->GetEntry(jentry);
        if(id != 0)
            continue;
        if(D_M >= blindMin && D_M <= blindMax)
            continue; // Blinding

        data_Bkg[0].push_back(D_IP);
        data_Bkg[1].push_back(D_FD);
        data_Bkg[2].push_back(D_FDt);
        data_Bkg[3].push_back(D_FDz);
        data_Bkg[4].push_back(D_time);
    }

    // =========================================================================
    // 3. CALCOLO E DISEGNO DELLE MATRICI
    // =========================================================================
    TH2D *hCorrSignal = new TH2D("hCorrSignal", "", nVars, 0, nVars, nVars, 0, nVars);
    TH2D *hCorrBkg = new TH2D("hCorrBkg", "", nVars, 0, nVars, nVars, 0, nVars);

    std::vector<TString> axisLabels
        = { "IP(D)", "FD(D)", "FD_{T}(D)", "FD_{z}(D)", "#tau(D) [ct]" };
    for(int i = 0; i < nVars; ++i)
    {
        hCorrSignal->GetXaxis()->SetBinLabel(i + 1, axisLabels[i]);
        hCorrSignal->GetYaxis()->SetBinLabel(i + 1, axisLabels[i]);
        hCorrBkg->GetXaxis()->SetBinLabel(i + 1, axisLabels[i]);
        hCorrBkg->GetYaxis()->SetBinLabel(i + 1, axisLabels[i]);
    }

    // Algoritmo di Pearson corretto (Scale-Invariant e immune a variabili con ordini di
    // grandezza minuscoli)
    auto GetPearsonCorrelation
        = [](const std::vector<double> &x, const std::vector<double> &y) -> double
    {
        if(x.empty() || x.size() != y.size())
            return 0.0;
        double sumX = 0, sumY = 0, sumXY = 0, sumX2 = 0, sumY2 = 0;
        size_t n = x.size();
        for(size_t i = 0; i < n; ++i)
        {
            sumX += x[i];
            sumY += y[i];
            sumXY += x[i] * y[i];
            sumX2 += x[i] * x[i];
            sumY2 += y[i] * y[i];
        }

        double varX = (double)n * sumX2 - sumX * sumX;
        double varY = (double)n * sumY2 - sumY * sumY;

        if(varX <= 0.0 || varY <= 0.0)
            return 0.0;

        double num = (double)n * sumXY - sumX * sumY;
        double den = std::sqrt(varX * varY);
        return num / den;
    };

    for(int i = 0; i < nVars; ++i)
    {
        for(int j = 0; j < nVars; ++j)
        {
            hCorrSignal->SetBinContent(
                i + 1, j + 1, GetPearsonCorrelation(data_Signal[i], data_Signal[j]));
            hCorrBkg->SetBinContent(i + 1, j + 1, GetPearsonCorrelation(data_Bkg[i], data_Bkg[j]));
        }
    }

    TCanvas *cCorr = new TCanvas("cCorr", "Correlations", 1500, 700);
    cCorr->Divide(2, 1);

    hCorrSignal->GetZaxis()->SetRangeUser(-1.0, 1.0);
    hCorrSignal->SetMarkerSize(1.5);
    hCorrSignal->SetMarkerColor(kBlack);

    hCorrBkg->GetZaxis()->SetRangeUser(-1.0, 1.0);
    hCorrBkg->SetMarkerSize(1.5);
    hCorrBkg->SetMarkerColor(kBlack);

    // Trucco temporaneo del font di testo a precisione-2 per disegnare correttamente i numeri
    gStyle->SetTextFont(42);

    cCorr->cd(1);
    gPad->SetLeftMargin(0.22);
    gPad->SetBottomMargin(0.18);
    gPad->SetRightMargin(0.15);
    hCorrSignal->Draw("COLZ TEXT");

    cCorr->cd(2);
    gPad->SetLeftMargin(0.22);
    gPad->SetBottomMargin(0.18);
    gPad->SetRightMargin(0.15);
    hCorrBkg->Draw("COLZ TEXT");

    cCorr->Update();

    if(savePlots)
    {
        cCorr->SaveAs("./_fig/ParentVertex_Correlations.pdf");
        cCorr->SaveAs("./_root/ParentVertex_Correlations.root");
    }

    // Ripristino del font in pixel originale
    gStyle->SetTextFont(43);

    cout << "--> Analisi delle correlazioni completata con successo!" << endl;
}

void analysis::DrawBlindedBkgFit(const TString &cutString, double a_param)
{
    SetLBStyle(); // Applica il tuo stile di pubblicazione
    gStyle->SetOptFit(0); // Evitiamo la stat box standard per non affollare il plot

    Double_t xMin = 1.777 - 12 * 0.0058;
    Double_t xMax = 1.777 + 12 * 0.0058;
    Double_t blindMin = 1.777 - 3 * 0.0058; // 1.7596 GeV
    Double_t blindMax = 1.777 + 3 * 0.0058; // 1.7944 GeV

    // =========================================================================
    // 1. CARICAMENTO DATI E RIEMPIMENTO ISTOGRAMMA
    // =========================================================================
    LoadDataset(0); // Dati reali
    if(!fChain)
    {
        cerr << "[ERROR] fChain nullo!" << endl;
        return;
    }

    TTreeFormula formulaData("cutData", cutString.Data(), fChain);
    TH1D *h_mass_raw
        = new TH1D("h_mass_raw", ";Invariant Mass [GeV/#it{c}^{2}];Entries", 50, xMin, xMax);
    h_mass_raw->SetDirectory(nullptr);
    AddBinSizeOnYTitle(h_mass_raw, "GeV/#it{c}^{2}");

    Int_t currentTreeNumber = -1;
    for(Long64_t jentry = 0; jentry < fChain->GetEntriesFast(); jentry++)
    {
        if(LoadTree(jentry) < 0)
            break;
        if(fChain->GetTreeNumber() != currentTreeNumber)
        {
            currentTreeNumber = fChain->GetTreeNumber();
            formulaData.UpdateFormulaLeaves();
        }
        fChain->GetEntry(jentry);
        if(id != 0)
            continue;

        if(formulaData.EvalInstance() > 0)
        {
            h_mass_raw->Fill(D_M);
        }
    }

    // =========================================================================
    // 2. ACCECAMENTO RIGOROSO DEI DATI DA DISEGNARE
    // =========================================================================
    TH1D *h_draw = (TH1D *)h_mass_raw->Clone("h_draw");
    for(int i = 1; i <= h_draw->GetNbinsX(); ++i)
    {
        double center = h_draw->GetBinCenter(i);
        if(center > blindMin && center < blindMax)
        {
            h_draw->SetBinContent(i, 0.0);
            h_draw->SetBinError(i, 0.0); // Nessun punto dati mostrato nella finestra del segnale!
        }
    }

    // =========================================================================
    // 3. FIT SULLE SIDEBAND
    // =========================================================================
    TF1 *f_bkg_fit = new TF1("f_bkg_fit", Pol2Blinded, xMin, xMax, 3);
    f_bkg_fit->SetParameters(h_draw->GetEntries() / 50.0, 0.0, 0.0);
    h_draw->Fit(f_bkg_fit, "R Q N"); // Fit silenzioso

    double p0 = f_bkg_fit->GetParameter(0);
    double p1 = f_bkg_fit->GetParameter(1);
    double p2 = f_bkg_fit->GetParameter(2);

    // Calcoliamo B_est
    double binWidth = h_draw->GetBinWidth(1);
    double integral = f_bkg_fit->Integral(blindMin, blindMax);
    double B_est = integral / binWidth;

    // =========================================================================
    // 4. CREAZIONE CURVE PER IL PLOT
    // =========================================================================
    // La funzione di fit originale (disegnata solida) mostrerà un vuoto nella regione blindata
    f_bkg_fit->SetLineColor(kBlue + 1);
    f_bkg_fit->SetLineWidth(3);
    f_bkg_fit->SetNpx(1000);

    // Creiamo una funzione continua che proietta il fit dentro la blinding region (linea
    // tratteggiata)
    TF1 *f_projection
        = new TF1("f_projection", "[0] + [1]*(x-1.85) + [2]*(x-1.85)*(x-1.85)", xMin, xMax);
    f_projection->SetParameters(p0, p1, p2);
    f_projection->SetLineColor(kRed);
    f_projection->SetLineStyle(2); // Tratteggiato
    f_projection->SetLineWidth(3);
    f_projection->SetNpx(1000);

    // =========================================================================
    // 5. DISEGNO SUL CANVAS
    // =========================================================================
    TCanvas *cDiag = new TCanvas("cDiag", "Diagnostic Sideband Fit", 800, 600);
    cDiag->cd();
    cDiag->SetGrid();

    double max_y = h_draw->GetMaximum() * 1.3;
    h_draw->SetMaximum(max_y);
    h_draw->SetMinimum(0.0);
    h_draw->Draw("E"); // Disegna i dati (esclusa la blinding region)

    // Disegniamo la proiezione tratteggiata rossa (sotto) e il fit solido blu (sopra)
    f_projection->Draw("SAME");
    f_bkg_fit->Draw("SAME");

    // Disegniamo le linee verticali di blinding
    TLine *lineL = new TLine(blindMin, 0, blindMin, max_y);
    TLine *lineH = new TLine(blindMax, 0, blindMax, max_y);
    lineL->SetLineStyle(2);
    lineL->SetLineColor(kGray + 2);
    lineL->SetLineWidth(2);
    lineH->SetLineStyle(2);
    lineH->SetLineColor(kGray + 2);
    lineH->SetLineWidth(2);
    lineL->Draw("SAME");
    lineH->Draw("SAME");

    // Scatola informativa con i risultati
    TPaveText *pt = new TPaveText(0.18, 0.65, 0.52, 0.88, "NDC");
    pt->SetBorderSize(1);
    pt->SetFillColor(kWhite);
    pt->SetTextFont(42);
    pt->SetTextAlign(12);
    pt->AddText(Form("Fit Range: [%.2f, %.2f] GeV", xMin, xMax));
    pt->AddText(Form("Blinded Window: [%.4f, %.4f] GeV", blindMin, blindMax));
    pt->AddText(Form("#chi^{2} / ndf = %.1f / %d", f_bkg_fit->GetChisquare(), f_bkg_fit->GetNDF()));
    pt->AddText(Form("Est. Background B = %.1f events", B_est));
    pt->Draw("SAME");

    TLegend *leg = new TLegend(0.55, 0.70, 0.88, 0.88);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextFont(42);
    leg->AddEntry(h_draw, "Data (Sidebands Only)", "ep");
    leg->AddEntry(f_bkg_fit, "Sideband Fit (Pol2)", "l");
    leg->AddEntry(f_projection, "B Projection in Blind Region", "l");
    leg->Draw("SAME");

    cDiag->Update();

    if(savePlots)
    {
        cDiag->SaveAs("./_fig/Diagnostic_SidebandFit.pdf");
    }

    // Cleanup
    delete h_mass_raw;
}

TString analysis::GetCutString(bool useOptimized)
{
    if(useOptimized)
    {
        // INTERRUTTORE ACCESO: Baseline + tagli ottimizzati del padre D
        return cut_baseline + " && D_IP <= 70e-6";
    }
    else
    {
        // INTERRUTTORE SPENTO: Solo i tagli di baseline originali (Tabella 2)
        return cut_baseline;
    }
}

// ----------- UNBLINDING -----------
void analysis::UnblindResults()
{
    // Applica lo stile globale di Lorenzo
    SetLBStyle();
    gStyle->SetOptFit(0);

    constexpr bool usePol1Bkg = false;
    Double_t xMin = 1.6;
    Double_t xMax = 2.10;

    cout << "\n=======================================================" << endl;
    cout << "          PERFORMING FINAL UNBLINDED ANALYSIS          " << endl;
    cout << "=======================================================" << endl;

    // =========================================================================
    // STEP 1: CARICAMENTO DATI REALI SBLINDATI
    // =========================================================================
    LoadDataset(0);
    if(!fChain)
    {
        cerr << "[ERROR] fChain is null!" << endl;
        return;
    }

    g_data_events.clear();
    auto h_mass
        = new TH1D("h_mass_unblinded", ";Invariant Mass [GeV/#it{c}^{2}];Entries", 200, xMin, xMax);
    h_mass->SetDirectory(nullptr);

    TString cutString = GetCutString(USE_OPTIMIZED_CUTS);
    TTreeFormula formula("cut", cutString.Data(), fChain);
    Int_t currentTreeNumber = -1;

    for(Long64_t jentry = 0; jentry < fChain->GetEntries(); jentry++)
    {
        if(LoadTree(jentry) < 0)
            break;
        if(fChain->GetTreeNumber() != currentTreeNumber)
        {
            currentTreeNumber = fChain->GetTreeNumber();
            formula.UpdateFormulaLeaves();
        }
        fChain->GetEntry(jentry);
        if(id != 0)
            continue;

        if(formula.EvalInstance() > 0)
        {
            if(D_M >= xMin && D_M <= xMax)
            {
                g_data_events.push_back(D_M);
                h_mass->Fill(D_M);
            }
        }
    }
    g_data_hist = h_mass;
    double nEntries = h_mass->GetEntries();
    double binWidth = h_mass->GetBinWidth(1);

    // Caricamento dei fit ausiliari MC per vincolare le forme dei background fisici
    LoadDataset(1);
    AuxFitResult res_sig = FitTemplateMass(44);
    AuxFitResult res_p1 = FitTemplateMass(34);
    AuxFitResult res_p2 = FitTemplateMass(41);
    AuxFitResult res_arg = FitTemplateMass(42);
    LoadDataset(0); // Ritorno ai dati reali

    g_pdf_unbinned = new FullUnbinnedPDF(xMin, xMax, res_sig, res_p1, res_p2, res_arg, usePol1Bkg);

    // =========================================================================
    // STEP 2: FIT MASSIMO SVINCOLATO (H1 - SEGNALE LIBERO)
    // =========================================================================
    ROOT::Math::Minimizer *minH1 = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
    minH1->SetMaxFunctionCalls(100000);
    minH1->SetTolerance(0.01);
    minH1->SetPrintLevel(0);

    int nPars = USE_EXTENDED ? 21 : 20;
    ROOT::Math::Functor fNLL(
        USE_EXTENDED ? &BinnedExtended2NLL : (USE_UNBINNED ? &Unbinned2NLL : &Binned2NLL), nPars);
    minH1->SetFunction(fNLL);

    minH1->SetVariable(0, "f_s", 0.002, 0.0005);
    minH1->SetVariable(1, "f_1", 0.021, 0.001);
    minH1->SetVariable(2, "f_2", 0.043, 0.001);
    minH1->SetVariable(3, "f_3", 0.656, 0.001);
    minH1->SetVariable(4, "slope", -0.93, 0.01);

    // Impostazione dei parametri di forma congelati ai fit ausiliari MC
    minH1->SetFixedVariable(5, "sig_mean", res_sig.params[1]);
    minH1->SetFixedVariable(6, "sig_sigma1", res_sig.params[2]);
    minH1->SetFixedVariable(7, "sig_sigma2", res_sig.params[3]);
    minH1->SetFixedVariable(8, "sig_frac1", res_sig.params[4]);

    minH1->SetFixedVariable(9, "p1_mean", res_p1.params[1]);
    minH1->SetFixedVariable(10, "p1_sigma1", res_p1.params[2]);
    minH1->SetFixedVariable(11, "p1_sigma2", res_p1.params[3]);
    minH1->SetFixedVariable(12, "p1_frac1", res_p1.params[4]);

    minH1->SetFixedVariable(13, "p2_mean", res_p2.params[1]);
    minH1->SetFixedVariable(14, "p2_sigma1", res_p2.params[2]);
    minH1->SetFixedVariable(15, "p2_sigma2", res_p2.params[3]);
    minH1->SetFixedVariable(16, "p2_frac1", res_p2.params[4]);

    minH1->SetFixedVariable(17, "argus_m0", res_arg.params[1]);
    minH1->SetFixedVariable(18, "argus_c", res_arg.params[2]);
    minH1->SetFixedVariable(19, "argus_p", res_arg.params[3]);

    if(USE_EXTENDED)
    {
        minH1->SetVariable(20, "N_tot", nEntries, std::sqrt(nEntries));
    }

    minH1->Minimize();
    minH1->Hesse();

    double nll_H1 = minH1->MinValue();
    const double *xs = minH1->X();
    const double *errs = minH1->Errors();

    double fit_fs = xs[0];
    double fit_fs_err = errs[0];
    double fit_f3 = xs[3];
    double fit_f3_err = errs[3];
    double cov_fs_f3 = minH1->CovMatrix(0, 3);

    // =========================================================================
    // STEP 3: FIT VINCOLATO (H0 - IPOTESI NULLA f_s = 0)
    // =========================================================================
    ROOT::Math::Minimizer *minH0 = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
    minH0->SetMaxFunctionCalls(100000);
    minH0->SetTolerance(0.01);
    minH0->SetPrintLevel(0);
    minH0->SetFunction(fNLL);

    minH0->SetFixedVariable(0, "f_s", 0.0);

    for(unsigned int p = 1; p < minH1->NDim(); ++p)
    {
        minH0->SetVariable(p, minH1->VariableName(p), xs[p], errs[p]);
        if(minH1->IsFixedVariable(p))
        {
            minH0->FixVariable(p);
        }
    }

    minH0->Minimize();
    double nll_H0 = minH0->MinValue();

    // =========================================================================
    // STEP 4: CALCOLO DEI PUNTI DI STIMA E DELLE INCERTEZZE
    // =========================================================================
    double r_eff = GetMCEfficiencyRatio(USE_OPTIMIZED_CUTS);
    const double br_taunu_nom = 5.39e-2;
    const double br_taunu_err = 0.09e-2;
    const double br_phimunu_nom = 2.24e-2;
    const double br_phimunu_err = 0.11e-2;
    double k_factor_nom = r_eff * (br_phimunu_nom / br_taunu_nom);

    double fit_Y = fit_fs / ((1.0 - fit_fs) * fit_f3);

    // Propagazione degli errori su Y (Delta Method con covarianza)
    double dF_dfs = 1.0 / ((1.0 - fit_fs) * (1.0 - fit_fs) * fit_f3);
    double dF_df3 = -fit_fs / ((1.0 - fit_fs) * fit_f3 * fit_f3);
    double var_Y = (dF_dfs * dF_dfs * fit_fs_err * fit_fs_err)
        + (dF_df3 * dF_df3 * fit_f3_err * fit_f3_err) + (2.0 * dF_dfs * dF_df3 * cov_fs_f3);
    double fit_Y_err = (var_Y > 0.0) ? std::sqrt(var_Y) : 0.0;

    double fit_BR = fit_Y * k_factor_nom;
    double fit_Rtau = fit_Y * r_eff;

    // Propagazione degli errori sistematici esterni per il BR
    double rel_err_k2 = (br_phimunu_err / br_phimunu_nom) * (br_phimunu_err / br_phimunu_nom)
        + (br_taunu_err / br_taunu_nom) * (br_taunu_err / br_taunu_nom);
    double err_k = k_factor_nom * std::sqrt(rel_err_k2);
    double fit_BR_err
        = std::sqrt((k_factor_nom * k_factor_nom * var_Y) + (fit_Y * fit_Y * err_k * err_k));
    double fit_Rtau_err = fit_Rtau * (fit_Y_err / fit_Y);

    // =========================================================================
    // STEP 5: CALCOLO DEL p-value E SIGNIFICATIVITÀ CON FATTORE 1/2
    // =========================================================================
    double delta_nll = nll_H0 - nll_H1;
    if(delta_nll < 0.0)
        delta_nll = 0.0;

    double p_value = 1.0;
    double significance = 0.0;

    if(fit_fs > 0.0)
    {
        // p-value ad una coda con barriera fisica (Teorema di Chernoff)
        p_value = 0.5 * TMath::Prob(delta_nll, 1);
        significance = TMath::NormQuantile(1.0 - p_value);
    }
    else
    {
        p_value = 0.5;
        significance = 0.0;
    }

    // =========================================================================
    // STEP 6: ESTRAZIONE DEI LIMITI DA ENTRAMBE LE BANDE (CONTINUE E DISCRETE)
    // =========================================================================
    auto getFCIntervalsAll
        = [](const TString &suffix, double observed_val, double &low_cont, double &up_cont,
              double &low_disc, double &up_disc) -> std::pair<bool, bool>
    {
        TString fName = Form("./_root/FC_Belt_Output_%s.root", suffix.Data());
        TFile *f = TFile::Open(fName, "READ");
        if(!f || f->IsZombie())
        {
            return { false, false };
        }

        bool ok_cont = false;
        TGraph *g_up_c = (TGraph *)f->Get("g_up_interpolated");
        TGraph *g_low_c = (TGraph *)f->Get("g_low_interpolated");
        if(g_up_c && g_low_c)
        {
            up_cont = g_up_c->Eval(observed_val);
            low_cont = g_low_c->Eval(observed_val);
            if(low_cont < 0.0 || std::isnan(low_cont))
                low_cont = 0.0;
            if(up_cont < 0.0 || std::isnan(up_cont))
                up_cont = 0.0;
            ok_cont = true;
        }

        bool ok_disc = false;
        TGraph *g_up_d = (TGraph *)f->Get("g_up_discrete_interpolated");
        TGraph *g_low_d = (TGraph *)f->Get("g_low_discrete_interpolated");
        if(g_up_d && g_low_d)
        {
            up_disc = g_up_d->Eval(observed_val);
            low_disc = g_low_d->Eval(observed_val);
            if(low_disc < 0.0 || std::isnan(low_disc))
                low_disc = 0.0;
            if(up_disc < 0.0 || std::isnan(up_disc))
                up_disc = 0.0;
            ok_disc = true;
        }

        f->Close();
        return { ok_cont, ok_disc };
    };

    double fs_low_c = 0.0, fs_up_c = 0.0, fs_low_d = 0.0, fs_up_d = 0.0;
    double BR_low_c = 0.0, BR_up_c = 0.0, BR_low_d = 0.0, BR_up_d = 0.0;
    double Rtau_low_c = 0.0, Rtau_up_c = 0.0, Rtau_low_d = 0.0, Rtau_up_d = 0.0;

    auto res_fs = getFCIntervalsAll("fs", fit_fs, fs_low_c, fs_up_c, fs_low_d, fs_up_d);
    auto res_BR = getFCIntervalsAll("BR", fit_BR, BR_low_c, BR_up_c, BR_low_d, BR_up_d);
    auto res_Rtau
        = getFCIntervalsAll("Rtau", fit_Rtau, Rtau_low_c, Rtau_up_c, Rtau_low_d, Rtau_up_d);

    // =========================================================================
    // STEP 7: LOG OUTPUT DETTAGLIATO (CON FORMATTAZIONE PDG ALLA SCALA 10^-4)
    // =========================================================================
    cout << "\n======================================================================" << endl;
    cout << "                    FINAL UNBLINDED PHYSICS RESULTS                    " << endl;
    cout << "======================================================================" << endl;
    cout << "  Fitted Signal Fraction (f_s) : " << FormatPDG(fit_fs, fit_fs_err) << endl;
    cout << "  Fitted Normalization (f_3)   : " << FormatPDG(fit_f3, fit_f3_err) << endl;
    cout << Form("  Profile Likelihood Ratio t_0 : %.4f (NLL_H0 = %.2f, NLL_H1 = %.2f)", delta_nll,
        nll_H0, nll_H1)
         << endl;
    cout << Form(
        "  Corrected One-Sided p-value  : %e (Chernoff Mixture 1/2 delta + 1/2 chi2)", p_value)
         << endl;
    cout << Form("  Signal Significance (Z)     : %.2f standard deviations (sigma)", significance)
         << endl;
    cout << "----------------------------------------------------------------------" << endl;

    // Risultati Branching Ratio formattati PDG alla scala naturale di 10^-4
    cout << "  8. RESULTS FOR B(tau+ -> phi mu+):" << endl;
    cout << "     - Point Estimate:  " << FormatPDG(fit_BR * 1e4, fit_BR_err * 1e4, "10^{-4}")
         << endl;

    if(res_BR.first)
    {
        cout << "     [Continuous Belt (Resolution Extrap.) 90% CL]:" << endl;
        if(BR_low_c <= 0.0)
        {
            cout << Form("       - Upper Limit: B < %.2f x 10^-4", BR_up_c * 1e4) << endl;
        }
        else
        {
            cout << Form(
                "       - Two-Sided Interval: [%.2f, %.2f] x 10^-4", BR_low_c * 1e4, BR_up_c * 1e4)
                 << endl;
        }
    }
    if(res_BR.second)
    {
        cout << "     [Discrete Smoothed Belt (Toy MC + Spline) 90% CL]:" << endl;
        if(BR_low_d <= 0.0)
        {
            cout << Form("       - Upper Limit: B < %.2f x 10^-4", BR_up_d * 1e4) << endl;
        }
        else
        {
            cout << Form(
                "       - Two-Sided Interval: [%.2f, %.2f] x 10^-4", BR_low_d * 1e4, BR_up_d * 1e4)
                 << endl;
        }
    }
    cout << "----------------------------------------------------------------------" << endl;

    // Risultati R_tau formattati PDG riscalando a 10^-4
    cout << "  9. RESULTS FOR R_tau:" << endl;
    cout << "     - Point Estimate:  " << FormatPDG(fit_Rtau * 1e4, fit_Rtau_err * 1e4, "10^{-4}")
         << endl;

    if(res_Rtau.first)
    {
        cout << "     [Continuous Belt (Resolution Extrap.) 90% CL]:" << endl;
        if(Rtau_low_c <= 0.0)
        {
            cout << Form("       - Upper Limit: R_tau < %.2f x 10^-4", Rtau_up_c * 1e4) << endl;
        }
        else
        {
            cout << Form("       - Two-Sided Interval: [%.2f, %.2f] x 10^-4", Rtau_low_c * 1e4,
                Rtau_up_c * 1e4)
                 << endl;
        }
    }
    if(res_Rtau.second)
    {
        cout << "     [Discrete Smoothed Belt (Toy MC + Spline) 90% CL]:" << endl;
        if(Rtau_low_d <= 0.0)
        {
            cout << Form("       - Upper Limit: R_tau < %.2f x 10^-4", Rtau_up_d * 1e4) << endl;
        }
        else
        {
            cout << Form("       - Two-Sided Interval: [%.2f, %.2f] x 10^-4", Rtau_low_d * 1e4,
                Rtau_up_d * 1e4)
                 << endl;
        }
    }
    cout << "======================================================================" << endl;

    // =========================================================================
    // STEP 8: PREPARAZIONE DELLE COPIE PROFONDE PER LE LAMBDA (RISOLVE IL CRASH)
    // =========================================================================
    FullUnbinnedPDF pdf_val = *g_pdf_unbinned;
    std::vector<double> xs_val(xs, xs + minH1->NDim());

    // =========================================================================
    // STEP 9: DISEGNO DEL PLOT SBLINDATO (STILE LORENZO)
    // =========================================================================
    // Nota: creiamo la canvas senza argomenti di dimensione per ereditare l'800x600 di lbStyle
    TCanvas *c_unblind = new TCanvas("c_unblind", "Unblinded Final Maximum Likelihood Fit");
    double splitPoint = 0.30;

    // Pad Superiore (Fit di Massa)
    TPad *pad1 = new TPad("pad1", "Main Fit Pad", 0.0, splitPoint, 1.0, 1.0);
    pad1->SetLeftMargin(0.16);
    pad1->SetRightMargin(0.08);
    pad1->SetTopMargin(0.08);
    pad1->SetBottomMargin(0.02); // Congiunge perfettamente con il pad inferiore
    pad1->Draw();
    pad1->cd();

    h_mass->GetYaxis()->SetTitle("Entries");
    AddBinSizeOnYTitle(h_mass, "GeV/#it{c}^{2}", "Entries");

    // Rimuove etichette asse X dal pad principale
    h_mass->GetXaxis()->SetLabelSize(0);
    h_mass->GetXaxis()->SetTitleSize(0);
    h_mass->SetMinimum(0.0);
    h_mass->SetMarkerStyle(20);
    h_mass->SetMarkerSize(0.9);
    h_mass->Draw("E");

    // Calcolo del fattore di scala per la normalizzazione
    double nTotFit = USE_EXTENDED ? xs_val[20] : nEntries;

    // Lambda che catturano per VALORE (Copia profonda sicura)
    auto scale_pdf_lambda = [binWidth, xs_val, nTotFit, pdf_val](double *x, double *par)
    {
        double xx = x[0];
        return nTotFit * binWidth * pdf_val(&xx, xs_val.data());
    };

    auto scale_sig_lambda = [binWidth, xs_val, nTotFit, pdf_val](double *x, double *par)
    {
        double xx = x[0];
        double p_sig = pdf_val.Eval2G(xx, xs_val[5], xs_val[6], xs_val[7], xs_val[8]);
        return nTotFit * binWidth * (xs_val[0] * p_sig);
    };

    auto scale_comb_lambda
        = [binWidth, xs_val, nTotFit, pdf_val, usePol1Bkg](double *x, double *par)
    {
        double xx = x[0];
        double p_4 = usePol1Bkg ? pdf_val.EvalPol1(xx, xs_val[4]) : pdf_val.EvalExpo(xx, xs_val[4]);
        double frac_comb = (1.0 - xs_val[0]) * (1.0 - xs_val[1] - xs_val[2] - xs_val[3]);
        return nTotFit * binWidth * (frac_comb * p_4);
    };

    auto scale_p1_lambda = [binWidth, xs_val, nTotFit, pdf_val](double *x, double *par)
    {
        double xx = x[0];
        double p_1 = pdf_val.Eval2G(xx, xs_val[9], xs_val[10], xs_val[11], xs_val[12]);
        return nTotFit * binWidth * ((1.0 - xs_val[0]) * xs_val[1] * p_1);
    };

    auto scale_p2_lambda = [binWidth, xs_val, nTotFit, pdf_val](double *x, double *par)
    {
        double xx = x[0];
        double p_2 = pdf_val.Eval2G(xx, xs_val[13], xs_val[14], xs_val[15], xs_val[16]);
        return nTotFit * binWidth * ((1.0 - xs_val[0]) * xs_val[2] * p_2);
    };

    auto scale_argus_lambda = [binWidth, xs_val, nTotFit, pdf_val](double *x, double *par)
    {
        double xx = x[0];
        double p_3 = pdf_val.EvalArgus(xx, xs_val[17], xs_val[18], xs_val[19]);
        return nTotFit * binWidth * ((1.0 - xs_val[0]) * xs_val[3] * p_3);
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

    TF1 *f_arg = new TF1("f_arg", scale_argus_lambda, xMin, xMax, 0);
    f_arg->SetNpx(10000);
    f_arg->SetLineColor(kCyan + 1);
    f_arg->SetLineStyle(5);
    f_arg->SetLineWidth(2);

    f_draw->Draw("SAME");
    f_sig->Draw("SAME");
    f_comb->Draw("SAME");
    f_p1->Draw("SAME");
    f_p2->Draw("SAME");
    f_arg->Draw("SAME");

    // =========================================================================
    // CALCOLO DEL CHI2/NDOF SUI DATI SBLINDATI (USATO PER IL PAVETEXT)
    // =========================================================================
    double chi2 = 0.0;
    int nBinsUsed = 0;
    for(int i = 1; i <= h_mass->GetNbinsX(); i++)
    {
        double x = h_mass->GetBinCenter(i);
        if(x < xMin || x > xMax)
            continue;

        double obs = h_mass->GetBinContent(i);
        double err = h_mass->GetBinError(i);
        if(err > 0)
        {
            double val = f_draw->Eval(h_mass->GetBinCenter(i));
            chi2 += (obs - val) * (obs - val) / (err * err);
            nBinsUsed++;
        }
    }
    int ndf = nBinsUsed - minH1->NFree();

    // Legenda trasparente con font assoluto 43 a 20 pixel (senza override manuali non necessari)
    TLegend *leg = new TLegend(0.48, 0.40, 0.92, 0.88);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextFont(43);
    leg->SetTextSize(20);
    leg->AddEntry(h_mass, "Data (Unblinded)", "ep");
    leg->AddEntry(f_draw, "Total Fit", "l");
    leg->AddEntry(f_sig, "Signal (#tau^{+}#rightarrow#phi#mu^{+})", "l");
    leg->AddEntry(f_p1, "D^{+}#rightarrow#phi#pi^{+} Bkg", "l");
    leg->AddEntry(f_p2, "D_{s}^{+}#rightarrow#phi#pi^{+} Bkg", "l");
    leg->AddEntry(f_arg, "D_{s}^{+}#rightarrow#phi#mu^{+}#nu_{#mu} (Argus) Bkg", "l");
    leg->AddEntry(f_comb, "Combinatorial Bkg", "l");
    leg->Draw("SAME");

    // Box informativa con chi2/ndf e il BR fittato (Stile pulito, senza loghi o bordi)
    TPaveText *paveStats = new TPaveText(0.18, 0.68, 0.46, 0.86, "NDC");
    paveStats->SetBorderSize(0);
    paveStats->SetFillStyle(0);
    paveStats->SetTextFont(43);
    paveStats->SetTextSize(22);
    paveStats->SetTextAlign(12);
    paveStats->AddText(Form("#chi^{2} / ndf = %.1f / %d", chi2, ndf));
    paveStats->AddText("#font[12]{B}(#tau^{+}#rightarrow#phi#mu^{+}) = "
        + FormatPDG(fit_BR * 1e4, fit_BR_err * 1e4, "10^{-4}"));
    paveStats->Draw();

    // Pad Inferiore (Residui / Pulls)
    c_unblind->cd();
    TPad *pad2 = new TPad("pad2", "Pull Pad", 0.0, 0.0, 1.0, splitPoint);
    pad2->SetLeftMargin(0.16);
    pad2->SetRightMargin(0.08);
    pad2->SetTopMargin(0.02);
    pad2->SetBottomMargin(0.38); // Margine per il titolo asse X
    pad2->SetGridy();
    pad2->Draw();
    pad2->cd();

    TH1D *hPull = (TH1D *)h_mass->Clone("hPull_unblinded");
    hPull->Reset();
    hPull->SetStats(0);

    // Configurazione degli assi per il pad dei residui (Default di lbStyle ereditati)
    hPull->GetYaxis()->SetTitle("Pull");
    hPull->GetYaxis()->SetRangeUser(-5.0, 5.0);

    hPull->GetXaxis()->SetTitle("M(D_{s}^{+}) [GeV/#it{c}^{2}]");

    for(int i = 1; i <= h_mass->GetNbinsX(); i++)
    {
        double obs = h_mass->GetBinContent(i);
        double err = h_mass->GetBinError(i);
        if(err > 0)
        {
            double val = f_draw->Eval(h_mass->GetBinCenter(i));
            hPull->SetBinContent(i, (obs - val) / err);
        }
        else
        {
            hPull->SetBinContent(i, -999);
        }
    }
    hPull->SetMarkerStyle(20);
    hPull->SetMarkerSize(0.8);
    hPull->SetLineColor(kBlack);
    hPull->Draw("P");

    TLine *line0 = new TLine(xMin, 0.0, xMax, 0.0);
    line0->SetLineColor(kRed);
    line0->SetLineWidth(2);
    line0->Draw();

    c_unblind->Update();

    // Salvataggio del plot dei residui
    if(savePlots)
    {
        SavePlot(c_unblind, "./_fig/Unblinded_FinalFit_Pulls", false);
        c_unblind->SaveAs("./_root/Unblinded_FinalFit_Pulls.root");
    }

    // =========================================================================
    // STEP 10: DISEGNO DELLA BANDA FELDMAN-COUSINS CONTINUA CON MISURA OVERLAY
    // =========================================================================
    TFile *fFC = TFile::Open("./_root/FC_Belt_Output_BR.root", "READ");
    if(fFC && !fFC->IsZombie())
    {
        // Carichiamo le curve continue anziché quelle discrete
        TGraph *g_belt = (TGraph *)fFC->Get("g_belt_interpolated");
        TGraph *g_up = (TGraph *)fFC->Get("g_up_interpolated");
        TGraph *g_low = (TGraph *)fFC->Get("g_low_interpolated");

        if(g_belt && g_up && g_low)
        {
            // Canvas quadrata (800x800) per preservare l'aspect ratio diagonale (Y = X)
            TCanvas *cFC_meas
                = new TCanvas("cFC_meas", "Feldman-Cousins Unblinded Measurement", 800, 800);
            cFC_meas->cd();
            cFC_meas->SetGrid();

            // Configurazione dei margini per garantire spazio adeguato al testo degli assi
            cFC_meas->SetLeftMargin(0.16);
            cFC_meas->SetRightMargin(0.08);
            cFC_meas->SetTopMargin(0.08);
            cFC_meas->SetBottomMargin(0.14);

            // Clonazione e applicazione del fattore di scala fisico noto (10^-4)
            TGraph *g_belt_scaled = (TGraph *)g_belt->Clone("g_belt_scaled_unblind");
            TGraph *g_up_scaled = (TGraph *)g_up->Clone("g_up_scaled_unblind");
            TGraph *g_low_scaled = (TGraph *)g_low->Clone("g_low_scaled_unblind");

            const double target_scale = 1e4;

            for(int i = 0; i < g_belt_scaled->GetN(); ++i)
            {
                double tx, ty;
                g_belt_scaled->GetPoint(i, tx, ty);
                g_belt_scaled->SetPoint(i, tx * target_scale, ty * target_scale);
            }
            for(int i = 0; i < g_up_scaled->GetN(); ++i)
            {
                double tx, ty;
                g_up_scaled->GetPoint(i, tx, ty);
                g_up_scaled->SetPoint(i, tx * target_scale, ty * target_scale);
            }
            for(int i = 0; i < g_low_scaled->GetN(); ++i)
            {
                double tx, ty;
                g_low_scaled->GetPoint(i, tx, ty);
                g_low_scaled->SetPoint(i, tx * target_scale, ty * target_scale);
            }

            g_belt_scaled->SetFillColorAlpha(kGreen - 9, 0.45);
            g_belt_scaled->SetLineWidth(0);

            g_up_scaled->SetLineColor(kGreen + 2);
            g_up_scaled->SetLineWidth(3);

            g_low_scaled->SetLineColor(kGreen + 2);
            g_low_scaled->SetLineWidth(3);

            // Determinazione deterministica dei limiti basata sul valore massimo (allineato
            // all'originale)
            double max_x = -1e9;
            double max_y = -1e9;
            for(int i = 0; i < g_up_scaled->GetN(); ++i)
            {
                double tx, ty;
                g_up_scaled->GetPoint(i, tx, ty);
                if(tx > max_x)
                    max_x = tx;
                if(ty > max_y)
                    max_y = ty;
            }

            double plot_x_min = -0.05 * max_x;
            double plot_x_max = 1.10 * max_x;
            double plot_y_min = 0.0;
            double plot_y_max = 1.10 * max_y;

            // HARD RESET DI gStyle: Forza l'annullamento di qualsiasi offset estremo o dimensione
            // zero ereditata
            gStyle->SetLabelFont(42, "XYZ");
            gStyle->SetLabelSize(0.04, "XYZ");
            gStyle->SetLabelColor(kBlack, "XYZ");
            gStyle->SetLabelOffset(
                0.007, "XYZ"); // Ripristina l'offset dei numeri a valori visibili standard

            gStyle->SetTitleFont(42, "XYZ");
            gStyle->SetTitleSize(0.045, "XYZ");
            gStyle->SetTitleColor(kBlack, "XYZ");
            gStyle->SetTitleOffset(1.2, "X");
            gStyle->SetTitleOffset(1.5, "Y");

            gStyle->SetAxisColor(kBlack, "XYZ");
            gStyle->SetTickLength(0.03, "XYZ");

            // Creazione cornice TH2F (Semplificato il TLatex usando #it{B} per evitare conflitti di
            // font)
            TH2F *hFrameFC = new TH2F("hFrameFC_unblind",
                ";Measured #hat{#it{B}}(#tau^{+}#rightarrow#phi#mu^{+}) [10^{-4}];True "
                "#it{B}(#tau^{+}#rightarrow#phi#mu^{+}) [10^{-4}]",
                100, plot_x_min, plot_x_max, 100, plot_y_min, plot_y_max);
            hFrameFC->SetStats(0);

            // Sovrascrittura locale di sicurezza su hFrameFC per ignorare lo stile globale
            hFrameFC->GetXaxis()->SetLabelFont(42);
            hFrameFC->GetXaxis()->SetLabelSize(0.04);
            hFrameFC->GetXaxis()->SetLabelColor(kBlack);
            hFrameFC->GetXaxis()->SetLabelOffset(0.007); // Riporta i numeri vicino all'asse
            hFrameFC->GetXaxis()->SetTitleFont(42);
            hFrameFC->GetXaxis()->SetTitleSize(0.045);
            hFrameFC->GetXaxis()->SetTitleColor(kBlack);
            hFrameFC->GetXaxis()->SetTitleOffset(1.2);
            hFrameFC->GetXaxis()->SetAxisColor(kBlack);

            hFrameFC->GetYaxis()->SetLabelFont(42);
            hFrameFC->GetYaxis()->SetLabelSize(0.04);
            hFrameFC->GetYaxis()->SetLabelColor(kBlack);
            hFrameFC->GetYaxis()->SetLabelOffset(0.007);
            hFrameFC->GetYaxis()->SetTitleFont(42);
            hFrameFC->GetYaxis()->SetTitleSize(0.045);
            hFrameFC->GetYaxis()->SetTitleColor(kBlack);
            hFrameFC->GetYaxis()->SetTitleOffset(1.5);
            hFrameFC->GetYaxis()->SetAxisColor(kBlack);

            hFrameFC->Draw();

            // Disegno elementi della banda continua riscalati
            g_belt_scaled->Draw("F SAME");
            g_up_scaled->Draw("L SAME");
            g_low_scaled->Draw("L SAME");

            // Diagonale di misura ideale (Y = X)
            TLine *diag = new TLine(0.0, 0.0, plot_y_max, plot_y_max);
            diag->SetLineStyle(2);
            diag->SetLineColor(kGray + 2);
            diag->SetLineWidth(2);
            diag->Draw("SAME");

            // Coordinate riscalate della misura reale (intersezione con la banda continua)
            double meas_x = fit_BR * target_scale;
            double limit_low = BR_low_c * target_scale;
            double limit_up = BR_up_c * target_scale;

            // 1. Linea verticale tratteggiata per la misura sperimentale X (fino al limite
            // superiore)
            TLine *l_meas_vert = new TLine(meas_x, plot_y_min, meas_x, limit_up);
            l_meas_vert->SetLineStyle(2); // Tratteggio classico
            l_meas_vert->SetLineColor(kBlack);
            l_meas_vert->SetLineWidth(2);
            l_meas_vert->Draw("SAME");

            // 2. Linee orizzontali tratteggiate di proiezione verso l'asse Y
            TLine *l_proj_up = new TLine(plot_x_min, limit_up, meas_x, limit_up);
            l_proj_up->SetLineStyle(3); // Linea punteggiata
            l_proj_up->SetLineColor(kRed);
            l_proj_up->SetLineWidth(2);
            l_proj_up->Draw("SAME");

            TLine *l_proj_low = nullptr;
            if(limit_low > 0.0)
            {
                l_proj_low = new TLine(plot_x_min, limit_low, meas_x, limit_low);
                l_proj_low->SetLineStyle(3);
                l_proj_low->SetLineColor(kRed);
                l_proj_low->SetLineWidth(2);
                l_proj_low->Draw("SAME");
            }

            // 3. Intervallo sblindato al 90% CL finale disegnato direttamente sull'asse Y (True B)
            // Lo scostiamo leggermente a destra dell'asse Y (2% del range X) per non sovrapporlo
            // alla linea dell'asse
            double y_axis_offset = plot_x_min + 0.02 * (plot_x_max - plot_x_min);
            TLine *l_err_y = new TLine(y_axis_offset, limit_low, y_axis_offset, limit_up);
            l_err_y->SetLineColor(kRed);
            l_err_y->SetLineWidth(4);
            l_err_y->Draw("SAME");

            // Trattini orizzontali di chiusura (Caps) sull'asse Y
            double cap_w = 0.015 * (plot_x_max - plot_x_min);
            TLine *l_cap_low
                = new TLine(y_axis_offset - cap_w, limit_low, y_axis_offset + cap_w, limit_low);
            TLine *l_cap_up
                = new TLine(y_axis_offset - cap_w, limit_up, y_axis_offset + cap_w, limit_up);
            l_cap_low->SetLineColor(kRed);
            l_cap_low->SetLineWidth(3);
            l_cap_up->SetLineColor(kRed);
            l_cap_up->SetLineWidth(3);
            l_cap_low->Draw("SAME");
            l_cap_up->Draw("SAME");

            // Legenda del plot di Feldman-Cousins (Font 43 a 24 pixel ereditato da stile)
            TLegend *legFC = new TLegend(0.18, 0.68, 0.58, 0.86);
            legFC->SetBorderSize(0);
            legFC->SetFillStyle(0);
            legFC->SetTextFont(43);
            legFC->SetTextSize(24);
            legFC->AddEntry(g_belt_scaled, "Feldman-Cousins Belt (90% CL)", "f");
            legFC->AddEntry(l_err_y, "Unblinded Interval (90% CL)", "l");
            legFC->AddEntry(l_meas_vert, "Point Estimate (#hat{B})", "l");
            legFC->Draw("SAME");

            gPad->RedrawAxis();
            cFC_meas->Update();

            if(savePlots)
            {
                SavePlot(cFC_meas, "./_fig/FeldmanCousins_FinalMeasurement", false);
                cFC_meas->SaveAs("./_root/FeldmanCousins_FinalMeasurement.root");
            }
        }
        fFC->Close();
    }

    // =========================================================================
    // STEP 11: PULIZIA DELLA MEMORIA
    // =========================================================================
    // delete minH1;
    // delete minH0;
    // delete g_pdf_unbinned;

    //    //g_pdf_unbinned = nullptr;
    // g_data_hist = nullptr;
}
