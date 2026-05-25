#include <RtypesCore.h>

#include "GuiTypes.h"
#define lifetimeANA_cxx
#include <iostream>

#include <TCanvas.h>
#include <TColor.h>
#include <TF1.h>
#include <TH1.h>
#include <TH2D.h>
#include <TLegend.h>
#include <TLine.h>
#include <TMath.h>
#include <TPaveStats.h>
#include <TProfile.h>
#include <TRandom3.h>
#include <TStyle.h>

#include "lbrootstyle.hh"
#include "lifetimeANA.h"

using namespace std;
using namespace TMath;

constexpr Double_t MC_LIFE = 410.3e-15; // s
constexpr Bool_t savePlots = false;
constexpr Double_t TIME_CUT = 1; // Taglio minimo su t/tau_MC

// --- FUNZIONI DI FIT ---
Double_t f_2G_Shared(Double_t *x, Double_t *par)
{
    return par[1] * exp(-0.5 * pow((x[0] - par[0]) / par[2], 2))
        + par[3] * exp(-0.5 * pow((x[0] - par[0]) / par[4], 2));
}
Double_t f_2G_Indep(Double_t *x, Double_t *par)
{
    return par[1] * exp(-0.5 * pow((x[0] - par[0]) / par[2], 2))
        + par[4] * exp(-0.5 * pow((x[0] - par[3]) / par[5], 2));
}
Double_t f_2G_Frac(Double_t *x, Double_t *par)
{
    // par[0] = Yield (Eventi * bin width)
    // par[1] = Mu_1
    // par[2] = Sigma_1
    // par[3] = Mu_2
    // par[4] = Sigma_2
    // par[5] = Frac_1 (Frazione della prima gaussiana)

    double dx1 = (x[0] - par[1]) / par[2];
    double dx2 = (x[0] - par[3]) / par[4];

    double g1 = TMath::Exp(-0.5 * dx1 * dx1) / (par[2] * TMath::Sqrt(TMath::TwoPi()));
    double g2 = TMath::Exp(-0.5 * dx2 * dx2) / (par[4] * TMath::Sqrt(TMath::TwoPi()));

    return par[0] * (par[5] * g1 + (1.0 - par[5]) * g2);
}
Double_t f_3G_Shared(Double_t *x, Double_t *par)
{
    return par[1] * exp(-0.5 * pow((x[0] - par[0]) / par[2], 2))
        + par[3] * exp(-0.5 * pow((x[0] - par[0]) / par[4], 2))
        + par[5] * exp(-0.5 * pow((x[0] - par[0]) / par[6], 2));
}
Double_t f_3G_Indep(Double_t *x, Double_t *par)
{
    return par[1] * exp(-0.5 * pow((x[0] - par[0]) / par[2], 2))
        + par[4] * exp(-0.5 * pow((x[0] - par[3]) / par[5], 2))
        + par[7] * exp(-0.5 * pow((x[0] - par[6]) / par[8], 2));
}

class ExpoMultiGaussConv
{
  private:
    int fNGauss;

    Double_t SingleConv(Double_t t, Double_t tau, Double_t mu, Double_t sigma) const
    {
        if(tau <= 0 || sigma <= 0)
            return 0;
        Double_t invTau = 1.0 / tau;
        Double_t argExp = (sigma * sigma) / (2.0 * tau * tau) - (t - mu) / tau;
        Double_t argErfc = sigma / (sqrt(2.0) * tau) - (t - mu) / (sqrt(2.0) * sigma);
        if(argExp > 700)
            return 0;
        return (invTau / 2.0) * TMath::Exp(argExp) * TMath::Erfc(argErfc);
    }

  public:
    ExpoMultiGaussConv(int n)
        : fNGauss(n)
    {
    }

    // Mapping dei parametri par[]:
    // [0] : Tau
    // [1] : Mu (punto zero)
    // [2] : Yield (Area totale segnale)
    // [3] a [3+N-1] : Le Sigma (sigma_1, sigma_2, ..., sigma_N)
    // [3+N] a [3+N+(N-2)] : Le Frazioni (f_1, f_2, ..., f_{N-1})
    // [Ultimo] : Background

    Double_t operator()(Double_t *x, Double_t *par)
    {
        Double_t t = x[0];
        Double_t tau = par[0];
        Double_t mu = par[1];
        Double_t yield = par[2];

        Double_t totalPDF = 0;
        Double_t sumFrac = 0;

        for(int i = 0; i < fNGauss; i++)
        {
            Double_t sigma = par[3 + i]; // Prende le sigma in sequenza
            Double_t frac = 0;

            if(i < fNGauss - 1)
            {
                // Frazioni libere
                frac = par[3 + fNGauss + i];
                sumFrac += frac;
            }
            else
            {
                // L'ultima frazione è vincolata per avere somma = 1
                frac = 1.0 - sumFrac;
            }

            totalPDF += frac * SingleConv(t, tau, mu, sigma);
        }

        // Il background è l'ultimo parametro del vettore
        Double_t bkg = par[GetNPar() - 1];

        return yield * totalPDF + bkg;
    }

    int GetNPar() const
    {
        // 3 (Tau, Mu, Yield) + N (Sigmas) + (N-1) (Fractions) + 1 (Bkg)
        return 2 * fNGauss + 3;
    }
};

class FullPDF_ConvAcc
{
  private:
    int fNGauss;
    bool fUseAcc;

    // Convoluzione Esponenziale x Gaussiana numericamente stabile
    Double_t SingleConv(Double_t t, Double_t tau, Double_t mu, Double_t sigma) const
    {
        if(tau <= 0 || sigma <= 0)
            return 0;

        Double_t argExp = (sigma * sigma) / (2.0 * tau * tau) - (t - mu) / tau;
        Double_t argErfc = sigma / (sqrt(2.0) * tau) - (t - mu) / (sqrt(2.0) * sigma);

        // RIPRISTINATO: Protezione contro NaN. Per argExp molto grandi,
        // l'espansione asintotica di Erfc annulla l'esponenziale portando il limite a 0.
        if(argExp > 700.0)
            return 0.0;

        return (1.0 / (2.0 * tau)) * TMath::Exp(argExp) * TMath::Erfc(argErfc);
    }

    // Accettanza normalizzata al plateau = 1 (Sum of two Erfs)
    Double_t Acceptance(Double_t t, Double_t *par, int offset) const
    {
        if(!fUseAcc)
            return 1.0;

        Double_t frac = par[offset];
        Double_t mu1 = par[offset + 1];
        Double_t sig1 = par[offset + 2];
        Double_t mu2 = par[offset + 3];
        Double_t sig2 = par[offset + 4];

        Double_t erf1 = 0.5 * (1.0 + TMath::Erf((t - mu1) / sig1));
        Double_t erf2 = 0.5 * (1.0 + TMath::Erf((t - mu2) / sig2));

        return frac * erf1 + (1.0 - frac) * erf2;
    }

  public:
    FullPDF_ConvAcc(int n, bool useAcc = true)
        : fNGauss(n)
        , fUseAcc(useAcc)
    {
    }

    // Nuovo Mapping dei parametri par[] per medie INDIPENDENTI:
    // [0] : Tau
    // [1] a [N] : Mu_1, ..., Mu_N (medie indipendenti)
    // [N+1] : Yield
    // [N+2] a [2N+1] : Sigma_1, ..., Sigma_N
    // [2N+2] a [3N] : Frazioni (f_1, ..., f_{N-1})
    // --- SE fUseAcc == true (Offset = 3*N + 1) ---
    // [3N+1] a [3N+5] : Parametri Accettanza (frac, mu1, sig1, mu2, sig2)
    // [Ultimo] : Background
    Double_t operator()(Double_t *x, Double_t *par)
    {
        Double_t t = x[0];
        Double_t tau = par[0];
        Double_t yield = par[fNGauss + 1];

        Double_t totalConvPDF = 0;
        Double_t sumFrac = 0;

        for(int i = 0; i < fNGauss; i++)
        {
            Double_t mu = par[1 + i];
            Double_t sigma = par[fNGauss + 2 + i];
            Double_t frac = 0;

            if(i < fNGauss - 1)
            {
                frac = par[2 * fNGauss + 2 + i];
                sumFrac += frac;
            }
            else
            {
                frac = 1.0 - sumFrac;
            }
            totalConvPDF += frac * SingleConv(t, tau, mu, sigma);
        }

        int accOffset = 3 * fNGauss + 1;
        Double_t acc = Acceptance(t, par, accOffset);

        return yield * (totalConvPDF * acc);
    }

    int GetNPar() const
    {
        // 1 (Tau) + N (Mus) + 1 (Yield) + N (Sigmas) + (N-1) (Fractions) + [5 if Acc]
        return 3 * fNGauss + 1 + (fUseAcc ? 5 : 0);
    }
};

// --- Methods ---

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
    TH1D *histo_data_MKpi = new TH1D("histo_data_MKpi", "", 100, 1.8, 1.95);
    TH1D *histo_mc_MKpi = new TH1D("histo_mc_MKpi", "", 100, 1.8, 1.95);
    TH1D *histo_mc_time = new TH1D("histo_mc_time", "", 100, 0, 10);

    histo_data_MKpi->Sumw2();
    histo_mc_MKpi->Sumw2();
    histo_mc_time->Sumw2();

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

        // #### My code #########
        if(id == 1)
        { // DATA

            histo_data_MKpi->Fill(M0_MKpi);
            // Print out the Mpipi values
            // std::cout <<M0_MKpi << std::endl;

        }; // end DATA

        if(id == 13)
        { // MC
            histo_mc_MKpi->Fill(M0_MKpi);
            histo_mc_time->Fill(M0_time / (410.3e-15));
            // std::cout <<M0_MKpi << std::endl;

        }; // end MC

    }; // #### end loop over jentry

    // Create Canvas
    TCanvas *c_histo_data_MKpi = new TCanvas("c_histo_data_MKpi", "canvas histo", 500, 500);
    c_histo_data_MKpi->cd();
    histo_data_MKpi->Draw();
    if(savePlots)
        c_histo_data_MKpi->Print("./_fig/c_histo_data_MKpi.pdf");
    if(savePlots)
        c_histo_data_MKpi->Print("./_fig/c_histo_data_MKpi.eps");

    // Create Canvas
    TCanvas *c_histo_mc_MKpi = new TCanvas("c_histo_mc_MKpi", "canvas histo", 500, 500);
    c_histo_mc_MKpi->cd();
    histo_mc_MKpi->Draw();
    if(savePlots)
        c_histo_mc_MKpi->Print("./_fig/c_histo_mc_MKpi.pdf");
    if(savePlots)
        c_histo_mc_MKpi->Print("./_fig/c_histo_mc_MKpi.eps");

    // Create Canvas
    TCanvas *c_histo_mc_time = new TCanvas("c_histo_mc_time", "canvas histo", 500, 500);
    c_histo_mc_time->cd();
    histo_mc_time->Draw();
    if(savePlots)
        c_histo_mc_time->Print("./_fig/c_histo_mc_time.pdf");
    if(savePlots)
        c_histo_mc_time->Print("./_fig/c_histo_mc_time.eps");

    // Create a new file to store histograms
    TFile *histo_file = new TFile("./_root/histo_file.root", "RECREATE", "put a title");
    histo_file->cd();
    histo_data_MKpi->Write();
    histo_mc_MKpi->Write();
    histo_mc_time->Write();
    histo_file->Close();

    return;
}

void lifetimeANA::Resolution(Int_t nGauss, Bool_t useSharedMean)
{
    // Stile globale
    gStyle->SetOptStat(1110);
    gStyle->SetOptFit(1111);
    gStyle->SetStatX(0.99);
    gStyle->SetStatY(0.94);

    auto hRes = new TH1D("hRes", "", 60, -0.2, 0.2);

    if(fChain == 0)
        return;
    Long64_t nentries = fChain->GetEntriesFast();
    for(Long64_t jentry = 0; jentry < nentries; jentry++)
    {
        if(LoadTree(jentry) < 0)
            break;
        fChain->GetEntry(jentry);
        if(id == 13)
        {
            Double_t t_reco = M0_time / MC_LIFE;

            // --- 1. APPLICAZIONE DEL TAGLIO (LIFETIME CUT) ---
            // Riempiamo l'istogramma solo se siamo fuori dalla zona di forte bias (< 0.15)
            if(t_reco > 0.15)
            {
                hRes->Fill((M0_time - M0_time_true) * 1e12);
            }
        }
    }

    // --- LOGICA FIT ---
    Double_t m = hRes->GetMean(), s = hRes->GetRMS(), a = hRes->GetMaximum();
    TF1 *fRes = nullptr;
    if(nGauss == 2)
    {
        if(useSharedMean)
        {
            fRes = new TF1("fRes", f_2G_Shared, -0.2, 0.2, 5);
            fRes->SetParameters(m, a, s * 0.6, a * 0.2, s * 1.5);
            fRes->SetParNames("Mean", "N1", "S1", "N2", "S2");
        }
        else
        {
            fRes = new TF1("fRes", f_2G_Indep, -0.2, 0.2, 6);
            fRes->SetParameters(m, a, s * 0.6, m, a * 0.2, s * 1.5);
            fRes->SetParNames("Mean1", "N1", "S1", "Mean2", "N2", "S2");
        }
    }
    else
    {
        if(useSharedMean)
        {
            fRes = new TF1("fRes", f_3G_Shared, -0.2, 0.2, 7);
            fRes->SetParameters(m, a, s * 0.5, a * 0.3, s * 1.2, a * 0.1, s * 2.5);
            fRes->SetParNames("Mean", "N1", "S1", "N2", "S2", "N3", "S3");
        }
        else
        {
            fRes = new TF1("fRes", f_3G_Indep, -0.2, 0.2, 9);
            fRes->SetParameters(m, a, s * 0.5, m, a * 0.3, s * 1.2, m, a * 0.1, s * 2.5);
            fRes->SetParNames("Mean1", "N1", "S1", "Mean2", "N2", "S2", "Mean3", "N3", "S3");
        }
    }
    hRes->Fit(fRes, "L I R");

    // --- DISEGNO ---
    TCanvas *c1 = new TCanvas("c1", "Fit and Residuals", 900, 900);

    // Proporzioni Pad
    double splitPoint = 0.33;
    double margin = 0.15; // Margine sinistro generoso per i titoli Y

    // Pad 1: Plot principale
    TPad *pad1 = new TPad("pad1", "pad1", 0, splitPoint, 1, 1.0);
    pad1->SetBottomMargin(0.015);
    pad1->SetLeftMargin(margin);
    pad1->SetRightMargin(0.08);
    pad1->Draw();
    pad1->cd();

    hRes->GetYaxis()->SetTitle(Form("Entries / (%.4f ps)", hRes->GetBinWidth(1)));
    hRes->GetYaxis()->SetTitleSize(0.045);
    hRes->GetYaxis()->SetTitleOffset(1.5);
    hRes->GetYaxis()->SetLabelSize(0.04);
    hRes->GetXaxis()->SetLabelSize(0); // Nascondi X sopra
    hRes->Draw("HIST E");
    fRes->Draw("same");

    // Componenti
    for(int i = 0; i < nGauss; i++)
    {
        TF1 *gi = new TF1(Form("g%d", i), "gaus", -0.2, 0.2);
        if(useSharedMean)
            gi->SetParameters(fRes->GetParameter(1 + i * 2), fRes->GetParameter(0),
                fRes->GetParameter(2 + i * 2));
        else
            gi->SetParameters(fRes->GetParameter(1 + i * 3), fRes->GetParameter(0 + i * 3),
                fRes->GetParameter(2 + i * 3));
        gi->SetLineStyle(2);
        gi->SetLineColor(i == 0 ? kGreen + 2 : (i == 1 ? kCyan + 1 : kMagenta));
        gi->Draw("same");
    }

    // Pad 2: Residui
    c1->cd();
    TPad *pad2 = new TPad("pad2", "pad2", 0, 0, 1, splitPoint);
    pad2->SetTopMargin(0.015);
    pad2->SetBottomMargin(0.35);
    pad2->SetLeftMargin(margin);
    pad2->SetRightMargin(0.08);
    pad2->SetGridy();
    pad2->Draw();
    pad2->cd();

    TH1D *hDelta = (TH1D *)hRes->Clone("hDelta");
    hDelta->Reset();
    hDelta->SetStats(0);
    for(int i = 1; i <= hRes->GetNbinsX(); i++)
    {
        hDelta->SetBinContent(i, hRes->GetBinContent(i) - fRes->Eval(hRes->GetBinCenter(i)));
        hDelta->SetBinError(i, hRes->GetBinError(i));
    }

    // --- PAREGGIO GRANDEZZA TESTI ---
    // Rapporto altezze pad1/pad2 = (1-0.33) / (0.33-0) = 0.67 / 0.33 = ~2.03
    double ratio = (1.0 - splitPoint) / splitPoint;

    hDelta->GetYaxis()->SetTitle("Residuals");
    hDelta->GetYaxis()->SetTitleSize(0.045 * ratio); // Pareggiato a Entries
    hDelta->GetYaxis()->SetTitleOffset(1.5 / ratio);
    hDelta->GetYaxis()->SetLabelSize(0.04 * ratio);
    hDelta->GetYaxis()->SetNdivisions(505);

    hDelta->GetXaxis()->SetTitle("(M0_time - M0_time_true) [ps]");
    hDelta->GetXaxis()->SetTitleSize(0.045 * ratio);
    hDelta->GetXaxis()->SetTitleOffset(1.1);
    hDelta->GetXaxis()->SetLabelSize(0.04 * ratio);

    hDelta->SetMarkerStyle(20);
    hDelta->SetMarkerSize(0.8);
    hDelta->Draw("E P");

    TLine *line0 = new TLine(-0.2, 0, 0.2, 0);
    line0->SetLineColor(kRed);
    line0->SetLineWidth(2);
    line0->Draw();

    // Finally
    if(savePlots)
        c1->SaveAs(Form(
            "fit_resolution_%dGauss_%sMean.png", nGauss, useSharedMean ? "Shared" : "Independent"));
}

void lifetimeANA::Resolution2D()
{
    SetLBStyle();
    // Stile globale per i plot
    gStyle->SetOptStat(0);
    // gStyle->SetPalette(kBird);
    //  gStyle->SetTitleSize(0.04, "XYZ");
    //  gStyle->SetLabelSize(0.035, "XYZ");

    // 1. Definizione Istogrammi 2D e Supporto
    auto hRes2D
        = new TH2D("hRes2D", "Residuals vs True Time;M0_time_true [ps];M0_time - M0_time_true [ps]",
            60, 0, 2, 80, -0.2, 0.2);

    if(fChain == 0)
        return;
    Long64_t nentries = fChain->GetEntriesFast();

    for(Long64_t jentry = 0; jentry < nentries; jentry++)
    {
        if(LoadTree(jentry) < 0)
            break;
        fChain->GetEntry(jentry);
        if(id == 13)
        {
            Double_t true_ps = M0_time_true * 1e12;
            Double_t res_ps = (M0_time - M0_time_true) * 1e12;
            hRes2D->Fill(true_ps, res_ps);
        }
    }

    // --- CALCOLO BIAS (Profile) ---
    TProfile *pBias = hRes2D->ProfileX("pBias");
    pBias->SetTitle("2D Residuals vs True Time;M0_time_true [ps];Residuals [ps]");
    // pBias->SetMarkerStyle(20);
    // pBias->SetMarkerSize(1.0);
    pBias->SetMarkerColor(kBlue + 1);
    pBias->SetLineColor(kBlue + 1);

    // --- CALCOLO RESOLUTION (Std Dev) ---
    TH1D *hStdDev = new TH1D("hStdDev",
        "Resolution (RMS of Residuals);M0_time_true [ps];Standard Deviation [ps]", 60, 0, 2);

    for(int i = 1; i <= hRes2D->GetNbinsX(); i++)
    {
        TH1D *hTmp = hRes2D->ProjectionY("_tmp", i, i);
        if(hTmp->GetEntries() > 5)
        {
            hStdDev->SetBinContent(i, hTmp->GetRMS());
            hStdDev->SetBinError(i, hTmp->GetRMSError());
        }
        delete hTmp;
    }
    // hStdDev->SetMarkerStyle(20);
    // hStdDev->SetMarkerSize(1.0);
    hStdDev->SetMarkerColor(kRed + 1);
    hStdDev->SetLineColor(kRed + 1);
    hStdDev->SetMinimum(0);

    // --- CANVAS 1: TH2D COLZ ---
    TCanvas *c_2D = new TCanvas("c_2D", "Residuals 2D");
    Fix2DMargins(c_2D);
    hRes2D->Draw("COLZ");
    if(savePlots)
        c_2D->SaveAs("plot_2D_residuals.png");

    // --- CANVAS 2: BIAS (Mean) ---
    TCanvas *c_Bias = new TCanvas("c_Bias", "Bias Plot");
    // c_Bias->SetLeftMargin(0.13);
    // pBias->GetYaxis()->SetTitleOffset(1.3);
    pBias->Draw("E1");

    TLine *l0 = new TLine(0, 0, 2, 0);
    l0->SetLineColor(kBlack);
    l0->SetLineStyle(7); // Tratteggiato
    l0->Draw();
    if(savePlots)
        c_Bias->SaveAs("plot_2D_bias.pdf");

    // --- CANVAS 3: RESOLUTION (Std Dev) ---
    TCanvas *c_StdDev = new TCanvas("c_StdDev", "Resolution Plot");
    // c_StdDev->SetLeftMargin(0.13);
    // hStdDev->GetYaxis()->SetTitleOffset(1.3);
    hStdDev->Draw("E1");
    if(savePlots)
        c_StdDev->SaveAs("plot_2D_resolution.pdf");
}

void lifetimeANA::ConvolvedFit(Bool_t useAcceptance)
{
    SetLBStyle();

    // Stile globale
    gStyle->SetOptStat(1110);
    gStyle->SetOptFit(1111);
    gStyle->SetStatX(0.99);
    gStyle->SetStatY(0.94);

    auto hTime_reco = new TH1D("hTime_reco", "Fit Temporale;t /#tau_{MC};Events", 100, 0., 10.);
    hTime_reco->Sumw2();

    if(fChain == 0)
        return;
    Long64_t nentries = fChain->GetEntriesFast();
    for(Long64_t jentry = 0; jentry < nentries; jentry++)
    {
        if(LoadTree(jentry) < 0)
            break;
        fChain->GetEntry(jentry);
        if(id == 13)
        {
            Double_t t_reco = M0_time / MC_LIFE;
            hTime_reco->Fill(t_reco);
        }
    }

    // Configurazione modello a 2 Gaussiane con medie indipendenti
    Int_t nGaus = 2;
    FullPDF_ConvAcc fModel(nGaus, useAcceptance);
    Int_t nPar = fModel.GetNPar();

    Double_t minT = 1.;
    Double_t maxT = 10.;
    auto fFit = new TF1("fFit", fModel, minT, maxT, nPar);
    fFit->SetNpx(1000);

    // --- INIZIALIZZAZIONE PARAMETRI ---
    fFit->SetParameter(0, 1.0);
    fFit->SetParName(0, "Tau");
    fFit->FixParameter(1, 0.0029);
    fFit->SetParName(1, "Mu_1");
    fFit->FixParameter(2, 0.0063);
    fFit->SetParName(2, "Mu_2"); // Media della seconda Gaussiana
    fFit->SetParameter(3, hTime_reco->Integral());
    fFit->SetParName(3, "Yield");

    // Risoluzione (Sigma e Frazioni)
    fFit->FixParameter(4, 0.0337);
    fFit->SetParName(4, "Sigma_1");
    // fFit->SetParLimits(4, 0.001, 0.5);
    fFit->FixParameter(5, 0.0517);
    fFit->SetParName(5, "Sigma_2");
    // fFit->SetParLimits(5, 0.001, 0.5);

    fFit->FixParameter(6, 0.5327);
    fFit->SetParName(6, "G_Frac_1");
    // fFit->SetParLimits(6, 0.0, 1.0);

    // Accettanza (Offset di partenza = 3*nGaus + 1 = 7)
    if(useAcceptance)
    {
        int aIdx = 7;
        fFit->FixParameter(aIdx + 0, 0.640);
        fFit->SetParName(aIdx + 0, "Acc_Frac");
        fFit->FixParameter(aIdx + 1, 0.657);
        fFit->SetParName(aIdx + 1, "Acc_Mu1");
        fFit->FixParameter(aIdx + 2, 0.223);
        fFit->SetParName(aIdx + 2, "Acc_Sig1");
        fFit->FixParameter(aIdx + 3, 1.085);
        fFit->SetParName(aIdx + 3, "Acc_Mu2");
        fFit->FixParameter(aIdx + 4, 0.454);
        fFit->SetParName(aIdx + 4, "Acc_Sig2");
    }

    // Esegui il FIT
    hTime_reco->Fit(fFit, "L I 0", "", minT, maxT);

    // --- DISEGNO CON PAD SDOPPIATI ---
    TCanvas *c_fit = new TCanvas("c_fit", "Convolved & Acceptance Fit with Pulls", 1200, 900);

    double splitPoint = 0.30;
    double leftMargin = 0.15;
    double rightMargin = 0.05;

    // Pad 1: Plot principale (Fit)
    TPad *pad1 = new TPad("pad1", "Main Fit Pad", 0.0, splitPoint, 1.0, 1.0);
    pad1->SetBottomMargin(0.02); // Margine piccolo per toccare il pad inferiore
    pad1->SetLeftMargin(leftMargin);
    pad1->SetRightMargin(rightMargin);
    pad1->SetLogy();
    pad1->Draw();
    pad1->cd();

    hTime_reco->GetYaxis()->SetTitle("Events");
    hTime_reco->GetYaxis()->SetTitleSize(0.05);
    hTime_reco->GetYaxis()->SetTitleOffset(1.35);
    hTime_reco->GetYaxis()->SetLabelSize(0.04);
    hTime_reco->GetXaxis()->SetLabelSize(0); // Nasconde i numeri sull'asse X superiore
    hTime_reco->GetXaxis()->SetTitleSize(0);

    hTime_reco->SetMinimum(0.5); // Per la scala logaritmica
    hTime_reco->Draw("E");
    fFit->SetLineColor(kRed);
    fFit->Draw("same");

    // Pad 2: Residui Normalizzati (Pulls)
    c_fit->cd();
    TPad *pad2 = new TPad("pad2", "Pull Pad", 0.0, 0.0, 1.0, splitPoint);
    pad2->SetTopMargin(0.02);
    pad2->SetBottomMargin(0.35);
    pad2->SetLeftMargin(leftMargin);
    pad2->SetRightMargin(rightMargin);
    pad2->SetGridy();
    pad2->Draw();
    pad2->cd();

    // Calcolo dei residui normalizzati (Pull)
    TH1D *hPull = (TH1D *)hTime_reco->Clone("hPull");
    hPull->Reset();
    hPull->SetStats(0);

    for(int i = 1; i <= hTime_reco->GetNbinsX(); i++)
    {
        double x = hTime_reco->GetBinCenter(i);
        double obs = hTime_reco->GetBinContent(i);
        double err = hTime_reco->GetBinError(i);

        // Calcola il pull solo all'interno del range di fit
        if(x >= 0. && x <= 8.0)
        {
            if(err > 0)
            {
                double fitVal = fFit->Eval(x);
                double pull = (obs - fitVal) / err;
                hPull->SetBinContent(i, pull);
                // hPull->SetBinError(i, 1.0); // Errore unitario per la barra del pull
            }
        }
    }

    // Adattamento proporzioni testi per il pad inferiore
    double ratio = (1.0 - splitPoint) / splitPoint; // ~2.33

    hPull->GetYaxis()->SetTitle("Pull");
    hPull->GetYaxis()->SetTitleSize(0.05 * ratio);
    hPull->GetYaxis()->SetTitleOffset(1.35 / ratio);
    hPull->GetYaxis()->SetLabelSize(0.04 * ratio);
    hPull->GetYaxis()->SetNdivisions(505);
    hPull->GetYaxis()->SetRangeUser(-5.0, 5.0); // Range standard per i pull

    hPull->GetXaxis()->SetTitle("t / #tau_{MC}");
    hPull->GetXaxis()->SetTitleSize(0.05 * ratio);
    hPull->GetXaxis()->SetTitleOffset(1.1);
    hPull->GetXaxis()->SetLabelSize(0.04 * ratio);

    hPull->SetMarkerStyle(20);
    hPull->SetMarkerSize(0.8);
    hPull->SetMarkerColor(kBlack);
    hPull->SetLineColor(kBlack);
    hPull->Draw("P");

    // Linea orizzontale di riferimento a Pull = 0
    TLine *line0 = new TLine(0.0, 0.0, 8.0, 0.0);
    line0->SetLineColor(kRed);
    line0->SetLineWidth(2);
    line0->Draw();

    c_fit->Update();

    if(savePlots)
        c_fit->SaveAs("plot_full_fit_pulls.pdf");
}

/*
void lifetimeANA::ConvolvedFit()
{
    // Stile globale
    gStyle->SetOptStat(1110);
    gStyle->SetOptFit(1111);
    gStyle->SetStatX(0.99);
    gStyle->SetStatY(0.94);

    auto hTime_true = new TH1D("hTime_true", "", 100, 0., 0.);

    if(fChain == 0)
        return;
    Long64_t nentries = fChain->GetEntriesFast();
    for(Long64_t jentry = 0; jentry < nentries; jentry++)
    {
        if(LoadTree(jentry) < 0)
            break;
        fChain->GetEntry(jentry);
        if(id == 13)
        {
            Bool_t cut = true;
            if(cut)
                hTime_true->Fill(M0_time / MC_LIFE);
        }
    }

    // Functor instance
    Int_t nGaus = 2;
    ExpoMultiGaussConv FConv(nGaus);
    Int_t nPar = FConv.GetNPar();

    // Fit function
    auto fConv = new TF1("fConv", FConv, 2., 8., nPar);

    // Parametri fissi
    fConv->SetParameter(0, 1.0); // Tau
    fConv->SetParameter(1, 1.0); // Mu
    fConv->SetParameter(2, hTime_true->Integral()); // Yield

    for(int i = 0; i < nGaus; i++)
    {
        fConv->SetParameter(3 + i, 0.05 * (i + 1)); // Sigma_i
        fConv->SetParLimits(3 + i, 0.001, 0.5);

        if(i < nGaus - 1)
        {
            fConv->SetParameter(3 + nGaus + i, 1.0 / nGaus); // Frazione f_i
            fConv->SetParLimits(3 + nGaus + i, 0.0, 1.0);
        }
    }

    // Background (L'ultimo)
    fConv->SetParameter(nPar - 1, 1.0);

    hTime_true->Fit(fConv, "L I R 0");

    TCanvas *c_fit = new TCanvas("c_fit", "Convolved Fit", 800, 700);
    hTime_true->Draw("E HISTO");
    fConv->Draw("same");

    if(savePlots)
        c_fit->SaveAs("plot_conv_fit.png");
}
*/

void lifetimeANA::Acceptance(UInt_t seed)
{
    TRandom3 rnd = TRandom3(seed);

    TH1D *hTime_gen = new TH1D("hTime_Gen", "", 200, 0., 10.);
    TH1D *hTime_rec = new TH1D("hTime_Rec", "", 200, 0., 10.);
    hTime_gen->Sumw2();
    hTime_rec->Sumw2();

    for(int i = 0; i < fChain->GetEntriesFast(); i++)
    {
        if(LoadTree(i) < 0)
            break;
        fChain->GetEntry(i);
        if(id == 13) // MC
        {
            Double_t t_reco = M0_time / MC_LIFE;

            // --- 1. APPLICAZIONE DEL TAGLIO (LIFETIME CUT) ---
            // Riempiamo l'istogramma solo se siamo fuori dalla zona di forte bias (< 0.15)
            if(t_reco > 0.15)
            {
                hTime_rec->Fill(M0_time_true / MC_LIFE);
            }
        }
    }

    for(int i = 0; i < 1e9; i++)
    {
        double x = rnd.Exp(MC_LIFE);
        hTime_gen->Fill(x / MC_LIFE);

        if(i % int(1e6) == 0)
            std::cout << "\rGenerated " << i << " events" << std::flush;
    }
    cout << endl;

    // 1. Esegui la divisione diretta
    hTime_rec->Divide(hTime_rec, hTime_gen, 1.0, 1.0, "B");

    // 2. Definisci l'intervallo del plateau (da 2.5 fino al limite superiore dell'istogramma)
    double xMinPlateau = 2.5;
    double xMaxPlateau
        = hTime_rec->GetXaxis()->GetXmax(); // Prende il limite destro dell'istogramma (es. 10.0)

    // 3. Crea la funzione costante e fitta l'istogramma nel range definito
    TF1 *fPlateau = new TF1("fPlateau", "pol0", xMinPlateau, xMaxPlateau);

    // Fit plateau
    hTime_rec->Fit(fPlateau, "RQ0");
    double plateauValue = fPlateau->GetParameter(0);

    // 4. Riscala l'istogramma per questo valore
    if(plateauValue > 0)
    {
        hTime_rec->Scale(1.0 / plateauValue);
    }
    else
    {
        std::cout << "Attenzione: valore del plateau non valido (<= 0)." << std::endl;
    }

    // Simple erf
    TF1 *fErf = new TF1("fErf", "[0]/2.0 * (1.0 + TMath::Erf((x - [1])/[2]))", 0., 8.);
    // Impostazione dei parametri iniziali (fondamentale per la convergenza)
    fErf->SetParameters(1.0, 1.0, 0.5);

    // Acceptance function
    // [0] = Altezza plateau
    // [1] = Punto di transizione (flesso)
    // [2] = Larghezza della coda bassa (salita a sinistra)
    // [3] = Larghezza della coda alta (avvicinamento al plateau a destra)
    TF1 *fAsymErf = new TF1("fAsymErf",
        "[0]/2.0 * (1.0 + (x < [1] ? TMath::Erf((x-[1])/[2]) : TMath::Erf((x-[1])/[3])))", 0., 8.);
    // Impostazione dei parametri iniziali
    fAsymErf->SetParameters(1.0, 0.8, 0.3, 1.2);

    // [0] = Plateau
    // [1] = Parametro di posizione
    // [2] = Parametro di larghezza
    // [3] = Parametro di asimmetria (esponente)
    TF1 *fRichards
        = new TF1("fRichards", "[0] / TMath::Power(1.0 + TMath::Exp(-(x-[1])/[2]), [3])", 0., 10.);
    // Impostazione dei parametri iniziali
    fRichards->SetParameters(1.0, 0.5, 0.3, 0.5);

    // [0] = Altezza Plateau
    // [1] = Frazione della prima Erf (tra 0 e 1)
    // [2] = Punto di flesso 1
    // [3] = Larghezza salita 1
    // [4] = Punto di flesso 2
    // [5] = Larghezza salita 2
    TF1 *fSumErf = new TF1("fSumErf",
        "[0]/2.0 * ( [1]*(1.0 + TMath::Erf((x-[2])/[3])) + (1.0-[1])*(1.0 + "
        "TMath::Erf((x-[4])/[5])) )",
        0., 8);
    fSumErf->SetParameters(1.0, 0.5, 0.5, 0.2, 1.5, 1.0);

    // [0] = Plateau
    // [1] = mu (parametro di posizione del logaritmo)
    // [2] = sigma (parametro di larghezza/asimmetria)
    TF1 *fLogNormCDF = new TF1("fLogNormCDF",
        "x > 0 ? [0]/2.0 * (1.0 + TMath::Erf((TMath::Log(x) - [1]) / ([2] * TMath::Sqrt2()))) : 0",
        0., 4.);
    fLogNormCDF->SetParameters(1.0, 0.0, 0.5);

    TF1 *fitFunc = fSumErf;
    fitFunc->SetNpx(1e4);
    hTime_rec->Fit(fitFunc, "L I R 0");

    auto cDiv = new TCanvas("cDiv", "", 800, 700);
    cDiv->cd();
    hTime_rec->Draw("E HISTO");
    fitFunc->Draw("SAME");
}

void lifetimeANA::AcceptanceWeight(Float_t life, UInt_t seed)
{
    TRandom3 rnd = TRandom3(seed);

    TH1D *hTime_gen = new TH1D("hTime_Gen", "", 200, 0., 10.);
    TH1D *hTime_true = new TH1D("hTime_true", "", 200, 0., 10.);
    hTime_gen->Sumw2();
    hTime_true->Sumw2();

    Double_t simTau = life * MC_LIFE;

    for(int i = 0; i < fChain->GetEntriesFast(); i++)
    {
        if(LoadTree(i) < 0)
            break;
        fChain->GetEntry(i);
        if(id == 13) // MC
        {
            Double_t weight
                = (MC_LIFE / simTau) * exp(-M0_time_true * (1. / simTau - 1. / MC_LIFE));
            hTime_true->Fill(M0_time_true / MC_LIFE, weight);
        }
    }

    for(int i = 0; i < 1e6; i++)
    {
        double x = rnd.Exp(simTau);
        hTime_gen->Fill(x / MC_LIFE);

        if(i % int(1e6) == 0)
            std::cout << "\rGenerated " << i << " events" << std::flush;
    }
    cout << endl;

    auto c = new TCanvas();
    c->cd();
    // hTime_gen->Draw("E HISTO");
    hTime_true->Draw("E HISTO SAME");
}

void lifetimeANA::AcceptanceScan(
    Float_t tauFactMin, Float_t tauFactMax, Float_t tauStep, UInt_t seed)
{
    SetLBStyle();
    gStyle->SetCanvasPreferGL();

    TRandom3 rnd = TRandom3(seed);

    // Prepariamo i Canvas e la Legenda
    TCanvas *cScan = new TCanvas("cScan", "Acceptance Scan", 800, 700);
    cScan->cd();
    TLegend *leg = new TLegend(0.65, 0.2, 0.88, 0.45); // Posizione: in basso a destra
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);

    TCanvas *cDist = new TCanvas("cDist", "Distributions Before Division", 1200, 600);
    cDist->Divide(2, 1);
    cDist->cd(1);
    gPad->SetLogy(); // Scala logaritmica per i Reco
    cDist->cd(2);
    gPad->SetLogy(); // Scala logaritmica per i Gen

    // 1. Calcola il numero esatto di step
    int nSteps = FloorNint((tauFactMax - tauFactMin) / tauStep) + 1;
    if(nSteps < 1)
        nSteps = 1;

    // 2. Calcola colonne e righe ottimali usando la radice quadrata
    int nCols = Ceil(Sqrt(nSteps));
    int nRows = Ceil((double)nSteps / nCols);

    // 3. Crea e dividi il canvas
    TCanvas *cPlat = new TCanvas("cPlat", "Plateau fits debug", 1200, 800);
    cPlat->Divide(nCols, nRows);

    // Array di colori belli per distinguere le curve (Nero, Rosso, Blu, Verde scuro, Magenta,
    // Ciano, Arancione)
    int colors[] = { kBlack, kRed, kBlue, kGreen + 2, kMagenta, kCyan + 2, kOrange + 7, kViolet + 2,
        kTeal + 2, kGray + 2 };
    int colorIdx = 0;
    bool isFirst = true;

    // Loop sui fattori moltiplicativi (es. da 0.5 a 1.5 con step 0.2)
    // Aggiungo 1e-5 a tauFactMax per evitare problemi di precisione dei float nel loop
    for(Float_t fact = tauFactMin; fact <= tauFactMax + 1e-5; fact += tauStep)
    {
        std::cout << "Elaborating tau = " << fact << " * MC_LIFE..." << std::endl;
        Double_t simTau = fact * MC_LIFE;

        // Nomi unici per gli istogrammi di questa iterazione
        TString hRecName = Form("hTime_Rec_%.2f", fact);
        TString hGenName = Form("hTime_Gen_%.2f", fact);

        TH1D *hTime_gen = new TH1D(hGenName, "Gen (Toys);t /#tau_{MC};Events", 200, 0., 10.);
        TH1D *hTime_rec = new TH1D(hRecName, "Reco (Reweighted);t /#tau_{MC};Events", 200, 0., 10.);
        hTime_gen->Sumw2();
        hTime_rec->Sumw2();

        // 1. LETTURA ALBERO E REWEIGHTING (Numeratore)
        for(int i = 0; i < fChain->GetEntriesFast(); i++)
        {
            if(LoadTree(i) < 0)
                break;
            fChain->GetEntry(i);
            if(id == 13) // MC selezionato
            {
                Double_t weight
                    = (MC_LIFE / simTau) * exp(-M0_time_true * (1. / simTau - 1. / MC_LIFE));
                hTime_rec->Fill(M0_time_true / MC_LIFE, weight);
            }
        }

        // 2. GENERAZIONE TOYS (Denominatore)
        int nToys = 1e9;
        for(int i = 0; i < nToys; i++)
        {
            double x = rnd.Exp(simTau); // Generiamo con il NUOVO tau
            hTime_gen->Fill(x / MC_LIFE);
        }

        // =======================================================
        // 3. PLOT PRIMA DELLA DIVISIONE
        // =======================================================
        int currentColor = colors[colorIdx % 7];

        TH1D *hPlotRec = (TH1D *)hTime_rec->Clone(Form("plotRec_%.2f", fact));
        TH1D *hPlotGen = (TH1D *)hTime_gen->Clone(Form("plotGen_%.2f", fact));

        hPlotRec->SetLineColor(currentColor);
        hPlotGen->SetLineColor(currentColor);
        hPlotRec->SetStats(0);
        hPlotGen->SetStats(0);

        cDist->cd(1);
        if(isFirst)
        {
            hPlotRec->SetMaximum(1e5); // Massimo fisso enorme per non tagliare niente
            hPlotRec->SetMinimum(0.1); // Minimo per evitare log(0)
            hPlotRec->Draw("HIST");
        }
        else
        {
            hPlotRec->Draw("HIST SAME");
        }

        cDist->cd(2);
        if(isFirst)
        {
            hPlotGen->SetMaximum(1e9);
            hPlotGen->SetMinimum(0.1);
            hPlotGen->Draw("HIST");
        }
        else
        {
            hPlotGen->Draw("HIST SAME");
        }

        // =======================================================
        // 4. CALCOLO ACCETTANZA
        // =======================================================
        // Dividiamo subito, così cerchiamo il plateau sull'accettanza vera
        hTime_rec->Divide(hTime_rec, hTime_gen, 1.0, 1.0);

        // =======================================================
        // 5. NORMALIZZAZIONE PLATEAU (Metodo della Moda visiva)
        // =======================================================
        double xMinPlateau = 0;
        double xMaxPlateau = 4;

        int binMin = hTime_rec->FindBin(xMinPlateau);
        int binMax = hTime_rec->FindBin(xMaxPlateau);

        // Cerchiamo il minimo e massimo dei valori in quel range per definire l'asse X di hMode
        double minVal = 1e9, maxVal = -1e9;
        std::vector<double> validValues;

        for(int b = binMin; b <= binMax; b++)
        {
            double val = hTime_rec->GetBinContent(b);
            if(val > 1e-5)
            {
                if(val < minVal)
                    minVal = val;
                if(val > maxVal)
                    maxVal = val;
                validValues.push_back(val);
            }
        }

        double plateauValue = 1.0;

        if(validValues.size() > 2)
        {
            cPlat->cd(colorIdx + 1); // Andiamo nel pad di debug

            // Istogramma dei valori dei bin
            TH1D *hMode = new TH1D(Form("hMode_%.2f", fact),
                Form("#tau = %.2f;Acceptance value;N. Bin", fact), 50, minVal * 0.9, maxVal * 1.1);

            for(double v : validValues)
            {
                hMode->Fill(v);
            }

            hMode->SetLineColor(currentColor);
            hMode->SetFillColorAlpha(currentColor, 0.3);
            hMode->Draw("HIST");

            // Troviamo il picco
            int maxBin = hMode->GetMaximumBin();
            plateauValue = hMode->GetBinCenter(maxBin);

            // Disegniamo una riga rossa verticale sul valore scelto
            TLine *lMode = new TLine(plateauValue, 0, plateauValue, hMode->GetMaximum());
            lMode->SetLineColor(kRed);
            lMode->SetLineWidth(2);
            lMode->Draw("SAME");
        }
        else if(validValues.size() > 0)
        {
            plateauValue = validValues[0];
        }

        // Scaliamo l'accettanza per il valore trovato
        if(plateauValue > 0)
        {
            hTime_rec->Scale(1.0 / plateauValue);
        }

        // =======================================================
        // 6. DISEGNO ACCETTANZA
        // =======================================================
        cScan->cd(); // Torniamo sulla canvas principale
        hTime_rec->GetYaxis()->SetTitle("Acceptance");
        hTime_rec->GetXaxis()->SetRangeUser(0, 4.1);
        hTime_rec->SetLineColor(currentColor);
        hTime_rec->SetMarkerColor(currentColor);
        // hTime_rec->SetFillColor(0);
        hTime_rec->SetMarkerStyle(20);
        hTime_rec->SetMarkerSize(0.6);
        hTime_rec->SetStats(0);

        if(isFirst)
        {
            hTime_rec->SetMaximum(1.5);
            hTime_rec->SetMinimum(0.0);
            hTime_rec->Draw("E HISTO");
            isFirst = false;
        }
        else
        {
            hTime_rec->Draw("E HISTO SAME");
        }

        leg->AddEntry(hTime_rec, Form("#tau = %.2f #tau_{MC}", fact), "lep");

        colorIdx++;
    } // <--- FINE DEL CICLO FOR

    cScan->cd();
    leg->Draw();

    cScan->Update();
    cPlat->Update();
    cDist->Update();

    if(savePlots)
    {
        cScan->SaveAs(Form("cScan_%.2f_%.2f.pdf", tauFactMin, tauFactMax));
        cPlat->SaveAs(Form("cPlat_%.2f_%.2f.pdf", tauFactMin, tauFactMax));
        cDist->SaveAs(Form("cDist_%.2f_%.2f.pdf", tauFactMin, tauFactMax));
    }

    std::cout << "Scan done." << std::endl;
}

// 3. CALCOLO ACCETTANZA
// Non uso l'opzione "B" perché i bin hanno pesi non interi. La propagazione standard di
// ROOT gestisce i pesi.
// hTime_rec->Divide(hTime_rec, hTime_gen, 1.0, 1.0);

// 4. NORMALIZZAZIONE PLATEAU
/*
double xMinPlateau = 2.5;
double xMaxPlateau = hTime_rec->GetXaxis()->GetXmax();

TF1 *fPlateau = new TF1(Form("fPlat_%.2f", fact), "pol0", xMinPlateau, xMaxPlateau);
hTime_rec->Fit(fPlateau, "RQ0"); // R=Range, Q=Quiet, 0=Non disegnare
double plateauValue = fPlateau->GetParameter(0);
*/

void lifetimeANA::CalibrateTimeBias()
{
    // Applichiamo il tuo stile pulito
    SetLBStyle();
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0); // Disattivato per ora, visto che non fittiamo

    // 1. DEFINIZIONE ISTOGRAMMA 2D
    // CORREZIONE: Impostiamo dei limiti reali sull'asse Y (es. da -0.1 a 0.1 ps)
    // altrimenti il TProfile fallisce.
    auto hRes2D_Reco = new TH2D("hRes2D_Reco",
        "Residuals vs Reconstructed Time;M0_time / MC_LIFE;M0_time - M0_time_true [ps]", 80, 0.0,
        1.0, // Asse X: da 0 a 1 vita media
        100, -0.2, 0.2 // Asse Y: limiti reali per i residui (in ps)
    );

    if(fChain == 0)
        return;
    Long64_t nentries = fChain->GetEntriesFast();

    for(Long64_t jentry = 0; jentry < nentries; jentry++)
    {
        if(LoadTree(jentry) < 0)
            break;
        fChain->GetEntry(jentry);
        if(id == 13) // Solo Monte Carlo
        {
            Double_t reco_ps = M0_time / MC_LIFE;
            Double_t res_ps = (M0_time - M0_time_true) * 1e12; // Residuo in ps

            // Riempiamo solo se siamo dentro la prima vita media
            if(reco_ps >= TIME_CUT && reco_ps <= 1.0)
            {
                hRes2D_Reco->Fill(reco_ps, res_ps);
            }
        }
    }

    // 2. Creazione del profilo
    TProfile *pBiasReco = hRes2D_Reco->ProfileX("pBiasReco");
    pBiasReco->SetTitle(
        "Time Bias (First Lifetime Only);M0_time / MC_LIFE;Bias (M0_time - M0_time_true) [ps]");
    pBiasReco->SetMarkerColor(kBlue + 1);
    pBiasReco->SetLineColor(kBlue + 1);
    pBiasReco->SetMarkerStyle(20);
    pBiasReco->SetMarkerSize(0.8);

    // Zoomiamo l'asse Y del profilo per vedere bene il trend del bias (da -0.01 a 0.10 ps)
    pBiasReco->GetYaxis()->SetRangeUser(-0.5, 0.5);

    // ========================================================================
    // 3. FUNZIONE (COMMENTATA / NON ATTIVA)
    // ========================================================================
    /*
    TF1 *fBias = new TF1("fBias",
        "x > [3] ? [0] * TMath::Power(x - [3], [1]) * TMath::Exp(-[2] * (x - [3])) : 0.0",
    0.0, 1.0); fBias->SetParameter(3, 0.01); fBias->SetParLimits(3, 0.0, 0.05);
    fBias->SetParameter(1, 1.5);
    fBias->SetParLimits(1, 0.5, 3.0);
    fBias->SetParameter(2, 12.0);
    fBias->SetParLimits(2, 2.0, 30.0);
    fBias->SetParameter(0, 10.0);
    fBias->SetLineColor(kRed);
    fBias->SetLineWidth(3);
    // pBiasReco->Fit(fBias, "L I R Q");
    */

    // 4. DISEGNO
    TCanvas *cCalib = new TCanvas("cCalib", "Bias Calibration Curve");
    pBiasReco->Draw("E1");

    // Linea dello zero per riferimento (da 0 a 1)
    TLine *l0 = new TLine(-0.5, 0, 1, 0);
    l0->SetLineColor(kBlack);
    l0->SetLineStyle(7);
    l0->Draw("same");

    if(savePlots)
    {
        cCalib->SaveAs("time_bias_calibration_1tau.pdf");
    }
}

//  ============================= PIPELINE ======================================

// ==============================================================================
// 1. FIT DELLA RISOLUZIONE (in unità normalizzate t/Tau)
// ==============================================================================
std::vector<double> lifetimeANA::FitResolutionNormalized()
{
    SetLBStyle();
    gStyle->SetOptFit(1111);

    // Asse X in unità normalizzate (t/tau_MC). Range +/- 0.5 vite medie.
    auto hRes = new TH1D("hRes_norm",
        "Time Resolution (Normalized);(t_{reco} - t_{true}) / #tau_{MC};Events", 100, -0.5, 0.5);

    if(fChain == 0)
        return {};
    for(Long64_t jentry = 0; jentry < fChain->GetEntriesFast(); jentry++)
    {
        if(LoadTree(jentry) < 0)
            break;
        fChain->GetEntry(jentry);
        if(id == 13)
        {
            Double_t t_reco = M0_time / MC_LIFE;
            // Usiamo lo STESSO TAGLIO che useremo nel fit finale
            // if(t_reco > TIME_CUT)
            hRes->Fill((M0_time - M0_time_true) / MC_LIFE);
        }
    }

    TF1 *fRes = new TF1("fRes_norm", f_2G_Frac, -0.5, 0.5, 6);
    fRes->SetParameters(hRes->Integral() * hRes->GetBinWidth(1), 0.0, 0.05, 0.0, 0.15, 0.7);
    fRes->SetParNames("Yield", "Mu_1", "Sigma_1", "Mu_2", "Sigma_2", "Frac_1");
    fRes->SetParLimits(5, 0.0, 1.0); // Frazione tra 0 e 1

    hRes->Fit(fRes, "L I R 0");

    // --- DISEGNO CON PULL ---
    TCanvas *c1 = new TCanvas("c_Res_Norm", "Resolution Fit", 800, 800);
    TPad *pad1 = new TPad("pad1_res", "", 0, 0.3, 1, 1.0);
    TPad *pad2 = new TPad("pad2_res", "", 0, 0.0, 1, 0.3);
    pad1->SetBottomMargin(0.02);
    pad2->SetTopMargin(0.02);
    pad2->SetBottomMargin(0.3);
    pad1->Draw();
    pad2->Draw();

    pad1->cd();
    hRes->GetXaxis()->SetLabelSize(0);
    hRes->GetXaxis()->SetTitleSize(0);
    hRes->Draw("E");
    fRes->Draw("same");

    TF1 *g1 = new TF1("g1_norm", "gaus", -0.5, 0.5);
    g1->SetParameters(fRes->GetParameter(0) * fRes->GetParameter(5)
            / (fRes->GetParameter(2) * TMath::Sqrt(TMath::TwoPi())),
        fRes->GetParameter(1), fRes->GetParameter(2));
    g1->SetLineStyle(2);
    g1->SetLineColor(kGreen + 2);
    g1->Draw("same");

    TF1 *g2 = new TF1("g2_norm", "gaus", -0.5, 0.5);
    g2->SetParameters(fRes->GetParameter(0) * (1.0 - fRes->GetParameter(5))
            / (fRes->GetParameter(4) * TMath::Sqrt(TMath::TwoPi())),
        fRes->GetParameter(3), fRes->GetParameter(4));
    g2->SetLineStyle(2);
    g2->SetLineColor(kCyan + 1);
    g2->Draw("same");

    pad2->cd();
    pad2->SetGridy();
    TH1D *hPull = (TH1D *)hRes->Clone("hPull_res");
    hPull->Reset();
    for(int i = 1; i <= hRes->GetNbinsX(); i++)
    {
        double err = hRes->GetBinError(i);
        if(err > 0)
            hPull->SetBinContent(
                i, (hRes->GetBinContent(i) - fRes->Eval(hRes->GetBinCenter(i))) / err);
    }
    hPull->GetYaxis()->SetTitle("Pull");
    hPull->GetYaxis()->SetTitleSize(gStyle->GetTitleSize("Y"));
    hPull->GetYaxis()->SetLabelSize(gStyle->GetLabelSize("Y"));
    hPull->GetYaxis()->SetTitleOffset(gStyle->GetTitleOffset("Y"));
    hPull->GetYaxis()->SetRangeUser(-5, 5);

    hPull->GetXaxis()->SetTitle("(t_{reco} - t_{true}) / #tau_{MC}");
    hPull->GetXaxis()->SetTitleSize(gStyle->GetTitleSize("X"));
    hPull->GetXaxis()->SetLabelSize(gStyle->GetLabelSize("X"));
    hPull->GetXaxis()->SetTitleOffset(gStyle->GetTitleOffset("X"));

    hPull->Draw("P");

    c1->Update();

    if(savePlots)
        c1->SaveAs("plot_1_Resolution.pdf");

    // Ritorna i parametri estratti
    return { fRes->GetParameter(1), fRes->GetParameter(2), fRes->GetParameter(3),
        fRes->GetParameter(4), fRes->GetParameter(5) }; // {mu1, sig1, mu2, sig2, frac1}
}

// ==============================================================================
// 2. FIT DELL'ACCETTANZA (in unità normalizzate t/Tau)
// ==============================================================================
std::vector<double> lifetimeANA::FitAcceptanceNormalized(UInt_t seed)
{
    SetLBStyle();
    TRandom3 rnd(seed);

    TH1D *hTime_gen = new TH1D("hTime_Gen_Acc", "", 100, 0., 10.);
    TH1D *hTime_rec
        = new TH1D("hTime_Rec_Acc", "Acceptance Fit;t / #tau_{MC};Acceptance", 100, 0., 10.);
    hTime_gen->Sumw2();
    hTime_rec->Sumw2();

    for(int i = 0; i < fChain->GetEntriesFast(); i++)
    {
        if(LoadTree(i) < 0)
            break;
        fChain->GetEntry(i);
        if(id == 13)
        {
            Double_t t_reco = M0_time / MC_LIFE;
            // if(t_reco > TIME_CUT)
            hTime_rec->Fill(M0_time_true / MC_LIFE);
        }
    }

    for(int i = 0; i < 5e7; i++)
        hTime_gen->Fill(rnd.Exp(1.0)); // 1.0 perché t è normalizzato a tau_MC

    hTime_rec->Divide(hTime_rec, hTime_gen, 1.0, 1.0, "B");

    // Normalizzazione al plateau
    TF1 *fPlateau = new TF1("fPlateau", "pol0", 2.5, 9.0);
    hTime_rec->Fit(fPlateau, "LRIQ0");
    if(fPlateau->GetParameter(0) > 0)
        hTime_rec->Scale(1.0 / fPlateau->GetParameter(0));

    // Funzione Accettanza (compatibile al 100% con FullPDF_ConvAcc)
    // Cambia il 4.0 in 10.0 per coprire tutto il dominio
    TF1 *fAcc = new TF1("fAcc_Fit",
        "0.5 * ( [0]*(1.0 + TMath::Erf((x-[1])/[2])) + (1.0-[0])*(1.0 + TMath::Erf((x-[3])/[4])) )",
        0.0, 10.0);
    fAcc->SetParameters(0.5, 0.5, 0.2, 1.5, 1.0);
    fAcc->SetParNames("Frac", "Mu_1", "Sig_1", "Mu_2", "Sig_2");
    fAcc->SetParLimits(0, 0.0, 1.0);

    hTime_rec->Fit(fAcc, "LRI0");

    TCanvas *c2 = new TCanvas("c_Acc", "Acceptance Fit", 800, 600);
    hTime_rec->SetMarkerStyle(20);
    hTime_rec->GetYaxis()->SetRangeUser(0, 1.2);
    hTime_rec->Draw("E");
    fAcc->SetLineColor(kRed);
    fAcc->Draw("same");

    if(savePlots)
        c2->SaveAs("plot_2_Acceptance.pdf");

    return { fAcc->GetParameter(0), fAcc->GetParameter(1), fAcc->GetParameter(2),
        fAcc->GetParameter(3), fAcc->GetParameter(4) }; // {frac, mu1, sig1, mu2, sig2}
}

// ==============================================================================
// 3. FIT GLOBALE CONVOLUTO
// ==============================================================================
void lifetimeANA::RunGlobalFit()
{
    SetLBStyle();

    std::cout << "\n=== STEP 1: Fitting Resolution ===" << std::endl;
    std::vector<double> resPars = FitResolutionNormalized();

    std::cout << "\n=== STEP 2: Fitting Acceptance ===" << std::endl;
    std::vector<double> accPars = FitAcceptanceNormalized(42);

    std::cout << "\n=== STEP 3: Final Convolved Fit ===" << std::endl;

    gStyle->SetOptFit(1111);
    auto hTime_reco
        = new TH1D("hTime_reco_final", "Final Time Fit;t / #tau_{MC};Events", 100, 0., 10.);
    hTime_reco->Sumw2();

    for(Long64_t jentry = 0; jentry < fChain->GetEntriesFast(); jentry++)
    {
        if(LoadTree(jentry) < 0)
            break;
        fChain->GetEntry(jentry);
        if(id == 13)
        {
            Double_t t_reco = M0_time / MC_LIFE;
            // if(t_reco > TIME_CUT)
            hTime_reco->Fill(t_reco);
        }
    }

    Int_t nGaus = 2;
    Bool_t useAcceptance = true;
    FullPDF_ConvAcc fModel(nGaus, useAcceptance);

    Double_t minT = 0.5;
    Double_t maxT = 10.0;
    auto fFit = new TF1("fFitFinal", fModel, minT, maxT, fModel.GetNPar());
    fFit->SetNpx(1000);

    // --- MAPPATURA ESATTA DEI PARAMETRI DELLA CLASSE FullPDF_ConvAcc ---
    // [0] : Tau (Libero, parametro fisico)
    fFit->SetParameter(0, 1.0);
    fFit->SetParName(0, "Tau");

    // [1], [2] : Mu_1, Mu_2 (Fissati dalla risoluzione)
    fFit->FixParameter(1, resPars[0]);
    fFit->SetParName(1, "Res_Mu1");
    fFit->FixParameter(2, resPars[2]);
    fFit->SetParName(2, "Res_Mu2");

    // [3] : Yield (Libero)
    fFit->SetParameter(3, hTime_reco->Integral() * hTime_reco->GetBinWidth(1));
    fFit->SetParName(3, "Yield");

    // [4], [5] : Sigma_1, Sigma_2 (Fissati dalla risoluzione)
    fFit->FixParameter(4, resPars[1]);
    fFit->SetParName(4, "Res_Sig1");
    fFit->FixParameter(5, resPars[3]);
    fFit->SetParName(5, "Res_Sig2");

    // [6] : Frazione G1 (Fissata dalla risoluzione)
    fFit->FixParameter(6, resPars[4]);
    fFit->SetParName(6, "Res_Frac1");

    // [7] a [11] : Parametri Accettanza (Fissati dall'accettanza)
    fFit->FixParameter(7, accPars[0]);
    fFit->SetParName(7, "Acc_Frac");
    fFit->FixParameter(8, accPars[1]);
    fFit->SetParName(8, "Acc_Mu1");
    fFit->FixParameter(9, accPars[2]);
    fFit->SetParName(9, "Acc_Sig1");
    fFit->FixParameter(10, accPars[3]);
    fFit->SetParName(10, "Acc_Mu2");
    fFit->FixParameter(11, accPars[4]);
    fFit->SetParName(11, "Acc_Sig2");

    // Eseguiamo il fit! (L=Likelihood, I=Integral over bins, R=Range)
    hTime_reco->Fit(fFit, "L I R 0");

    // --- DISEGNO FINALE BELLISSIMO (LogY + Pulls) ---
    TCanvas *c_fit = new TCanvas("c_fit_final", "Final Convolved Fit", 1000, 1000);
    TPad *pad1 = new TPad("pad1", "", 0, 0.3, 1, 1.0);
    TPad *pad2 = new TPad("pad2", "", 0, 0.0, 1, 0.3);
    pad1->SetBottomMargin(0.03);
    pad1->SetLogy(); // Log scale per vederci chiaro
    pad2->SetTopMargin(0.02);
    pad2->SetBottomMargin(0.35); // necessary for x-axis
    pad1->Draw();
    pad2->Draw();

    pad1->cd();
    hTime_reco->GetYaxis()->SetTitle("Events");
    hTime_reco->GetXaxis()->SetLabelSize(0);
    hTime_reco->GetXaxis()->SetTitleSize(0);
    hTime_reco->SetMinimum(0.5);
    hTime_reco->SetMarkerStyle(20);
    hTime_reco->Draw("E");
    fFit->SetLineColor(kBlue + 1);
    fFit->SetLineWidth(3);
    fFit->Draw("same");

    pad2->cd();
    pad2->SetGridy();
    TH1D *hPull = (TH1D *)hTime_reco->Clone("hPull_final");
    hPull->Reset();
    hPull->SetStats(0);
    for(int i = 1; i <= hTime_reco->GetNbinsX(); i++)
    {
        double x = hTime_reco->GetBinCenter(i);
        if(x < minT || x > maxT)
            continue; // Calcola i pull solo nel range fittato
        double obs = hTime_reco->GetBinContent(i);
        double err = hTime_reco->GetBinError(i);
        if(err > 0)
            hPull->SetBinContent(i, (obs - fFit->Eval(x)) / err);
    }

    hPull->GetYaxis()->SetTitle("Pull");
    hPull->GetYaxis()->SetTitleSize(gStyle->GetTitleSize("Y"));
    hPull->GetYaxis()->SetLabelSize(gStyle->GetLabelSize("Y"));
    hPull->GetYaxis()->SetTitleOffset(gStyle->GetTitleOffset("Y"));

    hPull->GetXaxis()->SetTitle("t / #tau_{MC}");
    hPull->GetXaxis()->SetTitleSize(gStyle->GetTitleSize("X"));
    hPull->GetXaxis()->SetLabelSize(gStyle->GetLabelSize("X"));
    hPull->GetXaxis()->SetTitleOffset(gStyle->GetTitleOffset("X"));

    hPull->GetYaxis()->SetNdivisions(505);
    hPull->GetYaxis()->SetRangeUser(-5.0, 5.0);
    hPull->SetMarkerStyle(20);
    hPull->SetFillColor(kBlue - 7);
    hPull->Draw("HIST P");

    if(savePlots)
        c_fit->SaveAs("plot_3_FinalConvolvedFit.pdf");

    std::cout << "\n=============================================" << std::endl;
    std::cout << " FIT CONCLUSO CON SUCCESSO!" << std::endl;
    std::cout << " Tempo MC atteso: 1.000" << std::endl;
    std::cout << " Tempo Fit      : " << fFit->GetParameter(0) << " +/- " << fFit->GetParError(0)
              << std::endl;
    std::cout << "=============================================\n" << std::endl;
}
