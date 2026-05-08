#include "RtypesCore.h"
#define lifetimeANA_cxx
#include <iostream>

#include <TCanvas.h>
#include <TF1.h>
#include <TH1.h>
#include <TH2D.h>
#include <TLine.h>
#include <TMath.h>
#include <TPaveStats.h>
#include <TProfile.h>
#include <TStyle.h>

#include "lifetimeANA.h"

using namespace std;
using namespace TMath;

constexpr Double_t MC_LIFE = 410.3e-15; // s

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
    c_histo_data_MKpi->Print("./_fig/c_histo_data_MKpi.pdf");
    c_histo_data_MKpi->Print("./_fig/c_histo_data_MKpi.eps");

    // Create Canvas
    TCanvas *c_histo_mc_MKpi = new TCanvas("c_histo_mc_MKpi", "canvas histo", 500, 500);
    c_histo_mc_MKpi->cd();
    histo_mc_MKpi->Draw();
    c_histo_mc_MKpi->Print("./_fig/c_histo_mc_MKpi.pdf");
    c_histo_mc_MKpi->Print("./_fig/c_histo_mc_MKpi.eps");

    // Create Canvas
    TCanvas *c_histo_mc_time = new TCanvas("c_histo_mc_time", "canvas histo", 500, 500);
    c_histo_mc_time->cd();
    histo_mc_time->Draw();
    c_histo_mc_time->Print("./_fig/c_histo_mc_time.pdf");
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

void lifetimeANA::Resolution(Int_t nGauss, Bool_t useSharedMean, Bool_t saveGraphs)
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
            hRes->Fill((M0_time - M0_time_true) * 1e12);
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
    if(saveGraphs)
        c1->SaveAs(Form(
            "fit_resolution_%dGauss_%sMean.png", nGauss, useSharedMean ? "Shared" : "Independent"));
}

void lifetimeANA::Resolution2D(Bool_t saveGraphs)
{
    // Stile globale per i plot
    gStyle->SetOptStat(0);
    gStyle->SetPalette(kBird);
    gStyle->SetTitleSize(0.04, "XYZ");
    gStyle->SetLabelSize(0.035, "XYZ");

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
    pBias->SetMarkerStyle(20);
    pBias->SetMarkerSize(1.0);
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
    hStdDev->SetMarkerStyle(20);
    hStdDev->SetMarkerSize(1.0);
    hStdDev->SetMarkerColor(kRed + 1);
    hStdDev->SetLineColor(kRed + 1);
    hStdDev->SetMinimum(0);

    // --- CANVAS 1: TH2D COLZ ---
    TCanvas *c_2D = new TCanvas("c_2D", "Residuals 2D", 800, 700);
    c_2D->SetRightMargin(0.15); // Spazio per la barra colore
    c_2D->SetLeftMargin(0.13);
    hRes2D->Draw("COLZ");
    if(saveGraphs)
        c_2D->SaveAs("plot_2D_residuals.png");

    // --- CANVAS 2: BIAS (Mean) ---
    TCanvas *c_Bias = new TCanvas("c_Bias", "Bias Plot", 800, 700);
    c_Bias->SetLeftMargin(0.13);
    pBias->GetYaxis()->SetTitleOffset(1.3);
    pBias->Draw("E1");

    TLine *l0 = new TLine(0, 0, 2, 0);
    l0->SetLineColor(kBlack);
    l0->SetLineStyle(7); // Tratteggiato
    l0->Draw();
    if(saveGraphs)
        c_Bias->SaveAs("plot_2D_bias.png");

    // --- CANVAS 3: RESOLUTION (Std Dev) ---
    TCanvas *c_StdDev = new TCanvas("c_StdDev", "Resolution Plot", 800, 700);
    c_StdDev->SetLeftMargin(0.13);
    hStdDev->GetYaxis()->SetTitleOffset(1.3);
    hStdDev->Draw("E1");
    if(saveGraphs)
        c_StdDev->SaveAs("plot_2D_resolution.png");
}

void lifetimeANA::ConvolvedFit(Bool_t saveGraphs)
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
    Int_t nGaus = 1;
    ExpoMultiGaussConv FConv(nGaus);
    Int_t nPar = FConv.GetNPar();

    // Fit function
    auto fConv = new TF1("fConv", FConv, 0., 8., nPar);

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

    hTime_true->Fit(fConv, "L I R");

    if(saveGraphs)
    {
        TCanvas *c_fit = new TCanvas("c_fit", "Convolved Fit", 800, 700);
        hTime_true->Draw("E HISTO");
        fConv->Draw("same");
        c_fit->SaveAs("plot_conv_fit.png");
    }
}
