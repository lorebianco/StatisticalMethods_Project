void plot_smearing() {
  gStyle->SetCanvasPreferGL(true);
  // --- 1. Parametri della simulazione ---
  double tau = 1.5;      // Tempo di vita medio vero (es. 1.5 picosecondi)
  double sigma = 0.3;    // Risoluzione temporale del rivelatore (0.3 ps)
  int n_events = 500000; // Numero di eventi da generare

  // Stile generale per rendere il plot più pulito
  gStyle->SetOptStat(0);

  // --- 2. Creazione degli istogrammi ---
  // Estendiamo l'asse X da -2 a 10 per vedere i tempi negativi
  TH1F *h_true = new TH1F("h_true", ";Time [ps];Events", 120, -2.0, 10.0);
  TH1F *h_smeared = new TH1F("h_smeared", ";Time [ps];Events", 120, -2.0, 10.0);

  // --- 3. Generazione degli eventi (Il nostro mini Monte Carlo) ---
  // Usiamo gRandom che è il generatore di default di ROOT
  gRandom->SetSeed(0);

  for (int i = 0; i < n_events; ++i) {
    // A. Generiamo il tempo fisico (vero) che segue l'esponenziale
    double t_true = gRandom->Exp(tau);

    // B. Generiamo l'errore di misura causato dal rivelatore (Gaussiana
    // centrata in 0)
    double error = gRandom->Gaus(0, sigma);

    // C. Il tempo che noi misureremmo è la somma dei due
    double t_meas = t_true + error;

    // Riempiamo gli istogrammi
    h_true->Fill(t_true);
    h_smeared->Fill(t_meas);
  }

  // --- 4. Estetica e Colori ---
  h_true->SetLineColor(kBlue);
  h_true->SetLineWidth(2);

  h_smeared->SetLineColor(kRed);
  h_smeared->SetLineWidth(2);
  // Riempiamo l'istogramma smeared per renderlo più visibile
  h_smeared->SetFillColorAlpha(kRed, 0.3);

  // --- 5. Disegno sul Canvas ---
  TCanvas *c1 = new TCanvas("c1", "Smearing", 800, 600);

  // Disegniamo prima lo smeared perché il picco del true è più alto
  h_true->Draw("HIST");
  h_smeared->Draw("HIST SAME");

  // --- 6. Aggiunta della Legenda ---
  TLegend *leg = new TLegend(0.55, 0.7, 0.88, 0.88);
  leg->AddEntry(h_true, Form("True Time (Exp, #tau=%.1f ps)", tau), "l");
  leg->AddEntry(h_smeared, Form("Measured Time (#sigma=%.1f ps)", sigma), "f");
  leg->SetBorderSize(0);
  leg->Draw();
}
