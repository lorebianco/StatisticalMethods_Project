{
    // 1. Carica il file con l'implementazione (assicurati che il percorso sia corretto)
    // Puoi usare ".L lifetimeANA.C+" se desideri compilarlo con ACLiC, oppure solo ".L
    // lifetimeANA.C"
    gROOT->ProcessLine(".L lifetimeANA.C+");

    // 2. Dichiara l'oggetto globale 't' in modo che sia disponibile sulla riga di comando di ROOT
    gROOT->ProcessLine("lifetimeANA t;");

    // Opzionale: un messaggio di conferma
    std::cout << "[D0_Tau] File lifetimeANA.C loaded and instance 't' created." << std::endl;
}
