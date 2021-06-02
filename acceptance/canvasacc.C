
void canvasacc() {
    TFile *input = TFile::Open("Analysis_Fe1_2.root");

    // Creating canvases
    TCanvas *craw = new TCanvas("craw","Raw data PhiPQ, Fe target");
    TCanvas *ccorr = new TCanvas("ccorr","Corrected data PhiPQ, Fe target");
    craw->Divide(3,2);
    ccorr->Divide(3,2);

    float limits[4] = {0., 1.35, 1.82, 5.};
    int Nhist = 6;

    // Copying histograms
    for (int j=0; j<Nhist; j++){
        TH1F *histraw = (TH1F*)input->Get(Form("hraw1%i",j));
        craw->cd(j+1);
        histraw->SetTitle(Form("%.2f < Q^{2} < %.2f, %i #leq Nphe < %i;#phi_{PQ} [deg];Counts",limits[1],limits[2],5*j,5*(j+1)));
        if (j==5) {histraw->SetTitle(Form("%.2f < Q^{2} < %.2f, %i #leq Nphe;#phi_{PQ} [deg];Counts",limits[1],limits[2],5*j));}
        histraw->Draw();

        TH1F *histcorr = (TH1F*)input->Get(Form("hcorr1%i",j));
        ccorr->cd(j+1);
        histcorr->SetTitle(Form("%.2f < Q^{2} < %.2f, %i #leq Nphe < %i;#phi_{PQ} [deg];Counts",limits[1],limits[2],5*j,5*(j+1)));
        if (j==5) {histcorr->SetTitle(Form("%.2f < Q^{2} < %.2f, %i #leq Nphe;#phi_{PQ} [deg];Counts",limits[1],limits[2],5*j));}
        histcorr->Draw();
    }
}
