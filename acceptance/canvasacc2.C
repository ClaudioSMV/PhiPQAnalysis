
void canvasacc2() {

    std::string names[] = {"D", "C", "Fe", "Pb"};
    std::vector<TFile*> input_vec;

    for (int i=0; i<4; i++){
        TFile *input = TFile::Open(Form("Acc_%s1_cuts.root",names[i].c_str()),"READ");
        input_vec.push_back(input);
    }

    // Creating canvases
    TCanvas *creco = new TCanvas("creco","Reconstructed");
    TCanvas *ctrue = new TCanvas("ctrue","Thrown (Monte Carlo)");
    TCanvas *cacc  = new TCanvas("cacc" ,"Acceptance");
    // TCanvas *cacc  = new TCanvas("cacc" ,Form("Acceptance, %s target",target.c_str()));
    creco->Divide(2,2);
    ctrue->Divide(2,2);
    cacc->Divide(2,2);

    // Copying histograms
    for (int j=0; j<4; j++){
        creco->cd(j+1);
        TH1F *histreco = (TH1F*)input_vec[j]->Get("hreco");
        histreco->SetTitle(Form("Reconstructed, %s target;#phi_{PQ} [deg];Counts",names[j].c_str()));
        histreco->Draw();

        TH1F *histtrue = (TH1F*)input_vec[j]->Get("htrue");
        ctrue->cd(j+1);
        histtrue->SetTitle(Form("Thrown (Monte Carlo), %s target;#phi_{PQ} [deg];Counts",names[j].c_str()));
        histtrue->Draw();

        TH1F *histacc = (TH1F*)input_vec[j]->Get("hacc");
        cacc->cd(j+1);
        histacc->SetTitle(Form("Acceptance, %s target;#phi_{PQ} [deg];Counts",names[j].c_str()));
        histacc->Draw();
    }

    creco->Print("PlotAcceptance_all.pdf(","pdf");
    ctrue->Print("PlotAcceptance_all.pdf","pdf");
    cacc->Print("PlotAcceptance_all.pdf)","pdf");

}
