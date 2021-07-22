// Nphe vs kinematic variables, one run/file per target

void canvasnphe() {
    std::string target = "Fe";

    TFile *input_sim = TFile::Open(Form("Nphevs_%ssim.root",target.c_str()),"READ");
    TFile *input_data = TFile::Open(Form("Nphevs_%sdata.root",target.c_str()),"READ");

    // hvec: {hQ2, hNu, hZh, hPt2, hPhiPQ};
    std::string var[] = {"Q2", "Xb", "Zh", "Pt2", "PhiPQ"};
    std::string var_unit[] = {"Q^{2} [GeV^{2}]", "X_{b}", "Z_{h}", "P_{t}^{2} [GeV^{2}]", "#phi_{PQ} [deg]"};

    // Creating canvases
    std::vector<TCanvas*> canv;
    int Nbins = sizeof(var)/sizeof(var[0]);
    for (int i=0; i<Nbins; i++){
        TCanvas *c = new TCanvas(Form("c%i",i),Form("Nphe vs %i",i));
        c->Divide(2,1);
        canv.push_back(c);
    }

    // Copying histograms
    for (int j=0; j<Nbins; j++){
        TH2F *hist_sim = (TH2F*)input_sim->Get(Form("h%s",var[j].c_str()));
        canv[j]->cd(1);
        hist_sim->SetStats(0);
        hist_sim->GetYaxis()->SetTitle(Form("%s",var_unit[j].c_str()));
        hist_sim->SetTitle(Form("Simulation, %s target",target.c_str()));
        hist_sim->Draw("COLZ");

        TH2F *hist_data = (TH2F*)input_data->Get(Form("h%s",var[j].c_str()));
        canv[j]->cd(2);
        hist_data->SetStats(0);
        hist_data->GetYaxis()->SetTitle(Form("%s",var_unit[j].c_str()));
        hist_data->SetTitle(Form("Data, %s target",target.c_str()));
        hist_data->Draw("COLZ");

        if (j==0) canv[j]->Print(Form("Nphevs_%scuts.pdf(",target.c_str()),"pdf");
        else if (j!=(Nbins-1)) canv[j]->Print(Form("Nphevs_%scuts.pdf",target.c_str()),"pdf");
        else if (j==(Nbins-1)) canv[j]->Print(Form("Nphevs_%scuts.pdf)",target.c_str()),"pdf");
        else std::cout << "Error in loop filling with histograms!" << std::endl;
    }
}
