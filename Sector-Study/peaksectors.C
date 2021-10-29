
void peaksectors(){
    TFile *input = TFile::Open("../analysis/root-files/AnalysisCSMV_Fe2.root","READ");
    TH1F *hdata = (TH1F*) input->Get("hdata");
    TH1F *hcorr = (TH1F*) input->Get("hcorr");

    int Limits[] = {0, 20, 50, 80, 100, 130, 160, 180};
    double x0[2][7], y1[2][7];

    TH1F *hdata_max = new TH1F("hdata_max","Relative maximum",180,-180,180);
    TH1F *hcorr_max = new TH1F("hcorr_max","Relative maximum",180,-180,180);
    hdata_max->SetMarkerStyle(8);
    hdata_max->SetMarkerColor(kRed+2);
    hdata_max->SetLineColor(kRed+1);
    hcorr_max->SetMarkerStyle(8);
    hcorr_max->SetMarkerColor(kRed+2);
    hcorr_max->SetLineColor(kRed+1);

    TCanvas *c1 = new TCanvas("c1","Data");
    TCanvas *c2 = new TCanvas("c2","Corrected");

    c1->cd();
    hdata->Draw();
    for (int i=0; i<7; i++){
        hdata->GetXaxis()->SetRange(Limits[i], Limits[i+1]);
        int maxbin = hdata->GetMaximumBin();
        x0[0][i] = hdata->GetBinCenter(maxbin);
        y1[0][i] = hdata->GetBinContent(maxbin);
        std::cout << "Data: Bin of maximum number " << i+1 << ": " << maxbin << " (value: " << y1[0][i] << ")" << std::endl;
        hdata_max->Fill(x0[0][i],y1[0][i]);
    }
    hdata->GetXaxis()->SetRange();
    hdata_max->Draw("SAME");

    c2->cd();
    hcorr->Draw();
    for (int i=0; i<7; i++){
        hcorr->GetXaxis()->SetRange(Limits[i], Limits[i+1]);
        int maxbin = hcorr->GetMaximumBin();
        x0[1][i] = hcorr->GetBinCenter(maxbin);
        y1[1][i] = hcorr->GetBinContent(maxbin);
        std::cout << "Corrected: Bin of maximum number " << i+1 << ": " << maxbin << " (value: " << y1[1][i] << ")" << std::endl;
        hcorr_max->Fill(x0[1][i],y1[1][i]);
    }
    hcorr->GetXaxis()->SetRange();
    hcorr_max->Draw("SAME");
}
