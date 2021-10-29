// Four targets acceptance comparison side by side

void canvasSector(TString target = "Fe", TString number = "1"){

    TFile *input = TFile::Open("./root-files/AccSector_"+target+"1.root","READ");

    // Creating canvases
    TCanvas *creco1 = new TCanvas("creco1","Reconstructed1 Fe");
    TCanvas *creco2 = new TCanvas("creco2","Reconstructed2 Fe");
    TCanvas *ctrue1 = new TCanvas("ctrue1","Thrown1 (Monte Carlo) Fe");
    TCanvas *ctrue2 = new TCanvas("ctrue2","Thrown2 (Monte Carlo) Fe");
    creco1->Divide(3,2);
    creco2->Divide(3,2);
    ctrue1->Divide(3,2);
    ctrue2->Divide(3,2);


    TH1F *hrecoTot = new TH1F("hrecoTot","Total Reconstructed",180,-180,180);
    TH1F *htrueTot = new TH1F("htrueTot","Total Thrown",180,-180,180);
    // Copying histograms
    for (int j=0; j<6; j++){
        creco1->cd(j+1);
        TH1F *hreco1 = (TH1F*) input->Get(Form("hrecoSec%i",j));
        hreco1->GetXaxis()->SetTitle("#phi_{PQ}");
        hreco1->Draw();
        hrecoTot->Add(hreco1);
        creco2->cd(j+1);
        if (j!=5){
            TH1F *hreco2 = (TH1F*) input->Get(Form("hrecoSec%i",j+6));
            hreco2->GetXaxis()->SetTitle("#phi_{PQ}");
            hreco2->Draw();
            hrecoTot->Add(hreco2);
        }
        ctrue1->cd(j+1);
        TH1F *htrue1 = (TH1F*) input->Get(Form("htrueSec%i",j));
        htrue1->GetXaxis()->SetTitle("#phi_{PQ}");
        htrue1->Draw();
        htrueTot->Add(htrue1);
        ctrue2->cd(j+1);
        if (j!=5){
            TH1F *htrue2 = (TH1F*) input->Get(Form("htrueSec%i",j+6));
            htrue2->GetXaxis()->SetTitle("#phi_{PQ}");
            htrue2->Draw();
            htrueTot->Add(htrue2);
        }
    }
    creco2->cd(6);
    hrecoTot->GetXaxis()->SetTitle("#phi_{PQ}");
    hrecoTot->Draw();
    ctrue2->cd(6);
    htrueTot->GetXaxis()->SetTitle("#phi_{PQ}");
    htrueTot->Draw();

    TString title = "PlotSectors"+number+".pdf";
    creco1->Print(title+"(","pdf");
    creco2->Print(title,"pdf");
    ctrue1->Print(title,"pdf");
    ctrue2->Print(title+")","pdf");

}
