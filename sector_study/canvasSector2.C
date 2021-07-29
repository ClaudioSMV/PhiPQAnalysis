// Four targets acceptance comparison side by side

void canvasSector2(TString target = "Fe", TString number = "2"){

    TFile *input = TFile::Open("./root-files/AccSector_"+target+"1.root","READ");

    // Creating canvases
    TCanvas *creco = new TCanvas("creco","Reconstructed Fe");
    TCanvas *ctrue = new TCanvas("ctrue","Throw (Monte Carlo) Fe");
    TCanvas *ctots = new TCanvas("ctots","Totals Fe");
    creco->Divide(3,2);
    ctrue->Divide(3,2);
    ctots->Divide(2,1);

    Color_t color[6] = {kRed, kMagenta, kBlue, kCyan, kGreen, kOrange};

    TH1F *hrecoTot = new TH1F("hrecoTot","Total Reconstructed",180,-180,180);
    hrecoTot->GetXaxis()->SetTitle("#phi_{PQ}");
    TH1F *htrueTot = new TH1F("htrueTot","Total Thrown",180,-180,180);
    htrueTot->GetXaxis()->SetTitle("#phi_{PQ}");

    THStack *hreco_stack = new THStack("hreco_stack","");
    THStack *htrue_stack = new THStack("htrue_stack","");

    TString histTitle = "SectorEl - Sector = ";
    // Copying histograms
    for (int j=0; j<6; j++){
        creco->cd(j+1);
        TH1F *hreco_tmp = new TH1F(Form("hr%i",j),"",180,-180,180);
        TH1F *hreco1 = (TH1F*) input->Get(Form("hrecoSec%i",j));
        if (j!=5){
            TH1F *hreco2 = (TH1F*) input->Get(Form("hrecoSec%i",j+6));
            hreco_tmp->Add(hreco1,hreco2);
            hreco_tmp->SetTitle(Form("Reco "+histTitle+"%i, %i;#phi_{PQ}",j-5,j+1));
        }
        else {
            hreco_tmp->Add(hreco1);
            hreco_tmp->SetTitle(Form("Reco "+histTitle+"%i;#phi_{PQ}",0));
        }
        hreco_tmp->SetFillColorAlpha(color[j],1.);
        hreco_tmp->DrawCopy();
        hrecoTot->Add(hreco_tmp);

        // ctots->cd(1);
        hreco_tmp->SetFillColorAlpha(color[j],0.3);
        hreco_stack->Add(hreco_tmp);
        // hreco_tmp->Draw("SAME");
        
        ctrue->cd(j+1);
        TH1F *htrue_tmp = new TH1F(Form("ht%i",j),"",180,-180,180);
        TH1F *htrue1 = (TH1F*) input->Get(Form("htrueSec%i",j));
        if (j!=5){
            TH1F *htrue2 = (TH1F*) input->Get(Form("htrueSec%i",j+6));
            htrue_tmp->Add(htrue1,htrue2);
            htrue_tmp->SetTitle(Form("Thrown "+histTitle+"%i, %i;#phi_{PQ}",j-5,j+1));
        }
        else {
            htrue_tmp->Add(htrue1);
            htrue_tmp->SetTitle(Form("Thrown "+histTitle+"%i;#phi_{PQ}",0));
        }
        htrue_tmp->SetFillColorAlpha(color[j],1.);
        htrue_tmp->DrawCopy();
        htrueTot->Add(htrue_tmp);

        // ctots->cd(2);
        htrue_tmp->SetFillColorAlpha(color[j],0.3);
        htrue_stack->Add(htrue_tmp);
        // htrue_tmp->Draw("SAME");
    }

    ctots->cd(1);
    hreco_stack->Add(hrecoTot);
    hreco_stack->SetTitle("Reco Total;#phi_{PQ}");
    hreco_stack->Draw("NOSTACK");
    ctots->cd(2);
    htrue_stack->Add(htrueTot);
    htrue_stack->SetTitle("Thrown Total;#phi_{PQ}");
    htrue_stack->Draw("NOSTACK");

    TString title = "PlotSectors"+number+".pdf";
    creco->Print(title+"(","pdf");
    ctrue->Print(title,"pdf");
    ctots->Print(title+")","pdf");

}
