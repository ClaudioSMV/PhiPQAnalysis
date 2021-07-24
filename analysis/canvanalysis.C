
void canvanalysis(TString target = "Fe", TString method = "CSMV", TString version = "2"){
    TString file_name = "Analysis"+method+"_"+target+version;
    TString title = "./root-files/"+file_name+".root";
    TFile *input = TFile::Open(title,"READ");

    // Draw canvases
    TCanvas *cdata = new TCanvas("cdata","Data");
    TCanvas *ccorr = new TCanvas("ccorr","Corrected");

    // Saving histograms
    TString hdata_name = "hdata";
    TString hcorr_name = "hcorr";
    if (method=="SM"){
        hdata_name += "_fin";
        hcorr_name += "_fin";
    }

    cdata->cd();
    (input->Get(hdata_name))->Draw();
    ccorr->cd();
    (input->Get(hcorr_name))->Draw();

    cdata->Print(file_name+".pdf(","pdf");
    ccorr->Print(file_name+".pdf)","pdf");
}
