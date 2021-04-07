
void canvasPhiPQ() {
    TFile *input = TFile::Open("/home/claudio/work/clas-data/PhiPQAnalysis/histsPhiPQ.root");

    // Creating canvases
    TCanvas *cD = new TCanvas("cD","Pion PhiPQ, D target");
    TCanvas *cX = new TCanvas("cX","Pion PhiPQ, Fe target");
    cD->Divide(3,2);
    cX->Divide(3,2);

    int Nhist = 3;

    // Copying histograms
    for (int j=0; j<Nhist; j++){
        TH1F *histogramD0 = (TH1F*)input->Get(Form("histPhiPQ_D_%d0",j));
        cD->cd(j+1);
        histogramD0->Draw();
        TH1F *histogramD1 = (TH1F*)input->Get(Form("histPhiPQ_D_%d1",j));
        cD->cd(j+1+Nhist);
        histogramD1->Draw();

        TH1F *histogramX0 = (TH1F*)input->Get(Form("histPhiPQ_Fe_%d0",j));
        cX->cd(j+1);
        histogramX0->Draw();
        TH1F *histogramX1 = (TH1F*)input->Get(Form("histPhiPQ_Fe_%d1",j));
        cX->cd(j+1+Nhist);
        histogramX1->Draw();
    }
}
