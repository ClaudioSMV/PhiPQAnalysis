
void quantiles() {
    TFile *input = TFile::Open("../data_Fe_skimmed.root");
    TTree *tree = (TTree*) input->Get("ntuple_data");

    float Q2 = 99;
    float TargType = 99;
    float pid = 99;

    tree->SetBranchAddress("Q2",&Q2);
    tree->SetBranchAddress("TargType",&TargType);
    tree->SetBranchAddress("pid",&pid);

    TH1F *Q2hist = new TH1F("Q2hist","Q2 hist",100,0,5);
    TH1F *Q2hist_D = new TH1F("Q2hist_D","Q2 hist D",100,0,5);
    TH1F *Q2hist_X = new TH1F("Q2hist_X","Q2 hist X",100,0,5);

    int Nentries = tree->GetEntries();
    for (int i = 0; i<Nentries; i++){
        // if (i == 25) break;

        tree->GetEntry(i);

        if (pid != 211) continue;

        if (TargType == 1){
            Q2hist->Fill(Q2);
            Q2hist_D->Fill(Q2);
        }
        if (TargType == 2){
            Q2hist->Fill(Q2);
            Q2hist_X->Fill(Q2);
        }
    }

    Q2hist->Draw();

    double xq[3] = {1./3., 2./3., 3./3.};
    double yq[3];
    double yqD[3];
    double yqX[3];

    Q2hist->GetQuantiles(3,yq,xq);
    Q2hist_D->GetQuantiles(3,yqD,xq);
    Q2hist_X->GetQuantiles(3,yqX,xq);

    for (int i=0; i<3; i++){
        std::cout << "The quantile number " << i+1 << " is " << yq[i] << ", D: " << yqD[i] << ", Fe: " << yqX[i] << std::endl;
    }

}
