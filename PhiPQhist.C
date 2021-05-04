
void PhiPQhist() {
    TFile *input = TFile::Open("~/work/clas-data/data_Fe_skimmed.root");
    TTree *tree = (TTree*) input->Get("ntuple_data");

    TFile *output = TFile::Open("histsPhiPQ.root","RECREATE");

    //// Defining variables of interest (default value 99)
    float Q2 = 99; // range [0,5] (GeV)^2
    float Nu = 99; // range [0,5] (GeV)
    // float W = 99; // range [0,3.5] (GeV)
    float pid = 99;
    float PhiPQ = 99;
    float TargType = 99; //TargType = {1: Deuterium (liquid); 2: Carbon-Iron-Lead (solid); 0: ?}

    //// Setting branches reading
    tree->SetBranchAddress("Q2",&Q2);
    tree->SetBranchAddress("Nu",&Nu);
    // tree->SetBranchAddress("W",&W);
    tree->SetBranchAddress("pid",&pid);
    tree->SetBranchAddress("PhiPQ",&PhiPQ);
    tree->SetBranchAddress("TargType",&TargType);

    //// Creating histograms
    std::vector<TH1F*> histPhiPQ_D;
    std::vector<TH1F*> histPhiPQ_X; // X stands for solid targets (C, Fe, Pb)

    char target[] = "Fe"; // The material of the target is needed!
    string prename = Form("histPhiPQ_%s_",target);

    float limits[4] = {0., 1.35, 1.82, 5.}; // Q2 limits on each bin. Calculated using quantiles!

    int Nbins = 3; //
    for (int i=0; i<Nbins; i++){ // histPhiPQ_D_XY: {X: Q2 bin; Y: Nu bin}
        auto histnameD0 = new TH1F(Form("histPhiPQ_D_%d0",i),Form("D target, %.2f < Q^{2} < %.2f, 0 < #nu < 2.5;#phi_{PQ} [degree];Counts",limits[i],limits[i+1]),360,-180,180);
        histPhiPQ_D.push_back(histnameD0);
        auto histnameD1 = new TH1F(Form("histPhiPQ_D_%d1",i),Form("D target, %.2f < Q^{2} < %.2f, 2.5 < #nu < 5.0;#phi_{PQ} [degree];Counts",limits[i],limits[i+1]),360,-180,180);
        histPhiPQ_D.push_back(histnameD1);

        string num = to_string(i);
        string name0 = prename + num + "0";
        auto histnameX0 = new TH1F(name0.c_str(), Form("%s target, %.2f < Q^{2} < %.2f, 0 < #nu < 2.5;#phi_{PQ} [degree];Counts",target,limits[i],limits[i+1]),360,-180,180);
        histPhiPQ_X.push_back(histnameX0);
        string name1 = prename + num + "1";
        auto histnameX1 = new TH1F(name1.c_str(), Form("%s target, %.2f < Q^{2} < %.2f, 2.5 < #nu < 5.0;#phi_{PQ} [degree];Counts",target,limits[i],limits[i+1]),360,-180,180);
        histPhiPQ_X.push_back(histnameX1);
    }

    //// Reading branches
    int Nentries = tree->GetEntries();

    for (int i = 0; i<Nentries; i++) {
        // if (i == 25) break; // Disable this line for real run!!!
        tree->GetEntry(i);
        
        if (pid != 211 || TargType == 0) continue;

        // Liquid target
        if (TargType == 1){
            for (int j=0; j<Nbins; j++){ // j: loop on Q2 bins -> {0., 1.35, 1.82, 5.}
                if ( Q2<limits[j+1] && Nu<2.5 ){ // 2.5 stands for a cut on Nu (Should be enhanced to quantiles)
                    histPhiPQ_D[2*j]->Fill(PhiPQ);
                    break;
                }

                if ( Q2<limits[j+1] && Nu<5.){
                    histPhiPQ_D[2*j+1]->Fill(PhiPQ);
                    break;
                }
            }
        }

        // Solid target
        if (TargType == 2){
            for (int j=0; j<Nbins; j++){ // j: loop on Q2 bins
                if ( Q2<limits[j+1] && Nu<2.5 ){ // 2.5 stands for a cut on Nu (Should be enhanced to quantiles)
                    histPhiPQ_X[2*j]->Fill(PhiPQ);
                    break;
                }

                if ( Q2<limits[j+1] && Nu<5.){
                    histPhiPQ_X[2*j+1]->Fill(PhiPQ);
                    break;
                }
            }
        }
    }

    output->Write();
    output->Close();
}
