
void analysisPhiPQ(){
    // Target names
    std::string target = "Fe"; // Change to analyse <target> = {D, C, Fe, Pb}
    std::string target_n = target;
    int target_type = 2; // Liquid target (D): 1; Solid target [DEFAULT] (C, Fe, Pb): 2; Not-recognized (?): 0;

    if (target == "D"){
        // target_n = "D"; // NEEDS COMPLETE!!
        target_type = 1;
    }
    if (target == "Fe"){
        target_n = "Fe4"; // Necessary to enter in the correct data-file
    }

    // Data
    TFile *data = TFile::Open(Form("/home/claudio/work/clas-data/data_%s.root",target_n.c_str()));
    TTree *tree = (TTree*) data->Get("ntuple_data");

    // I/O
    TFile *accinput = TFile::Open(Form("Acceptance_%s1.root",target.c_str()));
    TFile *output = TFile::Open(Form("Analysis_%s1.root",target.c_str()),"RECREATE");

    // Create acceptance histograms
    TH1F *hacc = (TH1F*) accinput->Get("hacc");
    TH1F *hweight = new TH1F("hweight","Weights data",100,5,5);
    TH1F *hweight_error = new TH1F("hweight_error","Error weights data",100,1,1);

    // Create raw and corrected histograms
    // TH1F *hraw = new TH1F("hraw",Form("Non-corrected data %s, 1.35 < Q^{2} #leq 1.82, #nu #geq 2",target.c_str()),360,-180,180);
    // TH1F *hcorr = new TH1F("hcorr",Form("Corrected data %s, 1.35 < Q^{2} #leq 1.82, #nu #geq 2",target.c_str()),360,-180,180);

    // Definitions
    // Electron's variables (float)
    float Q2 = -99;
    float Nu = -99;
    int TargType = -99;

    // Particles' variables (vector<float>)
    std::vector<int> *pid = 0;
    std::vector<float> *PhiPQ = 0;
    std::vector<float> *Nphe = 0;

    // Retrieve data
    tree->SetBranchAddress("Q2",&Q2);
    tree->SetBranchAddress("Nu",&Nu);
    tree->SetBranchAddress("TargType",&TargType);
    tree->SetBranchAddress("pid",&pid);
    tree->SetBranchAddress("PhiPQ",&PhiPQ);
    tree->SetBranchAddress("Nphe",&Nphe);

    float limits[4] = {0., 1.35, 1.82, 5.}; // Q2 limits on each bin. Calculated using quantiles!

    // //// Creating from-tree histograms
    std::vector<TH1F*> hraw_v;
    std::vector<TH1F*> hcorr_v;

    int Nbins = 6; // 6 bins for Nphe -> {0-5, 5-10, 10-15, 15-20, 20-25, >25}
    for (int i=0; i<Nbins; i++){ // Format: hraw_XY: {X: Q2 bin; Y: Nphe bin} Nu >= 2.
        // Create target histograms
        std::string nameraw = "hraw1" + to_string(i);

        auto hraw = new TH1F(nameraw.c_str(), Form("Raw data on %s target, %.2f < Q^{2} < %.2f, %i #leq Nphe < %i;#phi_{PQ} [deg];Counts",target.c_str(),limits[1],limits[2],5*i,5*(i+1)),360,-180,180); // #nu #geq 2
        hraw_v.push_back(hraw);

        std::string namecorr = "hcorr1" + to_string(i);

        auto hcorr = new TH1F(namecorr.c_str(), Form("Corrected data on %s target, %.2f < Q^{2} < %.2f, %i #leq Nphe < %i;#phi_{PQ} [deg];Counts",target.c_str(),limits[1],limits[2],5*i,5*(i+1)),360,-180,180); // #nu #geq 2
        hcorr_v.push_back(hcorr);
    }

    // Loop
    int Nentries = tree->GetEntries();
    for (int row=0; row<Nentries; row++){
        // if (row==180) break;
        tree->GetEntry(row);
        if (TargType != target_type) continue;

        int ientries = PhiPQ->size();
        for (int i=0; i<ientries; i++){
            if ((*pid)[i] != 211 || (*Nphe)[i] < 0) continue;

            if (Nu >= 2 && limits[1] < Q2 && Q2 <= limits[2]){ // Fe target
                int bin = hacc->FindBin((*PhiPQ)[i]);
                float accvalue = hacc->GetBinContent(bin);
                if (accvalue == 0) {continue;} // Eliminate problematic bins
                float accvalue_error = hacc->GetBinError(bin);
                float weight = 1./accvalue;
                float weight_error = weight*accvalue_error/accvalue; // Error propagation formula
                hweight->Fill(weight);
                hweight_error->Fill(weight_error); // Should I add these errors to hweight or just leave them as separated hists?

                int iNphe = (*Nphe)[i]/5.;
                if (iNphe > 5) {iNphe = 5;}

                hraw_v[iNphe]->Fill((*PhiPQ)[i]);
                hcorr_v[iNphe]->Fill((*PhiPQ)[i],weight); // Should I calculate the error as err(raw) -> weight*err or sth?
            }
        }
    }

    output->Write();
    output->Close();
}
