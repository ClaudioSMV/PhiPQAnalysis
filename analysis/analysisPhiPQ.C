
void analysisPhiPQ(){
    // Target names
    std::string target = "Fe"; // Change to analyse <target> = {D, C, Fe, Pb}
    std::string n = "1";

    std::string target_n = target+n;
    int real_targ = 2; // Liquid target (D): 1; Solid target [DEFAULT] (C, Fe, Pb): 2; Not-recognized (?): 0;

    if (target == "D"){
        // target_n = "D"; // NEEDS COMPLETE!!
        real_targ = 1;
    }

    // Data
    TFile *data = TFile::Open(Form("/home/claudio/work/clas-data/data_%s_light.root",target_n.c_str()),"READ");
    TTree *tree = (TTree*) data->Get("ntuple_data");

    // I/O
    TFile *accinput = TFile::Open(Form("Acceptance_%s_l25.root",target_n.c_str()),"READ");
    TFile *output = TFile::Open(Form("Analysis_%s_l25.root",target_n.c_str()),"RECREATE");

    // Create acceptance histograms
    TH1F *hacc = (TH1F*) accinput->Get("hacc");
    TH1F *hweight = new TH1F("hweight","Weights data",100,5,5);
    TH1F *hweight_error = new TH1F("hweight_error","Error weights data",100,1,1);

    // Create raw and corrected histograms
    // TH1F *hraw = new TH1F("hraw",Form("Non-corrected data %s, 1.35 < Q^{2} #leq 1.82, #nu #geq 2",target.c_str()),360,-180,180);
    // TH1F *hcorr = new TH1F("hcorr",Form("Corrected data %s, 1.35 < Q^{2} #leq 1.82, #nu #geq 2",target.c_str()),360,-180,180);

    // Definitions
    // Electron's variables (float or int)
    float Q2 = -99;
    float Nu = -99;
    int TargType = -99;

    // Particles' variables (vector<float or int>)
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

    int Nbins = 1;
    for (int i=1; i<=Nbins; i++){ // Format: hraw_X: {X: Q2 bin} Nu >= 2.
        // Create target histograms
        std::string nameraw = "hraw" + to_string(i);

        auto hraw = new TH1F(nameraw.c_str(), Form("Raw data on %s target, Q^{2} > %.2f, #nu #geq 2., Nphe < 25;#phi_{PQ} [deg];Counts",target_n.c_str(),limits[i]),360,-180,180); // #nu #geq 2
        hraw_v.push_back(hraw);

        std::string namecorr = "hcorr" + to_string(i);

        auto hcorr = new TH1F(namecorr.c_str(), Form("Corrected data on %s target, Q^{2} > %.2f, #nu #geq 2., Nphe < 25;#phi_{PQ} [deg];Counts",target_n.c_str(),limits[i]),360,-180,180); // #nu #geq 2
        hcorr_v.push_back(hcorr);
    }

    // Loop
    int Nentries = tree->GetEntries();
    for (int row=0; row<Nentries; row++){
        // if (row==80) break;
        tree->GetEntry(row);
        if (TargType != real_targ) continue;

        int ientries = PhiPQ->size();
        for (int i=0; i<ientries; i++){
            if ((*pid)[i] != 211 || (*Nphe)[i] < 0) continue;

            if (Nu >= 2 && limits[1] < Q2 && (*Nphe)[i]<25){
                int bin = hacc->FindBin((*PhiPQ)[i]);
                float accvalue = hacc->GetBinContent(bin);
                if (accvalue == 0) {continue;} // Eliminate problematic bins
                float accvalue_error = hacc->GetBinError(bin);
                float weight = 1./accvalue;
                float weight_error = weight*accvalue_error/accvalue; // Error propagation formula
                hweight->Fill(weight);
                hweight_error->Fill(weight_error); // Should I add these errors to hweight or just leave them as separated hists?

                hraw_v[0]->Fill((*PhiPQ)[i]);
                hcorr_v[0]->Fill((*PhiPQ)[i],weight); // Should I calculate the error as err(raw) -> weight*err or sth?
            }
        }
    }

    output->Write();
    output->Close();
}
