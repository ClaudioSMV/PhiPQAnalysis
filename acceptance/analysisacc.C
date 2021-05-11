
void analysisacc(){
    TFile *accinput = TFile::Open("Acceptance.root");

    TFile *output = TFile::Open("Analysis.root","RECREATE");

    // Create acceptance histograms
    TH1F *hacc = (TH1F*) accinput->Get("hacc");
    TH1F *hweight = new TH1F("hweight","Weights data",100,5,55);
    TH1F *hweight_error = new TH1F("hweight_error","Error weights data",100,1,1);

    // Create raw and corrected histograms
    TH1F *hraw = new TH1F("hraw","Non-corrected data D, 1.35 < Q^{2}< 1.82, #nu > 2",360,-180,180);
    TH1F *hcorr = new TH1F("hcorr","Corrected data D, 1.35 < Q^{2}< 1.82, #nu > 2",360,-180,180);

    // Data
    TFile *data = TFile::Open("/home/claudio/work/clas-data/data_Fe_skimmed.root");
    TTree *tree = (TTree*) data->Get("ntuple_data");

    // Definitions
    float Q2 = -99;
    float Nu = -99;
    float pid = -99;
    float PhiPQ = -999;
    float TargType = -99;

    // Retrieve data
    tree->SetBranchAddress("Q2",&Q2);
    tree->SetBranchAddress("Nu",&Nu);
    tree->SetBranchAddress("pid",&pid);
    tree->SetBranchAddress("PhiPQ",&PhiPQ);
    tree->SetBranchAddress("TargType",&TargType);


    // //// Creating from-tree histograms
    // std::vector<TH1F*> hraw_D;
    // std::vector<TH1F*> hraw_X; // X stands for solid targets (C, Fe, Pb)

    // std::vector<TH1F*> hcorr_D;
    // std::vector<TH1F*> hcorr_X; // X stands for solid targets (C, Fe, Pb)

    // char target[] = "Fe"; // The material of the target is needed!
    // string prenamerawX = Form("hraw_%s_",target);
    // string prenamecorrX = Form("hcorr_%s_",target);

    float limits[4] = {0., 1.35, 1.82, 5.}; // Q2 limits on each bin. Calculated using quantiles!

    // int Nbins = 3;
    // for (int i=0; i<Nbins; i++){ // Format: hraw_D_XY: {X: Q2 bin; Y: Nu bin}
    //     // Deuterium (liquid) target histograms
    //     auto hrawnameD0 = new TH1F(Form("hraw_D_%d0",i),Form("Raw data on D target, %.2f < Q^{2} < %.2f, 0 < #nu < 2.5;#phi_{PQ} [degree];Counts",limits[i],limits[i+1]),360,-180,180);
    //     hraw_D.push_back(hrawnameD0);
    //     auto hrawnameD1 = new TH1F(Form("hraw_D_%d1",i),Form("Raw data on D target, %.2f < Q^{2} < %.2f, 2.5 < #nu < 5.0;#phi_{PQ} [degree];Counts",limits[i],limits[i+1]),360,-180,180);
    //     hraw_D.push_back(hrawnameD1);

    //     auto hcorrnameD0 = new TH1F(Form("hcorr_D_%d0",i),Form("Corrected data on D target, %.2f < Q^{2} < %.2f, 0 < #nu < 2.5;#phi_{PQ} [degree];Counts",limits[i],limits[i+1]),360,-180,180);
    //     hcorr_D.push_back(hcorrnameD0);
    //     auto hcorrnameD1 = new TH1F(Form("hcorr_D_%d1",i),Form("Corrected data on D target, %.2f < Q^{2} < %.2f, 2.5 < #nu < 5.0;#phi_{PQ} [degree];Counts",limits[i],limits[i+1]),360,-180,180);
    //     hcorr_D.push_back(hcorrnameD1);

    //     // X (solid) target histograms
    //     // string num = to_string(i);
    //     string nameraw0 = prenamerawX + to_string(i) + "0";
    //     string nameraw1 = prenamerawX + to_string(i) + "1";

    //     auto hrawnameX0 = new TH1F(nameraw0.c_str(), Form("Raw data on %s target, %.2f < Q^{2} < %.2f, 0 < #nu < 2.5;#phi_{PQ} [degree];Counts",target,limits[i],limits[i+1]),360,-180,180);
    //     hraw_X.push_back(hrawnameX0);
    //     auto hrawnameX1 = new TH1F(nameraw1.c_str(), Form("Raw data on %s target, %.2f < Q^{2} < %.2f, 2.5 < #nu < 5.0;#phi_{PQ} [degree];Counts",target,limits[i],limits[i+1]),360,-180,180);
    //     hraw_X.push_back(hrawnameX1);

    //     string namecorr0 = prenamecorrX + to_string(i) + "0";
    //     string namecorr1 = prenamecorrX + to_string(i) + "1";

    //     auto hcorrnameX0 = new TH1F(namecorr0.c_str(), Form("Corrected data on %s target, %.2f < Q^{2} < %.2f, 0 < #nu < 2.5;#phi_{PQ} [degree];Counts",target,limits[i],limits[i+1]),360,-180,180);
    //     hcorr_X.push_back(hcorrnameX0);
    //     auto hcorrnameX1 = new TH1F(namecorr1.c_str(), Form("Corrected data on %s target, %.2f < Q^{2} < %.2f, 2.5 < #nu < 5.0;#phi_{PQ} [degree];Counts",target,limits[i],limits[i+1]),360,-180,180);
    //     hcorr_X.push_back(hcorrnameX1);
    // }

    // Acceptance variables definition
    int bin = -99;
    float accvalue = -99;
    float accvalue_error = -99;
    float weight = -99;
    float weight_error = -99;

    // Loop
    int Nentries = tree->GetEntries();
    for (int i=0; i<Nentries; i++){
        // if (i==50) break;
        tree->GetEntry(i);
        if (pid != 211 || TargType == 0) continue;

        if (limits[1] < Q2 < limits[2] && TargType == 1 && Nu > 2){ // Deuterium
            hraw->Fill(PhiPQ);
            bin = hacc->FindBin(PhiPQ);
            accvalue = hacc->GetBinContent(bin);
            if (accvalue == 0) {continue;} // Eliminate problematic bins
            accvalue_error = hacc->GetBinError(bin);
            weight = 1./accvalue;
            weight_error = weight*accvalue_error/accvalue; // Error propagation formula

            // std::cout << "Bin content: " << accvalue << "; Error: " << accvalue_error << "; Error^2: " << accvalue_error*accvalue_error << ";" << std::endl;
            // std::cout << "Weight:" << weight << "; Weight error: " << weight_error << "; Acc: " << accvalue << "; Acc error: " << accvalue_error << std::endl;

            hweight->Fill(weight);
            hweight_error->Fill(weight_error); // Should I add these errors to hweight or just leave them as separated hists?
            hcorr->Fill(PhiPQ,weight); // Should I calculate the error as err(raw) -> weight*err or sth?
        }
        // if (limits[1] < Q2 < limits[2] && TargType == 2){ // Fe
        // }
    }

    output->Write();
    output->Close();
}
