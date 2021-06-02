
void accPhiPQ(){
    std::string target = "Fe"; // Change to calculate acceptance of <target> = {D, C, Fe, Pb}

    TFile *input = TFile::Open(Form("/home/claudio/work/clas-data/hsim_%s1.root",target.c_str()));
    TTree *tree = (TTree*) input->Get("ntuple_sim");

    TFile *output = TFile::Open(Form("Acceptance_%s1.root",target.c_str()),"RECREATE");

    std::vector<float> *PhiPQ = 0;
    std::vector<float> *mc_PhiPQ = 0;
    std::vector<int> *pid = 0;
    std::vector<int> *mc_pid = 0;
    int TargType;
    int mc_TargType;

    tree->SetBranchAddress("PhiPQ",&PhiPQ);
    tree->SetBranchAddress("mc_PhiPQ",&mc_PhiPQ);
    tree->SetBranchAddress("pid",&pid);
    tree->SetBranchAddress("mc_pid",&mc_pid);
    tree->SetBranchAddress("TargType",&TargType);
    tree->SetBranchAddress("mc_TargType",&mc_TargType);

    // Create histograms
    TH1F *hreco = new TH1F("hreco",Form("Reconstructed, %s target;#phi_{PQ} [deg];Counts",target.c_str()),360,-180,180);
    TH1F *htrue = new TH1F("htrue",Form("Thrown (Monte Carlo), %s target;#phi_{PQ} [deg];Counts",target.c_str()),360,-180,180);
    TH1F *hacc = new TH1F("hacc",Form("Acceptance, %s target;#phi_{PQ} [deg];Counts",target.c_str()),360,-180,180);

    int Nentries = tree->GetEntries();

    for (int row=0; row<Nentries; row++){
        // if (row==5) break;
        tree->GetEntry(row);
        int ientries = PhiPQ->size();
        for (int i=0; i<ientries; i++){
            if ((*PhiPQ)[i] != -9999 && TargType != -9999 && (*pid)[i]==211) hreco->Fill((*PhiPQ)[i]); // It's recommended to avoid restrictions using reconstructed TargType
            if ((*mc_PhiPQ)[i] != -9999 && mc_TargType == 2 && (*mc_pid)[i]==211) htrue->Fill((*mc_PhiPQ)[i]); // Should I restrict Targ!=-9999 or this?
        }
        // Should I add another condition on pid or sth?
    }

    hacc->Sumw2(); // IMPORTANT to attach correct error values
    hacc->Divide(hreco,htrue);

    output->Write();
    output->Close();
}
