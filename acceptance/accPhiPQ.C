
void accPhiPQ(){
    std::string target = "Fe"; // Change to calculate acceptance of <target> = {D, C, Fe, Pb}
    int target_n=2;

    if (target == "D") target_n=1;

    TFile *input = TFile::Open(Form("/home/claudio/work/clas-data/hsim_%s1.root",target.c_str()),"READ");
    TTree *tree = (TTree*) input->Get("ntuple_sim");

    TFile *output = TFile::Open(Form("Acc_%s1_cuts2.root",target.c_str()),"RECREATE");

    float Q2_limits[] = {1.0, 4.0};
    float Xb_limits[] = {0.12, 0.57};
    float Zh_limits[] = {0., 1.};
    float Pt2_limits[] = {0.0, 1.0};
    float PhiPQ_limits[] = {-180.0, 180.0};

    int TargType;
    int mc_TargType;
    float Q2, Xb;
    float mc_Q2, mc_Xb;
    std::vector<float> *Zh = 0;
    std::vector<float> *mc_Zh = 0;
    std::vector<float> *Pt2 = 0;
    std::vector<float> *mc_Pt2 = 0;
    std::vector<float> *PhiPQ = 0;
    std::vector<float> *mc_PhiPQ = 0;
    std::vector<int> *pid = 0;
    std::vector<int> *mc_pid = 0;
    std::vector<float> *Nphe = 0;

    tree->SetBranchAddress("TargType",&TargType);
    tree->SetBranchAddress("mc_TargType",&mc_TargType);
    tree->SetBranchAddress("Q2",&Q2);
    tree->SetBranchAddress("mc_Q2",&mc_Q2);
    tree->SetBranchAddress("Xb",&Xb);
    tree->SetBranchAddress("mc_Xb",&mc_Xb);
    tree->SetBranchAddress("Zh",&Zh);
    tree->SetBranchAddress("mc_Zh",&mc_Zh);
    tree->SetBranchAddress("Pt2",&Pt2);
    tree->SetBranchAddress("mc_Pt2",&mc_Pt2);
    tree->SetBranchAddress("PhiPQ",&PhiPQ);
    tree->SetBranchAddress("mc_PhiPQ",&mc_PhiPQ);
    tree->SetBranchAddress("pid",&pid);
    tree->SetBranchAddress("mc_pid",&mc_pid);
    tree->SetBranchAddress("Nphe",&Nphe);

    // Create histograms
    TH1F *hreco = new TH1F("hreco",Form("Reconstructed, %s target;#phi_{PQ} [deg];Counts",target.c_str()),360,-180,180);
    TH1F *htrue = new TH1F("htrue",Form("Thrown (Monte Carlo), %s target;#phi_{PQ} [deg];Counts",target.c_str()),360,-180,180);
    TH1F *hacc = new TH1F("hacc",Form("Acceptance, %s target;#phi_{PQ} [deg];Counts",target.c_str()),360,-180,180);

    int Nentries = tree->GetEntries();

    for (int row=0; row<Nentries; row++){
        // if (row==5) break;
        tree->GetEntry(row);

        bool cut_el, cut_had; // Cuts applied over electron/hadron variables
        bool cut_mc_el, cut_mc_had; // Cuts applied over electron/hadron mc_variables
        if (TargType == target_n && Q2>Q2_limits[0] && Q2<Q2_limits[1] &&
          Xb>Xb_limits[0] && Xb<Xb_limits[1]) cut_el=true;
        else cut_el=false;
        if (mc_TargType == target_n && mc_Q2>Q2_limits[0] && mc_Q2<Q2_limits[1] &&
          mc_Xb>Xb_limits[0] && mc_Xb<Xb_limits[1]) cut_mc_el=true;
        else cut_mc_el=false;

        int ientries = PhiPQ->size();
        for (int i=0; i<ientries; i++){
            if ((*pid)[i]==211 && (*Nphe)[i]<25 && (*Zh)[i]>=Zh_limits[0] && (*Zh)[i]<=Zh_limits[1] && (*Pt2)[i]>=Pt2_limits[0] &&
              (*Pt2)[i]<=Pt2_limits[1] && (*PhiPQ)[i]>=PhiPQ_limits[0] && (*PhiPQ)[i]<=PhiPQ_limits[1]) cut_had=true;
            else cut_had=false;

            if ((*mc_pid)[i]==211 && (*mc_Zh)[i]>=Zh_limits[0] && (*mc_Zh)[i]<=Zh_limits[1] && (*mc_Pt2)[i]>=Pt2_limits[0] &&
              (*mc_Pt2)[i]<=Pt2_limits[1] && (*mc_PhiPQ)[i]>=PhiPQ_limits[0] && (*mc_PhiPQ)[i]<=PhiPQ_limits[1]) cut_mc_had=true;
            else cut_mc_had=false;

            if (cut_el && cut_had) hreco->Fill((*PhiPQ)[i]);
            if (cut_mc_el && cut_mc_had) htrue->Fill((*mc_PhiPQ)[i]);
        }
    }

    hacc->Sumw2(); // IMPORTANT to attach correct error values
    hacc->Divide(hreco,htrue);

    output->Write();
    output->Close();
}
