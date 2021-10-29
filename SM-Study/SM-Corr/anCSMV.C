/* Analysis made using acceptance calculated in Claudio's way */

void anCSMV(){
////// Beginning (Edit according to target name) <target> = {D, C, Fe, Pb}
    std::string target = "Fe"; // Change to calculate acceptance of <target> = {D, C, Fe, Pb}

////// I/O Files
    TFile *indata = TFile::Open(Form("../../data_%s1V2_light.root",target.c_str()),"READ");
    TTree *tree = (TTree*) indata->Get("ntuple_data");
    TFile *inacc = TFile::Open(Form("../acceptance/root-files/AccCSMV_%s2.root",target.c_str()),"READ");

    TFile *output = TFile::Open(Form("AnalysisCSMV_%s2.root",target.c_str()),"RECREATE");

////// Limits and sizes of cuts
    float Q2_limits[] = {1.0, 4.0};
    float Xb_limits[] = {0.12, 0.57};
    float Zh_limits[] = {0.0, 1.0};
    float Pt2_limits[] = {0.0, 1.0};
    float PhiPQ_limits[] = {-180.0, 180.0};

////// Definition of tree-variables
    int target_n = 2;
    if (target == "D") target_n = 1;

    int TargType;
    float Q2, Xb;
    float Yb, W, vyec;
    std::vector<float> *Zh = 0;
    std::vector<float> *Pt2 = 0;
    std::vector<float> *PhiPQ = 0;
    std::vector<int> *pid = 0;
    std::vector<float> *Nphe = 0;

    tree->SetBranchAddress("TargType",&TargType);
    tree->SetBranchAddress("Q2",&Q2);
    tree->SetBranchAddress("Xb",&Xb);
    tree->SetBranchAddress("Yb",&Yb);
    tree->SetBranchAddress("W",&W);
    tree->SetBranchAddress("vyec",&vyec);
    tree->SetBranchAddress("Zh",&Zh);
    tree->SetBranchAddress("Pt2",&Pt2);
    tree->SetBranchAddress("PhiPQ",&PhiPQ);
    tree->SetBranchAddress("pid",&pid);
    tree->SetBranchAddress("Nphe",&Nphe);

////// Create histograms
    TH1F *hdata = new TH1F("hdata",Form("Data, %s target;#phi_{PQ} [deg];Counts",target.c_str()),180,-180,180);
    TH1F *hcorr = new TH1F("hcorr",Form("Corrected, %s target;#phi_{PQ} [deg];Counts",target.c_str()),180,-180,180);
    hcorr->Sumw2();

////// Tree loop
    int Nentries = tree->GetEntries();
    for (int row=0; row<Nentries; row++){
        // if (row==60) break;
        tree->GetEntry(row);

        bool ecut=false, hcut=false;

        if (TargType==target_n && Q2>Q2_limits[0] && Q2<Q2_limits[1] && Xb>Xb_limits[0] && Xb<Xb_limits[1] &&
            Yb<0.85 && W>2 && vyec>-1.4 && vyec<1.4){
            ecut = true;
        }

        if (ecut==false) continue;

        int ientries = PhiPQ->size();
        for (int i=0; i<ientries; i++){
            if ((*pid)[i]==211 && (*Zh)[i]>Zh_limits[0] && (*Zh)[i]<Zh_limits[1] && (*Pt2)[i]>Pt2_limits[0] &&
                (*Pt2)[i]<Pt2_limits[1] && (*PhiPQ)[i]>PhiPQ_limits[0] && (*PhiPQ)[i]<PhiPQ_limits[1]){
                hcut = true;
            }
            if (ecut && hcut) hdata->Fill((*PhiPQ)[i]);
        } // end of hadrons' loop
    } // end filling loop

////// Filling final histograms
    TH1F *hacc = (TH1F*) inacc->Get("hacc");
    hcorr->Divide(hdata,hacc);

    output->Write();
    output->Close();
} // End of the macro
