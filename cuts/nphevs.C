
void nphevs(){
    std::string kind = "data"; // sim or data
    std::string target = "Fe"; // C, Fe, Pb

    std::string file_type;

    if (kind == "sim"){
        file_type = "hsim_"+target+"1_V2.root";
    }
    else if (kind == "data"){
        file_type = "data_"+target+"1V2_light.root";
    }

    TFile *input = TFile::Open(Form("../../%s",file_type.c_str()),"READ");
    TTree *tree = (TTree*)input->Get(Form("ntuple_%s",kind.c_str()));

    TFile *output = TFile::Open(Form("Nphevs_%s%s.root",target.c_str(),kind.c_str()),"RECREATE");

    std::string var[] = {"Q2", "Xb", "Zh", "Pt2", "PhiPQ"};
    std::string var_units[] = {"Q^{2} [GeV^{2}]", "X_{b}", "Z_{h}", "P_{t}^{2} [GeV^{2}]", "#phi_{PQ} [deg]"};

    float Q2_limits[] = {1.0, 4.0};
    float Xb_limits[] = {0.12, 0.57};
    float Zh_limits[] = {0.0, 1.0};
    float Pt2_limits[] = {0.0, 1.0};
    float PhiPQ_limits[] = {-180.0, 180.0};

    int var_bins[] = {30, 45, 40, 40, 120}; // {Q2_bins, Xb_bins, Zh_bins, Pt2_bins, PhiPQ_bins};
    float *var_limits[] = {Q2_limits, Xb_limits, Zh_limits, Pt2_limits, PhiPQ_limits};

    std::vector<TH2F*> hvec;
    int Ndim = sizeof(var)/sizeof(var[0]);
    for (int i=0; i<Ndim; i++){
        TH2F *hist = new TH2F(Form("h%s",var[i].c_str()), Form("Nphe vs %s, in %s target; Nphe; %s",var[i].c_str(),target.c_str(),var_units[i].c_str()),100,0,200,var_bins[i],var_limits[i][0],var_limits[i][1]);
        hvec.push_back(hist);
    }

 /*   
    // Definition of variables' limits (different bin sizes)
    float Q2_limits[] = {1.0, 1.17, 1.33, 1.51, 1.75, 2.12, 4.0};
    float Xb_limits[] = {0.12, 0.19, 0.23, 0.27, 0.33, 0.57};
    float Pt2_limits[] = {0.0, 0.03, 0.06, 0.1, 0.18, 1.0};

    int Q2_bins     = sizeof(Q2_limits)/sizeof(Q2_limits[0])-1; // 6
    int Xb_bins     = sizeof(Xb_limits)/sizeof(Xb_limits[0])-1; // 5
    int Pt2_bins    = sizeof(Pt2_limits)/sizeof(Pt2_limits[0])-1; // 5

    // Definition of variables' limits (equal bin sizes)
    int Zh_bins     = 10;
    int PhiPQ_bins  = 30;
    int Nphe_bins   = 100;

    float Zh_limits[Zh_bins+1];
    float PhiPQ_limits[PhiPQ_bins+1];
    float Nphe_limits[Nphe_bins+1];

    //Filling equal-sized-bins
    for (int i=0; i<=Zh_bins; i++){
        Zh_limits[i] = i/Zh_bins;
    }
    for (int i=0; i<=PhiPQ_bins; i++){
        PhiPQ_limits[i] = -180.+360.*i/PhiPQ_bins;
    }
    for (int i=0; i<=Nphe_bins; i++){
        Nphe_limits[i] = 200*i/Nphe_bins;
    }

    int var_bins[] = {Q2_bins, Xb_bins, Zh_bins, Pt2_bins, PhiPQ_bins};
    float *var_limits[] = {Q2_limits, Xb_limits, Zh_limits, Pt2_limits, PhiPQ_limits};

    std::vector<TH2F*> hvec;
    int Ndim = sizeof(var)/sizeof(var[0]);
    for (int i=0; i<Ndim; i++){
        TH2F *hist = new TH2F(Form("h%s",var[i].c_str()), Form("Nphe vs %s; Nphe; %s",var[i].c_str(),var_units[i].c_str()),Nphe_bins,Nphe_limits,var_bins[i],var_limits[i]);
        hvec.push_back(hist);
    }
*/

    int TargType = 0;
    std::vector<int> *pid = 0;

    float Q2 = 0;
    float Xb = 0;
    float Yb = 0;
    float W = 0;
    float vyec = 0;
    std::vector<float> *Zh = 0;
    std::vector<float> *Pt2 = 0;
    std::vector<float> *PhiPQ = 0;
    std::vector<float> *Nphe = 0;

    tree->SetBranchAddress("TargType",&TargType);
    tree->SetBranchAddress("pid",&pid);
    tree->SetBranchAddress("Q2",&Q2);
    tree->SetBranchAddress("Xb",&Xb);
    tree->SetBranchAddress("Yb",&Yb);
    tree->SetBranchAddress("W",&W);
    tree->SetBranchAddress("vyec",&vyec);
    tree->SetBranchAddress("Zh",&Zh);
    tree->SetBranchAddress("Pt2",&Pt2);
    tree->SetBranchAddress("PhiPQ",&PhiPQ);
    tree->SetBranchAddress("Nphe",&Nphe);

    int Nentries = tree->GetEntries();
    for (int row=0; row<Nentries; row++){
        // if (row == 500) break;
        tree->GetEntry(row);
        if (TargType!=2) continue;

        bool ecut;
        if (Q2>Q2_limits[0] && Q2<Q2_limits[1] && Xb>Xb_limits[0] && Xb<Xb_limits[1] &&
            W>2. && Yb<0.85 && vyec>-1.4 && vyec<1.4) ecut = true;
        else ecut = false;

        // bool ecut = true;        
        if (ecut){
            int ientries = pid->size();
            for (int i=0; i<ientries; i++){
                // hadron variable cuts
                if ((*pid)[i]==211 && (*Nphe)[i]>0 && (*Nphe)[i]<200 &&
                   (*Zh)[i]>=Zh_limits[0] && (*Zh)[i]<=Zh_limits[1] &&
                   (*Pt2)[i]>=Pt2_limits[0] && (*Pt2)[i]<=Pt2_limits[1] &&
                   (*PhiPQ)[i]>=PhiPQ_limits[0] && (*PhiPQ)[i]<=PhiPQ_limits[1]){
                    // hvec: {hQ2, hXb, hZh, hPt2, hPhiPQ};
                    hvec[0]->Fill((*Nphe)[i],Q2);
                    hvec[1]->Fill((*Nphe)[i],Xb);
                    hvec[2]->Fill((*Nphe)[i],(*Zh)[i]);
                    hvec[3]->Fill((*Nphe)[i],(*Pt2)[i]);
                    hvec[4]->Fill((*Nphe)[i],(*PhiPQ)[i]);
                }
            } // end of hadron's loop
        }
    } // end of electron's loop

    output->Write();
    output->Close();
}
