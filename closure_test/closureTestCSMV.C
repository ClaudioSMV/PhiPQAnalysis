/* Acceptance calculated by imposing general cuts to all kinematic variables at the same time (Claudio's acceptance) */

void closureTestCSMV(TString target = "Fe"){
    // <target> = {D, C, Fe, Pb}
    int target_n=2;
    if (target == "D") target_n=1;

    TString in_file = "../../clas-HSim/hsim_"+target+"1.root";
    TFile *input = TFile::Open(in_file,"READ");
    TTree *tree = (TTree*) input->Get("ntuple_sim");

    TString out_file = "ClosTestCSMV_"+target+"1.root";
    TFile *output = TFile::Open(out_file,"RECREATE");

    float Q2_limits[] = {1.0, 4.0};
    float Xb_limits[] = {0.12, 0.57};
    float Zh_limits[] = {0., 1.};
    float Pt2_limits[] = {0.0, 1.0};
    float PhiPQ_limits[] = {-180.0, 180.0};

    int TargType;
    int mc_TargType;
    float Q2, Xb;
    float Yb, W, vyec;
    float mc_Q2, mc_Xb;
    float mc_Yb, mc_W;
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
    tree->SetBranchAddress("Yb",&Yb);
    tree->SetBranchAddress("mc_Yb",&mc_Yb);
    tree->SetBranchAddress("W",&W);
    tree->SetBranchAddress("mc_W",&mc_W);
    tree->SetBranchAddress("vyec",&vyec);
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
    TString htitle = " ("+target+");#phi_{PQ} [deg];Counts";
    TH1F *hreco = new TH1F("hreco","Reco"+htitle,180,-180,180);
    TH1F *htrue = new TH1F("htrue","Thrown"+htitle,180,-180,180);
    TH1F *hacc = new TH1F("hacc","Acce"+htitle,180,-180,180);

    TH1F *hreco_CT = new TH1F("hreco_CT","Reco CT"+htitle,180,-180,180);
    TH1F *htrue_CT = new TH1F("htrue_CT","Thrown CT"+htitle,180,-180,180);
    TH1F *hcorr_CT = new TH1F("hcorr_CT","Corr CT"+htitle,180,-180,180);

    TH1F *hratio = new TH1F("hratio","CT Ratio;#phi_{PQ} [deg];Corr/True",180,-180,180);
    TH1F *hratio_rebin = new TH1F("hratio_rebin","CT Ratio;#phi_{PQ} [deg];Corr/True",10,-180,180);

    hacc->Sumw2();
    hcorr_CT->Sumw2();
    hratio->Sumw2();
    hratio_rebin->Sumw2();

    // int Nentries = 30; // for testing
    int Nentries = tree->GetEntries();
    int Nmiddle = Nentries/2;

    for (int row=0; row<Nentries; row++){
        tree->GetEntry(row);

        if (row==Nmiddle){
            hacc->Divide(hreco,htrue,1,1,"B");
        }

        bool cut_el, cut_had; // Cuts applied over electron/hadron variables
        bool cut_mc_el, cut_mc_had; // Cuts applied over electron/hadron mc_variables

        if (TargType == target_n && Q2>Q2_limits[0] && Q2<Q2_limits[1] && Xb>Xb_limits[0] &&
            Xb<Xb_limits[1] && Yb<0.85 && W>2 && vyec>-1.4 && vyec<1.4) cut_el=true;
        else cut_el=false;
        if (mc_TargType == target_n && mc_Q2>Q2_limits[0] && mc_Q2<Q2_limits[1] && mc_Xb>Xb_limits[0] &&
            mc_Xb<Xb_limits[1] && mc_Yb<0.85 && mc_W>2) cut_mc_el=true;
        else cut_mc_el=false;

        if (cut_el==false && cut_mc_el==false) continue; // Avoid entering a loop that won't add anything
        int ientries = PhiPQ->size();
        for (int i=0; i<ientries; i++){
            if ((*pid)[i]==211  && (*Zh)[i]>=Zh_limits[0] && (*Zh)[i]<=Zh_limits[1] && (*Pt2)[i]>=Pt2_limits[0] &&
              (*Pt2)[i]<=Pt2_limits[1] && (*PhiPQ)[i]>=PhiPQ_limits[0] && (*PhiPQ)[i]<=PhiPQ_limits[1]) cut_had=true;
            else cut_had=false;

            if ((*mc_pid)[i]==211 && (*mc_Zh)[i]>=Zh_limits[0] && (*mc_Zh)[i]<=Zh_limits[1] && (*mc_Pt2)[i]>=Pt2_limits[0] &&
              (*mc_Pt2)[i]<=Pt2_limits[1] && (*mc_PhiPQ)[i]>=PhiPQ_limits[0] && (*mc_PhiPQ)[i]<=PhiPQ_limits[1]) cut_mc_had=true;
            else cut_mc_had=false;

            if (cut_el && cut_had){
                if (row<Nmiddle) hreco->Fill((*PhiPQ)[i]);
                else{
                    hreco_CT->Fill((*PhiPQ)[i]);
                    int bin = hacc->FindBin((*PhiPQ)[i]);
                    float accept = hacc->GetBinContent(bin);
                    float err_accept = hacc->GetBinError(bin);
                    if (accept!=0){
                        float weight = 1./accept;
                        float err_weight = weight*err_accept/accept;
                        hcorr_CT->Fill((*PhiPQ)[i],weight);
                        // hcorr_CT->Fill((*PhiPQ)[i],weight);
                    }
                }
            }
            if (cut_mc_el && cut_mc_had){
                if (row<Nmiddle) htrue->Fill((*mc_PhiPQ)[i]);
                else{
                    htrue_CT->Fill((*mc_PhiPQ)[i]);
                }
            }
        }
    }
    hcorr_CT->Divide(hreco_CT,hacc);
    hratio->Divide(hcorr_CT,htrue_CT);

    TH1F *hcorr_CT_rebin = dynamic_cast<TH1F*>(hcorr_CT->Rebin(18,"hcorr_CT_rebin")); // Merge 18 bins in one, leaving 10 markers!
    TH1F *htrue_CT_rebin = dynamic_cast<TH1F*>(htrue_CT->Rebin(18,"htrue_CT_rebin"));

    hratio_rebin->Divide(hcorr_CT_rebin,htrue_CT_rebin);

    output->Write();
    output->Close();
}
