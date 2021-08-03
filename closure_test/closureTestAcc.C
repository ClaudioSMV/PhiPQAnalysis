
void closureTestAcc(){
    std::string target = "Fe"; // Change to calculate acceptance of <target> = {D, C, Fe, Pb}

    TFile *input = TFile::Open(Form("/home/claudio/work/clas-data/hsim_%s1.root",target.c_str()),"READ");
    TTree *tree = (TTree*) input->Get("ntuple_sim");

    TFile *output = TFile::Open(Form("Acc_%s1_CTest_l25.root",target.c_str()),"RECREATE"); // Closure Test

    std::vector<float> *PhiPQ = 0;
    std::vector<float> *mc_PhiPQ = 0;
    std::vector<int> *pid = 0;
    std::vector<int> *mc_pid = 0;
    std::vector<float> *Nphe = 0;
    int TargType;
    int mc_TargType;

    tree->SetBranchAddress("PhiPQ",&PhiPQ);
    tree->SetBranchAddress("mc_PhiPQ",&mc_PhiPQ);
    tree->SetBranchAddress("pid",&pid);
    tree->SetBranchAddress("mc_pid",&mc_pid);
    tree->SetBranchAddress("Nphe",&Nphe);
    tree->SetBranchAddress("TargType",&TargType);
    tree->SetBranchAddress("mc_TargType",&mc_TargType);

    // Create histograms
    TH1F *hreco1 = new TH1F("hreco1",Form("Reconstructed 1, %s target;#phi_{PQ} [deg];Counts",target.c_str()),360,-180,180);
    TH1F *htrue1 = new TH1F("htrue1",Form("Thrown 1 (Monte Carlo), %s target;#phi_{PQ} [deg];Counts",target.c_str()),360,-180,180);
    TH1F *hacc1 = new TH1F("hacc1",Form("Acceptance 1, %s target;#phi_{PQ} [deg];Counts",target.c_str()),360,-180,180);

    TH1F *hreco2 = new TH1F("hreco2",Form("Reconstructed 2, %s target;#phi_{PQ} [deg];Counts",target.c_str()),360,-180,180);
    TH1F *htrue2 = new TH1F("htrue2",Form("Thrown 2 (Monte Carlo), %s target;#phi_{PQ} [deg];Counts",target.c_str()),360,-180,180);
    TH1F *hcorr2 = new TH1F("hcorr2",Form("Corrected 2, %s target;#phi_{PQ} [deg];Counts",target.c_str()),360,-180,180);

    TH1F *hratio2 = new TH1F("hratio2",Form("Ratio Corr2/Thrown2, %s target;#phi_{PQ} [deg];Counts",target.c_str()),360,-180,180);
    TH1F *hratio3 = new TH1F("hratio3",Form("Ratio Corr2/Thrown2, %s target;#phi_{PQ} [deg];Counts",target.c_str()),20,-180,180);
    TH1F *hratio4 = new TH1F("hratio4",Form("Ratio Corr2/Thrown2, %s target;#phi_{PQ} [deg];Counts",target.c_str()),10,-180,180);

    int Nentries = tree->GetEntries();

    for (int row=0; row<Nentries; row++){
        // if (row==5) break;
        tree->GetEntry(row);
        int ientries = PhiPQ->size();

        if (row < Nentries/2){
            for (int i=0; i<ientries; i++){
                if ((*PhiPQ)[i] != -9999 && TargType != -9999 && (*pid)[i]==211 && (*Nphe)[i]<25) hreco1->Fill((*PhiPQ)[i]);
                if ((*mc_PhiPQ)[i] != -9999 && mc_TargType == 2 && (*mc_pid)[i]==211) htrue1->Fill((*mc_PhiPQ)[i]);
            }
        }
        else if (row >= Nentries/2){
            if (row == Nentries/2){
                hacc1->Sumw2(); // IMPORTANT to attach correct error values
                hacc1->Divide(hreco1,htrue1);
            }
            for (int i=0; i<ientries; i++){
                if ((*mc_PhiPQ)[i] != -9999 && mc_TargType == 2 && (*mc_pid)[i]==211) htrue2->Fill((*mc_PhiPQ)[i]);
                if ((*PhiPQ)[i] != -9999 && TargType != -9999 && (*pid)[i]==211 && (*Nphe)[i]<25){
                    hreco2->Fill((*PhiPQ)[i]);
                    int bin = hacc1->FindBin((*PhiPQ)[i]);
                    float accept = hacc1->GetBinContent(bin);
                    if (accept == 0) continue;
                    float weight = 1./accept;
                    hcorr2->Fill((*PhiPQ)[i],weight);
                }
            }
        }
        // Should I add another condition on pid or sth?
    }

    hratio2->Sumw2(); // IMPORTANT to attach correct error values
    hratio2->Divide(hcorr2,htrue2);

    TH1F *hcorr3 = dynamic_cast<TH1F*>(hcorr2->Rebin(18,"hcorr3")); // Merge 18 bins in one, leaving 20 markers!
    TH1F *htrue3 = dynamic_cast<TH1F*>(htrue2->Rebin(18,"htrue3"));

    hratio3->Sumw2(); // IMPORTANT to attach correct error values
    hratio3->Divide(hcorr3,htrue3);

    TH1F *hcorr4 = dynamic_cast<TH1F*>(hcorr2->Rebin(36,"hcorr4")); // Merge 36 bins in one, leaving 10 markers!
    TH1F *htrue4 = dynamic_cast<TH1F*>(htrue2->Rebin(36,"htrue4"));

    hratio4->Sumw2(); // IMPORTANT to attach correct error values
    hratio4->Divide(hcorr4,htrue4);

    output->Write();
    output->Close();
}
