
void nphe(){
    TFile *input = TFile::Open("/home/claudio/work/clas-data/hsim_Fe1.root","READ");
    TTree *tree = (TTree*)input->Get("ntuple_sim");

    TFile *output = TFile::Open("/home/claudio/work/clas-data/acceptance/Nphe_sim.root","RECREATE");

    TH1F *hNphe_1 = new TH1F("hNphe_1","Nphe reconstructed, TargType!=-9999; Nphe; Counts",200,0,200);// reconstructed data TT != -9999
    TH1F *hNphe_2 = new TH1F("hNphe_2","Nphe reconstructed, TargType==2; Nphe; Counts",200,0,200);// reconstructed data TT == 2
    // TH1F *hNphe = new TH1F("hNphe","Nphe data; Nphe; Counts",200,0,200);// actual data

    int TargType = 0;
    std::vector<int> *pid = 0;
    std::vector<float> *Nphe = 0;

    tree->SetBranchAddress("TargType",&TargType);
    tree->SetBranchAddress("pid",&pid);
    tree->SetBranchAddress("Nphe",&Nphe);

    int Nentries = tree->GetEntries();
    for (int row=0; row<Nentries; row++){
        // if (row == 150) break;
        tree->GetEntry(row);

        if (TargType==-9999) continue; // sim: TT==-9999; data: TT!=2;
        int ientries = pid->size();
        for (int i=0; i<ientries; i++){
            if ((*pid)[i]==211 && (*Nphe)[i]>0 && (*Nphe)[i]<200){
                hNphe_1->Fill((*Nphe)[i]);
                if (TargType==2){
                    hNphe_2->Fill((*Nphe)[i]);
                }
            }
        }
    }

    output->Write();
    output->Close();
}
