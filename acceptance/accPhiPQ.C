
void accPhiPQ(){
    TFile *input = TFile::Open("/home/claudio/work/clas-data/hsim_Fe1_skimmed.root");
    TTree *tree = (TTree*) input->Get("ntuple_sim");

    TFile *output = TFile::Open("Acceptance.root","RECREATE");

    float PhiPQ = -999;
    float mc_PhiPQ = -999;

    tree->SetBranchAddress("PhiPQ",&PhiPQ);
    tree->SetBranchAddress("mc_PhiPQ",&mc_PhiPQ);

    // Create histograms
    TH1F *hreco = new TH1F("hreco","Reconstructed",360,-180,180);
    TH1F *htrue = new TH1F("htrue","Monte Carlo",360,-180,180);

    // TH1F *hrecopion = new TH1F("hrecopion","Reconstructed",100,-180,180);//Should I calculate acceptance considering only the pion response?
    // TH1F *htruepion = new TH1F("htruepion","Monte Carlo",100,-180,180);

    TH1F *hacc = new TH1F("hacc","Acceptance",360,-180,180);
    // TH1F *haccpion = new TH1F("haccpion","Acceptance pion",100,-180,180);


    int Nentries = tree->GetEntries();

    for (int i=0; i<Nentries; i++){
        // if (i==50) break;
        tree->GetEntry(i);

        if (PhiPQ != -9999) hreco->Fill(PhiPQ);
        if (mc_PhiPQ != -9999) htrue->Fill(mc_PhiPQ);
        
        // if (pid == 211){
        //     hrecopion->Fill(PhiPQ);
        //     htruepion->Fill(mc_PhiPQ);
        // }
    }

    hacc->Sumw2(); // IMPORTANT to attach correct error values
    hacc->Divide(hreco,htrue);

    output->Write();
    output->Close();
}
