/* Acceptance calculated by imposing general cuts to all kinematic variables at the same time (Claudio's acceptance) */
#include "functions.h"

void closureTestCSMV_EX(TString target = "Fe", TString xvar = "Zh"){
    // <target> = {D, C, Fe, Pb}
    int target_n=2;
    if (target == "D") target_n=1;

    TString in_file = "../../clas-HSim/hsim_"+target+"1.root";
    TFile *input = TFile::Open(in_file,"READ");
    TTree *tree = (TTree*) input->Get("ntuple_sim");

    TString out_file = "ClosTestCSMV_"+target+"1_EX_"+xvar+".root";
    TFile *output = TFile::Open(out_file,"RECREATE");

    float Q2_limits[] = {1.0, 1.3, 1.8, 4.1};
    float Nu_limits[] = {2.2, 3.2, 3.7, 4.2};
    float Xb_limits[] = {0.12, 0.57};
    float Zh_limits[] = {0.0, 1.0};
    float Pt2_limits[] = {0.0, 1.0};
    float PhiPQ_limits[] = {-180.0, 180.0};

    int NQ2 = sizeof(Q2_limits)/sizeof(Q2_limits[0]) - 1;
    int NNu = sizeof(Nu_limits)/sizeof(Nu_limits[0]) - 1;
    int Nbins = NQ2*NNu;

    // Set TTree variables
    int TargType;
    int mc_TargType;
    float Q2, Nu, Xb;
    float Yb, W, vyec;
    float mc_Q2, mc_Nu, mc_Xb;
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
    tree->SetBranchAddress("Nu",&Nu);
    tree->SetBranchAddress("mc_Nu",&mc_Nu);
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

    std::vector<TH1F*> hreco_Vec;
    std::vector<TH1F*> htrue_Vec;
    std::vector<TH1F*> hacc_Vec;

    std::vector<TH1F*> hreco_CT_Vec;
    std::vector<TH1F*> htrue_CT_Vec;
    std::vector<TH1F*> hcorr_CT_Vec;

    std::vector<TH1F*> hratio_Vec;

    TString xvarname;
    if (xvar=="Zh") xvarname = "Z_{h}";
    else if (xvar=="Pt2") xvarname = "P_{t}^{2} [GeV^{2}]";
    else xvarname = "ERROR";

    for (int i=0; i<Nbins; i++){
        TString htitle = Form(" ("+target+"), %.1f<Q2<%.1f, %.1f<Nu<%.1f;"+xvarname+ ";Counts",
                              Q2_limits[i/NNu],Q2_limits[i/NNu+1],Nu_limits[i%NNu],Nu_limits[(i%NNu)+1]);
        TH1F *hreco = new TH1F(Form("hreco%i",i),"Reco"+htitle,10,0.,1.);
        hreco_Vec.push_back(hreco);
        TH1F *htrue = new TH1F(Form("htrue%i",i),"Thrown"+htitle,10,0.,1.);
        htrue_Vec.push_back(htrue);
        TH1F *hacc = new TH1F(Form("hacc%i",i),"Acc"+htitle,10,0.,1.);
        hacc->Sumw2();
        hacc_Vec.push_back(hacc);

        TH1F *hreco_CT = new TH1F(Form("hreco_CT%i",i),"Reco_CT"+htitle,10,0.,1.);
        hreco_CT_Vec.push_back(hreco_CT);
        TH1F *htrue_CT = new TH1F(Form("htrue_CT%i",i),"Thrown_CT"+htitle,10,0.,1.);
        htrue_CT_Vec.push_back(htrue_CT);
        TH1F *hcorr_CT = new TH1F(Form("hcorr_CT%i",i),"Corr_CT"+htitle,10,0.,1.);
        hcorr_CT->Sumw2();
        hcorr_CT_Vec.push_back(hcorr_CT);

        TH1F *hratio = new TH1F(Form("hratio%i",i),"CT Ratio;"+xvarname+";Corr/True",10,0.,1.);
        hratio->Sumw2();
        hratio_Vec.push_back(hratio);
    }

    // int Nentries = 30; // for testing
    int Nentries = tree->GetEntries();
    int Nmiddle = Nentries/2;

    std::cout << "Nmiddle: " << Nmiddle << "; Nentries: " << Nentries << std::endl;

    for (int row=0; row<Nentries; row++){
        tree->GetEntry(row);

        if (row==Nmiddle){
            for (int iacc=0; iacc<Nbins; iacc++){
                hacc_Vec[iacc]->Divide(hreco_Vec[iacc],htrue_Vec[iacc],1,1,"B");

                delete hreco_Vec[iacc];
                delete htrue_Vec[iacc];
            }
        }

        bool ecut, hcut; // Cuts applied over electron/hadron variables
        bool mc_ecut, mc_hcut; // Cuts applied over electron/hadron mc_variables

        if (TargType == target_n && Q2>Q2_limits[0] && Q2<Q2_limits[NQ2] && Xb>Xb_limits[0] &&
            Xb<Xb_limits[1] && Yb<0.85 && W>2 && vyec>-1.4 && vyec<1.4 && Nu>Nu_limits[0] &&
            Nu<Nu_limits[NNu]) ecut=true;
        else ecut=false;
        if (mc_TargType == target_n && mc_Q2>Q2_limits[0] && mc_Q2<Q2_limits[NQ2] && mc_Xb>Xb_limits[0] &&
            mc_Xb<Xb_limits[1] && mc_Yb<0.85 && mc_W>2 && mc_Nu>Nu_limits[0] &&
            mc_Nu<Nu_limits[NNu]) mc_ecut=true;
        else mc_ecut=false;

        if (ecut==false && mc_ecut==false) continue; // Avoid entering a loop that won't add anything

        int binQ2, binNu;
        if (ecut){
            binQ2 = var_position(NQ2, Q2, Q2_limits);
            binNu = var_position(NNu, Nu, Nu_limits);
        }
        int binmc_Q2, binmc_Nu;
        if (mc_ecut){
            binmc_Q2 = var_position(NQ2, mc_Q2, Q2_limits);
            binmc_Nu = var_position(NNu, mc_Nu, Nu_limits);
        }

        int ientries = PhiPQ->size();
        for (int i=0; i<ientries; i++){
            if ((*pid)[i]==211  && (*Zh)[i]>=Zh_limits[0] && (*Zh)[i]<=Zh_limits[1] && (*Pt2)[i]>=Pt2_limits[0] &&
              (*Pt2)[i]<=Pt2_limits[1] && (*PhiPQ)[i]>=PhiPQ_limits[0] && (*PhiPQ)[i]<=PhiPQ_limits[1]) hcut=true;
            else hcut=false;

            if ((*mc_pid)[i]==211 && (*mc_Zh)[i]>=Zh_limits[0] && (*mc_Zh)[i]<=Zh_limits[1] && (*mc_Pt2)[i]>=Pt2_limits[0] &&
              (*mc_Pt2)[i]<=Pt2_limits[1] && (*mc_PhiPQ)[i]>=PhiPQ_limits[0] && (*mc_PhiPQ)[i]<=PhiPQ_limits[1]) mc_hcut=true;
            else mc_hcut=false;

            if (ecut && hcut){
                if (row<Nmiddle){
                    if (xvar=="Zh") hreco_Vec[binNu + binQ2*NNu]->Fill((*Zh)[i]);
                    else if (xvar=="Pt2") hreco_Vec[binNu + binQ2*NNu]->Fill((*Pt2)[i]);
                }
                else{
                    if (xvar=="Zh") hreco_CT_Vec[binNu + binQ2*NNu]->Fill((*Zh)[i]);
                    else if (xvar=="Pt2") hreco_CT_Vec[binNu + binQ2*NNu]->Fill((*Pt2)[i]);
                }
            }
            if (mc_ecut && mc_hcut){
                if (row<Nmiddle){
                    if (xvar=="Zh") htrue_Vec[binmc_Nu + binmc_Q2*NNu]->Fill((*mc_Zh)[i]);
                    else if (xvar=="Pt2") htrue_Vec[binmc_Nu + binmc_Q2*NNu]->Fill((*mc_Pt2)[i]);
                }
                else{
                    if (xvar=="Zh") htrue_CT_Vec[binmc_Nu + binmc_Q2*NNu]->Fill((*mc_Zh)[i]);
                    else if (xvar=="Pt2") htrue_CT_Vec[binmc_Nu + binmc_Q2*NNu]->Fill((*mc_Pt2)[i]);
                }
            }
        }
    }

    for (int ifin=0; ifin<Nbins; ifin++){
        hcorr_CT_Vec[ifin]->Divide(hreco_CT_Vec[ifin],hacc_Vec[ifin]);
        hratio_Vec[ifin]->Divide(hcorr_CT_Vec[ifin],htrue_CT_Vec[ifin]);

        delete hreco_CT_Vec[ifin];
        delete htrue_CT_Vec[ifin];
    }

    output->Write();
    output->Close();
}

////// DEFINITION OF FUNCTIONS

int var_position(int Nvar, float var, float var_limits[]){
	for (int ivar=0; ivar<Nvar; ivar++){
		if (var>=var_limits[ivar] && var<var_limits[ivar+1]){
			return ivar;
		}
	}
	return -1;
}

int var_position(int i, int Nvar, std::vector<float> *var, float var_limits[]){
	for (int ivar=0; ivar<Nvar; ivar++){
		if ((*var)[i]>var_limits[ivar] && (*var)[i]<var_limits[ivar+1]){
			return ivar;
		}
	}
	return -1;
}
