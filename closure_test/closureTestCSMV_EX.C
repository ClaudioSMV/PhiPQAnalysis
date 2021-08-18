/* Acceptance calculated by imposing general cuts to all kinematic variables at the same time (Claudio's acceptance) */
#include "functions.h"

void closureTestCSMV_EX(TString target = "Fe", TString nfolder = "1"){ //, TString xvar = "Zh"){
    // <target> = {D, C, Fe, Pb}
    // if (xvar!="Zh" && xvar!="Pt2" && xvar!="PhiPQ"){
	// 	std::cout << "ERROR: Variable name" << std::endl;
	// 	return 0;
	// }

    // I/O Files
    TString in_file = "../../clas-HSim/hsim_"+target+nfolder+".root";
    TFile *input = TFile::Open(in_file,"READ");
    TTree *tree = (TTree*) input->Get("ntuple_sim");

    TString out_file = "ClosTestCSMV_"+target+nfolder+"_EX_all.root";
    TFile *output = TFile::Open(out_file,"RECREATE");

    // Set TTree variables
    int target_n = 2;
    if (target == "D") target_n = 1;

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

    // Limits of acceptance ranges
    float Q2_limits[] = {1.0, 4.1};
    float Nu_limits[] = {2.2, 4.2};
    // float Xb_limits[] = {0.12, 0.57};
    float Zh_limits[] = {0.0, 1.0};
    float Pt2_limits[] = {0.0, 1.0};
    int NBinsPhiPQ = 180;
    float PhiPQ_limits[NBinsPhiPQ+1];
    for (int i=0; i<NBinsPhiPQ+1; i++){
        PhiPQ_limits[i] = -180. + i*360./NBinsPhiPQ;
    }

    // Define histogram's variables
    std::vector<TString> xvar = {"Zh", "Pt2", "PhiPQ"};
    std::vector<TString> xvarname = {"Z_{h}", "P_{t}^{2} [GeV^{2}]", "#phi_{PQ} [deg]"};
    int Nxvar = xvar.size();

    // Limits of final histograms
    float Q2_Hlimits[] = {1.0, 1.3, 1.8, 4.1};
    float Nu_Hlimits[] = {2.2, 3.2, 3.7, 4.2};
    float xvar_Hlimits[3][2] = {{0.0, 1.0}, {0.0, 1.0}, {-180.0, 180.0}};

    int NHQ2 = sizeof(Q2_Hlimits)/sizeof(Q2_Hlimits[0]) - 1;
    int NHNu = sizeof(Nu_Hlimits)/sizeof(Nu_Hlimits[0]) - 1;
    int NHbins = NHQ2*NHNu;

    // Create acceptance histograms
    TH1F *htrue_in = new TH1F("htrue_in","True",NBinsPhiPQ,PhiPQ_limits);
    // htrue_in->Sumw2();
    TH1F *hreco_in = new TH1F("hreco_in","Reco",NBinsPhiPQ,PhiPQ_limits);
    // hreco_in->Sumw2();
    TH1F *hacc_in = new TH1F("hacc_in","Acc",NBinsPhiPQ,PhiPQ_limits);
    hacc_in->Sumw2();

    // Create correction and ratio (final) histograms 
    std::vector<TH1F*> htrue_CT_Vec;
    std::vector<TH1F*> hcorr_CT_Vec;

    std::vector<TH1F*> hratio_Vec;

    for (int ivar=0; ivar<Nxvar; ivar++){
        for (int i=0; i<NHbins; i++){
            TString htitle = Form(" ("+target+"), %.1f<Q2<%.1f, %.1f<Nu<%.1f;"+xvarname[ivar]+";Counts",
                                Q2_Hlimits[i/NHNu],Q2_Hlimits[i/NHNu+1],Nu_Hlimits[i%NHNu],Nu_Hlimits[(i%NHNu)+1]);
            TH1F *htrue_CT = new TH1F(Form("htrue_CT"+xvar[ivar]+"%i",i),"Thrown"+htitle,10,
                                      xvar_Hlimits[ivar][0],xvar_Hlimits[ivar][1]);
            htrue_CT_Vec.push_back(htrue_CT);
            TH1F *hcorr_CT = new TH1F(Form("hcorr_CT"+xvar[ivar]+"%i",i),"Corr"+htitle,10,
                                      xvar_Hlimits[ivar][0],xvar_Hlimits[ivar][1]);
            hcorr_CT->Sumw2();
            hcorr_CT_Vec.push_back(hcorr_CT);

            TString ratiotitle = Form("%.1f<Q2<%.1f, %.1f<Nu<%.1f",Q2_Hlimits[i/NHNu],Q2_Hlimits[i/NHNu+1],
                                    Nu_Hlimits[i%NHNu],Nu_Hlimits[(i%NHNu)+1]);
            TH1F *hratio = new TH1F(Form("hratio"+xvar[ivar]+"%i",i),"CT Ratio"+ratiotitle+";"+xvarname[ivar]+";Corr/True",10,
                                    xvar_Hlimits[ivar][0],xvar_Hlimits[ivar][1]);
            hratio->Sumw2();
            hratio_Vec.push_back(hratio);
        }
    }

    // TH1F *htrue_CT_all = new TH1F("htrue_CT_all","Thrown all",10,xvar_Hlimits[0],xvar_Hlimits[1]);
    // TH1F *hcorr_CT_all = new TH1F("hcorr_CT_all","Corr all",10,xvar_Hlimits[0],xvar_Hlimits[1]);
    // TH1F *hratio_all = new TH1F("hratio_all","Ratio all",10,xvar_Hlimits[0],xvar_Hlimits[1]);
    // hratio_all->Sumw2();

    // int Nentries = 200; // for testing
    int Nentries = tree->GetEntries();
    int Nmiddle = Nentries/2;

    std::cout << "Nmiddle: " << Nmiddle << "; Nentries: " << Nentries << std::endl;

    for (int row=0; row<Nentries; row++){
        tree->GetEntry(row);

        // Calculate acceptance when first half of the data has been processed
        if (row==Nmiddle){
            hacc_in->Divide(hreco_in,htrue_in,1,1,"B");

            // delete hreco_in;
            // delete htrue_in;
        }

        bool ecut, hcut; // Cuts applied over electron/hadron variables
        bool mc_ecut, mc_hcut; // Cuts applied over electron/hadron mc_variables

        if (TargType == target_n && Q2>Q2_limits[0] && Q2<Q2_limits[1] &&
            Yb<0.85 && W>2 && vyec>-1.4 && vyec<1.4 && Nu>Nu_limits[0] &&
            Nu<Nu_limits[1]) ecut=true;
        else ecut=false;
        if (mc_TargType == target_n && mc_Q2>Q2_limits[0] && mc_Q2<Q2_limits[1] &&
            mc_Yb<0.85 && mc_W>2 && mc_Nu>Nu_limits[0] &&
            mc_Nu<Nu_limits[1]) mc_ecut=true;
        else mc_ecut=false;

        if (ecut==false && mc_ecut==false) continue; // Avoid entering a loop that won't add anything

        int binQ2, binNu;
        int binmc_Q2, binmc_Nu;
        if (row>=Nmiddle){
            if (ecut){
                binQ2 = var_position(NHQ2, Q2, Q2_Hlimits);
                binNu = var_position(NHNu, Nu, Nu_Hlimits);
            }
            if (mc_ecut){
                binmc_Q2 = var_position(NHQ2, mc_Q2, Q2_Hlimits);
                binmc_Nu = var_position(NHNu, mc_Nu, Nu_Hlimits);
            }
        }

        int ientries = PhiPQ->size();
        for (int i=0; i<ientries; i++){
            if ((*pid)[i]==211  && (*Zh)[i]>=Zh_limits[0] && (*Zh)[i]<=Zh_limits[1] && (*Pt2)[i]>=Pt2_limits[0] &&
               (*Pt2)[i]<=Pt2_limits[1] && (*PhiPQ)[i]>=PhiPQ_limits[0] && (*PhiPQ)[i]<=PhiPQ_limits[NBinsPhiPQ]) hcut=true;
            else hcut=false;

            if ((*mc_pid)[i]==211 && (*mc_Zh)[i]>=Zh_limits[0] && (*mc_Zh)[i]<=Zh_limits[1] && (*mc_Pt2)[i]>=Pt2_limits[0] &&
               (*mc_Pt2)[i]<=Pt2_limits[1] && (*mc_PhiPQ)[i]>=PhiPQ_limits[0] && (*mc_PhiPQ)[i]<=PhiPQ_limits[NBinsPhiPQ]) mc_hcut=true;
            else mc_hcut=false;

            if (ecut && hcut){
                if (row<Nmiddle) hreco_in->Fill((*PhiPQ)[i]);
                else{
                    correctCSMV((*Zh)[i], (*PhiPQ)[i], hacc_in, hcorr_CT_Vec[binNu + binQ2*NHNu + 0*NHbins]);
                    correctCSMV((*Pt2)[i], (*PhiPQ)[i], hacc_in, hcorr_CT_Vec[binNu + binQ2*NHNu + 1*NHbins]);
                    correctCSMV((*PhiPQ)[i], (*PhiPQ)[i], hacc_in, hcorr_CT_Vec[binNu + binQ2*NHNu + 2*NHbins]);
                    
                    // if (xvar=="Zh"){correctCSMV((*Zh)[i], (*PhiPQ)[i], hacc_in, hcorr_CT_Vec[binNu + binQ2*NHNu]);
                    // correctCSMV((*Zh)[i], (*PhiPQ)[i], hacc_in, hcorr_CT_all);}
                    // else if (xvar=="Pt2"){correctCSMV((*Pt2)[i], (*PhiPQ)[i], hacc_in, hcorr_CT_Vec[binNu + binQ2*NHNu]);
                    // correctCSMV((*Pt2)[i], (*PhiPQ)[i], hacc_in, hcorr_CT_all);}
                    // else if (xvar=="PhiPQ"){correctCSMV((*PhiPQ)[i], (*PhiPQ)[i], hacc_in, hcorr_CT_Vec[binNu + binQ2*NHNu]);
                    // correctCSMV((*PhiPQ)[i], (*PhiPQ)[i], hacc_in, hcorr_CT_all);}
                }
            }
            if (mc_ecut && mc_hcut){
                if (row<Nmiddle) htrue_in->Fill((*mc_PhiPQ)[i]);
                else{
                    htrue_CT_Vec[binmc_Nu + binmc_Q2*NHNu + 0*NHbins]->Fill((*mc_Zh)[i]);
                    htrue_CT_Vec[binmc_Nu + binmc_Q2*NHNu + 1*NHbins]->Fill((*mc_Pt2)[i]);
                    htrue_CT_Vec[binmc_Nu + binmc_Q2*NHNu + 2*NHbins]->Fill((*mc_PhiPQ)[i]);

                    // if (xvar=="Zh"){htrue_CT_Vec[binmc_Nu + binmc_Q2*NHNu]->Fill((*mc_Zh)[i]);
                    //     htrue_CT_all->Fill((*mc_Zh)[i]);}
                    // else if (xvar=="Pt2"){htrue_CT_Vec[binmc_Nu + binmc_Q2*NHNu]->Fill((*mc_Pt2)[i]);
                    //     htrue_CT_all->Fill((*mc_Pt2)[i]);}
                    // else if (xvar=="PhiPQ"){htrue_CT_Vec[binmc_Nu + binmc_Q2*NHNu]->Fill((*mc_PhiPQ)[i]);
                    //     htrue_CT_all->Fill((*mc_PhiPQ)[i]);}
                }
            }
        }
    }

    for (int ifin=0; ifin<NHbins*Nxvar; ifin++){
        hratio_Vec[ifin]->Divide(hcorr_CT_Vec[ifin],htrue_CT_Vec[ifin]);

        delete htrue_CT_Vec[ifin];
        delete hcorr_CT_Vec[ifin];
    }
    // hratio_all->Divide(hcorr_CT_all,htrue_CT_all);
    // delete hacc_in;

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
	return -9999;
}

void correctCSMV(float var_to_corr, float var_acc, TH1F *hacc, TH1F *hcorr){
    int bin = hacc->FindBin(var_acc);
    float accept = hacc->GetBinContent(bin);
    if (accept!=0){
        float erracc = hacc->GetBinError(bin);
        hcorr->Fill(var_to_corr,1./accept);
    }
    else if (accept==0){
        std::cout << var_to_corr << " " << hcorr->GetXaxis()->GetTitle() << "; PhiPQ: " << var_acc << "; " << std::endl;
    }
    return 0;
}