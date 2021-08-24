/* Closure Test using a 5-d acceptance with THnSparse (S.Morán's acceptance) */

#include "functions.h"

void closureTestSM_Sparse(TString target = "Fe", TString nfold = "*"){
	if (nfold!="*" && nfold!="1" && nfold!="2" && nfold!="3"){
        std::cerr << "[ERROR] file doesn't exist" << std::endl;
        return -1;
    }

    // I/O Files
    TChain *tree = new TChain("ntuple_sim");
    TString in_file = "../../clas-HSim/hsim_"+target+nfold+".root";
    tree->Add(in_file);

    TString out_file;
    if (nfold=="*") out_file = "ClosTestSM_"+target+"FULL_SP.root";
    else out_file = "ClosTestSM_"+target+nfold+"_SP.root";
    TFile *output = TFile::Open(out_file,"RECREATE");

	// Definition of tree-variables
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

	// Limits of acceptance ranges
	// float Q2_limits[] = {1.0, 1.3, 1.8, 4.1}; // {1.0, 1.2, 1.3, 1.5, 1.8, 2.8, 4.1};
	// float Nu_limits[] = {2.2, 3.2, 3.7, 4.2}; // {2.2, 2.8, 3.2, 3.5, 3.7, 3.9, 4.2};
	// float Zh_limits[] = {0.0, 0.15, 0.25, 0.4, 0.7, 1.0};
	// float Pt2_limits[] = {0.0, 0.03, 0.06, 0.1, 0.18, 1.0};
	// int NBinsPhiPQ = 12; // Should be 180 but it's too much... 
    // float PhiPQ_limits[NBinsPhiPQ+1];
    // for (int i=0; i<NBinsPhiPQ+1; i++){
    //     PhiPQ_limits[i] = -180. + i*360./NBinsPhiPQ;
    // }

	// Create acceptance histograms
	const Int_t Ndim = 5;
	Int_t nbins[Ndim] = {3, 3, 5, 5, 12}; // One can use 120 instead of 12 for PhiPQ
	Double_t minbins[Ndim] = {1.0, 2.2, 0.0, 0.0, -180.0};
	Double_t maxbins[Ndim] = {4.1, 4.2, 1.0, 1.0, 180.0};

	THnSparse *hreco = new THnSparseD("hreco","Reco",Ndim,nbins,minbins,maxbins);
	THnSparse *htrue = new THnSparseD("htrue","True",Ndim,nbins,minbins,maxbins);
	THnSparse *hacc = new THnSparseD("hacc","Acc",Ndim,nbins,minbins,maxbins);

	// Set variable width bins
	Double_t Q2_limits[] = {1.0, 1.3, 1.8, 4.1};
	Double_t Nu_limits[] = {2.2, 3.2, 3.7, 4.2};
	Double_t Zh_limits[] = {0.0, 0.15, 0.25, 0.4, 0.7, 1.0};
	Double_t Pt2_limits[] = {0.0, 0.03, 0.06, 0.1, 0.18, 1.0};
    Double_t PhiPQ_limits[nbins[4]+1];
	for (int i=0; i<nbins[4]+1; i++){
        PhiPQ_limits[i] = -180. + i*360./nbins[4];
    }

	setBinVarSize(hreco,nbins,Q2_limits,Nu_limits,Zh_limits,Pt2_limits,PhiPQ_limits);
	setBinVarSize(htrue,nbins,Q2_limits,Nu_limits,Zh_limits,Pt2_limits,PhiPQ_limits);
	setBinVarSize(hacc,nbins,Q2_limits,Nu_limits,Zh_limits,Pt2_limits,PhiPQ_limits);
	hreco->Sumw2();
	htrue->Sumw2();
	hacc->Sumw2();

	// Define histogram's variables
	const Int_t Nxvar = 3;
    TString xvar[Nxvar] = {"Zh", "Pt2", "PhiPQ"};
    TString xvarname[Nxvar] = {"Z_{h}", "P_{t}^{2} [GeV^{2}]", "#phi_{PQ} [deg]"};

	// Limits of final histograms
	float Q2_Hlimits[4] = {1.0, 1.3, 1.8, 4.1};
	float Nu_Hlimits[4] = {2.2, 3.2, 3.7, 4.2};
	float xvar_Hlimits[3][2] = {{0.0, 1.0}, {0.0, 1.0}, {-180., 180.}};

	// N<var> is the number of bins (array's elements - 1), N<var>+1 is the number of limits (array's elements)
	int NQ2  = sizeof(Q2_limits)/sizeof(Q2_limits[0]) - 1;
	int NNu  = sizeof(Nu_limits)/sizeof(Nu_limits[0]) - 1;
	int NZh  = sizeof(Zh_limits)/sizeof(Zh_limits[0]) - 1;
	int NPt2 = sizeof(Pt2_limits)/sizeof(Pt2_limits[0]) - 1;
	int NPhiPQ  = sizeof(PhiPQ_limits)/sizeof(PhiPQ_limits[0]) - 1;

	int NHQ2 = sizeof(Q2_Hlimits)/sizeof(Q2_Hlimits[0]) - 1;
	int NHNu = sizeof(Nu_Hlimits)/sizeof(Nu_Hlimits[0]) - 1;

	int NTot = NQ2*NNu*NZh*NPt2*NPhiPQ; // Total number of acceptance-cut-bins
	int NHbins = NHQ2*NHNu; // Total number of plot bins

	// Acceptance-tree loop

	// int Nentries = 50; // for testing
	int Nentries = tree->GetEntries();
	int Nmiddle = Nentries/2;

	std::cout << "Nmiddle: " << Nmiddle << "; Nentries: " << Nentries << std::endl;

	for (int row=0; row<Nmiddle; row++){
		tree->GetEntry(row);

		bool ecut=false, mc_ecut=false; // hcut=false, mc_hcut=false;

		if (TargType==target_n && Q2>Q2_limits[0] && Q2<Q2_limits[NQ2] &&
			Yb<0.85 && W>2 && vyec>-1.4 && vyec<1.4 && Nu>Nu_limits[0] && Nu<Nu_limits[NNu]){
			ecut = true;
		}

		if (mc_TargType==target_n && mc_Q2>Q2_limits[0] && mc_Q2<Q2_limits[NQ2] &&
			mc_Yb<0.85 && mc_W>2 && mc_Nu>Nu_limits[0] && mc_Nu<Nu_limits[NNu]){
			mc_ecut = true;
		}

		if (ecut==false && mc_ecut==false) continue;

		int ientries = PhiPQ->size();
		for (int i=0; i<ientries; i++){
			if (ecut==true && (*pid)[i]==211 && (*Zh)[i]>Zh_limits[0] && (*Zh)[i]<Zh_limits[NZh] &&
			    (*Pt2)[i]>Pt2_limits[0] && (*Pt2)[i]<Pt2_limits[NPt2] && (*PhiPQ)[i]>PhiPQ_limits[0] &&
				(*PhiPQ)[i]<PhiPQ_limits[NPhiPQ]){
				// hcut = true;
				Double_t new_entry[] = {Q2, Nu, (*Zh)[i], (*Pt2)[i], (*PhiPQ)[i]};
				hreco->Fill(new_entry);
			}

			if (mc_ecut==true && (*mc_pid)[i]==211 && (*mc_Zh)[i]>Zh_limits[0] && (*mc_Zh)[i]<Zh_limits[NZh] &&
				(*mc_Pt2)[i]>Pt2_limits[0] && (*mc_Pt2)[i]<Pt2_limits[NPt2] && (*mc_PhiPQ)[i]>PhiPQ_limits[0] &&
				(*mc_PhiPQ)[i]<PhiPQ_limits[NPhiPQ]){
				// mc_hcut = true;
				Double_t new_mc_entry[] = {mc_Q2, mc_Nu, (*mc_Zh)[i], (*mc_Pt2)[i], (*mc_PhiPQ)[i]};
				htrue->Fill(new_mc_entry);
			}
		} // end of hadrons' loop
  	} // end filling loop

	// Calculating acceptance (5-dimensional)
	hacc->Divide(hreco,htrue,1,1,"B");
	delete hreco;
	delete htrue;

	// Create correction and ratio histograms
    std::vector<TH1F*> htrue_CT_Vec;
    std::vector<TH1F*> hcorr_CT_Vec;
    std::vector<TH1F*> hratio_Vec;

	for (int ivar=0; ivar<Nxvar; ivar++){
		for (int i=0; i<NHbins; i++){
			TString htitle = Form(" ("+target+"), %.1f<Q2<%.1f, %.1f<Nu<%.1f",Q2_Hlimits[i/NHNu],
								Q2_Hlimits[i/NHNu+1],Nu_Hlimits[i%NHNu],Nu_Hlimits[(i%NHNu)+1]);
			TString haxistitle = htitle+";"+xvarname[ivar]+ ";Counts";

			TH1F *htrue_CT = new TH1F(Form("htrue_CT"+xvar[ivar]+"%i",i),"Thrown_CT"+haxistitle,10,
										xvar_Hlimits[ivar][0],xvar_Hlimits[ivar][1]);
			htrue_CT_Vec.push_back(htrue_CT);
			TH1F *hcorr_CT = new TH1F(Form("hcorr_CT"+xvar[ivar]+"%i",i),"Corr_CT"+haxistitle,10,
										xvar_Hlimits[ivar][0],xvar_Hlimits[ivar][1]);
			hcorr_CT->Sumw2();
			hcorr_CT_Vec.push_back(hcorr_CT);
			TH1F *hratio = new TH1F(Form("hratio"+xvar[ivar]+"%i",i),"CT Ratio"+htitle+";"+xvarname[ivar]+";Corr/True",
										10,xvar_Hlimits[ivar][0],xvar_Hlimits[ivar][1]);
			hratio->Sumw2();
			hratio_Vec.push_back(hratio);
		}
	}

	// Correction-tree loop (Closure Test)
	for (int row=Nmiddle; row<Nentries; row++){
		tree->GetEntry(row);

		bool ecut=false, mc_ecut=false; // hcut=false, mc_hcut=false;

		int binHQ2, binHNu;
		if (TargType==target_n && Q2>Q2_Hlimits[0] && Q2<Q2_Hlimits[NHQ2] &&
			Yb<0.85 && W>2 && vyec>-1.4 && vyec<1.4 && Nu>Nu_Hlimits[0] && Nu<Nu_Hlimits[NHNu]){
			ecut = true;
			binHQ2 = var_position(NHQ2, Q2, Q2_Hlimits);
			binHNu = var_position(NHNu, Nu, Nu_Hlimits);
		}

		int binHmc_Q2, binHmc_Nu;
		if (mc_TargType==target_n && mc_Q2>Q2_Hlimits[0] && mc_Q2<Q2_Hlimits[NHQ2] &&
			mc_Yb<0.85 && mc_W>2 && mc_Nu>Nu_Hlimits[0] && mc_Nu<Nu_Hlimits[NHNu]){
			mc_ecut = true;
			binHmc_Q2 = var_position(NHQ2, mc_Q2, Q2_Hlimits);
			binHmc_Nu = var_position(NHNu, mc_Nu, Nu_Hlimits);
		}

		if (ecut==false && mc_ecut==false) continue;

		int ientries = PhiPQ->size();
		for (int i=0; i<ientries; i++){
			if (ecut==true && (*pid)[i]==211 && (*Zh)[i]>Zh_limits[0] && (*Zh)[i]<Zh_limits[NZh] &&
				(*Pt2)[i]>Pt2_limits[0] && (*Pt2)[i]<Pt2_limits[NPt2] && (*PhiPQ)[i]>PhiPQ_limits[0] &&
				(*PhiPQ)[i]<PhiPQ_limits[NPhiPQ]){

				Double_t accbin[] = {Q2, Nu, (*Zh)[i], (*Pt2)[i], (*PhiPQ)[i]};
				int corbin = binHNu + binHQ2*NHNu;

				Int_t bin = hacc->GetBin(accbin);
				Double_t acc_value = hacc->GetBinContent(bin);
				if (acc_value!=0){
					hcorr_CT_Vec[corbin + 0*NHbins]->Fill((*Zh)[i],1./acc_value);
					hcorr_CT_Vec[corbin + 1*NHbins]->Fill((*Pt2)[i],1./acc_value);
					hcorr_CT_Vec[corbin + 2*NHbins]->Fill((*PhiPQ)[i],1./acc_value);
				}
			}

			if (mc_ecut==true && (*mc_pid)[i]==211 && (*mc_Zh)[i]>Zh_limits[0] && (*mc_Zh)[i]<Zh_limits[NZh] &&
				(*mc_Pt2)[i]>Pt2_limits[0] && (*mc_Pt2)[i]<Pt2_limits[NPt2] && (*mc_PhiPQ)[i]>PhiPQ_limits[0] &&
				(*mc_PhiPQ)[i]<PhiPQ_limits[NPhiPQ]){

				int truebin = binHmc_Nu + binHmc_Q2*NHNu;
				htrue_CT_Vec[truebin + 0*NHbins]->Fill((*mc_Zh)[i]);
				htrue_CT_Vec[truebin + 1*NHbins]->Fill((*mc_Pt2)[i]);
				htrue_CT_Vec[truebin + 2*NHbins]->Fill((*mc_PhiPQ)[i]);
			}
		} // end of hadrons' loop
  	} // end filling loop

	// Apply Closure Test
	for (int i=0; i<NHbins*Nxvar; i++){
		hratio_Vec[i]->Divide(hcorr_CT_Vec[i],htrue_CT_Vec[i]);

		delete htrue_CT_Vec[i];
		delete hcorr_CT_Vec[i];
	}

	// Finishing (Save and close files)
	output->Write();
	output->Close();
} // End of the macro


////// DEFINITION OF FUNCTIONS

int var_position(int Nvar, float var, float var_limits[]){
	for (int ivar=0; ivar<Nvar; ivar++){
		if (var>=var_limits[ivar] && var<var_limits[ivar+1]){
			return ivar;
		}
	}
	return -9999;
}

void correct(float var_to_corr, float var_acc, TH1F *hacc, TH1F *hcorr){
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

void setBinVarSize(THnSparse *hist, Int_t nbins[], Double_t Q2_limits[], Double_t Nu_limits[],
					Double_t Zh_limits[], Double_t Pt2_limits[], Double_t PhiPQ_limits[]){
	hist->GetAxis(0)->Set(nbins[0], Q2_limits);
	hist->GetAxis(1)->Set(nbins[1], Nu_limits);
	hist->GetAxis(2)->Set(nbins[2], Zh_limits);
	hist->GetAxis(3)->Set(nbins[3], Pt2_limits);
	hist->GetAxis(4)->Set(nbins[4], PhiPQ_limits);
	return 0;
}