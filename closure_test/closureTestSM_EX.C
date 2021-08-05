/* Closure Test applied to S.Morán's acceptance */
#include "functions.h"

void closureTestSM_EX(TString target = "Fe", TString xvar = "Zh"){
	if (xvar!="Zh" && xvar!="Pt2"){
		std::cout << "ERROR: Variable name" << std::endl;
		return 0;
	}

	// I/O Files
	TString in_file = "../../clas-HSim/hsim_"+target+"1.root";
    TFile *input = TFile::Open(in_file,"READ");
    TTree *tree = (TTree*) input->Get("ntuple_sim");

    TString out_file = "ClosTestSM_"+target+"1_EX_"+xvar+".root";
    TFile *output = TFile::Open(out_file,"RECREATE");

	// Limits and sizes of cuts
	std::vector<float> Q2_Hlimits = {1.0, 1.3, 1.8, 4.1};
	std::vector<float> Nu_Hlimits = {2.2, 3.2, 3.7, 4.2};

	std::vector<float> Q2_limits = {1.0, 1.2, 1.3, 1.5, 1.8, 2.8, 4.1};
	std::vector<float> Nu_limits = {2.2, 2.8, 3.2, 3.5, 3.7, 3.9, 4.2};
	// std::vector<float> Xb_limits = {0.12, 0.57}; // {0.12, 0.19, 0.23, 0.27, 0.33, 0.57}; -> Not needed (Nu cut is stronger)
	std::vector<float> Zh_limits;
	std::vector<float> Pt2_limits;
	std::vector<float> PhiPQ_limits = {-180.0, 180.0};

	if (xvar=="Zh"){
		Zh_limits = {0.0,1.0};
		Pt2_limits = {0.0, 0.03, 0.06, 0.1, 0.18, 1.0};
	}
	else if (xvar=="Pt2"){
		Zh_limits = {0.0, 0.15, 0.25, 0.4, 0.7, 1.0};
		Pt2_limits = {0.0, 1.0};
	}

	// N<var> is the number of bins (array's elements - 1), N<var>+1 is the number of limits (array's elements)
	int NHQ2 = Q2_Hlimits.size() -1;
	int NHNu = Nu_Hlimits.size() -1;
	int NQ2 = Q2_limits.size() -1;
	int NNu = Nu_limits.size() -1;
	// int NXb = Xb_limits.size() -1;
	int NZh = Zh_limits.size() -1;
	int NPt2 = Pt2_limits.size() -1;
	int NPhiPQ = PhiPQ_limits.size() -1;

	int NHbins = NHQ2*NHNu; // Total number of plot bins
	int Nhadvars = NZh*NPt2*NPhiPQ;
	int NTot = NQ2*NNu*NZh*NPt2*NPhiPQ; // Total number of cut bins

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

	// Create histograms
    // std::vector<TH1F*> hreco_fin_Vec;
    // std::vector<TH1F*> htrue_fin_Vec;
    // std::vector<TH1F*> hacc_fin_Vec;

    std::vector<TH1F*> hreco_CT_fin_Vec;
    std::vector<TH1F*> htrue_CT_fin_Vec;
    std::vector<TH1F*> hcorr_CT_fin_Vec;

    std::vector<TH1F*> hratio_fin_Vec;

	TString xvarname;
    if (xvar=="Zh") xvarname = "Z_{h}";
    else if (xvar=="Pt2") xvarname = "P_{t}^{2} [GeV^{2}]";

	for (int i=0; i<NHbins; i++){
		TString htitle = Form(" ("+target+"), %.1f<Q2<%.1f, %.1f<Nu<%.1f;"+xvarname+ ";Counts",
                              Q2_Hlimits[i/NHNu],Q2_Hlimits[i/NHNu+1],Nu_Hlimits[i%NHNu],Nu_Hlimits[(i%NHNu)+1]);
		// TH1F *hreco_fin = new TH1F(Form("hreco_fin%i",i),"Reco"+htitle,10,0.,1.);
		// hreco_fin->Sumw2();
		// hreco_fin_Vec.push_back(hreco_fin);
		// TH1F *htrue_fin = new TH1F(Form("htrue_fin%i",i),"Thrown"+htitle,10,0.,1.);
		// htrue_fin->Sumw2();
		// htrue_fin_Vec.push_back(htrue_fin);
		// TH1F *hacc_fin = new TH1F(Form("hacc_fin%i",i),"Acc"+htitle,10,0.,1.);
		// hacc_fin->Sumw2();
		// hacc_fin_Vec.push_back(hacc_fin);

		TH1F *hreco_CT_fin = new TH1F(Form("hreco_CT_fin%i",i),"Reco_CT"+htitle,10,0.,1.);
		hreco_CT_fin->Sumw2();
		hreco_CT_fin_Vec.push_back(hreco_CT_fin);
		TH1F *htrue_CT_fin = new TH1F(Form("htrue_CT_fin%i",i),"Thrown_CT"+htitle,10,0.,1.);
		htrue_CT_fin->Sumw2();
		htrue_CT_fin_Vec.push_back(htrue_CT_fin);
		TH1F *hcorr_CT_fin = new TH1F(Form("hcorr_CT_fin%i",i),"Corr_CT"+htitle,10,0.,1.);
		hcorr_CT_fin->Sumw2();
		hcorr_CT_fin_Vec.push_back(hcorr_CT_fin);

		TH1F *hratio_fin = new TH1F(Form("hratio%i",i),"CT Ratio;"+xvarname+";Corr/True",10,0.,1.);
		hratio_fin->Sumw2();
		hratio_fin_Vec.push_back(hratio_fin);
	}

	std::vector<TH1F*> hreco_Vec;
	std::vector<TH1F*> htrue_Vec;
	std::vector<TH1F*> hacc_Vec;

	std::vector<TH1F*> hreco_CT_Vec;
	std::vector<TH1F*> htrue_CT_Vec;
	std::vector<TH1F*> hcorr_CT_Vec;

	for (int iQ2=0; iQ2<NQ2; iQ2++){
		for (int iNu=0; iNu<NNu; iNu++){
			for (int iZh=0; iZh<NZh; iZh++){
				for (int iPt2=0; iPt2<NPt2; iPt2++){
					int iTot = iPt2 + iZh*NPt2 + iNu*NZh*NPt2 + iQ2*NNu*NZh*NPt2;
					TH1F *hreco_tmp = new TH1F(Form("hreco%d",iTot),Form("Reco %d (%d, %d, %d, %d)",
											   iTot,iQ2,iNu,iZh,iPt2),10,0.,1.);
					hreco_Vec.push_back(hreco_tmp);
					TH1F *htrue_tmp = new TH1F(Form("htrue%d",iTot),Form("Thrown %d (%d, %d, %d, %d)",
										       iTot,iQ2,iNu,iZh,iPt2),10,0.,1.);
					htrue_Vec.push_back(htrue_tmp);
					TH1F *hacc_tmp = new TH1F(Form("hacc%d",iTot),Form("Acc %d (%d, %d, %d, %d)",
											  iTot,iQ2,iNu,iZh,iPt2),10,0.,1.);
					hacc_tmp->Sumw2();
					hacc_Vec.push_back(hacc_tmp);

					TH1F *hreco_CT_tmp = new TH1F(Form("hrecoCT%d",iTot),Form("Reco %d (%d, %d, %d, %d)",
											   iTot,iQ2,iNu,iZh,iPt2),10,0.,1.);
					hreco_CT_Vec.push_back(hreco_CT_tmp);
					TH1F *htrue_CT_tmp = new TH1F(Form("htrueCT%d",iTot),Form("Thrown %d (%d, %d, %d, %d)",
										       iTot,iQ2,iNu,iZh,iPt2),10,0.,1.);
					htrue_CT_Vec.push_back(htrue_CT_tmp);
					TH1F *hcorr_CT_tmp = new TH1F(Form("hcorrCT%d",iTot),Form("Acc %d (%d, %d, %d, %d)",
											  iTot,iQ2,iNu,iZh,iPt2),10,0.,1.);
					hcorr_CT_tmp->Sumw2();
					hcorr_CT_Vec.push_back(hcorr_CT_tmp);
				}
			}
		}
	}

	// Acceptance-tree loop

	// int Nentries = 100; // for testing
	int Nentries = tree->GetEntries();
	int Nmiddle = Nentries/2;

	std::cout << "Nmiddle: " << Nmiddle << "; Nentries: " << Nentries << std::endl;

	for (int row=0; row<Nmiddle; row++){
		tree->GetEntry(row);

		bool ecut=false, mc_ecut=false; // hcut=false, mc_hcut=false;

		int binQ2=-1, binNu=-1;
		if (TargType==target_n && Q2>Q2_limits[0] && Q2<Q2_limits[NQ2] && // Xb>Xb_limits[0] && Xb<Xb_limits[NXb] &&
			Yb<0.85 && W>2 && vyec>-1.4 && vyec<1.4 && Nu>Nu_limits[0] && Nu<Nu_limits[NNu]){
			ecut = true;
			binQ2 = var_position(NQ2, Q2, Q2_limits);
			binNu = var_position(NNu, Nu, Nu_limits);
		}

		int binmc_Q2=-1, binmc_Nu=-1;
		if (mc_TargType==target_n && mc_Q2>Q2_limits[0] && mc_Q2<Q2_limits[NQ2] && // mc_Xb>Xb_limits[0] && mc_Xb<Xb_limits[NXb] &&
			mc_Yb<0.85 && mc_W>2 && mc_Nu>Nu_limits[0] && mc_Nu<Nu_limits[NNu]){
			mc_ecut = true;
			binmc_Q2 = var_position(NQ2, mc_Q2, Q2_limits);
			binmc_Nu = var_position(NNu, mc_Nu, Nu_limits);
		}

		if (ecut==false && mc_ecut==false) continue;

		int ientries = PhiPQ->size();
		for (int i=0; i<ientries; i++){
			int binZh=-1, binPt2=-1;
			if (ecut==true && (*pid)[i]==211 && (*Zh)[i]>Zh_limits[0] && (*Zh)[i]<Zh_limits[NZh] && (*Pt2)[i]>Pt2_limits[0] &&
				(*Pt2)[i]<Pt2_limits[NPt2] && (*PhiPQ)[i]>PhiPQ_limits[0] && (*PhiPQ)[i]<PhiPQ_limits[NPhiPQ]){
				// hcut = true;
				binZh = var_position(i, NZh, Zh, Zh_limits);
				binPt2 = var_position(i, NPt2, Pt2, Pt2_limits);

				int fillbin = binPt2 + binZh*NPt2 + binNu*NZh*NPt2 + binQ2*NNu*NZh*NPt2;
				if (xvar=="Zh") hreco_Vec[fillbin]->Fill((*Zh)[i]);
				else if (xvar=="Pt2") hreco_Vec[fillbin]->Fill((*Pt2)[i]);
			}

			int binmc_Zh=-1, binmc_Pt2=-1;
			if (mc_ecut==true && (*mc_pid)[i]==211 && (*mc_Zh)[i]>Zh_limits[0] && (*mc_Zh)[i]<Zh_limits[NZh] &&
				(*mc_Pt2)[i]>Pt2_limits[0] && (*mc_Pt2)[i]<Pt2_limits[NPt2] && (*mc_PhiPQ)[i]>PhiPQ_limits[0] &&
				(*mc_PhiPQ)[i]<PhiPQ_limits[NPhiPQ]){
				// mc_hcut = true;
				binmc_Zh = var_position(i, NZh, mc_Zh, Zh_limits);
				binmc_Pt2 = var_position(i, NPt2, mc_Pt2, Pt2_limits);

				int fillmc_bin = binmc_Pt2 + binmc_Zh*NPt2 + binmc_Nu*NZh*NPt2 + binmc_Q2*NNu*NZh*NPt2;
				if (xvar=="Zh") htrue_Vec[fillmc_bin]->Fill((*mc_Zh)[i]);
				else if (xvar=="Pt2") htrue_Vec[fillmc_bin]->Fill((*mc_Pt2)[i]);
			}
		} // end of hadrons' loop
  	} // end filling loop

	// Calculating acceptance
	for (int i=0; i<NTot; i++){
		// int digitQ2 = i/(Nhadvars*NNu);
		// int digitNu = i/(Nhadvars)-digitQ2*NNu;

		// 2-factor comes from the fact that a bin in <var>_Hlimits is two bins of <var>_limits merged
		// (a,b,c,d) => (a,b)->(a+b=A) & (c,d)->(c+d=B) <=> (a,b,c,d)->(A,B)
		// hreco_fin_Vec[digitNu/2 + (digitQ2/2)*NHNu]->Add(hreco_Vec[i]); // I think that this is not necessary
		// htrue_fin_Vec[digitNu/2 + (digitQ2/2)*NHNu]->Add(htrue_Vec[i]); // I think that this is not necessary

		hacc_Vec[i]->Divide(hreco_Vec[i],htrue_Vec[i],1,1,"B");
		delete hreco_Vec[i];
		delete htrue_Vec[i];
	}

	// Correction-tree loop (Closure Test)
	for (int row=Nmiddle; row<Nentries; row++){
		tree->GetEntry(row);

		bool ecut=false, mc_ecut=false; // hcut=false, mc_hcut=false;

		int binQ2=-1, binNu=-1;
		if (TargType==target_n && Q2>Q2_limits[0] && Q2<Q2_limits[NQ2] && // Xb>Xb_limits[0] && Xb<Xb_limits[NXb] &&
			Yb<0.85 && W>2 && vyec>-1.4 && vyec<1.4 && Nu>Nu_limits[0] && Nu<Nu_limits[NNu]){
			ecut = true;
			binQ2 = var_position(NQ2, Q2, Q2_limits);
			binNu = var_position(NNu, Nu, Nu_limits);
		}

		int binmc_Q2=-1, binmc_Nu=-1;
		if (mc_TargType==target_n && mc_Q2>Q2_limits[0] && mc_Q2<Q2_limits[NQ2] && // mc_Xb>Xb_limits[0] && mc_Xb<Xb_limits[NXb] &&
			mc_Yb<0.85 && mc_W>2 && mc_Nu>Nu_limits[0] && mc_Nu<Nu_limits[NNu]){
			mc_ecut = true;
			binmc_Q2 = var_position(NQ2, mc_Q2, Q2_limits);
			binmc_Nu = var_position(NNu, mc_Nu, Nu_limits);
		}

		if (ecut==false && mc_ecut==false) continue;

		int ientries = PhiPQ->size();
		for (int i=0; i<ientries; i++){
			int binZh=-1, binPt2=-1;
			if (ecut==true && (*pid)[i]==211 && (*Zh)[i]>Zh_limits[0] && (*Zh)[i]<Zh_limits[NZh] && (*Pt2)[i]>Pt2_limits[0] &&
				(*Pt2)[i]<Pt2_limits[NPt2] && (*PhiPQ)[i]>PhiPQ_limits[0] && (*PhiPQ)[i]<PhiPQ_limits[NPhiPQ]){
				// hcut = true;
				binZh = var_position(i, NZh, Zh, Zh_limits);
				binPt2 = var_position(i, NPt2, Pt2, Pt2_limits);

				int fillbin = binPt2 + binZh*NPt2 + binNu*NZh*NPt2 + binQ2*NNu*NZh*NPt2;
				if (xvar=="Zh") hreco_CT_Vec[fillbin]->Fill((*Zh)[i]);
				else if (xvar=="Pt2") hreco_CT_Vec[fillbin]->Fill((*Pt2)[i]);
			}

			int binmc_Zh=-1, binmc_Pt2=-1;
			if (mc_ecut==true && (*mc_pid)[i]==211 && (*mc_Zh)[i]>Zh_limits[0] && (*mc_Zh)[i]<Zh_limits[NZh] &&
				(*mc_Pt2)[i]>Pt2_limits[0] && (*mc_Pt2)[i]<Pt2_limits[NPt2] && (*mc_PhiPQ)[i]>PhiPQ_limits[0] &&
				(*mc_PhiPQ)[i]<PhiPQ_limits[NPhiPQ]){
				// mc_hcut = true;
				binmc_Zh = var_position(i, NZh, mc_Zh, Zh_limits);
				binmc_Pt2 = var_position(i, NPt2, mc_Pt2, Pt2_limits);
				
				int fillmc_bin = binmc_Pt2 + binmc_Zh*NPt2 + binmc_Nu*NZh*NPt2 + binmc_Q2*NNu*NZh*NPt2;
				if (xvar=="Zh") htrue_CT_Vec[fillmc_bin]->Fill((*mc_Zh)[i]);
				else if (xvar=="Pt2") htrue_CT_Vec[fillmc_bin]->Fill((*mc_Pt2)[i]);
			}
		} // end of hadrons' loop
  	} // end filling loop

	// Calculating correction and Closure Test
	for (int i=0; i<NTot; i++){
		int digitQ2 = i/(Nhadvars*NNu);
		int digitNu = i/(Nhadvars)-digitQ2*NNu;

		// 2-factor comes from the fact that a bin in <var>_Hlimits is two bins of <var>_limits merged
		// (a,b,c,d) => (a,b)->(a+b=A) & (c,d)->(c+d=B) <=> (a,b,c,d)->(A,B)
		hreco_CT_fin_Vec[digitNu/2 + (digitQ2/2)*NHNu]->Add(hreco_CT_Vec[i]);
		htrue_CT_fin_Vec[digitNu/2 + (digitQ2/2)*NHNu]->Add(htrue_CT_Vec[i]);

		hcorr_CT_Vec[i]->Divide(hreco_CT_Vec[i],hacc_Vec[i]);
		hcorr_CT_fin_Vec[digitNu/2 + (digitQ2/2)*NHNu]->Add(hcorr_CT_Vec[i]);

		delete hreco_CT_Vec[i];
		delete htrue_CT_Vec[i];
		delete hcorr_CT_Vec[i];
		// delete hacc_Vec[i];
	}

	for (int i=0; i<NHbins; i++){
		hratio_fin_Vec[i]->Divide(hcorr_CT_fin_Vec[i],htrue_CT_fin_Vec[i]);

		delete hreco_CT_fin_Vec[i];
		delete htrue_CT_fin_Vec[i];
	}

	// Finishing (Save and close files)
	output->Write();
	output->Close();
} // End of the macro


////// DEFINITION OF FUNCTIONS

int var_position(int Nvar, float var, std::vector<float> var_limits){
	for (int ivar=0; ivar<Nvar; ivar++){
		if (var>=var_limits[ivar] && var<var_limits[ivar+1]){
			return ivar;
		}
	}
	return -1;
}

int var_position(int i, int Nvar, std::vector<float> *var, std::vector<float> var_limits){
	for (int ivar=0; ivar<Nvar; ivar++){
		if ((*var)[i]>=var_limits[ivar] && (*var)[i]<var_limits[ivar+1]){
			return ivar;
		}
	}
	return -1;
}
