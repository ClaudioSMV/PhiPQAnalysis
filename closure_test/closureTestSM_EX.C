/* Closure Test applied to S.Morán's acceptance */
#include "functions.h"

void closureTestSM_EX(TString target = "Fe", TString xvar = "Zh"){
	if (xvar!="Zh" && xvar!="Pt2" && xvar!="PhiPQ"){
		std::cout << "ERROR: Variable name" << std::endl;
		return 0;
	}

	// I/O Files
	TString in_file = "../../clas-HSim/hsim_"+target+"1.root";
    TFile *input = TFile::Open(in_file,"READ");
    TTree *tree = (TTree*) input->Get("ntuple_sim");

    TString out_file = "ClosTestSM_"+target+"1_EX_"+xvar+".root";
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
	float Q2_limits[] = {1.0, 1.3, 1.8, 4.1}; // {1.0, 1.2, 1.3, 1.5, 1.8, 2.8, 4.1};
	float Nu_limits[] = {2.2, 3.2, 3.7, 4.2}; // {2.2, 2.8, 3.2, 3.5, 3.7, 3.9, 4.2};
	float Zh_limits[] = {0.0, 0.15, 0.25, 0.4, 0.7, 1.0};
	float Pt2_limits[] = {0.0, 0.03, 0.06, 0.1, 0.18, 1.0};
	int NBinsPhiPQ = 12; // Should be 180 but it's too much... 
    float PhiPQ_limits[NBinsPhiPQ+1];
    for (int i=0; i<NBinsPhiPQ+1; i++){
        PhiPQ_limits[i] = -180. + i*360./NBinsPhiPQ;
    }

	// Limits of final histograms
	std::vector<float> Q2_Hlimits = {1.0, 1.3, 1.8, 4.1};
	std::vector<float> Nu_Hlimits = {2.2, 3.2, 3.7, 4.2};
	std::vector<float> xvar_Hlimits;
	if (xvar=="Zh" || xvar=="Pt2") xvar_Hlimits = {0.0, 1.0};
	else if (xvar=="PhiPQ") xvar_Hlimits = {-180., 180.};

	// N<var> is the number of bins (array's elements - 1), N<var>+1 is the number of limits (array's elements)
	int NQ2  = sizeof(Q2_limits)/sizeof(Q2_limits[0]) - 1;
	int NNu  = sizeof(Nu_limits)/sizeof(Nu_limits[0]) - 1;
	int NZh  = sizeof(Zh_limits)/sizeof(Zh_limits[0]) - 1;
	int NPt2 = sizeof(Pt2_limits)/sizeof(Pt2_limits[0]) - 1;
	int NPhiPQ  = sizeof(PhiPQ_limits)/sizeof(PhiPQ_limits[0]) - 1;

	int NHQ2 = Q2_Hlimits.size() -1;
	int NHNu = Nu_Hlimits.size() -1;

	int Nhadvars = NZh*NPt2*NPhiPQ;
	int NTot = NQ2*NNu*NZh*NPt2*NPhiPQ; // Total number of acceptance-cut-bins
	int NHbins = NHQ2*NHNu; // Total number of plot bins

	// Define histogram's variables
	TString xvarname;
	int xvarhbin;
    if (xvar=="Zh"){ xvarname = "Z_{h}"; xvarhbin = 0;}
    else if (xvar=="Pt2"){ xvarname = "P_{t}^{2} [GeV^{2}]"; xvarhbin = 0;}
	else if (xvar=="PhiPQ"){ xvarname = "#phi_{PQ} [deg]"; xvarhbin = 1;}

	// Create acceptance histograms
	std::vector<TH1F*> hreco_Vec;
	std::vector<TH1F*> htrue_Vec;
	std::vector<TH1F*> hacc_Vec;

	for (int iQ2=0; iQ2<NQ2; iQ2++){
		for (int iNu=0; iNu<NNu; iNu++){
			for (int iZh=0; iZh<NZh; iZh++){
				for (int iPt2=0; iPt2<NPt2; iPt2++){
					for (int iPhiPQ=0; iPhiPQ<NPhiPQ; iPhiPQ++){
						int iTot = iPhiPQ + iPt2*NPhiPQ + iZh*NPt2*NPhiPQ + iNu*NZh*NPt2*NPhiPQ + iQ2*NNu*NZh*NPt2*NPhiPQ;
						TH1F *hreco_tmp, *htrue_tmp, *hacc_tmp;
						if (xvar=="Zh"){
							hreco_tmp = new TH1F(Form("hreco%d",iTot),Form("Reco %d (%d, %d, %d, %d, %d)",
												iTot,iQ2,iNu,iZh,iPt2,iPhiPQ),NZh,Zh_limits);
							htrue_tmp = new TH1F(Form("htrue%d",iTot),Form("Thrown %d (%d, %d, %d, %d, %d)",
												iTot,iQ2,iNu,iZh,iPt2,iPhiPQ),NZh,Zh_limits);
							hacc_tmp = new TH1F(Form("hacc%d",iTot),Form("Acc %d (%d, %d, %d, %d, %d)",
												iTot,iQ2,iNu,iZh,iPt2,iPhiPQ),NZh,Zh_limits);
						}
						else if (xvar=="Pt2"){
							hreco_tmp = new TH1F(Form("hreco%d",iTot),Form("Reco %d (%d, %d, %d, %d, %d)",
												iTot,iQ2,iNu,iZh,iPt2,iPhiPQ),NPt2,Pt2_limits);
							htrue_tmp = new TH1F(Form("htrue%d",iTot),Form("Thrown %d (%d, %d, %d, %d, %d)",
												iTot,iQ2,iNu,iZh,iPt2,iPhiPQ),NPt2,Pt2_limits);
							hacc_tmp = new TH1F(Form("hacc%d",iTot),Form("Acc %d (%d, %d, %d, %d, %d)",
												iTot,iQ2,iNu,iZh,iPt2,iPhiPQ),NPt2,Pt2_limits);
						}
						else if (xvar=="PhiPQ"){
							hreco_tmp = new TH1F(Form("hreco%d",iTot),Form("Reco %d (%d, %d, %d, %d, %d)",
												iTot,iQ2,iNu,iZh,iPt2,iPhiPQ),NPhiPQ,PhiPQ_limits);
							htrue_tmp = new TH1F(Form("htrue%d",iTot),Form("Thrown %d (%d, %d, %d, %d, %d)",
												iTot,iQ2,iNu,iZh,iPt2,iPhiPQ),NPhiPQ,PhiPQ_limits);
							hacc_tmp = new TH1F(Form("hacc%d",iTot),Form("Acc %d (%d, %d, %d, %d, %d)",
												iTot,iQ2,iNu,iZh,iPt2,iPhiPQ),NPhiPQ,PhiPQ_limits);
						}
						hreco_Vec.push_back(hreco_tmp);
						htrue_Vec.push_back(htrue_tmp);
						hacc_tmp->Sumw2();
						hacc_Vec.push_back(hacc_tmp);
					}
				}
			}
		}
	}

	// Acceptance-tree loop

	// int Nentries = 30; // for testing
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
			int binZh=-1, binPt2=-1, binPhiPQ=-1;
			if (ecut==true && (*pid)[i]==211 && (*Zh)[i]>Zh_limits[0] && (*Zh)[i]<Zh_limits[NZh] && (*Pt2)[i]>Pt2_limits[0] &&
				(*Pt2)[i]<Pt2_limits[NPt2] && (*PhiPQ)[i]>PhiPQ_limits[0] && (*PhiPQ)[i]<PhiPQ_limits[NPhiPQ]){
				// hcut = true;
				binZh = var_position(i, NZh, Zh, Zh_limits);
				binPt2 = var_position(i, NPt2, Pt2, Pt2_limits);
				binPhiPQ = var_position(i, NPhiPQ, PhiPQ, PhiPQ_limits);

				int fillbin = binPhiPQ + binPt2*NPhiPQ + binZh*NPt2*NPhiPQ + binNu*NZh*NPt2*NPhiPQ +
							  binQ2*NNu*NZh*NPt2*NPhiPQ;
				if (xvar=="Zh") hreco_Vec[fillbin]->Fill((*Zh)[i]);
				else if (xvar=="Pt2") hreco_Vec[fillbin]->Fill((*Pt2)[i]);
				else if (xvar=="PhiPQ") hreco_Vec[fillbin]->Fill((*PhiPQ)[i]);
			}

			int binmc_Zh=-1, binmc_Pt2=-1, binmc_PhiPQ;
			if (mc_ecut==true && (*mc_pid)[i]==211 && (*mc_Zh)[i]>Zh_limits[0] && (*mc_Zh)[i]<Zh_limits[NZh] &&
				(*mc_Pt2)[i]>Pt2_limits[0] && (*mc_Pt2)[i]<Pt2_limits[NPt2] && (*mc_PhiPQ)[i]>PhiPQ_limits[0] &&
				(*mc_PhiPQ)[i]<PhiPQ_limits[NPhiPQ]){
				// mc_hcut = true;
				binmc_Zh = var_position(i, NZh, mc_Zh, Zh_limits);
				binmc_Pt2 = var_position(i, NPt2, mc_Pt2, Pt2_limits);
				binmc_PhiPQ = var_position(i, NPhiPQ, mc_PhiPQ, PhiPQ_limits);

				int fillmc_bin = binmc_PhiPQ + binmc_Pt2*NPhiPQ + binmc_Zh*NPt2*NPhiPQ + 
								 binmc_Nu*NZh*NPt2*NPhiPQ + binmc_Q2*NNu*NZh*NPt2*NPhiPQ;
				if (xvar=="Zh") htrue_Vec[fillmc_bin]->Fill((*mc_Zh)[i]);
				else if (xvar=="Pt2") htrue_Vec[fillmc_bin]->Fill((*mc_Pt2)[i]);
				else if (xvar=="PhiPQ") htrue_Vec[fillmc_bin]->Fill((*mc_PhiPQ)[i]);
			}
		} // end of hadrons' loop
  	} // end filling loop

	// Calculating acceptance
	for (int i=0; i<NTot; i++){
		hacc_Vec[i]->Divide(hreco_Vec[i],htrue_Vec[i],1,1,"B");
		delete hreco_Vec[i];
		delete htrue_Vec[i];
	}

	// Create correction and ratio histograms
    std::vector<TH1F*> htrue_CT_Vec;
    std::vector<TH1F*> hcorr_CT_Vec;
    std::vector<TH1F*> hratio_Vec;

	for (int i=0; i<NHbins; i++){
		TString htitle = Form(" ("+target+"), %.1f<Q2<%.1f, %.1f<Nu<%.1f",Q2_Hlimits[i/NHNu],
							  Q2_Hlimits[i/NHNu+1],Nu_Hlimits[i%NHNu],Nu_Hlimits[(i%NHNu)+1]);
		TString haxistitle = htitle+";"+xvarname+ ";Counts";

		// TH1F *hreco_CT = new TH1F(Form("hreco_CT%i",i),"Reco_CT"+haxistitle,10,
		// 							   xvar_Hlimits[0],xvar_Hlimits[1]);
		// hreco_CT->Sumw2();
		// hreco_CT_Vec.push_back(hreco_CT);
		TH1F *htrue_CT = new TH1F(Form("htrue_CT%i",i),"Thrown_CT"+haxistitle,10,
									   xvar_Hlimits[0],xvar_Hlimits[1]);
		// htrue_CT->Sumw2();
		htrue_CT_Vec.push_back(htrue_CT);
		TH1F *hcorr_CT = new TH1F(Form("hcorr_CT%i",i),"Corr_CT"+haxistitle,10,
								       xvar_Hlimits[0],xvar_Hlimits[1]);
		hcorr_CT->Sumw2();
		hcorr_CT_Vec.push_back(hcorr_CT);
		TH1F *hratio = new TH1F(Form("hratio%i",i),"CT Ratio"+htitle+";"+xvarname+";Corr/True",10,
						             xvar_Hlimits[0],xvar_Hlimits[1]);
		hratio->Sumw2();
		hratio_Vec.push_back(hratio);
	}
	TH1F *htrue_CT_all = new TH1F("htrue_CT_all","Thrown_CT All;"+xvarname+ ";Counts",10,xvar_Hlimits[0],xvar_Hlimits[1]);
	TH1F *hcorr_CT_all = new TH1F("hcorr_CT_all","Corr_CT All;"+xvarname+ ";Counts",10,xvar_Hlimits[0],xvar_Hlimits[1]);
	TH1F *hratio_all = new TH1F("hratio_all","CT Ratio All;"+xvarname+";Corr/True",10,xvar_Hlimits[0],xvar_Hlimits[1]);
	hratio_all->Sumw2();

	// Correction-tree loop (Closure Test)
	for (int row=Nmiddle; row<Nentries; row++){
		tree->GetEntry(row);

		bool ecut=false, mc_ecut=false; // hcut=false, mc_hcut=false;

		int binQ2=-1, binNu=-1;
		if (TargType==target_n && Q2>Q2_Hlimits[0] && Q2<Q2_Hlimits[NHQ2] &&
			Yb<0.85 && W>2 && vyec>-1.4 && vyec<1.4 && Nu>Nu_Hlimits[0] && Nu<Nu_Hlimits[NHNu]){
			ecut = true;
			binQ2 = var_position(NHQ2, Q2, Q2_Hlimits);
			binNu = var_position(NHNu, Nu, Nu_Hlimits);
		}

		int binmc_Q2=-1, binmc_Nu=-1;
		if (mc_TargType==target_n && mc_Q2>Q2_Hlimits[0] && mc_Q2<Q2_Hlimits[NHQ2] &&
			mc_Yb<0.85 && mc_W>2 && mc_Nu>Nu_Hlimits[0] && mc_Nu<Nu_Hlimits[NHNu]){
			mc_ecut = true;
			binmc_Q2 = var_position(NHQ2, mc_Q2, Q2_Hlimits);
			binmc_Nu = var_position(NHNu, mc_Nu, Nu_Hlimits);
		}

		if (ecut==false && mc_ecut==false) continue;

		int ientries = PhiPQ->size();
		for (int i=0; i<ientries; i++){
			int binZh=-1, binPt2=-1, binPhiPQ=-1;
			if (ecut==true && (*pid)[i]==211 && (*Zh)[i]>Zh_limits[0] && (*Zh)[i]<Zh_limits[NZh] && (*Pt2)[i]>Pt2_limits[0] &&
				(*Pt2)[i]<Pt2_limits[NPt2] && (*PhiPQ)[i]>PhiPQ_limits[0] && (*PhiPQ)[i]<PhiPQ_limits[NPhiPQ]){
				// hcut = true;
				binZh = var_position(i, NZh, Zh, Zh_limits);
				binPt2 = var_position(i, NPt2, Pt2, Pt2_limits);
				binPhiPQ = var_position(i, NPhiPQ, PhiPQ, PhiPQ_limits);

				int accbin = binPhiPQ + binPt2*NPhiPQ + binZh*NPt2*NPhiPQ +
							 binNu*NZh*NPt2*NPhiPQ + binQ2*NNu*NZh*NPt2*NPhiPQ;
				int corbin = binNu + binQ2*NHNu;
				if (xvar=="Zh"){
					// hreco_CT_Vec[corbin]->Fill((*Zh)[i]);
					correct((*Zh)[i], hacc_Vec[accbin], hcorr_CT_Vec[corbin]);
					correct((*Zh)[i], hacc_Vec[accbin], hcorr_CT_all);
				}
				else if (xvar=="Pt2"){
					// hreco_CT_Vec[corbin]->Fill((*Pt2)[i]);
					correct((*Pt2)[i], hacc_Vec[accbin], hcorr_CT_Vec[corbin]);
					correct((*Pt2)[i], hacc_Vec[accbin], hcorr_CT_all);
				}
				else if (xvar=="PhiPQ"){
					// hreco_CT_Vec[corbin]->Fill((*PhiPQ)[i]);
					correct((*PhiPQ)[i], hacc_Vec[accbin], hcorr_CT_Vec[corbin]);
					correct((*PhiPQ)[i], hacc_Vec[accbin], hcorr_CT_all);
				}
			}

			int binmc_Zh=-1, binmc_Pt2=-1, binmc_PhiPQ=-1;
			if (mc_ecut==true && (*mc_pid)[i]==211 && (*mc_Zh)[i]>Zh_limits[0] && (*mc_Zh)[i]<Zh_limits[NZh] &&
				(*mc_Pt2)[i]>Pt2_limits[0] && (*mc_Pt2)[i]<Pt2_limits[NPt2] && (*mc_PhiPQ)[i]>PhiPQ_limits[0] &&
				(*mc_PhiPQ)[i]<PhiPQ_limits[NPhiPQ]){
				// mc_hcut = true;
				// binmc_Zh = var_position(i, NZh, mc_Zh, Zh_limits);
				// binmc_Pt2 = var_position(i, NPt2, mc_Pt2, Pt2_limits);
				// int fillmc_bin = binmc_Pt2 + binmc_Zh*NPt2 + binmc_Nu*NZh*NPt2 + binmc_Q2*NNu*NZh*NPt2;

				int truebin = binmc_Nu + binmc_Q2*NHNu;
				if (xvar=="Zh"){htrue_CT_Vec[truebin]->Fill((*mc_Zh)[i]);
					htrue_CT_all->Fill((*mc_Zh)[i]);}
				else if (xvar=="Pt2"){htrue_CT_Vec[truebin]->Fill((*mc_Pt2)[i]);
					htrue_CT_all->Fill((*mc_Pt2)[i]);}
				else if (xvar=="PhiPQ"){htrue_CT_Vec[truebin]->Fill((*mc_PhiPQ)[i]);
					htrue_CT_all->Fill((*mc_PhiPQ)[i]);}
			}
		} // end of hadrons' loop
  	} // end filling loop

	for (int i=0; i<NTot; i++){
		delete hacc_Vec[i];
	}

	// Apply Closure Test
	for (int i=0; i<NHbins; i++){
		hratio_Vec[i]->Divide(hcorr_CT_Vec[i],htrue_CT_Vec[i]);

		// delete hreco_CT_Vec[i];
		delete htrue_CT_Vec[i];
		delete hcorr_CT_Vec[i];
	}
	hratio_all->Divide(hcorr_CT_all,htrue_CT_all);

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
	return -1;
}

int var_position(int Nvar, float var, std::vector<float> var_limits){
	for (int ivar=0; ivar<Nvar; ivar++){
		if (var>=var_limits[ivar] && var<var_limits[ivar+1]){
			return ivar;
		}
	}
	return -1;
}

int var_position(int i, int Nvar, std::vector<float> *var, float var_limits[]){
	for (int ivar=0; ivar<Nvar; ivar++){
		if ((*var)[i]>=var_limits[ivar] && (*var)[i]<var_limits[ivar+1]){
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

void correct(float var,TH1F *hacc,TH1F *hcorr){
    int bin = hacc->FindBin(var);
    float accept = hacc->GetBinContent(bin);
    if (accept!=0){
        float erracc = hacc->GetBinError(bin);
        hcorr->Fill(var,1./accept);
    }
    else if (accept==0){
        std::cout << "Zh: " << var << "; " << hacc->GetTitle() << std::endl;
    }
    return 0;
}