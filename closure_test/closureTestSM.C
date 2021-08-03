/* Closure Test applied to S.Morán's acceptance */
#include "functions.h"

void closureTestSM(TString target = "Fe"){
	// I/O Files
	TString in_file = "../../clas-HSim/hsim_"+target+"1.root";
    TFile *input = TFile::Open(in_file,"READ");
    TTree *tree = (TTree*) input->Get("ntuple_sim");

    TString out_file = "ClosTestSM_"+target+"1.root";
    TFile *output = TFile::Open(out_file,"RECREATE");

	// Limits and sizes of cuts
	// float Q2_limits[] = {1.0, 4.0};
	// float Xb_limits[] = {0.12, 0.57};
	// float Zh_limits[] = {0.0, 1.0};
	// float Pt2_limits[] = {0.0, 1.0};
	// float PhiPQ_limits[] = {-180.0, 180.0};

	float Q2_limits[] = {1.0, 1.17, 1.33, 1.51, 1.75, 2.12, 4.0};
	float Xb_limits[] = {0.12, 0.19, 0.23, 0.27, 0.33, 0.57};
	float Zh_limits[] = {0.0, 1.0};
	float Pt2_limits[] = {0.0, 0.03, 0.06, 0.1, 0.18, 1.0};
	float PhiPQ_limits[] = {-180.0, 180.0};

	// N<var> is the number of bins (array's elements - 1), N<var>+1 is the number of limits (array's elements)
	int NQ2  = sizeof(Q2_limits)/sizeof(Q2_limits[0]) - 1;
	int NXb  = sizeof(Xb_limits)/sizeof(Xb_limits[0]) - 1;
	int NZh  = sizeof(Zh_limits)/sizeof(Zh_limits[0]) - 1;
	int NPt2 = sizeof(Pt2_limits)/sizeof(Pt2_limits[0]) - 1;
	int NPhiPQ  = sizeof(PhiPQ_limits)/sizeof(PhiPQ_limits[0]) - 1;

	int NTot = NQ2*NXb*NZh*NPt2*NPhiPQ; // Total number of bins

	// Definition of tree-variables
	int target_n = 2;
	if (target == "D") target_n = 1;

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
	TH1F *hreco_fin = new TH1F("hreco_fin","Reconstructed, "+target+" target;#phi_{PQ} [deg];Counts",180,-180,180);
	TH1F *htrue_fin = new TH1F("htrue_fin","Thrown (Monte Carlo), "+target+" target;#phi_{PQ} [deg];Counts",180,-180,180);
	TH1F *hacc_fin = new TH1F("hacc_fin","Acceptance, "+target+" target;#phi_{PQ} [deg];Counts",180,-180,180);

	TH1F *hreco_CT_fin = new TH1F("hreco_CT_fin","Reconstructed, "+target+" target;#phi_{PQ} [deg];Counts",180,-180,180);
	TH1F *htrue_CT_fin = new TH1F("htrue_CT_fin","Thrown (Monte Carlo), "+target+" target;#phi_{PQ} [deg];Counts",180,-180,180);
	TH1F *hcorr_CT_fin = new TH1F("hcorr_CT_fin","Corrected, "+target+" target;#phi_{PQ} [deg];Counts",180,-180,180);

	TH1F *hratio = new TH1F("hratio","Closure Test Ratio;#phi_{PQ} [deg];Corr/True",180,-180,180);
	TH1F *hratio_rebin = new TH1F("hratio_rebin","Closure Test Ratio;#phi_{PQ} [deg];Corr/True",10,-180,180);

	hreco_fin->Sumw2();
	htrue_fin->Sumw2();
	hacc_fin->Sumw2();
	hreco_CT_fin->Sumw2();
	htrue_CT_fin->Sumw2();
	hcorr_CT_fin->Sumw2();
	hratio->Sumw2();
	hratio_rebin->Sumw2();

	std::vector<TH1F*> hreco_Vec;
	std::vector<TH1F*> htrue_Vec;
	std::vector<TH1F*> hacc_Vec;
	std::vector<TH1F*> hreco_CT_Vec;
	std::vector<TH1F*> htrue_CT_Vec;
	std::vector<TH1F*> hcorr_CT_Vec;

	for (int iQ2=0; iQ2<NQ2; iQ2++){
		for (int iXb=0; iXb<NXb; iXb++){
			for (int iZh=0; iZh<NZh; iZh++){
				for (int iPt2=0; iPt2<NPt2; iPt2++){
					int iTot = iPt2 + iZh*NPt2 + iXb*NZh*NPt2 + iQ2*NXb*NZh*NPt2;
					TH1F *hreco_tmp = new TH1F(Form("hreco%d",iTot),Form("Reco %d (%d, %d, %d, %d)",
											   iTot,iQ2,iXb,iZh,iPt2),180,-180,180);
					hreco_Vec.push_back(hreco_tmp);
					TH1F *htrue_tmp = new TH1F(Form("htrue%d",iTot),Form("Thrown %d (%d, %d, %d, %d)",
										       iTot,iQ2,iXb,iZh,iPt2),180,-180,180);
					htrue_Vec.push_back(htrue_tmp);
					TH1F *hacc_tmp = new TH1F(Form("hacc%d",iTot),Form("Acc %d (%d, %d, %d, %d)",
											  iTot,iQ2,iXb,iZh,iPt2),180,-180,180);
					hacc_tmp->Sumw2();
					hacc_Vec.push_back(hacc_tmp);

					TH1F *hreco_CT_tmp = new TH1F(Form("hrecoCT%d",iTot),Form("Reco %d (%d, %d, %d, %d)",
											   iTot,iQ2,iXb,iZh,iPt2),180,-180,180);
					hreco_CT_Vec.push_back(hreco_CT_tmp);
					TH1F *htrue_CT_tmp = new TH1F(Form("htrueCT%d",iTot),Form("Thrown %d (%d, %d, %d, %d)",
										       iTot,iQ2,iXb,iZh,iPt2),180,-180,180);
					htrue_CT_Vec.push_back(htrue_CT_tmp);
					TH1F *hcorr_CT_tmp = new TH1F(Form("hcorrCT%d",iTot),Form("Acc %d (%d, %d, %d, %d)",
											  iTot,iQ2,iXb,iZh,iPt2),180,-180,180);
					hcorr_CT_tmp->Sumw2();
					hcorr_CT_Vec.push_back(hcorr_CT_tmp);
				}
			}
		}
	}

	// Acceptance-tree loop

	// int Nentries = 30; // for testing
	int Nentries = tree->GetEntries();
	int Nmiddle = Nentries/2;

	for (int row=0; row<Nmiddle; row++){
		tree->GetEntry(row);

		bool ecut=false, mc_ecut=false; // hcut=false, mc_hcut=false;

		int binQ2=-1, binXb=-1;
		if (TargType==target_n && Q2>Q2_limits[0] && Q2<Q2_limits[NQ2] && Xb>Xb_limits[0] && Xb<Xb_limits[NXb] &&
			Yb<0.85 && W>2 && vyec>-1.4 && vyec<1.4){
			ecut = true;
			binQ2 = var_position(NQ2, Q2, Q2_limits);
			binXb = var_position(NXb, Xb, Xb_limits);
		}

		int binmc_Q2=-1, binmc_Xb=-1;
		if (mc_TargType==target_n && mc_Q2>Q2_limits[0] && mc_Q2<Q2_limits[NQ2] && mc_Xb>Xb_limits[0] && mc_Xb<Xb_limits[NXb] &&
			mc_Yb<0.85 && mc_W>2){
			mc_ecut = true;
			binmc_Q2 = var_position(NQ2, mc_Q2, Q2_limits);
			binmc_Xb = var_position(NXb, mc_Xb, Xb_limits);
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

				hreco_Vec[binPt2 + binZh*NPt2 + binXb*NZh*NPt2 + binQ2*NXb*NZh*NPt2]->Fill((*PhiPQ)[i]);
			}

			int binmc_Zh=-1, binmc_Pt2=-1;
			if (mc_ecut==true && (*mc_pid)[i]==211 && (*mc_Zh)[i]>Zh_limits[0] && (*mc_Zh)[i]<Zh_limits[NZh] &&
				(*mc_Pt2)[i]>Pt2_limits[0] && (*mc_Pt2)[i]<Pt2_limits[NPt2] && (*mc_PhiPQ)[i]>PhiPQ_limits[0] &&
				(*mc_PhiPQ)[i]<PhiPQ_limits[NPhiPQ]){
				// mc_hcut = true;
				binmc_Zh = var_position(i, NZh, mc_Zh, Zh_limits);
				binmc_Pt2 = var_position(i, NPt2, mc_Pt2, Pt2_limits);

				htrue_Vec[binmc_Pt2 + binmc_Zh*NPt2 + binmc_Xb*NZh*NPt2 + binmc_Q2*NXb*NZh*NPt2]->Fill((*mc_PhiPQ)[i]);
			}
		} // end of hadrons' loop
  	} // end filling loop

	// Calculating acceptance
	for (int i=0; i<NTot; i++){
		hreco_fin->Add(hreco_Vec[i]);
		htrue_fin->Add(htrue_Vec[i]);

		hacc_Vec[i]->Divide(hreco_Vec[i],htrue_Vec[i],1,1,"B");
		delete hreco_Vec[i];
		delete htrue_Vec[i];
	}

	hacc_fin->Divide(hreco_fin,htrue_fin);

	// Correction-tree loop (Closure Test)
	for (int row=Nmiddle; row<Nentries; row++){
		tree->GetEntry(row);

		bool ecut=false, mc_ecut=false; // hcut=false, mc_hcut=false;

		int binQ2=-1, binXb=-1;
		if (TargType==target_n && Q2>Q2_limits[0] && Q2<Q2_limits[NQ2] && Xb>Xb_limits[0] && Xb<Xb_limits[NXb] &&
			Yb<0.85 && W>2 && vyec>-1.4 && vyec<1.4){
			ecut = true;
			binQ2 = var_position(NQ2, Q2, Q2_limits);
			binXb = var_position(NXb, Xb, Xb_limits);
		}

		int binmc_Q2=-1, binmc_Xb=-1;
		if (mc_TargType==target_n && mc_Q2>Q2_limits[0] && mc_Q2<Q2_limits[NQ2] && mc_Xb>Xb_limits[0] && mc_Xb<Xb_limits[NXb] &&
			mc_Yb<0.85 && mc_W>2){
			mc_ecut = true;
			binmc_Q2 = var_position(NQ2, mc_Q2, Q2_limits);
			binmc_Xb = var_position(NXb, mc_Xb, Xb_limits);
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

				hreco_CT_Vec[binPt2 + binZh*NPt2 + binXb*NZh*NPt2 + binQ2*NXb*NZh*NPt2]->Fill((*PhiPQ)[i]);
			}

			int binmc_Zh=-1, binmc_Pt2=-1;
			if (mc_ecut==true && (*mc_pid)[i]==211 && (*mc_Zh)[i]>Zh_limits[0] && (*mc_Zh)[i]<Zh_limits[NZh] &&
				(*mc_Pt2)[i]>Pt2_limits[0] && (*mc_Pt2)[i]<Pt2_limits[NPt2] && (*mc_PhiPQ)[i]>PhiPQ_limits[0] &&
				(*mc_PhiPQ)[i]<PhiPQ_limits[NPhiPQ]){
				// mc_hcut = true;
				binmc_Zh = var_position(i, NZh, mc_Zh, Zh_limits);
				binmc_Pt2 = var_position(i, NPt2, mc_Pt2, Pt2_limits);

				htrue_CT_Vec[binmc_Pt2 + binmc_Zh*NPt2 + binmc_Xb*NZh*NPt2 + binmc_Q2*NXb*NZh*NPt2]->Fill((*mc_PhiPQ)[i]);
			}
		} // end of hadrons' loop
  	} // end filling loop

	// Calculating correction and Closure Test
	for (int i=0; i<NTot; i++){
		hreco_CT_fin->Add(hreco_CT_Vec[i]);
		htrue_CT_fin->Add(htrue_CT_Vec[i]);

		hcorr_CT_Vec[i]->Divide(hreco_CT_Vec[i],hacc_Vec[i]);
		hcorr_CT_fin->Add(hcorr_CT_Vec[i]);

		delete hreco_CT_Vec[i];
		delete htrue_CT_Vec[i];
		// delete hacc_Vec[i];
	}

	hratio->Divide(hcorr_CT_fin,htrue_CT_fin);

	TH1F *hcorr_CT_rebin = dynamic_cast<TH1F*>(hcorr_CT_fin->Rebin(18,"hcorr_CT_rebin")); // Merge 18 bins in one, leaving 10 markers!
    TH1F *htrue_CT_rebin = dynamic_cast<TH1F*>(htrue_CT_fin->Rebin(18,"htrue_CT_rebin"));

    hratio_rebin->Divide(hcorr_CT_rebin,htrue_CT_rebin);

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

int var_position(int i, int Nvar, std::vector<float> *var, float var_limits[]){
	for (int ivar=0; ivar<Nvar; ivar++){
		if ((*var)[i]>=var_limits[ivar] && (*var)[i]<var_limits[ivar+1]){
			return ivar;
		}
	}
	return -1;
}
