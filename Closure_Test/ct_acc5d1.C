/* Closure Test of 5-d acceptance with THnSparse, taking leading pion (higher Zh) and non-detector-biased filling */

#include "functions.h"

void ct_acc5d1(TString target = "Fe", TString nfold = "*", TString binning_name = ""){
	if (nfold!="*" && nfold!="1" && nfold!="2" && nfold!="3"){
        std::cerr << "[ERROR] file doesn't exist" << std::endl;
        return -1;
    }

	if (binning_name==""){
		std::cerr << "[ERROR] write a correct name for the output" << std::endl;
        return -1;
	}

    // I/O Files
    TChain *tree = new TChain("ntuple_sim");
    TString in_file = "../../clas-HSim/hsim_"+target+nfold+".root";
    tree->Add(in_file);

    TString out_file;
    if (nfold=="*") out_file = "./CT_Acc5d1_A"+target+"_"+binning_name+".root";
    else out_file = "./CT_Acc5d1_"+target+nfold+"_"+binning_name+".root";
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

	tree->SetBranchStatus("*",0);
	tree->SetBranchStatus("mc_TargType",1);
	tree->SetBranchStatus("Q2",1);
	tree->SetBranchStatus("mc_Q2",1);
	tree->SetBranchStatus("Nu",1);
    tree->SetBranchStatus("mc_Nu",1);
	tree->SetBranchStatus("Xb",1);
	tree->SetBranchStatus("mc_Xb",1);
	tree->SetBranchStatus("Yb",1);
	tree->SetBranchStatus("mc_Yb",1);
	tree->SetBranchStatus("W",1);
	tree->SetBranchStatus("mc_W",1);
	tree->SetBranchStatus("vyec",1);
	tree->SetBranchStatus("Zh",1);
	tree->SetBranchStatus("mc_Zh",1);
	tree->SetBranchStatus("Pt2",1);
	tree->SetBranchStatus("mc_Pt2",1);
	tree->SetBranchStatus("PhiPQ",1);
	tree->SetBranchStatus("mc_PhiPQ",1);
	tree->SetBranchStatus("pid",1);
	tree->SetBranchStatus("mc_pid",1);
	tree->SetBranchStatus("Nphe",1);

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

	// Create acceptance histograms
	const Int_t Ndim = 5;
	// OR : Original: {3, 3, 5, 5, 12} = 2700
	// CP : PhiPQ central peak: {3, 3, 5, 5, 40} = 9000
	// Int_t nbins[Ndim] = {3, 3, 5, 5, 12};
	Int_t nbins[Ndim] = {3, 3, 5, 5, 40}; // PhiPQ binning is really important due to the features seen!
	Double_t minbins[Ndim] = {1.0, 2.2, 0.0, 0.0, -180.0};
	Double_t maxbins[Ndim] = {4.1, 4.2, 1.0, 1.0, 180.0};

	// 'Classic' fills the numerator of acceptance with the values of the reconstructed particle, however, 'New' fills the
	// numerator with the values of the thrown particle that matches with the reconstructed one.
	THnSparse *hreco = new THnSparseD("hreco","Reco Classic",Ndim,nbins,minbins,maxbins);
	THnSparse *hreco_new = new THnSparseD("hreco_new","Reco New",Ndim,nbins,minbins,maxbins);
	THnSparse *htrue = new THnSparseD("htrue","True",Ndim,nbins,minbins,maxbins);
	THnSparse *hacc = new THnSparseD("hacc","Acc Classic",Ndim,nbins,minbins,maxbins);
	THnSparse *hacc_new = new THnSparseD("hacc_new","Acc New",Ndim,nbins,minbins,maxbins);

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
	setBinVarSize(hreco_new,nbins,Q2_limits,Nu_limits,Zh_limits,Pt2_limits,PhiPQ_limits);
	setBinVarSize(htrue,nbins,Q2_limits,Nu_limits,Zh_limits,Pt2_limits,PhiPQ_limits);
	setBinVarSize(hacc,nbins,Q2_limits,Nu_limits,Zh_limits,Pt2_limits,PhiPQ_limits);
	setBinVarSize(hacc_new,nbins,Q2_limits,Nu_limits,Zh_limits,Pt2_limits,PhiPQ_limits);
	hreco->Sumw2();
	hreco_new->Sumw2();
	htrue->Sumw2();
	hacc->Sumw2();
	hacc_new->Sumw2();

	// N<var> is the number of bins (array's elements - 1), N<var>+1 is the number of limits (array's elements)
	int NQ2  = sizeof(Q2_limits)/sizeof(Q2_limits[0]) - 1;
	int NNu  = sizeof(Nu_limits)/sizeof(Nu_limits[0]) - 1;
	int NZh  = sizeof(Zh_limits)/sizeof(Zh_limits[0]) - 1;
	int NPt2 = sizeof(Pt2_limits)/sizeof(Pt2_limits[0]) - 1;
	int NPhiPQ  = sizeof(PhiPQ_limits)/sizeof(PhiPQ_limits[0]) - 1;

	int NTotBins = NQ2*NNu*NZh*NPt2*NPhiPQ; // Total number of acceptance-cut-bins

	// Acceptance-tree loop

	// int Nentries = 200; // for testing
	int Nentries = tree->GetEntries();
	int Nmiddle = Nentries/2;
	std::cout << "Nentries: " << Nentries << ", Nmiddle: " << Nmiddle << std::endl;
	std::cout << "Calculating acceptance..." << std::endl;

	float perc = 1.;
	for (int row=0; row<Nmiddle; row++){
		tree->GetEntry(row);
		if (row > (perc/4.)*Nmiddle){
			std::cout << "\t" << perc*25 << "%" << std::endl;
			perc++;
		}

		bool ecut=false, mc_ecut=false, hcut=false, mc_hcut=false;

		if (TargType==target_n && Q2>Q2_limits[0] && Q2<Q2_limits[NQ2] &&
			Yb<0.85 && W>2 && vyec>-1.4 && vyec<1.4 && Nu>Nu_limits[0] && Nu<Nu_limits[NNu]){
			ecut = true;
		}

		if (mc_TargType==target_n && mc_Q2>Q2_limits[0] && mc_Q2<Q2_limits[NQ2] &&
			mc_Yb<0.85 && mc_W>2 && mc_Nu>Nu_limits[0] && mc_Nu<Nu_limits[NNu]){
			mc_ecut = true;
		}

		if (ecut==false && mc_ecut==false) continue;

		float lead_Zh=-999, lead_Pt2=-999, lead_PhiPQ=-999;
		float match_Zh=-999, match_Pt2=-999, match_PhiPQ=-999;
		float lead_mc_Zh=-999, lead_mc_Pt2=-999, lead_mc_PhiPQ=-999;
		int ientries = PhiPQ->size();
		for (int i=0; i<ientries; i++){
			if ((*pid)[i]==211 && (*Zh)[i]>Zh_limits[0] && (*Zh)[i]<Zh_limits[NZh] &&
			    (*Pt2)[i]>Pt2_limits[0] && (*Pt2)[i]<Pt2_limits[NPt2] && (*PhiPQ)[i]>PhiPQ_limits[0] &&
				(*PhiPQ)[i]<PhiPQ_limits[NPhiPQ]){
				hcut = true;
				if ((*Zh)[i] > lead_Zh){
					lead_Zh = (*Zh)[i];
					lead_Pt2 = (*Pt2)[i];
					lead_PhiPQ = (*PhiPQ)[i];
					if ((*mc_Zh)[i]>0){
						match_Zh = (*mc_Zh)[i];
						match_Pt2 = (*mc_Pt2)[i];
						match_PhiPQ = (*mc_PhiPQ)[i];
					}
				}
			}

			if ((*mc_pid)[i]==211 && (*mc_Zh)[i]>Zh_limits[0] && (*mc_Zh)[i]<Zh_limits[NZh] &&
				(*mc_Pt2)[i]>Pt2_limits[0] && (*mc_Pt2)[i]<Pt2_limits[NPt2] && (*mc_PhiPQ)[i]>PhiPQ_limits[0] &&
				(*mc_PhiPQ)[i]<PhiPQ_limits[NPhiPQ]){
				mc_hcut = true;
				if ((*mc_Zh)[i] > lead_mc_Zh){
					lead_mc_Zh = (*mc_Zh)[i];
					lead_mc_Pt2 = (*mc_Pt2)[i];
					lead_mc_PhiPQ = (*mc_PhiPQ)[i];
				}
			}
		} // end of hadrons' loop

		if (ecut==true && hcut==true){
			Double_t new_entry[] = {Q2, Nu, lead_Zh, lead_Pt2, lead_PhiPQ};
			hreco->Fill(new_entry);

			if (match_Zh!=-999){
				Double_t new_entry_2[] = {mc_Q2, mc_Nu, match_Zh, match_Pt2, match_PhiPQ};
				hreco_new->Fill(new_entry_2);
			}
		}
		if (mc_ecut==true && mc_hcut == true){
			Double_t new_mc_entry[] = {mc_Q2, mc_Nu, lead_mc_Zh, lead_mc_Pt2, lead_mc_PhiPQ};
			htrue->Fill(new_mc_entry);
		}
  	} // end filling loop

	std::cout << "\t100%" << " Acceptance calculated!\n" << std::endl;

	std::cout << "Total bins: " << NTotBins << std::endl;
	std::cout << "Non empty bins:" << std::endl;
	std::cout << "\tThrown - " << htrue->GetNbins() << " (" << 100.*htrue->GetNbins()/NTotBins << "%)" << std::endl;
	std::cout << "\tReconstructed - " << hreco->GetNbins() << " (" << 100.*hreco->GetNbins()/NTotBins << "%)" << std::endl;
	std::cout << "\tReconstructed New - " << hreco_new->GetNbins() << " (" << 100.*hreco_new->GetNbins()/NTotBins << "%)\n" << std::endl;

	// Calculating acceptance (5-dimensional)
	hacc->Divide(hreco,htrue,1,1,"B");
	hacc_new->Divide(hreco_new,htrue,1,1,"B");
	delete hreco;
	delete hreco_new;
	delete htrue;

	hacc->Write();
	hacc_new->Write();

	// Correction of data
	// Set final bin widths
	float Q2_Flimits[] = {1.0, 1.3, 1.8, 4.1};
	float Nu_Flimits[] = {2.2, 3.2, 3.7, 4.2};

	Int_t NQ2F  = sizeof(Q2_Flimits)/sizeof(Q2_Flimits[0]) - 1;
	Int_t NNuF  = sizeof(Nu_Flimits)/sizeof(Nu_Flimits[0]) - 1;

	// DIS limits
	// Q2_limits[0], Q2_limits[NQ2], PhiPQ_limits[0], PhiPQ_limits[NPhiPQ],... // Q2, Nu, Zh, Pt2, PhiPQ

	// Create raw and corrected histograms
	std::vector<TH3D*> hrawd_Vec;
	std::vector<TH3D*> hcorr_Vec;
	std::vector<TH3D*> hcorr_new_Vec;

	TString axes_title = "; Z_{h}; P_{t}^{2} [GeV^{2}]; #phi_{PQ} [deg]";

	for (int iQ2=0; iQ2<NQ2F; iQ2++){
		for (int iNu=0; iNu<NNuF; iNu++){
			TString bin_title = Form("%.1f<Q2<%.1f, %.1f<Nu<%.1f",Q2_Flimits[iQ2],Q2_Flimits[iQ2+1],
								     Nu_Flimits[iNu],Nu_Flimits[iNu+1]);

			TH3D *hrawd_tmp = new TH3D(Form("hraw%i%i",iQ2,iNu),"Raw, "+bin_title+axes_title,
			                           10,Zh_limits[0],Zh_limits[NZh],
									   10,Pt2_limits[0],Pt2_limits[NPt2],
									   40,PhiPQ_limits[0],PhiPQ_limits[NPhiPQ]);
			hrawd_tmp->Sumw2();
			hrawd_Vec.push_back(hrawd_tmp);

			TH3D *hcorr_tmp = new TH3D(Form("hcorr%i%i",iQ2,iNu),"Corr, "+bin_title+axes_title,
			                           10,Zh_limits[0],Zh_limits[NZh],
									   10,Pt2_limits[0],Pt2_limits[NPt2],
									   40,PhiPQ_limits[0],PhiPQ_limits[NPhiPQ]);
			hcorr_tmp->Sumw2();
			hcorr_Vec.push_back(hcorr_tmp);

			TH3D *hcorr_new_tmp = new TH3D(Form("hcorrnew%i%i",iQ2,iNu),"Corr New, "+bin_title+axes_title,
			                           10,Zh_limits[0],Zh_limits[NZh],
									   10,Pt2_limits[0],Pt2_limits[NPt2],
									   40,PhiPQ_limits[0],PhiPQ_limits[NPhiPQ]);
			hcorr_new_tmp->Sumw2();
			hcorr_new_Vec.push_back(hcorr_new_tmp);
		}
	}
	
	perc = 1.;
	std::cout << "Starting correction loop!" << std::endl;
	// Correction-tree loop (Closure Test)
	for (int row=Nmiddle; row<Nentries; row++){
		tree->GetEntry(row);
		if ((row - Nmiddle) > (perc/4.)*(Nentries-Nmiddle)){
			std::cout << "\t" << perc*25 << "%" << std::endl;
			perc++;
		}

		bool ecut=false, hcut=false;
		bool mc_ecut=false, mc_hcut=false;
		int binQ2, binNu;
		int mc_binQ2, mc_binNu;

		if (TargType==target_n && Q2>Q2_Flimits[0] && Q2<Q2_Flimits[NQ2F] &&
			Yb<0.85 && W>2 && vyec>-1.4 && vyec<1.4 && Nu>Nu_Flimits[0] && Nu<Nu_Flimits[NNuF]){
			ecut = true;
			binQ2 = var_position(NQ2F, Q2, Q2_Flimits);
			binNu = var_position(NNuF, Nu, Nu_Flimits);
		}
		if (mc_TargType==target_n && mc_Q2>Q2_Flimits[0] && mc_Q2<Q2_Flimits[NQ2F] &&
			mc_Yb<0.85 && mc_W>2 && mc_Nu>Nu_Flimits[0] && mc_Nu<Nu_Flimits[NNuF]){
			mc_ecut = true;
			mc_binQ2 = var_position(NQ2F, mc_Q2, Q2_Flimits);
			mc_binNu = var_position(NNuF, mc_Nu, Nu_Flimits);
		}

		if (!ecut && !mc_ecut) continue;

		float lead_Zh=-999, lead_Pt2=-999, lead_PhiPQ=-999;
		float lead_mc_Zh=-999, lead_mc_Pt2=-999, lead_mc_PhiPQ=-999;
		int ientries = PhiPQ->size();
		for (int i=0; i<ientries; i++){
			if ((*pid)[i]==211 && (*Zh)[i]>Zh_limits[0] && (*Zh)[i]<Zh_limits[NZh] && (*Pt2)[i]>Pt2_limits[0] &&
				(*Pt2)[i]<Pt2_limits[NPt2] && (*PhiPQ)[i]>PhiPQ_limits[0] && (*PhiPQ)[i]<PhiPQ_limits[NPhiPQ]){

				hcut = true;
				if ((*Zh)[i] > lead_Zh){
					lead_Zh = (*Zh)[i];
					lead_Pt2 = (*Pt2)[i];
					lead_PhiPQ = (*PhiPQ)[i];
				}
			}
			if ((*mc_pid)[i]==211 && (*mc_Zh)[i]>Zh_limits[0] && (*mc_Zh)[i]<Zh_limits[NZh] && (*mc_Pt2)[i]>Pt2_limits[0] &&
				(*mc_Pt2)[i]<Pt2_limits[NPt2] && (*mc_PhiPQ)[i]>PhiPQ_limits[0] && (*mc_PhiPQ)[i]<PhiPQ_limits[NPhiPQ]){

				mc_hcut = true;
				if ((*mc_Zh)[i] > lead_mc_Zh){
					lead_mc_Zh = (*mc_Zh)[i];
					lead_mc_Pt2 = (*mc_Pt2)[i];
					lead_mc_PhiPQ = (*mc_PhiPQ)[i];
				}
			}
		} // end of hadrons' loop

		if (ecut && hcut){
			Double_t accbin[] = {Q2, Nu, lead_Zh, lead_Pt2, lead_PhiPQ};
			int corbin = binNu + binQ2*NNuF;
			
			Int_t bin = hacc->GetBin(accbin);
			Double_t acc_value = hacc->GetBinContent(bin);
			if (acc_value!=0) hcorr_Vec[corbin]->Fill(lead_Zh,lead_Pt2,lead_PhiPQ,1./acc_value);

			Int_t bin_new = hacc_new->GetBin(accbin);
			Double_t acc_new_value = hacc_new->GetBinContent(bin_new);
			if (acc_new_value!=0) hcorr_new_Vec[corbin]->Fill(lead_Zh,lead_Pt2,lead_PhiPQ,1./acc_new_value);
		}
		if (mc_ecut && mc_hcut){
			int mc_bin = mc_binNu + mc_binQ2*NNuF;

			hrawd_Vec[mc_bin]->Fill(lead_mc_Zh,lead_mc_Pt2,lead_mc_PhiPQ);
		}
  	} // end filling loop
	  
	std::cout << "\t100%" << " Correction loop finished!" << std::endl;

	// Closure Test!
	TString xvar[3] = {"Zh", "Pt2", "PhiPQ"};
	TString xvartitle[3] = {"Z_{h}", "P_{t}^{2} [GeV^{2}]", "#phi_{PQ} [deg]"};
	Int_t xvarbin[3] = {10,10,40};
	Double_t xvarlimits[3][2] = {{Zh_limits[0],Zh_limits[NZh]},
								{Pt2_limits[0],Pt2_limits[NPt2]},
								{PhiPQ_limits[0],PhiPQ_limits[NPhiPQ]}};

	// std::vector<TH1D*> ctfactorZh_Vec;
	// std::vector<TH1D*> ctfactorPt2_Vec;
	// std::vector<TH1D*> ctfactorPhiPQ_Vec;

	std::cout << "Applying Closure Test..." << std::endl;
	for (int iQ2=0; iQ2<NQ2F; iQ2++){
		for (int iNu=0; iNu<NNuF; iNu++){
			TString bin_title = Form("%.1f<Q2<%.1f, %.1f<Nu<%.1f",Q2_Flimits[iQ2],Q2_Flimits[iQ2+1],
								     Nu_Flimits[iNu],Nu_Flimits[iNu+1]);
			int ibin = iNu+iQ2*NNuF;
			for (int ivar=0; ivar<3; ivar++){
				TH1D *hct_tmp = new TH1D(Form("hCT_"+xvar[ivar]+"%i%i",iQ2,iNu),
										 "Closure Test, "+bin_title+";"+xvartitle[ivar]+";Corr/True",
			                        	 xvarbin[ivar],xvarlimits[ivar][0],xvarlimits[ivar][1]);
				hct_tmp->Sumw2();
				if (ivar==0) hct_tmp->Divide(hcorr_Vec[ibin]->ProjectionX(),hrawd_Vec[ibin]->ProjectionX());
				else if (ivar==1) hct_tmp->Divide(hcorr_Vec[ibin]->ProjectionY(),hrawd_Vec[ibin]->ProjectionY());
				else if (ivar==2) hct_tmp->Divide(hcorr_Vec[ibin]->ProjectionZ(),hrawd_Vec[ibin]->ProjectionZ());
				// hct_tmp->Write();

				TH1D *hctnew_tmp = new TH1D(Form("hCTnew_"+xvar[ivar]+"%i%i",iQ2,iNu),
										 "Closure Test new, "+bin_title+";"+xvartitle[ivar]+";Corr/True",
			                        	 xvarbin[ivar],xvarlimits[ivar][0],xvarlimits[ivar][1]);
				hctnew_tmp->Sumw2();
				if (ivar==0) hctnew_tmp->Divide(hcorr_new_Vec[ibin]->ProjectionX(),hrawd_Vec[ibin]->ProjectionX());
				else if (ivar==1) hctnew_tmp->Divide(hcorr_new_Vec[ibin]->ProjectionY(),hrawd_Vec[ibin]->ProjectionY());
				else if (ivar==2) hctnew_tmp->Divide(hcorr_new_Vec[ibin]->ProjectionZ(),hrawd_Vec[ibin]->ProjectionZ());
				// hctnew_tmp->Write();
			}
		}
	}

	std::cout << "Saving and closing." << std::endl;

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

void setBinVarSize(THnSparse *hist, Int_t nbins[], Double_t Q2_limits[], Double_t Nu_limits[],
					Double_t Zh_limits[], Double_t Pt2_limits[], Double_t PhiPQ_limits[]){
	hist->GetAxis(0)->Set(nbins[0], Q2_limits);
	hist->GetAxis(1)->Set(nbins[1], Nu_limits);
	hist->GetAxis(2)->Set(nbins[2], Zh_limits);
	hist->GetAxis(3)->Set(nbins[3], Pt2_limits);
	hist->GetAxis(4)->Set(nbins[4], PhiPQ_limits);
	return 0;
}