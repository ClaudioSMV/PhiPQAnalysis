/* 5-d acceptance with THnSparse, taking leading pion (higher Zh) and non-detector-biased filling */

#include "functions.h"

void acc5d(TString target = "Fe", TString nfold = "*", TString binning_name = ""){
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
    if (nfold=="*") out_file = "./Acc5d/Acc5d_A"+target+"_"+binning_name+".root";
    else out_file = "./Acc5d/Acc5d_"+target+nfold+"_"+binning_name+".root";
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

	int NTot = NQ2*NNu*NZh*NPt2*NPhiPQ; // Total number of acceptance-cut-bins

	// Acceptance-tree loop

	// int Nentries = 200; // for testing
	int Nentries = tree->GetEntries();
	std::cout << "Nentries: " << Nentries << std::endl;

	int n_lead_pi = 0, n_pi = 0;
	int n_lead_pi_new = 0, n_pi_new = 0;

	float perc = 1.;
	for (int row=0; row<Nentries; row++){
		tree->GetEntry(row);
		if (row > (perc/4.)*Nentries){
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
				n_pi++;
				if (match_Zh!=-999) n_pi_new++;
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
			n_lead_pi++;

			if (match_Zh!=-999){
				Double_t new_entry_2[] = {mc_Q2, mc_Nu, match_Zh, match_Pt2, match_PhiPQ};
				hreco_new->Fill(new_entry_2);
				n_lead_pi_new++;
			}
		}
		if (mc_ecut==true && mc_hcut == true){
			Double_t new_mc_entry[] = {mc_Q2, mc_Nu, lead_mc_Zh, lead_mc_Pt2, lead_mc_PhiPQ};
			htrue->Fill(new_mc_entry);
		}
  	} // end filling loop

	std::cout << "\t100%" << " Loop finished!" << std::endl;
	std::cout << "Number of pions for acceptance: " << n_lead_pi << " out of " << n_pi;
	std::cout << " (" << Form("%.1f",100*(float)n_lead_pi/n_pi) << "%)" << std::endl;
	std::cout << "Number of pions for acceptance_new: " << n_lead_pi_new << " out of " << n_pi_new;
	std::cout << " (" << Form("%.1f",100*(float)n_lead_pi_new/n_pi_new) << "%)" << std::endl;

	// Calculating acceptance (5-dimensional)
	hacc->Divide(hreco,htrue,1,1,"B");
	hacc_new->Divide(hreco_new,htrue,1,1,"B");
	delete hreco;
	delete hreco_new;
	delete htrue;

	hacc->Write();
	hacc_new->Write();
	
	std::cout << "DONE!" << std::endl;
	// Finishing (Save and close files)
	output->Write();
	output->Close();
} // End of the macro

////// DEFINITION OF FUNCTIONS

void setBinVarSize(THnSparse *hist, Int_t nbins[], Double_t Q2_limits[], Double_t Nu_limits[],
					Double_t Zh_limits[], Double_t Pt2_limits[], Double_t PhiPQ_limits[]){
	hist->GetAxis(0)->Set(nbins[0], Q2_limits);
	hist->GetAxis(1)->Set(nbins[1], Nu_limits);
	hist->GetAxis(2)->Set(nbins[2], Zh_limits);
	hist->GetAxis(3)->Set(nbins[3], Pt2_limits);
	hist->GetAxis(4)->Set(nbins[4], PhiPQ_limits);
	return 0;
}