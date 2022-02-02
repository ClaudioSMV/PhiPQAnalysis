/* 5-d acceptance with THnSparse-SectorDependent, NLP (No Leading Pion (higher Zh)), good matching, and Xf cut */

#include "functions.h"

void acc5d_Sector(TString target = "Fe", TString nfold = "*", TString binning_name = ""){
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
    if (nfold=="*") out_file = "./Acc5d_Sector/Acc5dSector_A"+target+"_"+binning_name+".root";
    else out_file = "./Acc5d_Sector/Acc5dSector_"+target+nfold+"_"+binning_name+".root";
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
	std::vector<int> *Sector = 0; // 0-5
	std::vector<int> *mc_Sector = 0;
	std::vector<float> *Xf = 0;
	std::vector<float> *mc_Xf = 0;

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
	tree->SetBranchStatus("Sector",1);
	tree->SetBranchStatus("mc_Sector",1);
	tree->SetBranchStatus("Xf",1);
	tree->SetBranchStatus("mc_Xf",1);

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
	tree->SetBranchAddress("Sector",&Sector);
	tree->SetBranchAddress("mc_Sector",&mc_Sector);
	tree->SetBranchAddress("Xf",&Xf);
	tree->SetBranchAddress("mc_Xf",&mc_Xf);

	// Create acceptance histograms
	const Int_t Ndim = 5;
	// OR : Original: {3, 3, 5, 5, 12} = 2700
	// CP : PhiPQ central peak: {3, 3, 5, 5, 40} = 9000 // PhiPQ binning is really important due to the features seen!
	Int_t nbins[Ndim] = {3, 3, 5, 5, 12};
	if (binning_name=="CP") nbins[4] = 40;
	Double_t minbins[Ndim] = {1.0, 2.2, 0.0, 0.0, -180.0};
	Double_t maxbins[Ndim] = {4.1, 4.2, 1.0, 1.0, 180.0};

	// Set variable width bins
	Double_t Q2_limits[] = {1.0, 1.3, 1.8, 4.1};
	Double_t Nu_limits[] = {2.2, 3.2, 3.7, 4.2};
	Double_t Zh_limits[] = {0.0, 0.15, 0.25, 0.4, 0.7, 1.0};
	Double_t Pt2_limits[] = {0.0, 0.03, 0.06, 0.1, 0.18, 1.0};
    Double_t PhiPQ_limits[nbins[4]+1];
	for (int i=0; i<nbins[4]+1; i++){
        PhiPQ_limits[i] = -180. + i*360./nbins[4];
    }

	std::vector<THnSparse*> hreco_Vec;
	std::vector<THnSparse*> htrue_Vec;
	std::vector<THnSparse*> hacc_Vec;

	for (int j=0; j<6; j++){
		THnSparse *hreco_tmp = new THnSparseD(Form("hreco%i",j),Form("Reco Sector %i",j),Ndim,nbins,minbins,maxbins);
		THnSparse *htrue_tmp = new THnSparseD(Form("htrue%i",j),Form("True Sector %i",j),Ndim,nbins,minbins,maxbins);
		THnSparse *hacc_tmp  = new THnSparseD(Form("hacc%i",j), Form("Acc Sector %i",j), Ndim,nbins,minbins,maxbins);
		setBinVarSize(hreco_tmp,nbins,Q2_limits,Nu_limits,Zh_limits,Pt2_limits,PhiPQ_limits);
		setBinVarSize(htrue_tmp,nbins,Q2_limits,Nu_limits,Zh_limits,Pt2_limits,PhiPQ_limits);
		setBinVarSize(hacc_tmp,nbins,Q2_limits,Nu_limits,Zh_limits,Pt2_limits,PhiPQ_limits);
		hreco_tmp->Sumw2();
		htrue_tmp->Sumw2();
		hacc_tmp->Sumw2();
		hreco_Vec.push_back(hreco_tmp);
		htrue_Vec.push_back(htrue_tmp);
		hacc_Vec.push_back(hacc_tmp);
	}
	

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
	std::cout << "Nentries: " << Nentries << std::endl;
	std::cout << "Calculating acceptance..." << std::endl;

	std::vector<int> count_all_match_pi = {0,0,0,0,0,0}, count_reco_pi = {0,0,0,0,0,0};
	int bad_Sector_reco = 0;

	float perc = 1.;
	for (int row=0; row<Nentries; row++){
		tree->GetEntry(row);
		if (row > (perc/4.)*Nentries){
			std::cout << "\t" << perc*25 << "%" << std::endl;
			perc++;
		}

		bool good_mc_ecut=false, good_reco_ecut=false;

		if (mc_TargType==target_n && mc_Q2>Q2_limits[0] && mc_Q2<Q2_limits[NQ2] &&
			mc_Yb<0.85 && mc_W>2 && mc_Nu>Nu_limits[0] && mc_Nu<Nu_limits[NNu]){
			good_mc_ecut = true;
		}

		if (!good_mc_ecut) continue;

		if (TargType==target_n && Q2>Q2_limits[0] && Q2<Q2_limits[NQ2] &&
			Yb<0.85 && W>2 && vyec>-1.4 && vyec<1.4 && Nu>Nu_limits[0] && Nu<Nu_limits[NNu]){
			good_reco_ecut = true;
		}

		int ientries = PhiPQ->size();
		for (int i=0; i<ientries; i++){
			if ((*mc_pid)[i]!=211) continue;

			bool good_mc_hcut = false;

			if ((*mc_Xf)[i]>0 && (*mc_Zh)[i]>Zh_limits[0] && (*mc_Zh)[i]<Zh_limits[NZh] && (*mc_Pt2)[i]>Pt2_limits[0] &&
				(*mc_Pt2)[i]<Pt2_limits[NPt2] && (*mc_PhiPQ)[i]>PhiPQ_limits[0] &&
				(*mc_PhiPQ)[i]<PhiPQ_limits[NPhiPQ] && (*mc_Sector)[i]!=-9999){
				good_mc_hcut = true;
				Double_t new_mc_entry[] = {mc_Q2, mc_Nu, (*mc_Zh)[i], (*mc_Pt2)[i], (*mc_PhiPQ)[i]};
				htrue_Vec[(*mc_Sector)[i]]->Fill(new_mc_entry);
				count_all_match_pi[(*mc_Sector)[i]]++;
			}

			if (!good_mc_hcut) continue;
			if (good_reco_ecut && (*pid)[i]==211 && (*Xf)[i]>0 && (*Zh)[i]>Zh_limits[0] && (*Zh)[i]<Zh_limits[NZh] &&
			(*Pt2)[i]>Pt2_limits[0] && (*Pt2)[i]<Pt2_limits[NPt2] && (*PhiPQ)[i]>PhiPQ_limits[0] &&
			(*PhiPQ)[i]<PhiPQ_limits[NPhiPQ] && (*Sector)[i]!=-9999){
				// I'm only filling with generated values, so the unique use of the reconstructed variables
				// is to check if the particle is well reconstructed, but it doesn't matter the actual value
				Double_t new_entry[] = {mc_Q2, mc_Nu, (*mc_Zh)[i], (*mc_Pt2)[i], (*mc_PhiPQ)[i]};
				hreco_Vec[(*mc_Sector)[i]]->Fill(new_entry);
				if ((*mc_Sector)[i]!=(*Sector)[i]){
					std::cout << "Sector != mc_Sector in row " << row << ", entry " << i << std::endl;
					bad_Sector_reco++;
				}
				count_reco_pi[(*mc_Sector)[i]]++;
			}
		} // end of hadrons' loop
  	} // end filling loop

	std::cout << "\t100%" << " Acceptance calculated!\n" << std::endl;
	std::cout << "Number of well reconstructed pions per Sector: " << std::endl;
	for (int j=0; j<6; j++){
		std::cout << Form("\tS%i ",j) << count_reco_pi[j] << " out of " << count_all_match_pi[j];
		std::cout << " (" << Form("%.3f",100*(float)count_reco_pi[j]/count_all_match_pi[j]) << "%)" << std::endl;
	}

	std::cout << "\nTotal bins: " << NTotBins << std::endl;
	std::cout << "Filled bins per Sector:" << std::endl;
	for (int j=0; j<6; j++){
		std::cout << Form("\tReconstructed \tS%i ",j) << hreco_Vec[j]->GetNbins() << " (" << Form("%.3f",100.*hreco_Vec[j]->GetNbins()/NTotBins) << "%)" << std::endl;
		std::cout << Form("\tThrown \t\tS%i ",j) << htrue_Vec[j]->GetNbins() << " (" << Form("%.3f",100.*htrue_Vec[j]->GetNbins()/NTotBins) << "%)" << std::endl;
	}

	// Calculating acceptance (5-dimensional)
	for (int j=0; j<6; j++){
		hacc_Vec[j]->Divide(hreco_Vec[j],htrue_Vec[j],1,1,"B");
		delete hreco_Vec[j];
		delete htrue_Vec[j];
		hacc_Vec[j]->Write();
	}

	std::cout << "\nSaving and closing." << std::endl;

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