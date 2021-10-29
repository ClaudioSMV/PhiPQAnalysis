/* Correction of data using a 5-d acceptance with THnSparse */

#include "functions.h"

void corr5d(TString target = "Fe",TString file_n = "*", TString binning_name = "", bool clust = false){

	if (binning_name==""){
		std::cerr << "[ERROR] write a correct extension for the acceptance/output file" << std::endl;
        return -1;
	}

    // I/O Files
    TChain *tree = new TChain("ntuple_data");
	TString file_type = "_light";
	if (clust) file_type = "";
    TString data_file = "../../clas-data/data_"+target+"1"+file_type+".root";
    tree->Add(data_file);

	TString acc_file_title;
	if (file_n!="*") acc_file_title = "../acceptance/Acc5d/Acc5d_"+target+file_n+"_"+binning_name+".root";
	else acc_file_title = "../acceptance/Acc5d/Acc5d_A"+target+"_"+binning_name+".root";
	TFile *acc_file = TFile::Open(acc_file_title,"READ");
	THnSparseD *hacc = (THnSparseD*) acc_file->Get("hacc");
	THnSparseD *hacc_new = (THnSparseD*) acc_file->Get("hacc_new");

	TString out_file;
    if (file_n!="*") out_file = "Corr5d_"+target+file_n+"_"+binning_name+".root";
	else out_file = "Corr5d_A"+target+"_"+binning_name+".root";
    TFile *output = TFile::Open(out_file,"RECREATE");

	// Definition of tree-variables
	int target_n = 2;
	if (target == "D") target_n = 1;

	int TargType;
	float Q2, Nu, Xb;
	float Yb, W, vyec;
	std::vector<float> *Zh = 0;
	std::vector<float> *Pt2 = 0;
	std::vector<float> *PhiPQ = 0;
	std::vector<int> *pid = 0;
	std::vector<float> *Nphe = 0;

	tree->SetBranchStatus("*",0);
	tree->SetBranchStatus("Q2",1);
	tree->SetBranchStatus("Nu",1);
	tree->SetBranchStatus("Xb",1);
	tree->SetBranchStatus("Yb",1);
	tree->SetBranchStatus("W",1);
	tree->SetBranchStatus("vyec",1);
	tree->SetBranchStatus("Zh",1);
	tree->SetBranchStatus("Pt2",1);
	tree->SetBranchStatus("PhiPQ",1);
	tree->SetBranchStatus("pid",1);
	tree->SetBranchStatus("Nphe",1);

	tree->SetBranchAddress("TargType",&TargType);
	tree->SetBranchAddress("Q2",&Q2);
	tree->SetBranchAddress("Nu",&Nu);
	tree->SetBranchAddress("Xb",&Xb);
	tree->SetBranchAddress("Yb",&Yb);
	tree->SetBranchAddress("W",&W);
	tree->SetBranchAddress("vyec",&vyec);
	tree->SetBranchAddress("Zh",&Zh);
	tree->SetBranchAddress("Pt2",&Pt2);
	tree->SetBranchAddress("PhiPQ",&PhiPQ);
	tree->SetBranchAddress("pid",&pid);
	tree->SetBranchAddress("Nphe",&Nphe);

	// Set final bin widths
	float Q2_limits[] = {1.0, 1.3, 1.8, 4.1};
	float Nu_limits[] = {2.2, 3.2, 3.7, 4.2};

	Int_t NQ2  = sizeof(Q2_limits)/sizeof(Q2_limits[0]) - 1;
	Int_t NNu  = sizeof(Nu_limits)/sizeof(Nu_limits[0]) - 1;

	// DIS limits
	float Limits[5][2] = {{1.0, 4.1}, {2.2, 4.2}, {0.0, 1.0}, {0.0, 1.0}, {-180.0, 180.0}}; // Q2, Nu, Zh, Pt2, PhiPQ

	// Create raw and corrected histograms
	std::vector<TH3D*> hrawd_Vec;
	std::vector<TH3D*> hcorr_Vec;
	std::vector<TH3D*> hcorr_new_Vec;

	TString axes_title = "; Z_{h}; P_{t}^{2} [GeV^{2}]; #phi_{PQ} [deg]";

	for (int iQ2=0; iQ2<NQ2; iQ2++){
		for (int iNu=0; iNu<NNu; iNu++){
			TString bin_title = Form("%.1f<Q2<%.1f, %.1f<Nu<%.1f",Q2_limits[iQ2],Q2_limits[iQ2+1],
								     Nu_limits[iNu],Nu_limits[iNu+1]);

			TH3D *hrawd_tmp = new TH3D(Form("hraw%i%i",iQ2,iNu),"Raw, "+bin_title+axes_title,
			                           10,Limits[2][0],Limits[2][1],
									   10,Limits[3][0],Limits[3][1],
									   40,Limits[4][0],Limits[4][1]);
			hrawd_tmp->Sumw2();
			hrawd_Vec.push_back(hrawd_tmp);

			TH3D *hcorr_tmp = new TH3D(Form("hcorr%i%i",iQ2,iNu),"Corr, "+bin_title+axes_title,
			                           10,Limits[2][0],Limits[2][1],
									   10,Limits[3][0],Limits[3][1],
									   40,Limits[4][0],Limits[4][1]);
			hcorr_tmp->Sumw2();
			hcorr_Vec.push_back(hcorr_tmp);

			TH3D *hcorr_new_tmp = new TH3D(Form("hcorrnew%i%i",iQ2,iNu),"Corr New, "+bin_title+axes_title,
			                           10,Limits[2][0],Limits[2][1],
									   10,Limits[3][0],Limits[3][1],
									   40,Limits[4][0],Limits[4][1]);
			hcorr_new_tmp->Sumw2();
			hcorr_new_Vec.push_back(hcorr_new_tmp);
		}
	}

	int Nentries = tree->GetEntries();
	
	float perc = 1.;
	std::cout << "Starting correction loop!" << std::endl;
	// Correction-tree loop (Closure Test)
	for (int row=0; row<Nentries; row++){
		tree->GetEntry(row);
		if (row > (perc/4.)*Nentries){
			std::cout << "\t" << perc*25 << "%" << std::endl;
			perc++;
		}

		bool ecut=false, hcut=false;
		int binQ2, binNu;

		if (TargType==target_n && Q2>Limits[0][0] && Q2<Limits[0][1] &&
			Yb<0.85 && W>2 && vyec>-1.4 && vyec<1.4 && Nu>Limits[1][0] && Nu<Limits[1][1]){
			ecut = true;
			binQ2 = var_position(NQ2, Q2, Q2_limits);
			binNu = var_position(NNu, Nu, Nu_limits);
		}

		if (!ecut) continue;

		float lead_Zh=-999, lead_Pt2=-999, lead_PhiPQ=-999;
		int ientries = PhiPQ->size();
		for (int i=0; i<ientries; i++){
			if ((*pid)[i]==211 && (*Zh)[i]>Limits[2][0] && (*Zh)[i]<Limits[2][1] && (*Pt2)[i]>Limits[3][0] &&
				(*Pt2)[i]<Limits[3][1] && (*PhiPQ)[i]>Limits[4][0] && (*PhiPQ)[i]<Limits[4][1]){

				hcut = true;
				if ((*Zh)[i] > lead_Zh){
					lead_Zh = (*Zh)[i];
					lead_Pt2 = (*Pt2)[i];
					lead_PhiPQ = (*PhiPQ)[i];
				}
			}
		} // end of hadrons' loop

		if (ecut && hcut){
			Double_t accbin[] = {Q2, Nu, lead_Zh, lead_Pt2, lead_PhiPQ};
			int corbin = binNu + binQ2*NNu;

			hrawd_Vec[corbin]->Fill(lead_Zh,lead_Pt2,lead_PhiPQ);
			
			Int_t bin = hacc->GetBin(accbin);
			Double_t acc_value = hacc->GetBinContent(bin);
			if (acc_value!=0) hcorr_Vec[corbin]->Fill(lead_Zh,lead_Pt2,lead_PhiPQ,1./acc_value);

			Int_t bin_new = hacc_new->GetBin(accbin);
			Double_t acc_new_value = hacc_new->GetBinContent(bin_new);
			if (acc_new_value!=0) hcorr_new_Vec[corbin]->Fill(lead_Zh,lead_Pt2,lead_PhiPQ,1./acc_new_value);
		}
  	} // end filling loop
	  
	std::cout << "\t100%" << " Loop finished!" << std::endl;
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
