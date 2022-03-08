/* Difference in correction with a Full Acceptance or Sector dependent Acceptance */

#include "functions.h"

void analysisFullVsSectAcc(TString target = "Fe", TString nfold = "*", TString binning_name = "", bool cluster = false){
	if (nfold!="*" && nfold!="1" && nfold!="2" && nfold!="3"){
        std::cerr << "[ERROR] file doesn't exist" << std::endl;
        return -1;
    }

	if (binning_name==""){
		std::cerr << "[ERROR] write a correct name for the output" << std::endl;
        return -1;
	}

    // I/O Files
	TString SA_path;
    if (nfold=="*") SA_path = "./Corr5d_Sector/Corr5dSector_A"+target+"_"+binning_name+".root";
    else SA_path = "./Corr5d_Sector/Corr5dSector_"+target+nfold+"_"+binning_name+".root";
	TFile *SA_file = TFile::Open(SA_path,"READ");

	TString FA_path;
    if (nfold=="*") FA_path = "./Corr5d_Xf/Corr5dXf_A"+target+"_"+binning_name+".root";
    else FA_path = "./Corr5d_Xf/Corr5dXf_"+target+nfold+"_"+binning_name+".root";
	TFile *FA_file = TFile::Open(FA_path,"READ");

	TString out_path;
    if (nfold=="*") out_path = "Diff_A"+target+"_"+binning_name+".root";
    else out_path = "Diff_"+target+nfold+"_"+binning_name+".root";
	TFile *out_file = TFile::Open(out_path,"RECREATE");

	std::vector<TH1D*> hDiff_Vec;
	std::vector<TH1D*> hRate_Vec;

	int NQ2=3, NNu=3, NZh=5;

	for (int iQ2=0; iQ2<NQ2; iQ2++){
		for (int iNu=0; iNu<NNu; iNu++){
			TH3D *hcorr_SA = (TH3D*) SA_file->Get(Form("hcorrNLP%i%i",iQ2,iNu));
			TH3D *hcorr_FA = (TH3D*) FA_file->Get(Form("hcorrNLP%i%i",iQ2,iNu));

			for (int iZh=0; iZh<NZh; iZh++){
				TH1D *hPPQ_SA = (TH1D*) hcorr_SA->ProjectionZ(Form("hPPQ%i%i_SA%i",iQ2,iNu,iZh),2*iZh+1,2*iZh+2);
				TH1D *hPPQ_FA = (TH1D*) hcorr_FA->ProjectionZ(Form("hPPQ%i%i_FA%i",iQ2,iNu,iZh),2*iZh+1,2*iZh+2);
				// NOTE: Change the y-label range, so it looks pretty when getting them in the final .root file
				// hPPQ_SA->GetYaxis()->SetRangeUser(0,);

				TH1D *hDiff_tmp = new TH1D(Form("hDiff%i%i_%i",iQ2,iNu,iZh),Form("SectorAcc - FullAcc, %i%i, Zh%i",iQ2,iNu,iZh),40,-180.0, 180.0);
				hDiff_tmp->Add(hPPQ_SA, hPPQ_FA, 1, -1);
				hDiff_Vec.push_back(hDiff_tmp);

				TH1D *hRate_tmp = new TH1D(Form("hRate%i%i_%i",iQ2,iNu,iZh),Form("SectorAcc / FullAcc, %i%i, Zh%i",iQ2,iNu,iZh),40,-180.0, 180.0);
				hRate_tmp->Sumw2();
				hRate_tmp->Divide(hPPQ_SA, hPPQ_FA);
				hRate_tmp->GetYaxis()->SetRangeUser(0.4,1.6);
				hRate_Vec.push_back(hRate_tmp);
			}
		}
	}

	out_file->Write();
	out_file->Close();
}
