/* Analysis made using acceptance calculated as shown in SMorán's thesis */
#include "functions.h"

void anSMoran(){
////// Beginning (Edit according to target name) <target> = {D, C, Fe, Pb}
    std::string target = "Fe"; // Change to calculate acceptance of <target> = {D, C, Fe, Pb}

////// I/O Files
    TFile *indata = TFile::Open(Form("../../data_%s1V2_light.root",target.c_str()),"READ");
    TTree *tree = (TTree*) indata->Get("ntuple_data");
    TFile *inacc = TFile::Open(Form("../acceptance/root-files/AccSM_%s1.root",target.c_str()),"READ");

    TFile *output = TFile::Open(Form("AnalysisSM_%s1.root",target.c_str()),"RECREATE");

////// Limits and sizes of cuts
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

////// Definition of tree-variables
    int target_n = 2;
    if (target == "D") target_n = 1;

    int TargType;
    float Q2, Xb;
    // float Yb, W, vyec;
    std::vector<float> *Zh = 0;
    std::vector<float> *Pt2 = 0;
    std::vector<float> *PhiPQ = 0;
    std::vector<int> *pid = 0;
    std::vector<float> *Nphe = 0;

    tree->SetBranchAddress("TargType",&TargType);
    tree->SetBranchAddress("Q2",&Q2);
    tree->SetBranchAddress("Xb",&Xb);
    tree->SetBranchAddress("Zh",&Zh);
    tree->SetBranchAddress("Pt2",&Pt2);
    tree->SetBranchAddress("PhiPQ",&PhiPQ);
    tree->SetBranchAddress("pid",&pid);
    tree->SetBranchAddress("Nphe",&Nphe);

////// Create histograms
    TH1F *hdata_fin = new TH1F("hdata_fin",Form("Data, %s target;#phi_{PQ} [deg];Counts",target.c_str()),180,-180,180);
    TH1F *hcorr_fin = new TH1F("hcorr_fin",Form("Corrected, %s target;#phi_{PQ} [deg];Counts",target.c_str()),180,-180,180);
    hdata_fin->Sumw2();
    hcorr_fin->Sumw2();

    std::vector<TH1F*> hdata_Vec;
    // std::vector<TH1F*> hweight_Vec;
    std::vector<TH1F*> hcorr_Vec;

    // TF1 *unit = new TF1("unit","1",0,0);
    for (int iQ2=0; iQ2<NQ2; iQ2++){
        for (int iXb=0; iXb<NXb; iXb++){
            for (int iZh=0; iZh<NZh; iZh++){
                for (int iPt2=0; iPt2<NPt2; iPt2++){
                    int iTot = iPt2 + iZh*NPt2 + iXb*NZh*NPt2 + iQ2*NXb*NZh*NPt2;
                    std::string hdata_title = Form("Data %d (%d, %d, %d, %d)",iTot,iQ2,iXb,iZh,iPt2);
                    TH1F *hdata_tmp = new TH1F(Form("hdata%d",iTot),hdata_title.c_str(),180,-180,180);
                    hdata_tmp->Sumw2();
                    hdata_Vec.push_back(hdata_tmp);

                    // char hweight_title[] = Form("Weight %d (%d, %d, %d, %d)",iTot,iQ2,iXb,iZh,iPt2);
                    // TH1F *hweight_tmp = new TH1F(Form("hweight%d",iTot),hweight_title,180,-180,180);
                    // hweight_tmp->Eval(unit);
                    // hweight_Vec.push_back(hweight_tmp);

                    std::string hcorr_title = Form("Corr %d (%d, %d, %d, %d)",iTot,iQ2,iXb,iZh,iPt2);
                    TH1F *hcorr_tmp = new TH1F(Form("hcorr%d",iTot),hcorr_title.c_str(),180,-180,180);
                    hcorr_tmp->Sumw2();
                    hcorr_Vec.push_back(hcorr_tmp);
                }
            }
        }
    }

////// Tree loop
    int Nentries = tree->GetEntries();
    for (int row=0; row<Nentries; row++){
        // if (row==60) break;
        tree->GetEntry(row);

        bool ecut=false;

        int binQ2=-1, binXb=-1;
        if (TargType==target_n && Q2>Q2_limits[0] && Q2<Q2_limits[NQ2] && Xb>Xb_limits[0] && Xb<Xb_limits[NXb]){
            ecut = true;
            binQ2 = var_position(NQ2, Q2, Q2_limits);
            binXb = var_position(NXb, Xb, Xb_limits);
        }

        if (ecut==false) continue;

        int ientries = PhiPQ->size();
        for (int i=0; i<ientries; i++){
            int binZh=-1, binPt2=-1;
            if ((*pid)[i]==211 && (*Zh)[i]>Zh_limits[0] && (*Zh)[i]<Zh_limits[NZh] && (*Pt2)[i]>Pt2_limits[0] &&
                (*Pt2)[i]<Pt2_limits[NPt2] && (*PhiPQ)[i]>PhiPQ_limits[0] && (*PhiPQ)[i]<PhiPQ_limits[NPhiPQ]){
                binZh = var_position(i, NZh, Zh, Zh_limits);
                binPt2 = var_position(i, NPt2, Pt2, Pt2_limits);

                hdata_Vec[binPt2 + binZh*NPt2 + binXb*NZh*NPt2 + binQ2*NXb*NZh*NPt2]->Fill((*PhiPQ)[i]);
            }
        } // end of hadrons' loop
    } // end filling loop

////// Filling final histograms


    for (int i=0; i<NTot; i++){
        hdata_fin->Add(hdata_Vec[i]);

        TH1F *hacc_tmp = (TH1F*) inacc->Get(Form("hacc%d",i));

        // hweight_Vec[i]->Divide(hacc_tmp);
        hcorr_Vec[i]->Divide(hdata_Vec[i],hacc_tmp);

        hcorr_fin->Add(hcorr_Vec[i]);
        delete hacc_tmp;
        if (i%40!=0){
            delete hdata_Vec[i];
            // delete hweight_Vec[i];
            // delete hcorr_Vec[i];
        }
    }
    // delete unit;

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
