
#include <stdio.h>

int var_position(int Nvar, float var, float var_limits[]);
int var_position(int Nvar, float var, std::vector<float> var_limits);

int var_position(int i, int Nvar, std::vector<float> *var, float var_limits[]);
int var_position(int i, int Nvar, std::vector<float> *var, std::vector<float> var_limits);

// void correct(float var,TH1F *hacc,TH1F *hcorr);
void correct(float var_to_corr, float var_acc,TH1F *hacc,TH1F *hcorr);

void setBinVarSize(THnSparse *hist, Int_t nbins[], Double_t Q2_limits[], Double_t Nu_limits[],
					Double_t Zh_limits[], Double_t Pt2_limits[], Double_t PhiPQ_limits[]);
