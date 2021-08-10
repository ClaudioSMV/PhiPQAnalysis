
#include <stdio.h>

int var_position(int Nvar, float var, float var_limits[]);
int var_position(int Nvar, float var, std::vector<float> var_limits);

int var_position(int i, int Nvar, std::vector<float> *var, float var_limits[]);
int var_position(int i, int Nvar, std::vector<float> *var, std::vector<float> var_limits);

void correct(float var,TH1F *hacc,TH1F *hcorr);
void correctCSMV(float var_to_corr, float var_acc,TH1F *hacc,TH1F *hcorr);
