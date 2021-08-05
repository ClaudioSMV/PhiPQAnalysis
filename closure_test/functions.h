
#include <stdio.h>

int var_position(int Nvar, float var, float var_limits[]);
int var_position(int Nvar, float var, std::vector<float> var_limits);

int var_position(int i, int Nvar, std::vector<float> *var, float var_limits[]);
int var_position(int i, int Nvar, std::vector<float> *var, std::vector<float> var_limits);
