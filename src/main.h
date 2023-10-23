#pragma once
#include <vector>
#include <system.h>

using namespace std;

int load_systems(vector<System>&, Parameters&, vector<string>&);
void apply_parameters(vector<System>&, Parameters&);
double calc_error(vector<System>&,Parameters&);
void clear_forces(System&);
void calc_forces(vector<System>&, Parameters&);
double train(vector<System>&, Parameters&, bool);
void print_parameters(Parameters&);
int main(int, char*[]);

