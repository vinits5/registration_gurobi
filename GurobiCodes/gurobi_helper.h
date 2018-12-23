#ifndef gurobi_helper_H
#define gurobi_helper_H
#include <iostream>
#include <string>
#include "/home/vinit/gurobi/gurobi752/linux64/include/gurobi_c++.h"

using namespace std;

void define_GRBVar(GRBModel *, GRBVar *, int, int, long double, long double, char, string);

void define_3D_GRBVar(GRBModel *, GRBVar *, int, int, int, long double, long double, char, string);

#endif