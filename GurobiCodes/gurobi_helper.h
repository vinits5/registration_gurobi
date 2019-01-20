#ifndef gurobi_helper_H
#define gurobi_helper_H
#include <iostream>
#include <string>
#include "/home/vinit/gurobi/gurobi752/linux64/include/gurobi_c++.h"
#include <eigen3/Eigen/Dense>
#include "../helper.h"

using namespace Eigen;
using namespace std;

void define_AlphaM2SConstr(GRBModel *, vector<Matrix<double,3,3>> *, MatrixXd *, GRBVar *, GRBVar *, GRBVar *, GRBVar *, MatrixXd *,int *, int *);

void define_CorrespondenceConstr(GRBModel *, GRBVar *, int *, int *);

void define_GRBVar(GRBModel *, GRBVar *, int, int, long double, long double, char, string);

void define_phiConstr(GRBModel *, GRBVar *, GRBVar *, int *);

void define_RConstr(GRBModel *, GRBVar *, MatrixXd *);

void define_SOSConstr(GRBModel *, GRBVar *, GRBVar *, GRBVar *, double *, int *);

void define_SOS2Constr(GRBModel *, GRBVar *);

void define_TConstr(GRBModel *, GRBVar *, MatrixXd *);

void define_WConstr(GRBModel *, GRBVar *);

void define_3D_GRBVar(GRBModel *, GRBVar *, int, int, int, long double, long double, char, string);

void provide_initialSol(GRBVar *, GRBVar *, GRBVar *, GRBVar *, GRBVar *, GRBVar *, GRBVar *, MatrixXd *, OptVariables *, int *, int *, int *);

#endif