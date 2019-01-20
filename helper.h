#ifndef helper_H
#define helper_H
#include <iostream>
#include <eigen3/Eigen/Dense>
#include <string>
#include "/home/vinit/gurobi/gurobi752/linux64/include/gurobi_c++.h"
#include "KDTree/kd_tree.h"
#include "KDTree/type_defs.h"
#include "ICP/ICP.h"

using namespace std;
using namespace Eigen;

struct OptVariables{
	MatrixXd Cb;						// (Ns x Nm)
	vector<Matrix<double,3,3>> lam;		// (num_partitions_SOS2 x 3 x 3)
	Matrix<double,3,3> w;				
	MatrixXd alpha;						// (3 x Ns)
	MatrixXd phi;						// (1 x Ns)
};

MatrixXd bucket_refinement(MatrixXd *, Matrix<double,3,1> *, int *, MatrixXd *, MatrixXd *, MatrixXd *, KDTree *, MatrixXd *, KDTree *, MatrixXd *, int *, MatrixXd *, MatrixXd *, int *, int *, int *);

vector<Matrix<double,3,3>> B_Nsx6_to_3x3(MatrixXd *);

void data_shape(string, int *, int *);

void find_all_opt_variables(OptVariables *, MatrixXd *, MatrixXd *, KDTree *, MatrixXd *, int *, vector<Matrix<double,3,3>> *);

void find_alpha(MatrixXd *, MatrixXd *, MatrixXd *,MatrixXd *, MatrixXd *, vector<Matrix<double,3,3>> *);

void find_Cb(MatrixXd *, MatrixXd *, KDTree *, int);

void find_lam(vector<Matrix<double,3,3>> *, MatrixXd *, int *);

void find_phi(MatrixXd *, MatrixXd *);

void find_w(Matrix<double,3,3> *, vector<Matrix<double,3,3>> *, int *);

void flatten_matrix(MatrixXd *, double *);

void flatten_vector(vector<Matrix<double,3,3>> *, double *);

void flatten_w(Matrix<double,3,3> *, double *);

void linspace(double *,int, int, int);

void printMatrixSize(MatrixXd *);

MatrixXd read_file(string);

void sampleModelPoints(MatrixXd *, MatrixXd *, int *);


#endif
