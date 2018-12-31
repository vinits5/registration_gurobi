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
	MatrixXf Cb;						// (Ns x Nm)
	vector<Matrix<float,3,3>> lam;		// (num_partitions_SOS2 x 3 x 3)
	Matrix<float,3,3> w;				
	MatrixXf alpha;						// (3 x Ns)
	MatrixXf phi;						// (1 x Ns)
};

MatrixXf bucket_refinement(MatrixXf *, Matrix<float,3,1> *, int *, MatrixXf *, MatrixXf *, MatrixXf *, KDTree *, MatrixXf *, KDTree *, MatrixXf *, int *, MatrixXf *, MatrixXf *, int *, int *, int *);

vector<Matrix<float,3,3>> B_Nsx6_to_3x3(MatrixXf *);

void data_shape(string, int *, int *);

void find_all_opt_variables(OptVariables *, MatrixXf *, MatrixXf *, KDTree *, MatrixXf *, int *, vector<Matrix<float,3,3>> *);

void find_alpha(MatrixXf *, MatrixXf *, MatrixXf *,MatrixXf *, MatrixXf *, vector<Matrix<float,3,3>> *);

void find_Cb(MatrixXf *, MatrixXf *, KDTree *, int);

void find_lam(vector<Matrix<float,3,3>> *, MatrixXf *, int *);

void find_phi(MatrixXf *, MatrixXf *);

void find_w(Matrix<float,3,3> *, vector<Matrix<float,3,3>> *, int *);

void flatten_matrix(MatrixXf *, double *);

void flatten_vector(vector<Matrix<float,3,3>> *, double *);

void flatten_w(Matrix<float,3,3> *, double *);

void linspace(float *,int, int, int);

void printMatrixSize(MatrixXf *);

MatrixXf read_file(string);

void sampleModelPoints(MatrixXf *, MatrixXf *, int *);


#endif
