#ifndef ICP_H
#define ICP_H
#include <eigen3/Eigen/Dense>
#include "../KDTree/kd_tree.h"
using namespace Eigen;

struct Transformation{
	Matrix<double,3,3> rotation_matrix;
	Matrix<double,3,1> translation;
};

struct ICP_Result{
	Matrix<double,4,4> g;
	// MatrixXd p;
	double error;
	// MatrixXd target;
};

MatrixXd icp_test(MatrixXd *, MatrixXd *, MatrixXd *, KDTree *, MatrixXd *, KDTree *, MatrixXd *, 
	MatrixXd *, int *, int *, int *);

class ICP{
public:
	MatrixXd *p0_ICP,*p1_ICP,*SigmaS_ICP;
	KDTree *tree_M_sampled_ICP;

	ICP(MatrixXd*, MatrixXd*, KDTree*, MatrixXd*);

	ICP_Result compute(int);

	Transformation find_rigid_transform(MatrixXd*, MatrixXd*);
	Matrix<double,4,4> icp_Rt_to_matrix(Transformation);
};
#endif