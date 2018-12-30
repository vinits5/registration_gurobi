#ifndef CALLBACKMTOS_H
#define CALLBACKMTOS_H

#include <iostream>
#include <eigen3/Eigen/Dense>
#include <string>
#include "/home/vinit/gurobi/gurobi752/linux64/include/gurobi_c++.h"
#include "../KDTree/kd_tree.h"
#include "../KDTree/type_defs.h"
#include "../ICP/ICP.h"
#include "../helper.h"
#include "../GurobiCodes/gurobi_helper.h"

using namespace std;
using namespace Eigen;

// Structure for parameters need for callback function of Gurobi.
struct Params{
	int *Ns_sampled;						// Points in Model Data
	MatrixXf *S;							// 3 x Ns
	int *num_sampled_sens_points_ICP;		
	MatrixXf *V;							// 3 x Nv
	MatrixXf *M;							// 3 x Nm
	KDTree *tree_M;							// KD Tree of Model Points
	MatrixXf *M_sampled;					// 3 x Ns_sampled 
	KDTree *tree_M_sampled;					// KD tree of sampled points
	MatrixXf *SigmaS;						// n x 6
	MatrixXf *F;							// Nf x Nm
	vector<Matrix<float,3,3>> *B;
	MatrixXf *M_global;
	int *num_partitions_SOS2;
	int *points_per_face;					// Nf / Nm
	int *ICP_triangle_proj_switch_callback;	
	int *ICP_or_GICP_switch_callback;
};

class callbackMtoS: public GRBCallback{
public:
	// Variables
	int sol_repeat_counter = 0;
	int solflag = 0;
	Params *callback_params;
	GRBModel *model;

	// Methods
	callbackMtoS(GRBModel *, Params *);

protected:
	// Method
	void callback();
};


#endif