#include "/home/vinit/gurobi/gurobi752/linux64/include/gurobi_c++.h"
#include <iostream>
#include <eigen3/Eigen/Dense>
#include "helper.h"
#include "KDTree/kd_tree.h"
#include "GurobiCodes/gurobi_helper.h"
#include "ICP/ICP.h"
#include "Callback/callbackMtoS.h"

using namespace std;
using namespace Eigen;

#define PI 3.14159265358979323846

int main(){
	cout<<"Welcome to Gurobi Optimization: Part 2"<<endl;

	// Define Switches.
	// int outlier_switch = 0;
	int M_sampled_switch_global_optimiser = 1; 			// if 1, uses M_sampled in global optimizer 
	int set_start_sol_switch = 2; 						// 1 sets the initial seed value using ICP/GICP 
	int callback_switch = 1; 							// on or off
	int ICP_or_GICP_switch_callback = 1;				// 1 for ICP, 2 for GICP
	int ICP_triangle_proj_switch_callback = 0;			// on or off
	int ICP_or_GICP_switch_bucket = 1;					// 1 for ICP, 2 for GICP
	int ICP_triangle_proj_switch_bucket = 0;			// on or off
	int alpha_switch = 1;								// for alpha switch with sigma.

	// Define Constants
	int Ns_sampled = 20;								// fewer sensor points for optimization 
	int num_sampled_sens_points_ICP = 500;
	int approx_sampled_model_points = 20;

	// Read the data from text files.
	MatrixXf S = read_file("datasets/dragon_6000_noisy/S.txt"); 			// file is nx3  (S_given)
	MatrixXf St = read_file("datasets/dragon_6000_noisy/St.txt");			// file is nx3 	(St_given)
	MatrixXf V = read_file("datasets/dragon_6000_noisy/M.txt"); 			// file is nx3 	(V_given)
	MatrixXf M = read_file("datasets/dragon_6000_noisy/Msampled.txt"); 		// file is nx3 	(M_given)
	MatrixXf gt = read_file("datasets/dragon_6000_noisy/gt.txt");			// file is 4x4  (gt_given)
	MatrixXf B_given = read_file("datasets/dragon_6000_noisy/B.txt"); 			// file is nx6
	MatrixXf SigmaS = read_file("datasets/dragon_6000_noisy/SigmaS.txt"); // file is nx6 	(SigmaS_given)
	MatrixXf F = read_file("datasets/dragon_6000_noisy/F.txt"); 			// file is NfxNm (F_given)

	// Change Shape of Matrices.
	S.transposeInPlace();			  			// sensor points (3 x Ns)
	St.transposeInPlace();						// sensor points (3 x Ns)
	V.transposeInPlace(); 						// Vertices	(3 x Nv)
	M.transposeInPlace();						// Model points (3 x Nm)
	vector<Matrix<float,3,3>> B = B_Nsx6_to_3x3(&B_given);  		// (Ns x 3 x 3) 	(B_Nsx6 = B_given)

	// Display Matrix Shapes. (function defined in helper)
	printMatrixSize(&S);
	printMatrixSize(&F);
	printMatrixSize(&M);
	printMatrixSize(&V);

	// Important sizes.
	int Ns = S.cols();				// Number of Points in Sensor Data
	int Nm = M.cols();				// Number of Points in Model Data
	// int Nv = V.cols();				// Number of Vertices.
	int Nf = F.rows();				// rows of Faces.
	int num_partitions_SOS2 = 10;	// for SOS2 constraint.
	int points_per_face = Nm/Nf;	// Points per Face.

	// Sample Model Points from the Model Data.
	cout<<"Approx Sampled Model Points: "<<approx_sampled_model_points<<endl;
	int interval_sampled_model_points = Nm/approx_sampled_model_points;
	// MatrixXf M_sampled(M.rows(),(Nm/interval_sampled_model_points)+1);
	MatrixXf M_sampled(M.rows(),(Nm/interval_sampled_model_points));
	sampleModelPoints(&M, &M_sampled, &interval_sampled_model_points);

	int Nm_sampled = M_sampled.cols();			// Points in Sampled Model data.

	// Create KD-Tree of Model Points.
	KDTree tree_M=NULL, tree_M_sampled=NULL;
	for(int i=0; i<M.cols(); i++){					// M (3xNi)
		// 1. Find ith column from M & cast it to Matrix<long double,3,1>
		// 2. Call insert function from kd_tree.h for each point.
		// 3. Arguments for insert-> i] 3D point as Matrix<long double,3,1>
								//	ii] index of the point in point cloud.
								// iii] Memory location of KDTree.
		insert(M.col(i).cast<long double>(),i,&tree_M);
	}
	// Create KD-Tree of Sampled Model Points.
	for(int i=0; i<M_sampled.cols(); i++){					// M_sampled (3xNi)
		insert(M_sampled.col(i).cast<long double>(),i,&tree_M_sampled);
	}

	// ####################### Gurobi Codes #######################
	try{
		// Create Gurobi Environment and Gurobi Model.
		GRBEnv env = GRBEnv();
		GRBModel m = GRBModel(env);

		// Set Parameters.
		m.set(GRB_DoubleParam_TimeLimit,2*3600);
		// m.set("PoolSolutions","100");

		cout<<"Error Before Optimization: "<<endl;					// Error Before Optimization.
		// Find Angle Error.
		Matrix<float,3,3> errorMat = gt.block(0,0,3,3);				// Extract Rotation Matrix from Transformation Matrix.
		MatrixXf ea = errorMat.eulerAngles(2, 1, 0)*(180/PI);		// actual angles [ea(2,0), ea(1,0), ea(0,0)]
		// (Ignored the swap of angles. Used just to calculate angle norm.)
		cout<<"Angle Error: "<<ea.norm()<<endl;

		// Find Position Error.
		Matrix<float,3,1> positionError = gt.block(0,3,3,1);
		cout<<"Position Error: "<<positionError.norm()<<endl;

		// Define Model Data for optimization.
		MatrixXf *M_global;
		KDTree *tree_M_global;
		if(M_sampled_switch_global_optimiser==1){
			M_global = &M_sampled;						// Store the pointer in M_global.
			tree_M_global = &tree_M_sampled;			// Store the pointer of KDTree in global.
		}
		else{
			M_global = &M;								// Store the pointer in M.
			tree_M_global = &tree_M;
		}
		
		int Nm_global = M_global->cols();				// Store the Number of Points in Model data.

		GRBVar alpha[3][Ns_sampled];
		GRBVar phi[1][Ns_sampled];
		GRBVar T[3][1];
		GRBVar R[3][3];
		GRBVar W[3][3];
		GRBVar Cb_sampled[Ns_sampled][Nm_global];
		GRBVar lam[num_partitions_SOS2][3][3];

		// Define Gurobi Variables.
		define_GRBVar(&m, &alpha[0][0], 3, Ns_sampled, 0, GRB_INFINITY, GRB_CONTINUOUS, "alpha");
		define_GRBVar(&m, &phi[0][0], 1, Ns_sampled, 0, GRB_INFINITY, GRB_CONTINUOUS, "phi");
		define_GRBVar(&m, &T[0][0], 3, 1, -GRB_INFINITY, GRB_INFINITY, GRB_CONTINUOUS, "T");
		define_GRBVar(&m, &R[0][0], 3, 3, -1, 1, GRB_CONTINUOUS, "rotationMatrix");
		define_GRBVar(&m, &W[0][0], 3, 3, 0, 1, GRB_CONTINUOUS, "W");
		define_GRBVar(&m, &Cb_sampled[0][0], Ns_sampled, Nm_global, 0, 1, GRB_BINARY, "Cb_sampled");
		define_3D_GRBVar(&m, &lam[0][0][0], num_partitions_SOS2, 3, 3, 0, GRB_INFINITY, GRB_CONTINUOUS, "lambda");

		float q[num_partitions_SOS2];
		linspace(&q[0],-1,1,num_partitions_SOS2);

		// GRBVar o[1][Ns];
		// define_GRBVar(&m, &o[0][0], 1, Ns, 0, 1, GRB_BINARY, "o");
		// int bigM = 5000, phiMax = 500;

		m.update();							// Update the model.

		// Starting to compute minimum objective value.
		MatrixXf gt_inv = gt.inverse();		// Inverse of ground truth transformation. (Transformation model to sensor)

		OptVariables opt_vars;				// Store Opt Variables.
		find_all_opt_variables(&opt_vars, &S, &M, &tree_M, &gt, &num_partitions_SOS2, &B);	// Complete this function.

		cout<<"Best Objective Stop: "<<opt_vars.phi.sum()/Ns<<endl;
		m.set(GRB_DoubleParam_BestObjStop,opt_vars.phi.sum()/Ns);		// Best Objective Stop Parameter.

		if(set_start_sol_switch==1){
			m.set(GRB_IntParam_StartNodeLimit,2000000000-2);
			MatrixXf S_Sampled = S.block(0,0,3,num_sampled_sens_points_ICP);	// (3 x 500)
			// Perform ICP to find Transformation.
			MatrixXf transformation = icp_test(&V, &S_Sampled, &M, &tree_M, &M_sampled, &tree_M_sampled, &SigmaS, &F, &points_per_face, &ICP_triangle_proj_switch_bucket, &ICP_or_GICP_switch_bucket);
			MatrixXf S_NsSampled = S.block(0,0,3,Ns_sampled);				// (3 x 20)
			transformation = transformation.inverse();

			// Find new opt_variables.
			OptVariables opt_vars1;
			find_all_opt_variables(&opt_vars1, &S_NsSampled, M_global, tree_M_global, &transformation, &num_partitions_SOS2, &B);

			// Provide Initial Solution.
			provide_initialSol(&T[0][0], &R[0][0], &Cb_sampled[0][0], &lam[0][0][0], &W[0][0], &alpha[0][0], &phi[0][0], &transformation, &opt_vars1, &Ns_sampled, &Nm_global, &num_partitions_SOS2);
		}

		// Define All Constraints Here.
		// define_TConstr(&m, &T[0][0], &gt);
		define_RConstr(&m, &R[0][0], &gt);
		define_WConstr(&m, &W[0][0]);
		define_SOSConstr(&m, &lam[0][0][0], &R[0][0], &W[0][0], q, &num_partitions_SOS2);
		define_SOS2Constr(&m, &lam[0][0][0]);
		define_CorrespondenceConstr(&m, &Cb_sampled[0][0], &Ns_sampled, &Nm_global);

		// Alpha for model to sensor.
		if(alpha_switch == 1){
			// Stopped Here
			// Complete this function.
			define_AlphaM2SConstr(&m, &B, &S, &T[0][0], &R[0][0], &Cb_sampled[0][0], &alpha[0][0], M_global, &Ns_sampled, &Nm_global);
		}

		define_phiConstr(&m, &phi[0][0], &alpha[0][0], &Ns_sampled);

		// Define Objective function for minimization.
		GRBLinExpr obj_sum = 0;
		for(int i=0; i<Ns_sampled; i++){
			obj_sum = obj_sum + phi[0][i];
		}
		m.setObjective(obj_sum/Ns_sampled, GRB_MINIMIZE);

		Params params;
		params.Ns_sampled = &Ns_sampled;
		params.S = &S;
		params.num_sampled_sens_points_ICP = &num_sampled_sens_points_ICP;
		params.V = &V;
		params.M = &M;
		params.tree_M = &tree_M;
		params.M_sampled = &M_sampled;
		params.tree_M_sampled = &tree_M_sampled;
		params.SigmaS = &SigmaS;
		params.F = &F;
		params.B = &B;
		params.M_global = M_global;
		params.num_partitions_SOS2 = &num_partitions_SOS2;
		params.points_per_face = &points_per_face;
		params.ICP_triangle_proj_switch_callback = &ICP_triangle_proj_switch_callback;
		params.ICP_or_GICP_switch_callback = &ICP_or_GICP_switch_callback;

		m.update();
		// GRBConstr *c=0;
		// c = m.getConstrs();  
		// for(int i=0;i<m.get(GRB_IntAttr_NumConstrs);++i){
		// 	cout << c[i].get(GRB_StringAttr_ConstrName) << endl;
		// }

		m.write("debug.lp");
		if(callback_switch==1){
			callbackMtoS cb = callbackMtoS(&m, &params);
			m.setCallback(&cb);
			m.optimize();
		}
		else{
			m.optimize();
		}
	}
	catch(GRBException e){
		cout<<"Error code = "<<e.getErrorCode()<<endl;
		cout<<e.getMessage()<<endl;
	}
	catch(...){
		cout<<"Exception during optimization"<<endl;
	}


	return 0;
}