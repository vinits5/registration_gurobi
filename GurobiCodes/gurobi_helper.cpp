#include "gurobi_helper.h"
#include "/home/vinit/gurobi/gurobi752/linux64/include/gurobi_c++.h"
#include <iostream>
#include <string>

using namespace std;

void define_AlphaM2SConstr(GRBModel *model, vector<Matrix<float,3,3>> *B, MatrixXf *S, GRBVar *T, GRBVar *R, GRBVar *Cb_sampled, GRBVar *alpha, MatrixXf *M_global, int *Ns_sampled, int *Nm_global){
	// Arguments:
		// model:		Gurobi Model
		// B:			Ns x 3 x 3
		// S:			Sensor Data (3 x Ns)
		// T:			Translation variable in gurobi (3 x 1)
		// R:			Rotation variable in gurobi (3 x 3)
		// Cb_sampled:	Correspondence variable in gurobi (Ns_sampled x Nm_global)
		// M_global:	Model data (3 x N)

	for(int i=0; i<*Ns_sampled; i++){
		for(int j=0; j<3; j++){				// Corresponding to jth row of B.
			float sum=0;
			for(int k=0; k<3; k++){
				// sum = sum + B[i][j][k] * S[k][i]
				sum = sum + (*B)[i](j,k) * S->coeff(k,i);
			}
			float term_BSi = sum;

			GRBLinExpr sum1 = 0;
			for(int k=0; k<3; k++){
				// sum1 = sum1 + B[i][j][k] * T[k][0]
				sum1 = sum1 + (*B)[i](j,k) * (*(T+1*k+0));
			}
			GRBLinExpr term_BT = sum1;

			GRBQuadExpr sum3 = 0;
			for(int m=0; m<*Nm_global; m++){
				GRBLinExpr sum4 = 0;
				for(int l=0; l<3; l++){
					GRBLinExpr sum5 = 0;
					for(int k=0; k<3; k++){
						// sum5 = sum5 + B[i][j][k] * R[k][l]
						sum5 = sum5 + (*B)[i](j,k) * (*(R+3*k+l));
					}
					sum4 = sum4 + sum5 * M_global->coeff(l,m);
				}
				// sum3 = sum3 + sum4 * Cb_sampled[i][m]
				sum3 = sum3 + sum4 * (*(Cb_sampled+*Nm_global*i+m));
			}
			GRBQuadExpr term_BRMCi = sum3;

			model->addQConstr(*(alpha+*Ns_sampled*j+i) >= (term_BSi - term_BT - term_BRMCi), "alpha plus");
			model->addQConstr(*(alpha+*Ns_sampled*j+i) >= -(term_BSi - term_BT - term_BRMCi), "alpha plus");
		}
	}
}

void define_CorrespondenceConstr(GRBModel *model, GRBVar *var, int *rows, int *cols){
	// Arguments:
		// model:		Gurobi Model
		// var:			Cb_sampled (Ns_sampled x Nm_global)
		// rows:		Ns_sampled
		// cols:		Nm_global

	// Summation of Cij + oi (Correspondence matrix with outliers)
	for(int i=0; i<*rows; i++){
		GRBLinExpr sum = 0;
		for(int j=0; j<*cols; j++){
			// sum = sum + var[i][j]
			sum = sum + (*(var+(*cols)*i+j));
		}
		model->addConstr(sum==1,"Correspondence Matrix");
	}
}

// Define Gurobi Variables in Array.
void define_GRBVar(GRBModel *model, GRBVar *var, int row_size, int col_size, long double lb, long double ub, char dtype, string name){
	// Arguments
		// model: 				Pointer of the Gurobi Model.
		// var:					Pointer of the first element of the array of Gurobi Variables.
		// row_size, col_size:	Size of Array.
		// lb:					lower bound
		// ub:					upper bound
		// dtype:				Data type of variable (GRB_CONTINUOUS or GRB_BINARY)
		// name:				Name of the variable as string.

	// Loop to allocate a Gurobi Variable to each element of given array.
	for(int i=0; i<row_size; i++){
		for(int j=0; j<col_size; j++){
			// var[i][j]
			*(var+col_size*i+j)=model->addVar(lb,ub,0,dtype,name+"{"+to_string(i)+","+to_string(j)+"}");
		}
	}
}

void define_phiConstr(GRBModel *model, GRBVar *var, GRBVar *alpha, int *cols){
	// Arguments:
		// model:		Gurobi Model
		// var:			phi (1 x Ns_sampled)
		// alpha:		(3 x Ns_sampled)
		// cols:		Ns_sampled

	// for(int i=0; i<*cols; i++){		// Try later with quick sum.
	// 	model->addConstr(*(var+(*cols)*0+i) >= 0,"phi");
	// }

	for(int i=0; i<*cols; i++){
		GRBLinExpr sum = 0;
		for(int l=0; l<3; l++){
			// sum = sum + alpha[l][j]
			sum = sum + *(alpha+(*cols)*l+i);
		}
		model->addConstr(*(var+(*cols)*0+i) >= sum);
	}
}

void define_RConstr(GRBModel *model, GRBVar *var, MatrixXf *gt){
	// Arguments:
		// model:		Gurobi Model
		// var:			R (3 x 3)
		// gt:			Ground Truth Matrix (4 x 4)

	int var_cols = 3;

	// for(int i=0; i<3; i++){
	// 	for(int j=0; j<3; j++){
	// 		if(i==j){
	// 			model->addConstr(*(var+var_cols*i+j)==1,"rotConstr2u");
	// 		}
	// 		else{
	// 			model->addConstr(*(var+var_cols*i+j)==0,"rotConstr2u");
	// 		}
	// 	}
	// }

	// for(int i=0; i<3; i++){
	// 	for(int j=0; j<3; j++){
	// 		model->addConstr(*(var+var_cols*i+j)==gt->coeffRef(i,j),"rotConstr2u");
	// 	}
	// }

	for(int i=0; i<3; i++){
		GRBQuadExpr sum = 0;
		for(int j=0; j<3; j++){
			// sum = sum + var[j][i]
			sum = sum + (*(var+var_cols*j+i))*(*(var+var_cols*j+i));
		}
		model->addQConstr(sum<=1, "rotConstr");
	}

	for(int i=0; i<3; i++){
		for(int j=i+1; j<3; j++){
			GRBQuadExpr sum = 0;
			for(int k=0; k<3; k++){
				// sum = sum + (var[k][i]+var[k][j]) * (var[k][i]+var[k][j])
				sum = sum + (*(var+var_cols*k+i) + *(var+var_cols*k+j)) * (*(var+var_cols*k+i) + *(var+var_cols*k+j));
			}
			model->addQConstr(sum<=2, "rotConstr1");
		}
	}

	// mic check 123
	for(int i=0; i<3; i++){
		for(int j=i+1; j<3; j++){
			GRBQuadExpr sum = 0;
			for(int k=0; k<3; k++){
				// sum = sum + (var[k][i]-var[k][j]) * (var[k][i]-var[k][j])
				sum = sum + (*(var+var_cols*k+i) - *(var+var_cols*k+j)) * (*(var+var_cols*k+i) - *(var+var_cols*k+j));
			}
			model->addQConstr(sum<=2, "rotConstr1");
		}
	}

	int tmpConstr[2]={-1,1};
	for(int i=0; i<2; i++){
		for(int j=0; j<2; j++){
			GRBQuadExpr sum = 0;
			for(int k=0; k<3; k++){
				// sum = sum + (var[k][0] + C1 * var[k][1] + C2 * var[k][2]) * (var[k][0] + C1 * var[k][1] + C2 * var[k][2])
				sum = sum + (*(var+var_cols*k+0) + tmpConstr[i]*(*(var+var_cols*k+1))+tmpConstr[j]*(*(var+var_cols*k+2))) * (*(var+var_cols*k+0) + tmpConstr[i]*(*(var+var_cols*k+1))+tmpConstr[j]*(*(var+var_cols*k+2)));
			}
			model->addQConstr(sum<=3, "rotConstr1");
		}
	}
}

void define_SOSConstr(GRBModel *model, GRBVar *var, GRBVar *R, GRBVar *W, double *q, int *num_partitions_SOS2){
	// Arguments:
		// model:					Gurobi Model
		// var:						lambda (num_partitions_SOS2 x 3 x 3)
		// R:						Rotation Variable for Gurobi (3 x 3)
		// W:						W (3 x 3)
		// q:						linspace (1 x 10)
		// num_partitions_SOS2:		Size of var.

	// R & W Constraints
	int var_cols = 3;
	int var_rows = 3;
	for(int i=0; i<3; i++){
		for(int j=0; j<3; j++){
			GRBLinExpr sum1=0, sum2=0, sum3=0;
			for(int k=0; k<(*num_partitions_SOS2); k++){
				sum1 = sum1 + (*(var+var_cols*var_rows*k+var_cols*i+j));					// var[k][i][j]
				sum2 = sum2 + (*(var+var_cols*var_rows*k+var_cols*i+j)) * q[k];
				sum3 = sum3 + (*(var+var_cols*var_rows*k+var_cols*i+j)) * q[k] * q[k];
			}
			model->addConstr(sum1==1);
			model->addConstr((*(R+var_cols*i+j)) == sum2);
			model->addConstr((*(W+var_cols*i+j)) == sum3);
		}
	}
}

void define_SOS2Constr(GRBModel *model, GRBVar *var){
	// Arguments:
		// model:		Gurobi Model
		// var:			lambda (num_partitions_SOS2 x 3 x 3)

	// SOS2 Constraints
	int var_cols = 3;
	int var_rows = 3;
	for(int i=0; i<3; i++){
		for(int j=0; j<3; j++){
			GRBVar vars[10];
			for(int k=0; k<10; k++){
				vars[k]=(*(var+var_cols*var_rows*k+var_cols*i+j));		// var[k][i][j]
			}
			double weights[] = {1,2,3,4,5,6,7,8,9,10};
			model->addSOS(vars, weights, 10, GRB_SOS_TYPE2);
		}
	}
}

void define_TConstr(GRBModel *model, GRBVar *var, MatrixXf *gt){
	// Arguments:
		// model:		Gurobi Model
		// var:			Translation Var for Gurobi (3 x 1)
		// gt:			Ground Truth Transformation (4 x 4)

	// Setting Constraints for T values
	int var_cols = 1;
	for(int i=0; i<3; i++){
		
		model->addConstr(*(var+var_cols*i+0)==0,"transConstr2u");		// var[i][0]
	}
	for(int i=0; i<3; i++){
		model->addConstr(*(var+var_cols*i+0)==gt->coeffRef(i,3),"transConstr2u");	// var[i][0]
	}
}

void define_WConstr(GRBModel *model, GRBVar *var){
	// Arguments:
		// model:		Gurobi Model
		// var:			W (3 x 3)

	int var_cols = 3;
	// W Constraints
	for(int i=0; i<3; i++){
		GRBLinExpr sum = 0;
		for(int j=0; j<3; j++){
			sum = sum + (*(var+var_cols*i+j));		// var[i][j]
		}
		model->addConstr(sum>=1, "W");
	}
}

// Define Gurobi Variables in 3D Array.
void define_3D_GRBVar(GRBModel *model, GRBVar *var, int dim3_size, int col_size, int row_size, long double lb, long double ub, char dtype, string name){
	// Arguments
		// model: 							Pointer of the Gurobi Model.
		// var:								Pointer of the first element of the array of Gurobi Variables.
		// row_size, col_size, dim3_size:	Size of Array.
		// lb:								lower bound
		// ub:								upper bound
		// dtype:							Data type of variable (GRB_CONTINUOUS or GRB_BINARY)
		// name:							Name of the variable as string.

	// Loop to allocate a Gurobi Variable to each element of given array.
	for(int i=0; i<dim3_size; i++){
		for(int j=0; j<row_size; j++){
			for(int k=0; k<col_size; k++){
				*(var+row_size*col_size*i+col_size*j+k)=model->addVar(lb,ub,0,dtype,name+"{"+to_string(i)+","+to_string(j)+","+to_string(k)+"}");	// var[i][j][k]
			}
		}
	}
}

void provide_initialSol(GRBVar *T, GRBVar *R, GRBVar *Cb_sampled, GRBVar *lam, GRBVar *W, GRBVar *alpha, GRBVar *phi, MatrixXf *transformation, OptVariables *opt_vars, int *Ns_sampled, int *Nm_global, int *num_partitions_SOS2){
	// Provide a better start for R & T.
	for(int i=0; i<3; i++){
		(*(T+1*i+0)).set(GRB_DoubleAttr_Start,transformation->coeff(i,3));		// Set the Translation Part. T[i][0]
		for(int j=0; j<3; j++){
			(*(R+3*i+j)).set(GRB_DoubleAttr_Start, transformation->coeff(i,j));	// Set the Rotation Part. R[i][j]
		}
	}

	// Provide a better start for correspondences.
	for(int i=0; i<(opt_vars->Cb).rows(); i++){
		for(int j=0; j<(opt_vars->Cb).cols(); j++){
			// Ref: sudoku_c++.cpp in examples of gurobi. Cb_sampled[i][j]
			(*(Cb_sampled+*Nm_global*i+j)).set(GRB_DoubleAttr_Start, (opt_vars->Cb)(i,j));		// Set the correspondences to Cb_sampled Gurobi variable.
		}
	}	

	// Provide a better start for lambda
	for(int k=0; k<*num_partitions_SOS2; k++){
		for(int i=0; i<3; i++){
			for(int j=0; j<3; j++){
				(*(lam+9*k+3*i+j)).set(GRB_DoubleAttr_Start, (opt_vars->lam)[k](i,j));		// lam[k][i][j]
			}
		}
	}

	// Provide a better start for W.
	for(int i=0; i<3; i++){
		for(int j=0; j<3; j++){
			(*(W+3*i+j)).set(GRB_DoubleAttr_Start, (opt_vars->w)(i,j));			// W[i][j]
		}
	}

	// Provide a better start for alpha.
	for(int i=0; i<3; i++){
		for(int j=0; j<1; j++){

			(*(alpha+*Ns_sampled*i+j)).set(GRB_DoubleAttr_Start,(opt_vars->alpha)(i,j));		// alpha[i][j]
		}
	}

	// Provide a better start for phi.
	for(int i=0; i<*Ns_sampled; i++){
		(*(phi+*Ns_sampled*0+i)).set(GRB_DoubleAttr_Start,(opt_vars->phi)(0,i));			//phi[0][i]
	}
}
