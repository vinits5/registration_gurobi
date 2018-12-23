#include "gurobi_helper.h"
#include "/home/vinit/gurobi/gurobi752/linux64/include/gurobi_c++.h"
#include <iostream>
#include <string>

using namespace std;

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
			*(var+col_size*i+j)=model->addVar(lb,ub,0,dtype,name);
		}
	}
}

// ################################################################
// If error comes in future. Do check the below function. It could be potential reason of failure.
// ################################################################

// Define Gurobi Variables in 3D Array.
void define_3D_GRBVar(GRBModel *model, GRBVar *var, int row_size, int col_size, int dim3_size, long double lb, long double ub, char dtype, string name){
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
				*(var+row_size*col_size*i+col_size*j+k)=model->addVar(lb,ub,0,dtype,name);
			}
		}
	}
}