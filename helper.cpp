#include "helper.h"
#include "/home/vinit/gurobi/gurobi752/linux64/include/gurobi_c++.h"
#include <iostream>
#include <eigen3/Eigen/Dense>
#include <string>
#include <fstream>

using namespace std;
using namespace Eigen;

// Convert B matrix from (Nx6) to Nx3x3) matrix.
vector<Matrix<float,3,3>> B_Nsx6_to_3x3(MatrixXf *B_given){
	// Arguments:
		// B_given:		Matrix from B.txt file. (Nx6)
	// Output:
		// B_return:	Convert it to (Nx3x3)

	vector<Matrix<float,3,3>> B_return;										// Vector to store 3x3 matrices.
	for(int i=0; i<B_given->rows(); i++){
		Matrix<float,3,3> temp = MatrixXf::Zero(3,3);						// Temporary vector to store each 3x3 matrix.

		// B_given[i]:		[a,b,c,d,e,f]
		// B_return[i]:		[[a,b,c],
		//  				 [b,d,e],
		//  				 [c,e,f]]

		temp(0,0) = B_given->coeff(i,0);
		temp(0,1) = B_given->coeff(i,1); temp(1,0) = B_given->coeff(i,1);	// Same diagonal elements. 
		temp(0,2) = B_given->coeff(i,2); temp(2,0) = B_given->coeff(i,2);	// Same diagonal elements.
		temp(1,1) = B_given->coeff(i,3);
		temp(1,2) = B_given->coeff(i,4); temp(2,1) = B_given->coeff(i,4);	// Same diagonal elements.
		temp(2,2) = B_given->coeff(i,5);
		B_return.push_back(temp);
	}
	return B_return;
}

// Used to calculate the columns and rows in .txt file.
void data_shape(string file_name, int *row_size, int *col_size){
	// Args.
	// file_name: name of file in string.
	// row_size: memory location to store no. of rows. (values seperated by ',')
	// col_size: memory location to store no. of cols. (total no of lines in a file)

	ifstream file(file_name);						// Open a file
	string x;
	getline(file,x,'\n');							// Read first line and store in x.
	*col_size = count(x.begin(),x.end(),',');		// count number of commas appeared in first line.
	*col_size = *col_size+1;						// no of columns = no of commas in first line + 1
	*row_size = *row_size+1;						// add 1 to no of rows as first line is read.

	while(!file.eof()){								// till the end of file.
		*row_size= *row_size+1;						// keep adding 1 to no of rows.
		getline(file,x,'\n');						// read each line.
	}
	file.close();									// close the file.
}

// Define the use of this function by looking at code.
void find_all_opt_variables(OptVariables *opt_vars, MatrixXf *S, MatrixXf *M, KDTree *tree_M, MatrixXf *gt, int *num_partitions_SOS2, vector<Matrix<float,3,3>> *B){
	// Arguments:
		// S: 			Pointer of Sensor Data (3 x Ns)
		// M:			Pointer of Model Data (3 x Nm)
		// tree_M:		KDTree of Model Data
		// gt:			Pointer of transformation matrix from model to sensor (4x4)
		// B: 			Data from the text file. (Ns x 3 x 3)

	Matrix<float,3,1> translation = gt->block(0,3,3,1);
	MatrixXf transformed_sens_pts = (gt->block(0,0,3,3)*(*S)).colwise()+translation;	// Transform the sensor points (R*S+t) (3xNs)
	find_Cb(&(opt_vars->Cb), &transformed_sens_pts, tree_M, M->cols());		// Function to find correspondences.
	find_lam(&(opt_vars->lam), gt, num_partitions_SOS2);
	find_w(&(opt_vars->w), &(opt_vars->lam), num_partitions_SOS2);
	find_alpha(&(opt_vars->alpha), gt, S, M, &(opt_vars->Cb), B);
	// ##########################################################
	// Stopped Here.
	// Some Problem with sizes of matrix. Check it
	// ##########################################################
	// find_phi(&(opt_vars->phi), &(opt_vars->alpha));
}

// Fucntion to find alpha 
void find_alpha(MatrixXf *alpha, MatrixXf *gt, MatrixXf *S, MatrixXf *M, MatrixXf *Cb , vector<Matrix<float,3,3>> *B){
	// Arguments:
		// gt: 		Pointer of ground truth transformation (4 x 4)
		// S:		Pointer of Sensor Point Cloud (3 x Ns)
		// M:		Pointer of Model Point Cloud (3 x Nm)
		// Cb:		Correspondence Matrix. (Ns x Nm)
		// B: 			Data from the text file. (Ns x 3 x 3)
	// Output:
		// alpha:	Alpha parameter (3 x Ns)

	MatrixXf gt_inv = gt->inverse();
	*alpha = MatrixXf::Zero(3,S->cols());
	for(int i=0; i<S->cols(); i++){
		// This part is very slow in computation (even slower than python in my PC. Discuss it!!!)
		alpha->col(i)=(*B)[i]*(S->col(i)-gt_inv.block(0,3,3,1)-((gt_inv.block(0,0,3,3)*(*M))*(Cb->row(i)).transpose()));
	}
	*alpha = alpha->cwiseAbs();
}

// Function to find Correspondences Matrix.
void find_Cb(MatrixXf *Cb, MatrixXf *transformed_sens_pts, KDTree *tree_M, int Cb_cols){
	// Arguments:
		// transfomed_sens_pts:		Pointer of Transformed Sensor Points (3xN)
		// tree_M:					KD-Tree of Model Points
	// Output:
		// Cb:						Binary Output matrix for correspondences. (Ns x Nm)

	*Cb = MatrixXf::Zero(transformed_sens_pts->cols(),Cb_cols);	// Define a matrix for binary correspondences. (Ns x Nm)
	PointCloud sensor_ptCloud = transformed_sens_pts->cast<long double>();		// Store the sensor points to PointCloud datatype.
	MatrixXld neighbours = kd_search_targets(&sensor_ptCloud, *tree_M);			// Find Neighbouring points and their indices.
	for(int i=0; i<neighbours.cols(); i++){
		Cb->coeffRef(i,neighbours.coeff(3,i))=1;		// Assign 1 to correspondences given by KNN.
	}
}

// Find lambda value
void find_lam(vector<Matrix<float,3,3>> *lam, MatrixXf *gt, int *num_partitions_SOS2){
	// Arguments:
		// gt:					Pointer of the groud truth tranformation. (4x4)
		// num_partitions_SOS2:	Size of the lambda array.
	// Output:
		// lam:					Name of the lamda array to store lambda matrices (num_partitions_SOS2 x 3 x 3)

	MatrixXf gt_inv = gt->inverse();				// Find the inverse of ground truth rotation matrix.
	for(int i=0; i<*num_partitions_SOS2; i++){
		lam->push_back(Matrix<float,3,3>::Zero());			// Initialize zero matrices in each element of array.
	}
	float q[*num_partitions_SOS2];
	linspace(&q[0], -1, 1, *num_partitions_SOS2);	// function similar to np.linspace() in python.
	// Check condition and set the values of lamda parameter.
	for(int i=0; i<3; i++){
		for(int j=0; j<3; j++){
			for(int k=0; k<(*num_partitions_SOS2-1); k++){
				if(q[k]<=gt_inv(i,j) && gt_inv(i,j)<=q[k+1]){
					(*lam)[k](i,j)=(gt_inv(i,j)-q[k+1])/(q[k]-q[k+1]);
					(*lam)[k+1](i,j)=1-((gt_inv(i,j)-q[k+1])/(q[k]-q[k+1]));
				}
			}
		}
	}
}

void find_phi(MatrixXf *phi, MatrixXf *alpha){
	*phi = (alpha->colwise()).sum();
}

// Function to find W parameter.
void find_w(Matrix<float,3,3> *w, vector<Matrix<float,3,3>> *lam, int *num_partitions_SOS2){
	// Arguments:
		// lam:					Pointer of the lamda array (only name of array from main function) (num_partitions_SOS2 x 3 x 3)
		// num_partitions_SOS2:	Size of the lambda array & q linspace.
	// Output:
		// w:					Pointer of the matrix w_start (3 x 3)

	float q[*num_partitions_SOS2];
	linspace(&q[0], -1, 1, *num_partitions_SOS2);	// function similar to np.linspace() in python.
	for(int i=0; i<3; i++){
		for(int j=0; j<3; j++){
			float sum = 0;
			for(int k=0; k<*num_partitions_SOS2; k++){
				sum = sum + (*lam)[k](i,j)*q[k]*q[k];
			}
			w->coeffRef(i,j)=sum;
		}
	}
}

// function similar to numpy.linspace
void linspace(float *result_array, int lower_value, int upper_value, int size){
	// Arguments:
		// result_array: 	Pointer of first element of result array.
		// lower_value:		Starting value.
		// upper_value:		Ending value.
		// size:			Size of the array.

	for(int i=0; i<size; i++){
		*(result_array+i)=lower_value + i*((float)(upper_value - lower_value)/(size-1));
	}
}

// Just to check the size of given matrix.
void printMatrixSize(MatrixXf *mat){
	// Just to save memory give the address of matrix.
	// mat->rows()		mat: Matrix class in Eigen library & rows() is a function in that class. 
	// calling function of a pointer needs arrow(->) operator.
	cout<<"Rows: "<<mat->rows()<<" & "<<"Cols: "<<mat->cols()<<endl;
}

// Function to read given file and return data as a Matrix in Eigen library.
MatrixXf read_file(string file_name){
	// Args.
	// file_name: name of file in string.
	// Output
	// MatrixXf: matrix data type in eigen library with the data in given file.

	int row_size=0,col_size=0;						// to define shape of matrix.
	data_shape(file_name,&row_size,&col_size);		// to find shape of matrix.
	MatrixXf data(row_size,col_size);				// define a matrix.
	ifstream file(file_name);						// read the file.
	string x;										// to store the string from the file.
	for(int i=0; i<row_size; i++){					// loop to read rows.
		for(int j=0; j<col_size; j++){				// loop to read cols.
			if(j!=(col_size-1)){					// delimiter is ',' if not the last column.	
				getline(file,x,',');				// store the string in x.
			}
			else{									// delimiter is '\n' if it is the last column.
				getline(file,x,'\n');				// store the string in x.
			}
			data(i,j) = stof(x);					// convert string to float value in matrix.
		}
	}
	file.close();									// close the file.
	return data;									// return the matrix.
}

// Sample model points after certain interval.
void sampleModelPoints(MatrixXf *M, MatrixXf *M_sampled, int *num_sampled_model_points){
	// Args.
	// M:				 			Model points (3xN)
	// M_sampled: 					Sampled Model Points (3xN1)
	// num_sampled_model_points: 	Interval to sample points from M.
	
	for(int i=0; i<M_sampled->cols(); i++){
		M_sampled->block(0,i,1,1) = M->block(0,i*(*num_sampled_model_points),1,1);
		M_sampled->block(1,i,1,1) = M->block(1,i*(*num_sampled_model_points),1,1);
		M_sampled->block(2,i,1,1) = M->block(2,i*(*num_sampled_model_points),1,1);
	}
}