#include <iostream>
#include <eigen3/Eigen/Dense>

using namespace Eigen;
using namespace std;

void print3DArray(int *array_name, int rows, int cols, int depth){	
	for(int i=0; i<rows; i++){
		for(int j=0; j<cols; j++){
			for(int k=0; k<depth; k++){
				cout<<*(array_name+rows*cols*i+cols*j+k)<<endl;
			}
		}
	}
}

void print3Darray(int *array_name){
	for(int i=0; i<3; i++){
		for(int j=0; j<2; j++){
			for(int k=0; k<2; k++){
				cout<<*(array_name+4*i+2*j+k)<<endl;
			}
		}
	}
}

void print_val(int *array_name, int i, int j){
	cout<<*(array_name+3*i+j)<<endl;
}

int main(){
	const int rows=2,cols=2,depth=2;
	// int array_name[rows][cols][depth]={{{1,2},{3,4}},{{5,6},{7,8}}};
	// print3DArray(&array_name[0][0][0], rows, cols, depth);

	int array3d[3][2][2]={{{1,2},{3,4}},{{5,6},{7,8}},{{9,10},{11,12}}};
	// for(int i=0; i<3; i++){
	// 	for(int j=0; j<2; j++){
	// 		for(int k=0; k<2; k++){
	// 			cout<<array3d[i][j][k]<<endl;
	// 		}
	// 	}
	// }
	// print3Darray(&array3d[0][0][0]);

	int array2d[2][3]={{1,2,3},{4,5,6}};
	int i=1,j=2;
	print_val(&array2d[0][0],i,j);
	return 0;
}