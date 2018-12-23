#include <iostream>
#include <eigen3/Eigen/Dense>

using namespace Eigen;
using namespace std;

void printArray(int *array_name, int rows, int cols, int depth){	
	for(int i=0; i<rows; i++){
		for(int j=0; j<cols; j++){
			for(int k=0; k<depth; k++){
				cout<<*(array_name+rows*cols*i+cols*j+k)<<endl;
			}
		}
	}
}

int main(){
	const int rows=2,cols=2,depth=2;
	int array_name[rows][cols][depth]={{{1,2},{3,4}},{{5,6},{7,8}}};
	printArray(&array_name[0][0][0], rows, cols, depth);
	return 0;
}