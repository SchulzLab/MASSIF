#include "Matrix_new.hpp"

//column -> row
// notice that the matrix index starts by 0 but number of rows and number of columns starts by one
int main(){

	cout << "HIHI" << endl;
	//create matrix with number of rows and number of columns
	Matrix<int> m_(2,2);
	cout << "........." << endl;
	cout << m_ << endl; 
	//read element on position  1,1 (index matrix 0,0)
	cout << m_(1,1)	<< endl;

	//write element
	m_(1,2) = 2;
	cout << m_(1,2)<< endl;
	cout << "......" << endl;
	return 0;
}
