#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <Eigen/Core>


using namespace std;
using namespace Eigen;

Matrix2d ConvertToEigenMatrix(vector< vector<double> > data) {
    Matrix2d eMatrix(data.size(), data[0].size());
    for (int i = 0; i < data.size(); ++i)
        eMatrix.row(i) = VectorXd::Map(&data[i][0], data[0].size());
    return eMatrix;
}

int main() {
	vector< vector<double> > cal_matrix(2, vector<double>(2, 0));

	cal_matrix[0][0] = 1;
	cal_matrix[0][1] = 2;
	cal_matrix[1][0] = 3;
	cal_matrix[1][1] = 4;

	double mat[2][2];
	mat[0][0] = 1;
	mat[0][1] = 2;
	mat[1][0] = 3;
	mat[1][1] = 4;

	Matrix2d mat2d = ConvertToEigenMatrix(cal_matrix);

	// Matrix2d mat2d(mat);

	/*
	Matrix2d mat(2, 2);
	mat[0][0] = cal_matrix[0][0];
	mat[0][1] = cal_matrix[0][1];
	mat[1][0] = cal_matrix[1][0];
	mat[1][1] = cal_matrix[1][1];
	*/
		//mat.row(i) = VectorXd::Map(cal_matrix[i][0], cal_matrix[i].size());
	// for (int i = 0; i < 2; i++) {
	// }
	cout << mat2d << endl;
}	

