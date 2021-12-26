#include <iostream>
#include <vector>
#include <cmath>
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

double compute_one_norm_of_matrix(MatrixXd mat2d) {
	double one_norm = 0;
	for (int j = 0; j < mat2d.cols(); j++) {
		double col_sum = 0;
		for (int i = 0; i < mat2d.rows(); i++) {
			col_sum += abs(mat2d(i, j));
		}
		if (one_norm < col_sum) {
			one_norm = col_sum;
		}
	}
	return one_norm;
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
	cout << mat2d.inverse() << endl;
	cout << mat2d.inverse().lpNorm<1>() << endl;
	cout << compute_one_norm_of_matrix(mat2d.inverse()) << endl;
	cout << mat2d.inverse().col(0) << endl;
	cout << mat2d.inverse().col(0).lpNorm<1>() << endl;
	cout << mat2d.inverse().col(0).lpNorm<Infinity>() << endl;
	cout << mat2d.inverse().col(1) << endl;
	cout << mat2d.inverse().col(1).lpNorm<1>() << endl;
	cout << mat2d.inverse().col(1).lpNorm<Infinity>() << endl;
}	

