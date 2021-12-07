#include <iostream>
#include <vector>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

int main() {
	vector< vector<double> > cal_matrix(2, vector<double>(2, 0));

	cal_matrix[0][0] = 1;
	cal_matrix[0][1] = 2;
	cal_matrix[1][0] = 3;
	cal_matrix[1][1] = 4;

	Matrix2d mat;
	for (int i = 0; i < 2; i++) {
		mat.row(i) = VectorXd::Map(cal_matrix[i][0], cal_matrix[i].size());
	}
	cout << mat << endl;
}	

