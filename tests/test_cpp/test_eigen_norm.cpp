#include <cmath>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/IterativeLinearSolvers>
#include <unsupported/Eigen/IterativeSolvers>
#include <iostream>

#include "../../cpp/eigen_utils.hpp"

using namespace std;
using namespace Eigen;

int main() {
    MatrixXd A(2,2);
    A(0,0) = 1;
    A(0,1) = 10;
    A(1,0) = 100;
    A(1,1) = 1000;
    cout << "A" << endl;
    cout << A << endl;
    cout << "1 norm " << A.lpNorm<1>() << endl;
    cout << "infty norm " << A.lpNorm<Infinity>() << endl;

    A(0,0) = 1;
    A(0,1) = 2;
    A(1,0) = 3;
    A(1,1) = 4;

    cout << "A" << endl;
    cout << A << endl;
    cout << "A.inverse()" << endl;
    cout << A.inverse() << endl;
    MatrixXd B = A.inverse();
    cout << "B" << endl;
    cout << B << endl;
    cout << "1 norm " << A.inverse().lpNorm<1>() << endl; //! wrong!!!!!
    cout << "1 norm " << B.lpNorm<1>() << endl; //! wrong!!!!!
    cout << "infty norm " << A.inverse().lpNorm<Infinity>() << endl;
    cout << "infty norm " << B.lpNorm<Infinity>() << endl;
    cout << "original " << "1 norm " << libs_qrem::compute_one_norm_of_MatrixXd(A.inverse()) << endl;
    cout << "original " << "1 norm " << libs_qrem::compute_one_norm_of_MatrixXd(B) << endl;
    // cout << "original " << "infty norm " << A.inverse().lpNorm<Infinity>() << endl;

    VectorXd v(4);
    v(0) = 1;
    v(1) = 10;
    v(2) = 100;
    v(3) = 1000;
    cout << v << endl;
    cout << "1 norm " << v.lpNorm<1>() << endl;
    cout << "infty norm " << v.lpNorm<Infinity>() << endl;
    cout << "original " << "1 norm " << libs_qrem::compute_one_norm_of_VectorXd(v) << endl;
    cout << "original " << "infty norm " << libs_qrem::compute_infty_norm_of_VectorXd(v) << endl;
}