#include <vector>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/SVD>

#include "eigen_utils.hpp"

using namespace std;
using namespace Eigen;

namespace libs_qrem {

    Matrix2d stdvec2d_to_matrixXd(vector< vector<double> > data) {
        MatrixXd eMatrix(data.size(), data[0].size());
        for (size_t i = 0; i < data.size(); ++i) {
            eMatrix.row(i) = VectorXd::Map(&data[i][0], data[0].size());
        }
        return eMatrix;
    }

    vector<double> stdvec2d_stdvec_prod(vector< vector<double> >& A, vector<double>& y) {
        // compute x = Ay
        vector<double> x(y.size(), 0);
        for (size_t i = 0; i < x.size(); i++) {
            for (size_t j = 0; j < x.size(); j++) {
                x[i] += A[i][j] * y[j];
            }
        }        
        return x;
    }

    vector<double> matrixXd_stdvec_prod(MatrixXd& A, vector<double>& y) {
        // compute x = Ay
        vector<double> x(y.size(), 0);
        for (size_t i = 0; i < x.size(); i++) {
            for (size_t j = 0; j < x.size(); j++) {
                x[i] += A(i, j) * y[j];
            }
        }        
        return x;
    }

}