#include <vector>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/SVD>

#include "eigen_utils.hpp"

using namespace std;
using namespace Eigen;

namespace libs_qrem {

    vector<double> VectorXd_to_stdvec1d(VectorXd data) {
        return vector<double>(&data[0], data.data() + data.cols() * data.rows()); //? really ?
    }

    VectorXd stdvec1d_to_VectorXd(vector< double > data) {
        return Map<VectorXd, Unaligned>(data.data(), data.size());
    }

    MatrixXd stdvec2d_to_MatrixXd(vector< vector<double> > data) {
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

    vector<double> MatrixXd_stdvec_prod(MatrixXd& A, vector<double>& y) {
        // compute x = Ay
        vector<double> x(y.size(), 0);
        for (size_t i = 0; i < x.size(); i++) {
            for (size_t j = 0; j < x.size(); j++) {
                x[i] += A(i, j) * y[j];
            }
        }
        return x;
    }

    void normalize_cols(vector< vector<double> >& A) {
        vector<double> sum_of_cols(A.size());
        for (size_t i = 0; i < A.size(); i++) {
            for (size_t j = 0; j < A.size(); j++) {
                sum_of_cols[j] += A[i][j];
            }
        }
        for (size_t i = 0; i < A.size(); i++) {
            for (size_t j = 0; j < A.size(); j++) {
                A[i][j] /= sum_of_cols[j];
            }
        }
    }

}