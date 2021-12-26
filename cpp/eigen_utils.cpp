#include <vector>
#include <cmath>
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


    double compute_one_norm_of_MatrixXd(MatrixXd matrix) {
        double one_norm = 0;
        for (int j = 0; j < matrix.cols(); j++) {
            double col_sum = 0;
            for (int i = 0; i < matrix.rows(); i++) {
                col_sum += abs(matrix(i, j));
            }
            if (one_norm < col_sum) {
                one_norm = col_sum;
            }
        }
        return one_norm;
    }

    double compute_one_norm_of_stdvec2d(vector< vector<double> > matrix) {
        double one_norm = 0;
        for (int j = 0; j < matrix.cols(); j++) {
            double col_sum = 0;
            for (int i = 0; i < matrix.rows(); i++) {
                col_sum += abs(matrix[i][j]));
            }
            if (one_norm < col_sum) {
                one_norm = col_sum;
            }
        }
        return one_norm;
    }

    double compute_one_norm_of_VectorXd(VectorXd vec) {
        double col_sum = 0;
        for (int i = 0; i < vec.size(); i++) {
            col_sum += abs(vec(i));
        }
        return col_sum;
    }

    double compute_one_norm_of_stdvec1d(vector<double> vec) {
        double col_sum = 0;
        for (int i = 0; i < vec.size(); i++) {
            col_sum += abs(vec[i]);
        }
        return col_sum;
    }

    double compute_infty_norm_of_VectorXd(VectorXd vec) {
        double col_max = 0;
        for (int i = 0; i < vec.size(); i++) {
            if (col_max < abs(vec(i))) {
                col_max = abs(vec(i));
            }
        }
        return col_max;
    }

    double compute_infty_norm_of_stdvec1d(vector<double> vec) {
        double col_max = 0;
        for (int i = 0; i < vec.size(); i++) {
            if (col_max < abs(vec[i])) {
                col_max = abs(vec[i]);
            }
        }
        return col_max;
    }
}