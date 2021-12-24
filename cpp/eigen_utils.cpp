#include <vector>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/SVD>

#include "eigen_utils.hpp"

using namespace std;
using namespace Eigen;

namespace libs_qrem {

    Matrix2d vector2d_to_matrix2d(vector< vector<double> > data) {
        Matrix2d eMatrix(data.size(), data[0].size());
        for (size_t i = 0; i < data.size(); ++i)
            eMatrix.row(i) = VectorXd::Map(&data[i][0], data[0].size());
        return eMatrix;
    }

}