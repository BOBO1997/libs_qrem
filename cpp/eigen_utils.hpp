#include <vector>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/SVD>

using namespace std;
using namespace Eigen;

namespace libs_qrem {

    Matrix2d stdvec2d_to_matrixXd(vector< vector<double> > data);
    vector<double> stdvec2d_stdvec_prod(vector< vector<double> >& A, vector<double>& y);
    vector<double> matrixXd_stdvec_prod(MatrixXd& A, vector<double>& y);

}