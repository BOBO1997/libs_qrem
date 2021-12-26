#include <vector>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/SVD>

using namespace std;
using namespace Eigen;

namespace libs_qrem {

    vector<double> VectorXd_to_stdvec1d(VectorXd data);
    VectorXd stdvec1d_to_VectorXd(vector< double > data);
    MatrixXd stdvec2d_to_MatrixXd(vector< vector<double> > data);
    vector<double> stdvec2d_stdvec_prod(vector< vector<double> >& A, vector<double>& y);
    vector<double> MatrixXd_stdvec_prod(MatrixXd& A, vector<double>& y);
    void normalize_cols(vector< vector<double> >& A);
    double compute_one_norm_of_MatrixXd(MatrixXd matrix);
    double compute_one_norm_of_stdvec2d(vector<double> matrix);
    double compute_one_norm_of_VectorXd(VectorXd vec);
    double compute_one_norm_of_stdvec1d(vector<double> vec);
    double compute_infty_norm_of_VectorXd(VectorXd vec);
    double compute_infty_norm_of_stdvec1d(vector<double> vec);

}