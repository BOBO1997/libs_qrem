#include <Eigen/Core>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

namespace libs_qrem {
    double sign(double x);
    pair<double, int> max_and_argmax(VectorXd x);
    double harger_higham_lu(Matrix2d A);
    double harger_higham_bicgstab(Matrix2d A);
}