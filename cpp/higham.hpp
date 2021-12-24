#include <Eaigen/Dense>

using namespace std;
using namespace Eigen;

namespace libs_qrem {
    double sign(double x);
    pair<double, int> max_and_argmax(VectorXd x);
    double higham_direct(Matrix2d A);
    double higham_bicgstab(Matrix2d A);
}