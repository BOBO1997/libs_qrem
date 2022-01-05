#include <iostream>
#include <set>
#include <map>
#include <vector>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/SVD>
#include <Eigen/LU>
#include <Eigen/IterativeLinearSolvers>
#include <unsupported/Eigen/IterativeSolvers>
#include <algorithm>
#include <ctime>

using namespace Eigen;
using namespace std;

int main() {
    int n = 10000;
    VectorXd x(n), b(n);
    MatrixXd A(n,n);
    for (int i = 0; i < n; i++) {
        b(i) = i;
        for (int j = 0; j < n; j++) {
            A(i, j) = i - j;
        }
    }
    /* ... fill A and b ... */ 
    BiCGSTAB<MatrixXd> solver;
    solver.compute(A);
    x = solver.solve(b);
    std::cout << "#iterations:     " << solver.iterations() << std::endl;
    std::cout << "estimated error: " << solver.error()      << std::endl;
    /* ... update b ... */
    // x = solver.solve(b); // solve again
}