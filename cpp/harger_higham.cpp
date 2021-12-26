#include <cmath>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/IterativeLinearSolvers>
#include <unsupported/Eigen/IterativeSolvers>
#include <iostream>

#include "harger_higham.hpp"

using namespace std;
using namespace Eigen;

namespace libs_qrem {

    double sign(double x) {
        if (x >= 0) {
            return 1.0;
        }
        else {
            return -1.0;
        }
    }

    pair<double, int> max_and_argmax(VectorXd x) {
        double max_val = x(0);
        int arg_max = 0;
        for (int i = 1; i < x.size(); i++) {
            if (max_val < x(i)) {
                max_val = x(i);
                arg_max = i;
            }
        }
        return make_pair(max_val, arg_max);
    }

    double harger_higham_lu(MatrixXd A) {
        int n = A.rows();
        VectorXd e = VectorXd::Constant(n, 1 / (double)n);
        PartialPivLU<MatrixXd> solver(A), solver_T(A.transpose());
        VectorXd v = solver.solve(e);
        if (n == 1) {
            return v.lpNorm<1>(); // v(0); is also ok
        }
        double gamma = v.lpNorm<1>();
        VectorXd gsi = v.unaryExpr(ptr_fun(sign));
        VectorXd x = solver_T.solve(gsi);
        int k = 2;
        while (k < 6) {
            pair<double, int> x_infty_norm = max_and_argmax(x);
            int j = x_infty_norm.second;
            VectorXd e_j = VectorXd::Unit(n, j);
            v = solver.solve(e_j);
            double gamma_bar = gamma;
            double gamma = v.lpNorm<1>();
            if ((v - gsi).isMuchSmallerThan(1 / (double)n, 1e-3) | (gamma <= gamma_bar)) {
                break;
            }
            gsi = v.unaryExpr(ptr_fun(sign));
            x = solver_T.solve(gsi);
            k++;
            if (x_infty_norm.first == x.lpNorm<Infinity>()) {
                break;
            }
        }
        for (int i = 1; i < x.size() + 1; i++) {
            x(i - 1) = (double)pow(-1, i + 1) * (1 + (double)(i - 1) / (double)(n - 1));
        }
        x = solver.solve(x);
        if (2 * x.lpNorm<1>() / (double)(3 * n) > gamma) {
            v = x;
            gamma = 2 * x.lpNorm<1>() / (double)(2 * n);
        }
        return gamma;
    }
    
    double harger_higham_bicgstab(MatrixXd A) {
        int n = A.rows();
        VectorXd e = VectorXd::Constant(n, 1 / (double)n);
        BiCGSTAB<MatrixXd> solver(A), solver_T(A.transpose());
        VectorXd v = solver.solve(e);
        if (n == 1) {
            return v.lpNorm<1>(); // v(0); is also ok
        }
        double gamma = v.lpNorm<1>();
        VectorXd gsi = v.unaryExpr(ptr_fun(sign));
        VectorXd x = solver_T.solve(gsi);
        int k = 2;
        while (k < 6) {
            pair<double, int> x_infty_norm = max_and_argmax(x);
            int j = x_infty_norm.second;
            VectorXd e_j = VectorXd::Unit(n, j);
            v = solver.solve(e_j);
            double gamma_bar = gamma;
            double gamma = v.lpNorm<1>();
            if ((v - gsi).isMuchSmallerThan(1 / (double)n, 1e-3) | (gamma <= gamma_bar)) {
                break;
            }
            gsi = v.unaryExpr(ptr_fun(sign));
            x = solver_T.solve(gsi);
            k++;
            if (x_infty_norm.first == x.lpNorm<Infinity>()) {
                break;
            }
        }
        for (int i = 1; i < x.size() + 1; i++) {
            x(i - 1) = (double)pow(-1, i + 1) * (1 + (double)(i - 1) / (double)(n - 1));
        }
        x = solver.solve(x);
        if (2 * x.lpNorm<1>() / (double)(3 * n) > gamma) {
            v = x;
            gamma = 2 * x.lpNorm<1>() / (double)(2 * n);
        }
        return gamma;
    }

}