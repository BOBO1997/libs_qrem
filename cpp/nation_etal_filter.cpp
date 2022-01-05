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
#include <cassert>

#include "eigen_utils.hpp"
#include "hamming.hpp"
#include "sgs_algorithm.hpp"
#include "harger_higham.hpp"
#include "qrem_filter.hpp"
#include "nation_etal_filter.hpp"

using namespace std;
using namespace Eigen;

namespace libs_qrem {

     NationEtal_Filter::NationEtal_Filter(int num_clbits,
                              vector< vector< vector<double> > > cal_matrices,
                              vector< vector<int> > mit_pattern = vector< vector<int> >(0),
                              vector<int> meas_layout = vector<int>(0)) : 
                              QREM_Filter(num_clbits, cal_matrices, mit_pattern, meas_layout) {};

    void NationEtal_Filter::apply(Args args) {

        this->_exact_one_norm_of_inv_reduced_A = -1;
        this->_iterative_one_norm_of_inv_reduced_A = -1;
        
        map<string, int> hist = args.hist;
        int d = args.d;
        string method = args.method;

        clock_t t_start = clock();

        /*------------ preprocess ------------*/
        
        vector<double> extended_y = this->preprocess(hist, d); 
        clock_t t_prep = clock();

        /*------------ inverse operation ------------*/

        this->compute_reduced_A(this->_indices_to_keys_vector.size());
        normalize_cols(this->_reduced_A);
        clock_t t_mat = clock();

        /*------------ matrix free method ------------*/

        MatrixXd A = stdvec2d_to_MatrixXd(this->_reduced_A);
        VectorXd b = stdvec1d_to_VectorXd(extended_y);
        if (method == "iterative" | method == "bicgstab") {
            BiCGSTAB<MatrixXd> solver(A);
            VectorXd v = solver.solve(b);
            this->_iterations = solver.iterations();
            this->_error = solver.error();
            this->_x_s = VectorXd_to_stdvec1d(v);
        }
        else {
            PartialPivLU<MatrixXd> solver(A);
            VectorXd v = solver.solve(b);
            this->_x_s = VectorXd_to_stdvec1d(v);
        }
        this->_sum_of_x = 0;
        for (size_t i = 0; i < this->_x_s.size(); i++) {
            this->_sum_of_x += this->_x_s[i];
        }
        clock_t t_inv = clock();

        /*------------ sgs algorithm ------------*/

        this->_x_tilde = sgs_algorithm(this->_x_s, false);
        clock_t t_sgs = clock();

        /*------------ recovering histogram ------------*/

        this->recover_histogram();
        clock_t t_finish = clock();

        this->_durations.insert(make_pair("preprocess", (double)(t_prep - t_start) / CLOCKS_PER_SEC));
        this->_durations.insert(make_pair("reduction", (double)(t_mat - t_prep) / CLOCKS_PER_SEC));
        this->_durations.insert(make_pair("inverse", (double)(t_inv - t_mat) / CLOCKS_PER_SEC));
        this->_durations.insert(make_pair("sgs_algorithm", (double)(t_sgs - t_inv) / CLOCKS_PER_SEC));
        this->_durations.insert(make_pair("postprocess", (double)(t_finish - t_sgs) / CLOCKS_PER_SEC));
        this->_durations.insert(make_pair("total", (double)(t_finish - t_start) / CLOCKS_PER_SEC)); 

        return;
    }
}