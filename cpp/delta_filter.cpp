#include <iostream>
#include <set>
#include <map>
#include <vector>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <algorithm>
#include <ctime>

#include "eigen_utils.hpp"
#include "hamming.hpp"
#include "sgs_algorithm.hpp"
#include "qrem_filter.hpp"
#include "delta_filter.hpp"

using namespace std;
using namespace Eigen;

namespace libs_qrem {

     Delta_Filter::Delta_Filter(int num_clbits,
                                vector< vector< vector<double> > > cal_matrices,
                                vector< vector<int> > mit_pattern = vector< vector<int> >(0),
                                vector<int> meas_layout = vector<int>(0)) : 
                                QREM_Filter(num_clbits, cal_matrices, mit_pattern, meas_layout) {};
    
    void Delta_Filter::apply(Args args) {
        
        map<string, int> hist = args.hist;
        int d = args.d;
        clock_t t_start = clock();

        /*------------ preprocess ------------*/

        vector<double> extended_y = this->preprocess(hist, d);
        clock_t t_prep = clock();

        /*------------ inverse operation ------------*/
        
        this->compute_reduced_inv_A(this->_indices_to_keys_vector.size());
        this->_x_s = stdvec2d_stdvec_prod(this->_reduced_inv_A, extended_y);
        this->_sum_of_x = 0;
        for (size_t i = 0; i < this->_x_s.size(); i++) {
            this->_sum_of_x += this->_x_s[i];
        }
        clock_t t_inv = clock();

        /*------------ correction by delta ------------*/

        double sum_of_vi = this->sum_of_tensored_vector(this->choose_vecs(string(this->_num_clbits, '0'), this->_pinvVs));
        double lambda_i = this->sum_of_tensored_vector(this->choose_vecs(string(this->_num_clbits, '0'), this->_pinvSigmas));
        double delta_denom = (sum_of_vi * sum_of_vi) / (lambda_i * lambda_i);
        double delta_coeff = (1 - this->_sum_of_x) / delta_denom;
        double delta_col = sum_of_vi / (lambda_i * lambda_i);
        vector<double> v_col = this->col_basis(0, this->_pinvVs, this->_indices_to_keys_vector);

        this->_x_hat = vector<double>(this->_x_s.size(), 0);
        this->_sum_of_x_hat = 0;
        for (size_t state_idx = 0; state_idx < this->_indices_to_keys_vector.size(); state_idx++) {
            this->_x_hat[state_idx] = this->_x_s[state_idx] + delta_coeff * delta_col * v_col[state_idx];
            this->_sum_of_x_hat += this->_x_hat[state_idx];
        }
        clock_t t_delta = clock();

        /*------------ sgs algorithm ------------*/

        this->_x_tilde = sgs_algorithm(this->_x_hat, false);
        clock_t t_sgs = clock();

        /*------------ recovering histogram ------------*/

        this->recover_histogram();
        clock_t t_finish = clock();

        /*------------ add time durations ------------*/ 
        this->_durations.insert(make_pair("preprocess", (double)(t_prep - t_start) / CLOCKS_PER_SEC));
        this->_durations.insert(make_pair("inverse", (double)(t_inv - t_prep) / CLOCKS_PER_SEC));
        this->_durations.insert(make_pair("delta", (double)(t_delta - t_inv) / CLOCKS_PER_SEC));
        this->_durations.insert(make_pair("sgs_algorithm", (double)(t_sgs - t_delta) / CLOCKS_PER_SEC));
        this->_durations.insert(make_pair("postprocess", (double)(t_finish - t_sgs) / CLOCKS_PER_SEC));
        this->_durations.insert(make_pair("total", (double)(t_finish - t_start) / CLOCKS_PER_SEC));

        return;
    }
}