#include <iostream>
#include <set>
#include <map>
#include <vector>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/SVD>
#include <time.h>
#include <algorithm>
#include <chrono>
#include <ctime>

#include "eigen_utils.hpp"
#include "hamming.hpp"
#include "sgs_algorithm.hpp"
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

    void NationEtal_Filter::apply(map<string, int> hist,
                                  int d,
                                  double threshold) {

        tp_now t_start = chrono::system_clock::now();

        /*------------ preprocess ------------*/
        
        vector<double> extended_y = this->preprocess(hist, d); 

        // time for preprocess
        tp_now t_prep = chrono::system_clock::now();
        this->_durations.insert(make_pair("preprocess", chrono::duration_cast<chrono::milliseconds>(t_prep - t_start).count()));

        /*------------ inverse operation ------------*/

        compute_reduced_inv_A(this->_indices_to_keys_vector.size());
        this->_x_s = stdvec2d_stdvec_prod(this->_reduced_inv_A, extended_y);
        this->_sum_of_x = 0;
        for (size_t i = 0; i < this->_x_s.size(); i++) {
            this->_sum_of_x += this->_x_s[i];
        }

        // time for inverse operation
        tp_now t_inv = chrono::system_clock::now();
        this->_durations.insert(make_pair("inverse", chrono::duration_cast<chrono::milliseconds>(t_inv - t_prep).count()));

        /*------------ correction by least norm problem ------------*/

        this->_x_hat = vector<double>(this->_x_s.size(), 0);
        this->_sum_of_x_hat = 0;
        double diff = (1.0 - this->_sum_of_x) / (double)this->_x_s.size();
        for (size_t state_idx = 0; state_idx < this->_x_s.size(); state_idx++) {
            this->_x_hat[state_idx] = this->_x_s[state_idx] + diff;
            this->_sum_of_x_hat += this->_x_hat[state_idx];
        }

        // time for correction by delta
        tp_now t_lnp = chrono::system_clock::now();
        double dur_lnp = chrono::duration_cast<chrono::milliseconds>(t_lnp - t_inv).count();
        this->_durations.insert(make_pair("least_norm", dur_lnp));

        /*------------ sgs algorithm ------------*/

        this->_x_tilde = sgs_algorithm(this->_x_hat, false);

        // time for sgs algorithm
        tp_now t_sgs = chrono::system_clock::now();
        this->_durations.insert(make_pair("sgs_algorithm", chrono::duration_cast<chrono::milliseconds>(t_sgs - t_lnp).count()));

        /*------------ recovering histogram ------------*/

        this->_sum_of_x_tilde = 0;
        this->_mitigated_hist.clear();
        for (size_t i = 0; i < this->_indices_to_keys_vector.size(); i++) {
            if (this->_x_tilde[i] != 0) {
                this->_sum_of_x_tilde += this->_x_tilde[i];
                this->_mitigated_hist.insert(make_pair(this->_indices_to_keys_vector[i], this->_x_tilde[i] * this->_shots));
            }
        }

        // time for postprocess and total time
        tp_now t_finish = chrono::system_clock::now();
        this->_durations.insert(make_pair("postprocess", chrono::duration_cast<chrono::milliseconds>(t_finish - t_sgs).count()));
        this->_durations.insert(make_pair("total", chrono::duration_cast<chrono::milliseconds>(t_finish - t_start).count())); 

        return;
    }
}