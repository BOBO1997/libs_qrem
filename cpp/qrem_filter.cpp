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
#include "qrem_filter.hpp"
#include "hamming.hpp"
#include "sgs_algorithm.hpp"
#include "qrem_filter_base.hpp"

using namespace std;
using namespace Eigen;

namespace libs_qrem {

     QREM_Filter::QREM_Filter(int num_clbits,
                              vector< vector< vector<double> > > cal_matrices,
                              vector< vector<int> > mit_pattern = vector< vector<int> >(0),
                              vector<int> meas_layout = vector<int>(0)) : 
                              QREM_Filter_Base(num_clbits, cal_matrices, mit_pattern, meas_layout) {
        
    };

    void QREM_Filter::apply(map<string, int> hist,
                                    int d = 0,
                                    double threshold = 0.1) {
        int shots = 0;
        for (const auto& item: hist) {
            shots += item.second;
        }

        chrono::system_clock::time_point t_start = chrono::system_clock::now();

        /*------------ preprocess ------------*/

        set<string> keys;
        map<string, double> prob_dist;
        for (auto& item: hist) {
            keys.insert(item.first);
            prob_dist.insert(make_pair(item.first, (double)item.second / shots));
        }

        set<string> extended_keys = extend_keys(keys, d);
        vector<string> extended_keys_vector(extended_keys.begin(), extended_keys.end());

        this->_indices_to_keys_vector = vector<string>(extended_keys.size());
        int i = 0;
        for (const auto& key: extended_keys) {
            this->_indices_to_keys_vector[i] = key;
            i++;
        }

        vector<double> extended_y = extend_vectors(prob_dist, this->_indices_to_keys_vector);

        // set indices_of_matrices table
        this->_indices_of_matrices = vector< vector<int> >(this->_indices_to_keys_vector.size(), vector<int>(this->_pinv_matrices.size(), 0));
        for (size_t source_index = 0; source_index < this->_indices_to_keys_vector.size(); source_index++) {
            string source_state = this->_indices_to_keys_vector[source_index];
            for (size_t i = 0; i < this->_pinv_matrices.size(); i++) {
                int matrix_index = this->index_of_matrix(source_state, this->_poses_clbits[i]);
                this->_indices_of_matrices[source_index][i] = matrix_index;
            }
        }
        
        // time for preprocess
        chrono::system_clock::time_point t_prep = chrono::system_clock::now();
        double dur_prep = std::chrono::duration_cast<std::chrono::milliseconds>(t_prep - t_start).count();
        this->_durations.insert(make_pair("preprocess", dur_prep));

        /*------------ inverse operation ------------*/

        this->_x_s = vector<double>(extended_y.size(), 0);
        this->_sum_of_x = 0;
        for (size_t state_idx = 0; state_idx < extended_keys.size(); state_idx++) {
            double mitigated_value = mitigate_one_state(state_idx, extended_y, this->_indices_to_keys_vector);
            this->_x_s[state_idx] = mitigated_value;
            this->_sum_of_x += mitigated_value;
        }

        // time for inverse operation
        chrono::system_clock::time_point t_inv = chrono::system_clock::now();
        double dur_inv = std::chrono::duration_cast<std::chrono::milliseconds>(t_inv - t_prep).count();
        this->_durations.insert(make_pair("inverse", dur_inv));

        /*------------ correction by delta ------------*/

        double sum_of_vi = this->sum_of_tensored_vector(this->choose_vecs(string(this->_num_clbits, '0'), this->_pinvVs));
        double lambda_i = this->sum_of_tensored_vector(this->choose_vecs(string(this->_num_clbits, '0'), this->_pinvSigmas));
        double delta_denom = (sum_of_vi * sum_of_vi) / (lambda_i * lambda_i);
        double delta_coeff = (1 - this->_sum_of_x) / delta_denom;
        double delta_col = sum_of_vi / (lambda_i * lambda_i);
        vector<double> v_col = this->col_basis(0, this->_pinvVs, this->_indices_to_keys_vector);

        this->_x_hat = vector<double>(this->_x_s.size(), 0);
        this->_sum_of_x_hat = 0;
        for (size_t state_idx = 0; state_idx < extended_keys.size(); state_idx++) {
            this->_x_hat[state_idx] = this->_x_s[state_idx] + delta_coeff * delta_col * v_col[state_idx];
            this->_sum_of_x_hat += this->_x_hat[state_idx];
        }

        // time for correction by delta
        chrono::system_clock::time_point t_delta = chrono::system_clock::now();
        double dur_delta = std::chrono::duration_cast<std::chrono::milliseconds>(t_delta - t_inv).count();
        this->_durations.insert(make_pair("delta", dur_delta));

        /*------------ sgs algorithm ------------*/

        this->_x_tilde = sgs_algorithm(this->_x_hat);

        // time for sgs algorithm
        chrono::system_clock::time_point t_sgs = chrono::system_clock::now();
        double dur_sgs = std::chrono::duration_cast<std::chrono::milliseconds>(t_sgs - t_delta).count();
        this->_durations.insert(make_pair("sgs_algorithm", dur_sgs));

        /*------------ recovering histogram ------------*/
        this->_sum_of_x_tilde = 0;
        this->_mitigated_hist.clear();
        for (size_t i = 0; i < this->_indices_to_keys_vector.size(); i++) {
            if (this->_x_tilde[i] != 0) {
                this->_sum_of_x_tilde += this->_x_tilde[i];
                this->_mitigated_hist.insert(make_pair(this->_indices_to_keys_vector[i], this->_x_tilde[i] * shots));
            }
        }

        // time for postprocess
        chrono::system_clock::time_point t_finish = chrono::system_clock::now();
        double dur_postprocess = std::chrono::duration_cast<std::chrono::milliseconds>(t_finish - t_sgs).count();
        this->_durations.insert(make_pair("postprocess", dur_postprocess));
        
        // total time
        double dur_total = std::chrono::duration_cast<std::chrono::milliseconds>(t_finish - t_start).count();
        this->_durations.insert(make_pair("total", dur_total));
        
        return;
    }
}