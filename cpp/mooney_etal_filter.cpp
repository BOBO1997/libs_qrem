#include <iostream>
#include <set>
#include <map>
#include <vector>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <algorithm>
#include <ctime>
#include <cstdlib>
#include <bitset>

#include "eigen_utils.hpp"
#include "hamming.hpp"
#include "sgs_algorithm.hpp"
#include "qrem_filter.hpp"
#include "mooney_etal_filter.hpp"

using namespace std;
using namespace Eigen;

namespace libs_qrem {

     MooneyEtal_Filter::MooneyEtal_Filter(int num_clbits,
                              vector< vector< vector<double> > > cal_matrices,
                              vector< vector<int> > mit_pattern = vector< vector<int> >(0),
                              vector<int> meas_layout = vector<int>(0)) : 
                              QREM_Filter(num_clbits, cal_matrices, mit_pattern, meas_layout) {};

    int MooneyEtal_Filter::flip_state(int state_idx, int mat_idx, vector<int>& flip_poses) {
        for (size_t i = 0; i < flip_poses.size(); i++) { // if completely local, size is 1
            if ((mat_idx >> i) & 1) {
                state_idx = state_idx ^ (1 << (this->_num_clbits - 1 - flip_poses[i]));
            }
        }
        return state_idx;
    }

    string MooneyEtal_Filter::flip_state(string state, int mat_idx, vector<int>& flip_poses) {
        for (size_t i = 0; i < flip_poses.size(); i++) { // if completely local, size is 1
            if ((mat_idx >> i) & 1) { // if completely local, i = 0 (means "if (mat_idx == 1)" )
                if (state[flip_poses[i]] == '0') {
                    state[flip_poses[i]] = '1';
                }
                else {
                    state[flip_poses[i]] = '0';
                }
            }
        }
        return state;
    }

    void MooneyEtal_Filter::apply(Args args) {
        
        map<string, int> hist = args.hist;
        double threshold = args.threshold;
        clock_t t_start = clock();

        /*------------ preprocess ------------*/

        this->_shots = 0;
        for (const auto& item: hist) {
            this->_shots += item.second;
        }
        map<string, double> prob_dist;
        for (auto& item: hist) {
            prob_dist.insert(make_pair(item.first, (double)item.second / this->_shots));
        }
        clock_t t_prep = clock();

        /*------------ inverse operation ------------*/

        for (size_t i = 0; i < this->_pinv_matrices.size(); i++) { // O(n)
            map<string, double> x;
            for (const auto& item: prob_dist) { // O(1/t)
                int first_index = this->index_of_matrix(item.first, this->_poses_clbits[i]);
                double sum_of_count = 0;
                for (int k = 0; k < this->_pinv_matrices[i].rows(); k++) { // if completely local, k = 0 or 1
                    string source_state = this->flip_state(item.first, k, this->_poses_clbits[i]);
                    if (prob_dist.count(source_state) > 0) { // if the key exists
                        int second_index = this->index_of_matrix(source_state, this->_poses_clbits[i]);
                        sum_of_count += this->_pinv_matrices[i](first_index, second_index) * prob_dist[source_state];
                    }
                }
                if (abs(sum_of_count) >= threshold) {
                    x[item.first] = sum_of_count;
                }
            }
            map<string, double>().swap(prob_dist);
            prob_dist = x;
        }
        this->_dim = prob_dist.size();

        clock_t t_inv = clock();

        /*------------ sgs algorithm ------------*/

        this->_sum_of_x = 0;
        this->_indices_to_keys_vector = vector<string>(prob_dist.size());
        this->_x_s = vector<double>(prob_dist.size());
        int i = 0;
        for (auto& item: prob_dist) {
            this->_indices_to_keys_vector[i] = item.first;
            this->_x_s[i] = item.second;
            this->_sum_of_x += item.second;
            i++;
        }

        this->_x_tilde = sgs_algorithm(this->_x_s, false);

        clock_t t_sgs = clock();

        /*------------ recovering histogram ------------*/

        this->_sum_of_x_tilde = 0;
        this->_mitigated_hist.clear();
        for (size_t i = 0; i < this->_indices_to_keys_vector.size(); i++) {
            if (this->_x_tilde[i] != 0) {
                this->_sum_of_x_tilde += this->_x_tilde[i];
                this->_mitigated_hist.insert(make_pair(this->_indices_to_keys_vector[i], this->_x_tilde[i] * this->_shots));
            }
        }
        clock_t t_finish = clock();

        this->_durations.insert(make_pair("preprocess", (double)(t_prep - t_start) / CLOCKS_PER_SEC));
        this->_durations.insert(make_pair("inverse", (double)(t_inv - t_prep) / CLOCKS_PER_SEC));
        this->_durations.insert(make_pair("sgs_algorithm", (double)(t_sgs - t_inv) / CLOCKS_PER_SEC));
        this->_durations.insert(make_pair("postprocess", (double)(t_finish - t_sgs) / CLOCKS_PER_SEC));
        this->_durations.insert(make_pair("total", (double)(t_finish - t_start) / CLOCKS_PER_SEC));
        
        return;
    }
}