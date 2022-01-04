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
        for (size_t i = 0; i < flip_poses.size(); i++) {
            if ((mat_idx >> i) & 1) {
                state_idx = state_idx ^ ((1 << (this->_num_clbits - 1)) >> flip_poses[i]);
            }
        }
        return state_idx;
    }

    void MooneyEtal_Filter::apply(Args args) {
        
        map<string, int> hist = args.hist;
        double threshold = args.threshold;

        int shots = 0;
        for (const auto& item: hist) {
            shots += item.second;
        }

        clock_t t_start = clock();

        /*------------ preprocess ------------*/

        set<string> keys;
        map<string, double> prob_dist;
        for (auto& item: hist) {
            keys.insert(item.first);
            prob_dist.insert(make_pair(item.first, (double)item.second / shots));
        }

        clock_t t_prep = clock();

        /*------------ inverse operation ------------*/

        this->_sum_of_x = 0;
        for (size_t i = 0; i < this->_pinv_matrices.size(); i++) {
            map<string, double> x;
            for (const auto& key: keys) {
                int first_index = this->index_of_matrix(key, this->_poses_clbits[i]);
                double sum_of_count = 0;
                for (size_t k = 0; k < (size_t)this->_pinv_matrices[i].rows(); k++) {
                    string source_state = this->btos( this->flip_state(stoi(key, nullptr, 2), k, this->_poses_clbits[i]), this->_num_clbits);
                    if (prob_dist.count(source_state) > 0) {
                        int second_index = this->index_of_matrix(source_state, this->_poses_clbits[i]);
                        sum_of_count += this->_pinv_matrices[i](first_index, second_index) * prob_dist[source_state];
                    }
                }
                if (abs(sum_of_count) >= threshold) {
                    x[key] = sum_of_count;
                    this->_sum_of_x += sum_of_count;
                }
            }
            map<string, double>().swap(prob_dist);
            prob_dist = x;
        }

        this->_dim = prob_dist.size();
        if (this->_sum_of_x < 0) {
            cout << "negative counts" << endl;
        }

        clock_t t_inv = clock();

        /*------------ sgs algorithm ------------*/

        this->_indices_to_keys_vector = vector<string>(prob_dist.size());
        this->_x_s = vector<double>(prob_dist.size());
        int i = 0;
        for (auto& item: prob_dist) {
            this->_indices_to_keys_vector[i] = item.first;
            this->_x_s[i] = item.second;
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
                this->_mitigated_hist.insert(make_pair(this->_indices_to_keys_vector[i], this->_x_tilde[i] * shots));
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