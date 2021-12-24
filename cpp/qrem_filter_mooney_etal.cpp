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
#include <cstdlib>
#include <bitset>

#include "eigen_utils.hpp"
#include "qrem_filter_mooney_etal.hpp"
#include "hamming.hpp"
#include "sgs_algorithm.hpp"
#include "qrem_filter_base.hpp"

using namespace std;
using namespace Eigen;

namespace libs_qrem {

     QREM_Filter_MooneyEtal::QREM_Filter_MooneyEtal(int num_clbits,
                              vector< vector< vector<double> > > cal_matrices,
                              vector< vector<int> > mit_pattern = vector< vector<int> >(0),
                              vector<int> meas_layout = vector<int>(0)) : 
                              QREM_Filter_Base(num_clbits, cal_matrices, mit_pattern, meas_layout) {
        
    };

    string QREM_Filter_MooneyEtal::btos(int target_int, int n) {
        string s;
        for (int i = 0; i < n; i++) {
            s += to_string( (target_int >> (n - 1 - i)) & 1 );
        }
        return s;
    }

    int QREM_Filter_MooneyEtal::flip_state(int state_idx, int mat_idx, vector<int>& flip_poses) {
        for (size_t i = 0; i < flip_poses.size(); i++) {
            if ((mat_idx >> i) & 1) {
                state_idx = state_idx ^ ((1 << (this->_num_clbits - 1)) >> flip_poses[i]);
            }
        }
        return state_idx;
    }

    void QREM_Filter_MooneyEtal::apply(map<string, int> hist,
                                    int d = 0,
                                    double threshold = 0.1) {
        int shots = 0;
        for (const auto& item: hist) {
            shots += item.second;
        }

        tp_now t_start = chrono::system_clock::now();

        /*------------ preprocess ------------*/

        set<string> keys;
        map<string, double> prob_dist;
        for (auto& item: hist) {
            keys.insert(item.first);
            prob_dist.insert(make_pair(item.first, (double)item.second / shots));
        }

        // time for preprocess
        tp_now t_prep = chrono::system_clock::now();
        double dur_prep = chrono::duration_cast<chrono::milliseconds>(t_prep - t_start).count();
        this->_durations.insert(make_pair("preprocess", dur_prep));

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
            prob_dist = x;
        }
        if (this->_sum_of_x < 0) {
            cout << "negative counts" << endl;
        }

        // time for inverse operation
        tp_now t_inv = chrono::system_clock::now();
        this->_durations.insert(make_pair("inverse", chrono::duration_cast<chrono::milliseconds>(t_inv - t_prep).count()));

        /*------------ sgs algorithm ------------*/

        this->_indices_to_keys_vector = vector<string>(prob_dist.size());
        this->_x_s = vector<double>(prob_dist.size());
        int i = 0;
        for (auto& item: prob_dist) {
            this->_indices_to_keys_vector[i] = item.first;
            this->_x_s[i] = item.second;
            i++;
        }

        this->_x_tilde = sgs_algorithm(this->_x_s);

        // time for sgs algorithm
        tp_now t_sgs = chrono::system_clock::now();
        this->_durations.insert(make_pair("sgs_algorithm", chrono::duration_cast<chrono::milliseconds>(t_sgs - t_lnp).count()));

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
        tp_now t_finish = chrono::system_clock::now();
        this->_durations.insert(make_pair("postprocess", chrono::duration_cast<chrono::milliseconds>(t_finish - t_sgs).count()));
        this->_durations.insert(make_pair("total", chrono::duration_cast<chrono::milliseconds>(t_finish - t_start).count()));
        
        return;
    }
}