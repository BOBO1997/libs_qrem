#include <iostream>
#include <set>
#include <map>
#include <vector>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <algorithm>
#include <ctime>
#include <thread>

#include "eigen_utils.hpp"
#include "hamming.hpp"
#include "sgs_algorithm.hpp"
#include "qrem_filter.hpp"
#include "ignis_filter.hpp"

using namespace std;
using namespace Eigen;

namespace libs_qrem {

     Ignis_Filter::Ignis_Filter(int num_clbits,
                              vector< vector< vector<double> > > cal_matrices,
                              vector< vector<int> > mit_pattern = vector< vector<int> >(0),
                              vector<int> meas_layout = vector<int>(0)) : 
                              QREM_Filter(num_clbits, cal_matrices, mit_pattern, meas_layout) {};

    int Ignis_Filter::flip_state(int target_idx, int mat_index, vector<int>& flip_poses) {
        int source_idx = target_idx;
        for (size_t i = 0; i < flip_poses.size(); i++) {
            if ((mat_index >> i) & 1) {
                source_idx = source_idx ^ (1 << (this->_num_clbits - 1 - flip_poses[i]));
            }
        }
        return source_idx;
    }

    void Ignis_Filter::apply(Args args) {
        
        map<string, int> hist = args.hist;
        clock_t t_start = clock();

        /*------------ preprocess ------------*/

        this->_shots = 0;
        for (const auto& item: hist) {
            this->_shots += item.second;
        }
        this->_x_s = vector<double>(1 << this->_num_clbits, 0);
        for (const auto& item: hist) {
            this->_x_s[stoi(item.first, nullptr, 2)] = (double)item.second / this->_shots;
        }
        clock_t t_prep = clock();

        /*------------ inverse operation ------------*/

        for (size_t i = 0; i < this->_pinv_matrices.size(); i++) { // O(n)
            vector<double> temp_y(1 << this->_num_clbits, 0);
            // thread th([&]{
                for (int target_idx = 0; target_idx < (1 << this->_num_clbits); target_idx++) { // O(2^n)
                    int first_index = this->index_of_matrix(target_idx, this->_poses_clbits[i]);
                    for (int j = 0; j < this->_pinv_matrices[i].rows(); j++) { // for each index of inverse calibration matrix
                        int source_idx = this->flip_state(target_idx, j, this->_poses_clbits[i]);
                        int second_index = this->index_of_matrix(source_idx, this->_poses_clbits[i]);
                        temp_y[target_idx] += this->_pinv_matrices[i](first_index, second_index) * this->_x_s[source_idx];
                    }
                }
            // });
            // th.join();
            vector<double>().swap(this->_x_s); // release previous vector;
            this->_x_s = temp_y; // update the vector
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

        this->_sum_of_x_tilde = 0;
        this->_mitigated_hist.clear();
        for (size_t i = 0; i < (1 << this->_num_clbits); i++) {
            if (this->_x_tilde[i] != 0) {
                this->_sum_of_x_tilde += this->_x_tilde[i];
                this->_mitigated_hist.insert(make_pair(btos(i, this->_num_clbits), this->_x_tilde[i] * this->_shots));
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