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

#include "mitigation.hpp"
#include "hamming.hpp"
#include "sgs_algorithm.hpp"

using namespace std;
using namespace Eigen;

namespace libs_qrem {

    Matrix2d vector_to_matrix2d(vector< vector<double> > data) {
        Matrix2d eMatrix(data.size(), data[0].size());
        for (size_t i = 0; i < data.size(); ++i)
            eMatrix.row(i) = VectorXd::Map(&data[i][0], data[0].size());
        return eMatrix;
    }

     QREM_Filter::QREM_Filter(int num_clbits,
                              vector< vector< vector<double> > > cal_matrices,
                              vector< vector<int> > mit_pattern = vector< vector<int> >(0),
                              vector<int> meas_layout = vector<int>(0)) {
        
        this->_num_clbits = num_clbits;

        // convert vector obj to Matrix obj
        this->_cal_matrices = vector<Matrix2d>(cal_matrices.size());
        for (size_t i = 0; i < cal_matrices.size(); i++) {
            this->_cal_matrices[i] = vector_to_matrix2d(cal_matrices[i]);
        }
        
        // inverse of each matrix
        this->_pinv_matrices = vector<Matrix2d>(this->_cal_matrices.size());
        for (size_t i = 0; i < this->_cal_matrices.size(); i++) {
            this->_pinv_matrices[i] = this->_cal_matrices[i].inverse();
        }

        // svd of each matrix
        this->_svd_matrices = vector< JacobiSVD<Matrix2d> >(this->_cal_matrices.size());
        this->_Us = vector<Matrix2d>(this->_cal_matrices.size());
        this->_Sigmas = vector<Matrix2d>(this->_cal_matrices.size());
        this->_Vs = vector<Matrix2d>(this->_cal_matrices.size());
        for (size_t i = 0; i < this->_cal_matrices.size(); i++) {
            JacobiSVD<Matrix2d> svd(this->_cal_matrices[i], Eigen::ComputeFullU | Eigen::ComputeFullV);
            this->_svd_matrices[i] = svd;
            this->_Us[i] = svd.matrixU();
            this->_Sigmas[i] = svd.singularValues().asDiagonal();
            this->_Vs[i] = svd.matrixV();
        }

        // inverse of each svd matrices
        this->_pinvUs = vector<Matrix2d>(this->_cal_matrices.size());
        this->_pinvSigmas = vector<Matrix2d>(this->_cal_matrices.size()); 
        this->_pinvVs = vector<Matrix2d>(this->_cal_matrices.size());
        for (size_t i = 0; i < this->_svd_matrices.size(); i++) {
            this->_pinvUs[i] = this->_svd_matrices[i].matrixU().inverse();
            this->_pinvSigmas[i] = this->_svd_matrices[i].singularValues().asDiagonal().inverse();
            this->_pinvVs[i] = this->_svd_matrices[i].matrixV().inverse();
        }

        // set mit_pattern
        if (mit_pattern.size() == 0) {
            for (int i = 0; i < num_clbits; i++) {
                mit_pattern.push_back(vector<int>(1, i));
            }
        }
        this->_mit_pattern = mit_pattern;

        // set meas_layout
        if (meas_layout.size() == 0) {
            for (int i = 0; i < num_clbits; i++) {
                meas_layout.push_back(num_clbits - 1 - i);
            }
        }
        this->_meas_layout = meas_layout;
        
        // set qubits_to_clbits
        this->_qubits_to_clbits = vector<int>(1 + *max_element(this->_meas_layout.begin(), this->_meas_layout.end()), -1);
        for (size_t i = 0; i < this->_meas_layout.size(); i++) {
            this->_qubits_to_clbits[this->_meas_layout[i]] = i;
        }

        // set poses_clbits
        this->_poses_clbits = vector< vector<int> >(this->_pinv_matrices.size());
        for (size_t i = 0; i < this->_pinv_matrices.size(); i++) {
            vector<int> pos_clbits(this->_mit_pattern[i].size());
            for (size_t j = 0; j < this->_mit_pattern[i].size(); j++) {
                pos_clbits[j] = this->_qubits_to_clbits[ this->_mit_pattern[i][j] ];
            }
            this->_poses_clbits[i] = pos_clbits;
        }

    };

    int QREM_Filter::index_of_matrix(string state, vector<int>& pos_clbits) {
        int index = 0;
        int i = 0;
        for (const auto& pos: pos_clbits) {
            if (state[pos] == '1') {
                index += 1 << i;
            }
            i++;
        }
        return index;
    }

    double QREM_Filter::mitigate_one_state(int target_index, 
                                           vector<double>& extended_hist, 
                                           vector<string>& indices_to_keys_vector) {
        double new_count = 0;
        for (size_t source_index = 0; source_index < indices_to_keys_vector.size(); source_index++) {
            double tensor_elem = 1;
            for (size_t i = 0; i < this->_pinv_matrices.size(); i++) {
                int first_index = this->_indices_of_matrices[target_index][i];
                int second_index = this->_indices_of_matrices[source_index][i];
                tensor_elem *= this->_pinv_matrices[i](first_index, second_index);
            }
            new_count += tensor_elem * extended_hist[source_index]; 
        }
        return new_count;
    }

    vector<double> QREM_Filter::col_basis(int col_index, 
                                          vector<Matrix2d>& pinv_mats, 
                                          vector<string>& indices_to_keys_vector) {
        
        vector<double> col_i(indices_to_keys_vector.size());
        for (size_t source_index = 0; source_index < indices_to_keys_vector.size(); source_index++) {
            double tensor_elem = 1;
            for (size_t i = 0; i < pinv_mats.size(); i++) {
                int first_index = this->_indices_of_matrices[source_index][i];
                int second_index = this->_indices_of_matrices[col_index][i];
                tensor_elem *= pinv_mats[i](first_index, second_index);
            }
            col_i[source_index] = tensor_elem;
        }
        return col_i;
    }
    
    vector<Vector2d> QREM_Filter::choose_vecs(string state, vector<Matrix2d> matrices) {
        vector<Vector2d> vecs(matrices.size());
        for (size_t i = 0; i < matrices.size(); i++) {
            vecs[i] = matrices[i].row(this->index_of_matrix(state, this->_poses_clbits[i]));
        }
        return vecs;
    }

    double QREM_Filter::sum_of_tensored_vector(vector<Vector2d> vecs) {
        double sum_val = 1;
        for (const auto& vec: vecs) {
            sum_val *= vec.sum();
        }
        return sum_val;
    }

    map<string, double> QREM_Filter::apply(map<string, int> hist,
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

        map<string, int> keys_to_indices;
        map<int, string> indices_to_keys;
        vector<string> indices_to_keys_vector(extended_keys.size());
        int i = 0;
        for (const auto& key: extended_keys) {
            keys_to_indices.insert(make_pair(key, i));
            indices_to_keys.insert(make_pair(i, key));
            indices_to_keys_vector[i] = key;
            i++;
        }

        vector<double> extended_y = extend_vectors(prob_dist, keys_to_indices);

        this->_indices_of_matrices = vector< vector<int> >(indices_to_keys_vector.size(), vector<int>(this->_pinv_matrices.size(), 0));
        for (size_t source_index = 0; source_index < indices_to_keys_vector.size(); source_index++) {
            string source_state = indices_to_keys_vector[source_index];
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

        vector<double> x_s(extended_y.size(), 0);
        int sum_of_x = 0;
        for (size_t state_idx = 0; state_idx < extended_keys.size(); state_idx++) {
            double mitigated_value = mitigate_one_state(state_idx, extended_y, indices_to_keys_vector);
            x_s[state_idx] = mitigated_value;
            sum_of_x += mitigated_value;
        }

        // time for inverse operation
        chrono::system_clock::time_point t_inv = chrono::system_clock::now();
        double dur_inv = std::chrono::duration_cast<std::chrono::milliseconds>(t_inv - t_prep).count();
        this->_durations.insert(make_pair("inverse", dur_inv));

        /*------------ correction by delta ------------*/

        double sum_of_vi = this->sum_of_tensored_vector(this->choose_vecs(string(this->_num_clbits, '0'), this->_pinvVs));
        double lambda_i = this->sum_of_tensored_vector(this->choose_vecs(string(this->_num_clbits, '0'), this->_pinvSigmas));
        double delta_denom = (sum_of_vi * sum_of_vi) / (lambda_i * lambda_i);
        double delta_coeff = (1 - sum_of_x) / delta_denom;
        double delta_col = sum_of_vi / (lambda_i * lambda_i);
        vector<double> v_col = this->col_basis(keys_to_indices[string(this->_num_clbits, '0')], this->_pinvVs, indices_to_keys_vector);

        for (size_t state_idx = 0; state_idx < extended_keys.size(); state_idx++) {
            x_s[state_idx] += delta_coeff * delta_col * v_col[state_idx];
        }

        // time for correction by delta
        chrono::system_clock::time_point t_delta = chrono::system_clock::now();
        double dur_delta = std::chrono::duration_cast<std::chrono::milliseconds>(t_delta - t_inv).count();
        this->_durations.insert(make_pair("delta", dur_delta));

        /*------------ sgs algorithm ------------*/
        vector<double> x_tilde = sgs_algorithm(x_s);

        // time for sgs algorithm
        chrono::system_clock::time_point t_sgs = chrono::system_clock::now();
        double dur_sgs = std::chrono::duration_cast<std::chrono::milliseconds>(t_sgs - t_delta).count();
        this->_durations.insert(make_pair("sgs_algorithm", dur_sgs));

        /*------------ recovering histogram ------------*/
        i = 0;
        for (const auto& key: extended_keys) {
            this->_mitigated_hist.insert(make_pair(key, x_tilde[i] * shots));
            i++;
        }

        // time for postprocess
        chrono::system_clock::time_point t_finish = chrono::system_clock::now();
        double dur_postprocess = std::chrono::duration_cast<std::chrono::milliseconds>(t_finish - t_sgs).count();
        this->_durations.insert(make_pair("postprocess", dur_postprocess));
        
        // total time
        double dur_total = std::chrono::duration_cast<std::chrono::milliseconds>(t_finish - t_start).count();
        this->_durations.insert(make_pair("total", dur_total));
        
        return this->_mitigated_hist;
    }
}