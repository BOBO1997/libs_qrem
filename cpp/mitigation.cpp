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
        for (int i = 0; i < (int)data.size(); ++i)
            eMatrix.row(i) = VectorXd::Map(&data[i][0], data[0].size());
        return eMatrix;
    }

     QREM_Filter::QREM_Filter(int num_clbits,
                              vector< vector< vector<double> > > cal_matrices,
                              vector< vector<int> > mit_pattern = vector< vector<int> >(0),
                              vector<int> meas_layout = vector<int>(0)) {
        
        this->_num_clbits = num_clbits;

        // convert vector obj to Matrix obj
        for (auto& cal_matrix: cal_matrices) {
            this->_cal_matrices.push_back(vector_to_matrix2d(cal_matrix));
        }
        
        // inverse of each matrix
        for (auto& matrix: this->_cal_matrices) {
            this->_pinv_matrices.push_back(matrix.inverse());
        }

        // svd of each matrix
        for (auto& matrix: this->_cal_matrices) {
            JacobiSVD<Matrix2d> svd(matrix, Eigen::ComputeFullU | Eigen::ComputeFullV);
            this->_svd_matrices.push_back(svd);
            this->_Us.push_back(svd.matrixU());
            this->_Sigmas.push_back(svd.singularValues().asDiagonal());
            this->_Vs.push_back(svd.matrixV());
        }

        // inverse of each svd matrices
        for (auto& svd: this->_svd_matrices) {
            this->_Us.push_back(svd.matrixU().inverse());
            this->_Sigmas.push_back(svd.singularValues().asDiagonal().inverse());
            this->_Vs.push_back(svd.matrixV().inverse());
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
        for (int i = 0; i < (int)this->_meas_layout.size(); i++) {
            this->_qubits_to_clbits[this->_meas_layout[i]] = i;
        }
    };

    int QREM_Filter::index_of_matrix(string state, vector<int> pos_clbits) {
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

    double QREM_Filter::mitigate_one_state(string target_state, 
                                           vector<double> extended_hist, 
                                           map<string, int> keys_to_indices) {
        double new_count = 0;
        for (const auto& item: keys_to_indices) {
            string source_state = item.first;
            int source_index = item.second;
            double tensor_elem = 1;
            for (size_t i = 0; i < this->_pinv_matrices.size(); i++) {
                vector<int> pos_clbits;
                for (const auto& qubit: this->_mit_pattern[i]) {
                    pos_clbits.push_back(this->_qubits_to_clbits[qubit]);
                }
                int first_index = this->index_of_matrix(target_state, pos_clbits);
                int second_index = this->index_of_matrix(source_state, pos_clbits);
                tensor_elem *= this->_pinv_matrices[i](first_index, second_index);
            }
            new_count += tensor_elem * extended_hist[source_index]; 
        }
        return new_count;
    }

    vector<double> QREM_Filter::col_basis(string col_state, 
                                          set<string> labels, 
                                          vector<Matrix2d> pinv_mats, 
                                          map<string, int> keys_to_indices) {
        
        vector<double> col_i(labels.size());
        for (const auto& label: labels) {
            double tensor_elem = 1;
            for (size_t i = 0; i < pinv_mats.size(); i++) {
                vector<int> pos_clbits;
                for (const auto& qubit: this->_mit_pattern[i]) {
                    pos_clbits.push_back(this->_qubits_to_clbits[qubit]);
                }
                int first_index = this->index_of_matrix(label, pos_clbits);
                int second_index = this->index_of_matrix(col_state, pos_clbits);
                tensor_elem *= pinv_mats[i](first_index, second_index);
            }
            col_i[keys_to_indices[label]] = tensor_elem;
        }
        return col_i;
    }
    
    vector<Vector2d> QREM_Filter::choose_vecs(string state, vector<Matrix2d> matrices) {
        vector<Vector2d> vecs(matrices.size());
        for (size_t i = 0; i < matrices.size(); i++) {
            vector<int> pos_clbits;
            for (const auto& qubit: this->_mit_pattern[i]) {
                pos_clbits.push_back(this->_qubits_to_clbits[qubit]);
            }
            vecs[i] = matrices[i].row(this->index_of_matrix(state, pos_clbits));
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

        // preprocess

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
        int i = 0;
        for (const auto& key: extended_keys) {
            keys_to_indices.insert(make_pair(key, i));
            indices_to_keys.insert(make_pair(i, key));
            i++;
        }

        vector<double> extended_y = extend_vectors(prob_dist, keys_to_indices);
        
        // time for preprocess
        chrono::system_clock::time_point t_prep = chrono::system_clock::now();
        double dur_prep = std::chrono::duration_cast<std::chrono::milliseconds>(t_prep - t_start).count();
        this->_durations.insert(make_pair("preprocess", dur_prep));

        // proposed method
        vector<double> x_s = vector<double>(extended_y.size(), 0);
        int sum_of_x = 0;
        for (int state_idx = 0; state_idx < (int)extended_keys.size(); state_idx++) {
            double mitigated_value = mitigate_one_state(extended_keys_vector[state_idx], extended_y, keys_to_indices);
            x_s[state_idx] = mitigated_value;
            sum_of_x += mitigated_value;
        }

        // time for inverse operation
        chrono::system_clock::time_point t_inv = chrono::system_clock::now();
        double dur_inv = std::chrono::duration_cast<std::chrono::milliseconds>(t_inv - t_prep).count();
        this->_durations.insert(make_pair("inverse", dur_inv));

        double sum_of_vi = this->sum_of_tensored_vector(this->choose_vecs(string(this->_num_clbits, '0'), this->_pinvVs));
        double lambda_i = this->sum_of_tensored_vector(this->choose_vecs(string(this->_num_clbits, '0'), this->_pinvSigmas));
        double delta_denom = (sum_of_vi * sum_of_vi) / (lambda_i * lambda_i);
        double delta_coeff = (1 - sum_of_x) / delta_denom;
        double delta_col = sum_of_vi / (lambda_i * lambda_i);
        vector<double> v_col = this->col_basis(string(this->_num_clbits, '0'), extended_keys, this->_pinvVs, keys_to_indices);

        for (int state_idx = 0; state_idx < (int)extended_keys.size(); state_idx++) {
            x_s[state_idx] += delta_coeff * delta_col * v_col[state_idx];
        }

        // time for correction by delta
        chrono::system_clock::time_point t_delta = chrono::system_clock::now();
        double dur_delta = std::chrono::duration_cast<std::chrono::milliseconds>(t_delta - t_inv).count();
        this->_durations.insert(make_pair("delta", dur_delta));

        // sgs algorithm
        vector<double> x_tilde = sgs_algorithm(x_s);

        // time for sgs algorithm
        chrono::system_clock::time_point t_sgs = chrono::system_clock::now();
        double dur_sgs = std::chrono::duration_cast<std::chrono::milliseconds>(t_sgs - t_delta).count();
        this->_durations.insert(make_pair("sgs_algorithm", dur_sgs));

        // recover the shots
        map<string, double> mitigated_hist;
        i = 0;
        for (const auto& key: extended_keys) {
            mitigated_hist.insert(make_pair(key, x_tilde[i] * shots));
            i++;
        }

        // time for postprocess
        chrono::system_clock::time_point t_finish = chrono::system_clock::now();
        double dur_postprocess = std::chrono::duration_cast<std::chrono::milliseconds>(t_finish - t_sgs).count();
        this->_durations.insert(make_pair("postprocess", dur_postprocess));
        
        // total time
        double dur_total = std::chrono::duration_cast<std::chrono::milliseconds>(t_finish - t_start).count();
        this->_durations.insert(make_pair("total", dur_total));
        
        return mitigated_hist;
    }
}