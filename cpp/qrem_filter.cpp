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
#include <cmath>
#include <ctime>

#include "eigen_utils.hpp"
#include "qrem_filter_base.hpp"
#include "hamming.hpp"
#include "sgs_algorithm.hpp"

using namespace std;
using namespace Eigen;

namespace libs_qrem {

    QREM_Filter::QREM_Filter(int num_clbits,
                             vector< vector< vector<double> > > cal_matrices,
                             vector< vector<int> > mit_pattern = vector< vector<int> >(0),
                             vector<int> meas_layout = vector<int>(0)) {
        
        this->_num_clbits = num_clbits;

        // convert vector obj to Matrix obj
        this->_cal_matrices = vector<Matrix2d>(cal_matrices.size());
        for (size_t i = 0; i < cal_matrices.size(); i++) {
            this->_cal_matrices[i] = vector2d_to_matrix2d(cal_matrices[i]);
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

    void QREM_Filter::compute_reduced_A(vector<string>& indices_to_keys_vector) {
        this->_reduced_A = vector< vector<double> >(indices_to_keys_vector.size(), vector<double>(indices_to_keys_vector.size(), 0));
        for (size_t i = 0; i < indices_to_keys_vector.size(); i++) { // target
            for (size_t j = 0; j < indices_to_keys_vector.size(); j++) { // source
                double tensor_elem = 1;
                for (size_t k = 0; k < this->_cal_matrices.size(); k++) {
                    int first_index = this->_indices_of_matrices[i][k]; // target
                    int second_index = this->_indices_of_matrices[j][k]; // source
                    tensor_elem *= this->_cal_matrices[k](first_index, second_index);
                }
                this->_reduced_A[i][j] = tensor_elem; // original: y[target] = A[target][source] * x[source]
            }
        }
    }

    void QREM_Filter::compute_reduced_inv_A(vector<string>& indices_to_keys_vector) {
        this->_reduced_inv_A = vector< vector<double> >(indices_to_keys_vector.size(), vector<double>(indices_to_keys_vector.size(), 0));
        this->_max_element = 0;
        vector<double> abs_sum_of_rows(indices_to_keys_vector.size(), 0);
        for (size_t i = 0; i < indices_to_keys_vector.size(); i++) { // target
            for (size_t j = 0; j < indices_to_keys_vector.size(); j++) { // source
                double tensor_elem = 1;
                for (size_t k = 0; k < this->_pinv_matrices.size(); k++) {
                    int first_index = this->_indices_of_matrices[i][k]; // target
                    int second_index = this->_indices_of_matrices[j][k]; // source
                    tensor_elem *= this->_pinv_matrices[k](first_index, second_index);
                }
                this->_reduced_inv_A[i][j] = tensor_elem; // original: x[target] = A^-1[target][source] * y[source]
                abs_sum_of_rows[i] += abs(tensor_elem);
                if (tensor_elem > this->_max_element) {
                    this->_max_element = tensor_elem;
                }
            }
        }
        this->_one_norm = 0;
        for (size_t i = 0; i < abs_sum_of_rows.size(); i++) {
            if (this->_one_norm < abs_sum_of_rows[i]) {
                this->_one_norm = abs_sum_of_rows[i];
            }
        }
    }

    vector<double> QREM_Filter::mat_vec_prod(vector< vector<double> > A, vector<double> y) {
        // compute x = Ay
        vector<double> x(y.size(), 0);
        for (size_t i = 0; i < x.size(); i++) {
            for (size_t j = 0; j < x.size(); j++) {
                x[i] += A[i][j] * y[j];
            }
        }        
        return x;
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

    vector<double> QREM_Filter::preprocess(map<string, int> hist, int d) {
        this->_shots = 0;
        for (const auto& item: hist) {
            this->_shots += item.second;
        }

        set<string> keys;
        map<string, double> prob_dist;
        for (auto& item: hist) {
            keys.insert(item.first);
            prob_dist.insert(make_pair(item.first, (double)item.second / this->_shots));
        }

        set<string> extended_keys = extend_keys(keys, d);
        this->_indices_to_keys_vector = vector<string>(extended_keys.size());
        int i = 0;
        for (const auto& key: extended_keys) {
            this->_indices_to_keys_vector[i] = key;
            i++;
        }

        // set indices_of_matrices table
        this->_indices_of_matrices = vector< vector<int> >(this->_indices_to_keys_vector.size(), vector<int>(this->_pinv_matrices.size(), 0));
        for (size_t source_index = 0; source_index < this->_indices_to_keys_vector.size(); source_index++) {
            string source_state = this->_indices_to_keys_vector[source_index];
            for (size_t i = 0; i < this->_pinv_matrices.size(); i++) {
                int matrix_index = this->index_of_matrix(source_state, this->_poses_clbits[i]);
                this->_indices_of_matrices[source_index][i] = matrix_index;
            }
        }
        return extend_vectors(prob_dist, this->_indices_to_keys_vector); // extended_y
    }

    void QREM_Filter::recover_histogram() {
        this->_sum_of_x_tilde = 0;
        this->_mitigated_hist.clear();
        for (size_t i = 0; i < this->_indices_to_keys_vector.size(); i++) {
            if (this->_x_tilde[i] != 0) {
                this->_sum_of_x_tilde += this->_x_tilde[i];
                this->_mitigated_hist.insert(make_pair(this->_indices_to_keys_vector[i], this->_x_tilde[i] * this->_shots));
            }
        }
    }
}