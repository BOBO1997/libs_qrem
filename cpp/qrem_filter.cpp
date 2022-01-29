#include <iostream>
#include <set>
#include <map>
#include <vector>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/SVD>
#include <algorithm>
#include <cmath>
#include <ctime>
#include <thread>

#include "eigen_utils.hpp"
#include "hamming.hpp"
#include "harger_higham.hpp"
#include "qrem_filter.hpp"
#include "sgs_algorithm.hpp"

using namespace std;
using namespace Eigen;

namespace libs_qrem {

    QREM_Filter::QREM_Filter(int num_clbits,
                             vector< vector< vector<double> > > cal_matrices,
                             vector< vector<int> > mit_pattern,
                             vector<int> meas_layout) {
        
        this->_num_clbits = num_clbits;

        // convert vector obj to Matrix obj
        this->_cal_matrices = vector<MatrixXd>(cal_matrices.size());
        for (size_t i = 0; i < cal_matrices.size(); i++) {
            this->_cal_matrices[i] = stdvec2d_to_MatrixXd(cal_matrices[i]);
        }
        
        // inverse of each matrix
        this->_pinv_matrices = vector<MatrixXd>(this->_cal_matrices.size());
        for (size_t i = 0; i < this->_cal_matrices.size(); i++) {
            this->_pinv_matrices[i] = this->_cal_matrices[i].inverse();
        }

        // svd of each matrix
        this->_svd_matrices = vector< JacobiSVD<MatrixXd> >(this->_cal_matrices.size());
        this->_Us = vector<MatrixXd>(this->_cal_matrices.size());
        this->_Sigmas = vector<MatrixXd>(this->_cal_matrices.size());
        this->_Vs = vector<MatrixXd>(this->_cal_matrices.size());
        for (size_t i = 0; i < this->_cal_matrices.size(); i++) {
            JacobiSVD<MatrixXd> svd(this->_cal_matrices[i], Eigen::ComputeFullU | Eigen::ComputeFullV);
            this->_svd_matrices[i] = svd;
            this->_Us[i] = svd.matrixU();
            this->_Sigmas[i] = svd.singularValues().asDiagonal();
            this->_Vs[i] = svd.matrixV();
        }

        // inverse of each svd matrices
        this->_pinvUs = vector<MatrixXd>(this->_cal_matrices.size());
        this->_pinvSigmas = vector<MatrixXd>(this->_cal_matrices.size()); 
        this->_pinvVs = vector<MatrixXd>(this->_cal_matrices.size());
        for (size_t i = 0; i < this->_svd_matrices.size(); i++) {
            this->_pinvUs[i] = this->_svd_matrices[i].matrixU().inverse();
            this->_pinvSigmas[i] = this->_svd_matrices[i].singularValues().asDiagonal().inverse();
            this->_pinvVs[i] = this->_svd_matrices[i].matrixV().inverse();
        }

        // set mit_pattern: default = [[0], [1], [2], [3], ...]
        // related to the index of calibration matrix
        if (mit_pattern.size() == 0) {
            for (int i = 0; i < num_clbits; i++) {
                mit_pattern.push_back(vector<int>(1, i));
            }
        }
        this->_mit_pattern = mit_pattern;

        // set meas_layout: default = [n-1, n-2, ..., 2, 1, 0]
        // related to the endian of measurement results
        // since the measurement results in qiskit takes big-endian, the indices in meas_layout are sorted in the descending order
        // meas_layout = clbits_to_qubits
        if (meas_layout.size() == 0) {
            for (int i = 0; i < num_clbits; i++) {
                meas_layout.push_back(num_clbits - 1 - i);
            }
        }
        this->_meas_layout = meas_layout;
        
        // set qubits_to_clbits: [-1, -1, ..., -1]
        // clbit = qubits_to_clbits[qubit]
        this->_qubits_to_clbits = vector<int>(1 + *max_element(this->_meas_layout.begin(), this->_meas_layout.end()), -1);
        for (size_t i = 0; i < this->_meas_layout.size(); i++) {
            this->_qubits_to_clbits[this->_meas_layout[i]] = i;
        }

        // indicate the qubits block for each calibration matrix
        // set poses_clbits: 
        // [[0],     A_0
        //  [1],     A_1
        //  [2],     A_2
        // ...]      ...
        this->_poses_clbits = vector< vector<int> >(this->_pinv_matrices.size());
        for (size_t i = 0; i < this->_pinv_matrices.size(); i++) {
            vector<int> pos_clbits(this->_mit_pattern[i].size());
            for (size_t j = 0; j < this->_mit_pattern[i].size(); j++) {
                pos_clbits[j] = this->_qubits_to_clbits[ this->_mit_pattern[i][j] ];
            }
            this->_poses_clbits[i] = pos_clbits;
        }
    }

    QREM_Filter::~QREM_Filter() {}

    // Ready for general calibration matrix blocks
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

    int QREM_Filter::index_of_matrix(int state, vector<int>& pos_clbits) {
        int index = 0;
        int i = 0;
        for (const auto& pos: pos_clbits) {
            if (state >> (this->_num_clbits - 1 - pos) & 1) {
                index += 1 << i;
            }
            i++;
        }
        return index;
    }

    // generalized to the arbitrary qubit blocks: essence: this->_indices_of_matrices
    void QREM_Filter::compute_reduced_A(size_t size) {
        this->_reduced_A = vector< vector<double> >(size, vector<double>(size, 0));
        thread th([&] {
            for (size_t i = 0; i < size; i++) { // target
                for (size_t j = 0; j < size; j++) { // source
                    double tensor_elem = 1;
                    for (size_t k = 0; k < this->_cal_matrices.size(); k++) {
                        int first_index = this->_indices_of_matrices[i][k]; // target
                        int second_index = this->_indices_of_matrices[j][k]; // source
                        tensor_elem *= this->_cal_matrices[k](first_index, second_index);
                    }
                    this->_reduced_A[i][j] = tensor_elem; // original: y[target] = A[target][source] * x[source]
                }
            }
        });
        th.join();
    }

    // generalized to the arbitrary qubit blocks: essence: this->_indices_of_matrices
    void QREM_Filter::compute_reduced_inv_A(size_t size) {
        this->_reduced_inv_A = vector< vector<double> >(size, vector<double>(size, 0));
        vector<double> abs_sum_of_rows(size, 0);
        // thread th([&] {
            for (size_t i = 0; i < size; i++) { // target
                for (size_t j = 0; j < size; j++) { // source
                    double tensor_elem = 1;
                    for (size_t k = 0; k < this->_pinv_matrices.size(); k++) {
                        int first_index = this->_indices_of_matrices[i][k]; // target
                        int second_index = this->_indices_of_matrices[j][k]; // source
                        tensor_elem *= this->_pinv_matrices[k](first_index, second_index);
                    }
                    this->_reduced_inv_A[i][j] = tensor_elem; // original: x[target] = A^-1[target][source] * y[source]
                    abs_sum_of_rows[i] += abs(tensor_elem);
                }
            }
        // });
        // th.join();
        this->_exact_one_norm_of_reduced_inv_A = 0;
        for (size_t i = 0; i < size; i++) {
            if (this->_exact_one_norm_of_reduced_inv_A < abs_sum_of_rows[i]) {
                this->_exact_one_norm_of_reduced_inv_A = abs_sum_of_rows[i];
            }
        }
    }

    void QREM_Filter::exact_one_norm_of_inv_reduced_A() {
        this->_exact_one_norm_of_inv_reduced_A = compute_one_norm_of_MatrixXd(stdvec2d_to_MatrixXd(this->_reduced_A).inverse());
    }

    void QREM_Filter::exact_one_norm_of_reduced_inv_A() {
        this->compute_reduced_inv_A(this->_indices_to_keys_vector.size());
    }

    void QREM_Filter::iterative_one_norm_of_inv_reduced_A(string method) {
        if (method == "iterative" | method == "bicgstab") {
            this->_iterative_one_norm_of_inv_reduced_A = harger_higham_bicgstab(stdvec2d_to_MatrixXd(this->_reduced_A));
        }
        else {
            this->_iterative_one_norm_of_inv_reduced_A = harger_higham_lu(stdvec2d_to_MatrixXd(this->_reduced_A));
        }
    }

    // generalized to the arbitrary qubit blocks: essence: this->_indices_of_matrices
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

    // generalized to the arbitrary qubit blocks: essence: this->_indices_of_matrices
    vector<double> QREM_Filter::col_basis(int col_index, 
                                          vector<MatrixXd>& pinv_mats, 
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
    
    vector<VectorXd> QREM_Filter::choose_vecs(string state, vector<MatrixXd> matrices) {
        vector<VectorXd> vecs(matrices.size());
        for (size_t i = 0; i < matrices.size(); i++) {
            vecs[i] = matrices[i].row(this->index_of_matrix(state, this->_poses_clbits[i]));
        }
        return vecs;
    }

    double QREM_Filter::sum_of_tensored_vector(vector<VectorXd> vecs) {
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
        this->_dim = this->_indices_to_keys_vector.size();
        int i = 0;
        for (const auto& key: extended_keys) {
            this->_indices_to_keys_vector[i] = key;
            i++;
        }

        // set indices_of_matrices table
        // generalized to the arbitrary qubit blocks
        this->_indices_of_matrices = vector< vector<int> >(this->_indices_to_keys_vector.size(), vector<int>(this->_pinv_matrices.size(), 0));
        for (size_t source_index = 0; source_index < this->_indices_to_keys_vector.size(); source_index++) {
            string source_state = this->_indices_to_keys_vector[source_index];
            for (size_t i = 0; i < this->_pinv_matrices.size(); i++) {
                int matrix_index = this->index_of_matrix(source_state, this->_poses_clbits[i]);
                this->_indices_of_matrices[source_index][i] = matrix_index;
            }
        }
        return extend_vector(prob_dist, this->_indices_to_keys_vector); // extended_y
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

    //! unncecessary to be included in the class method
    //! only used in ignis_filter.cpp and mooney_etal_filter.cpp
    string QREM_Filter::btos(int target_int, int n) {
        string s;
        for (int i = 0; i < n; i++) {
            s += to_string( (target_int >> (n - 1 - i)) & 1 );
        }
        return s;
    }

}