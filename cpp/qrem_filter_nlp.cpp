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
#include "qrem_filter_nlp.hpp"
#include "hamming.hpp"
#include "sgs_algorithm.hpp"
#include "qrem_filter_base.hpp"

using namespace std;
using namespace Eigen;

namespace libs_qrem {

     QREM_Filter_Nlp::QREM_Filter_Nlp(int num_clbits,
                              vector< vector< vector<double> > > cal_matrices,
                              vector< vector<int> > mit_pattern = vector< vector<int> >(0),
                              vector<int> meas_layout = vector<int>(0)) :
                              QREM_Filter_Base(num_clbits, cal_matrices, mit_pattern, meas_layout) {
        
    };

    void QREM_Filter_Nlp::apply(map<string, int> hist,
                                int d = 0) {

        chrono::system_clock::time_point t_start = chrono::system_clock::now();

        /*------------ preprocess ------------*/

        vector<double> extended_y = this->preprocess(hist, d);

        // time for preprocess
        chrono::system_clock::time_point t_prep = chrono::system_clock::now();
        double dur_prep = std::chrono::duration_cast<std::chrono::milliseconds>(t_prep - t_start).count();
        this->_durations.insert(make_pair("preprocess", dur_prep));

        /*------------ inverse operation ------------*/

        compute_reduced_inv_A(this->_indices_to_keys_vector);
        this->_x_s = this->mat_vec_prod(this->_reduced_inv_A, extended_y);
        this->_sum_of_x = 0;
        for (size_t i = 0; i < this->_x_s.size(); i++) {
            this->_sum_of_x += this->_x_s[i];
        }

        // time for inverse operation
        chrono::system_clock::time_point t_inv = chrono::system_clock::now();
        double dur_inv = std::chrono::duration_cast<std::chrono::milliseconds>(t_inv - t_prep).count();
        this->_durations.insert(make_pair("inverse", dur_inv));

        return;
    }
}