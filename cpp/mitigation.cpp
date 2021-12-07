#include <iostream>
#include <map>
#include <vector>
#include <Eigen/Dense>
#include <time.h>
#include <algorithm>

#include "mitigation.hpp"

using namespace std;
using namespace Eigen;

namespace libs_qrem {

     QREM_Filter::QREM_Filter(int num_clbits,
                              vector< vector< vector<double> > > cal_matrices,
                              vector< vector<int> > mit_pattern = vector< vector<int> >(0),
                              vector<int> meas_layout = vector<int>(0)) {
        
        for (auto& cal_matrix: cal_matrices) {
            Matrix2d mat;
            mat[0][0] = cal_matrix[0][0];
            mat[0][1] = cal_matrix[0][1];
            mat[1][0] = cal_matrix[1][0];
            mat[1][1] = cal_matrix[1][1]];
            _cal_matrices.push_back(mat);
        }
        
        for (auto& cal_matrix: _cal_matrices) {
            _pinv_matrices.push_back(cal_matrix.inverse());
        }


        if (meas_layout.size() == 0) {
            for (int i = 0; i < num_clbits; i++) {
                meas_layout.push_back(num_clbits - 1 - i)
            }
        }
        else {
            _meas_layout = meas_layout;
        }

        _qubits_to_clbits = vector<int>(1 + *max_element(_meas_layout.begin(), _meas_layout.end()), -1);
        for (int i = 0; i < _meas_layout.size(); i++) {
            _qubits_to_clbits[_meas_layout[i]] = i
        }

        /*
        vector<int> v;

        const int N = 10000 * 10000;

        auto t1 = std::chrono::system_clock::now();
        for (int i = 0; i < N; i++){ 
            v.push_back(i);
        }
        auto t2 = std::chrono::system_clock::now();
        auto sec = std::chrono::duration_cast<std::chrono::seconds>(t2 - t1).count();
        auto msec = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
        cout << sec << " s " << endl;
        cout << msec << " s " << endl;
        */
    };

    void QREM_Filter::apply(map<string, int> hist) {
        int shots = 0;
        for (const auto& item: hist) {
            shots += item.second
        }

        return;
    }
}