#include <iostream>
#include <map>
#include <vector>
#include <Eigen/Dense>
#include <time.h>

#include "mitigation.hpp"

using namespace std;
using Eigen::MatrixXd;

namespace libs_qrem {

    // QREM_Filter::QREM_Filter(vector<MatrixXd> cal_matrices) {
    QREM_Filter::QREM_Filter(vector<double> cal_matrices) {
        // vector<MatrixXd> cal_matrices = cal_matrices;
        // cal_matrices = cal_matrices;

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
    };

    void QREM_Filter::apply(map<string, double> hist) {
        return;
    }
}