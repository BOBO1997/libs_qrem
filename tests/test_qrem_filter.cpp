#include <iostream>

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

#include "../cpp/qrem_filter.hpp"
// #include "../cpp/qrem_filter_mooney_etal.hpp"

using namespace std;
using namespace Eigen;
using namespace libs_qrem;

int main() {
    int n = 3;
    vector< vector< vector<double> > > cal_matrices(n, vector< vector<double> >(2, vector<double>(2, 1)));
    cal_matrices[0][0][0] = 0.9;
    cal_matrices[0][0][1] = 0.1;
    cal_matrices[0][1][0] = 0.2;
    cal_matrices[0][1][1] = 0.8;
    cal_matrices[1][0][0] = 0.9;
    cal_matrices[1][0][1] = 0.1;
    cal_matrices[1][1][0] = 0.2;
    cal_matrices[1][1][1] = 0.8;
    cal_matrices[2][0][0] = 0.9;
    cal_matrices[2][0][1] = 0.1;
    cal_matrices[2][1][0] = 0.2;
    cal_matrices[2][1][1] = 0.8;
    vector< vector<int> > mit_pattern(n, vector<int>(1, 0));
    mit_pattern[0][0] = 0;
    mit_pattern[1][0] = 1;
    mit_pattern[2][0] = 2;
    vector<int> meas_layout(n, 0);
    meas_layout[0] = 0;
    meas_layout[1] = 1;
    meas_layout[2] = 2;
    QREM_Filter qf(n, cal_matrices, mit_pattern, meas_layout);
    // QREM_Filter_MooneyEtal qf(n, cal_matrices, mit_pattern, meas_layout);
    // QREM_Filter qf;
    map<string, int> hist;
    hist.insert(make_pair("000", 40));
    hist.insert(make_pair("001", 10));
    hist.insert(make_pair("010", 10));
    hist.insert(make_pair("111", 40));
    vector<int> pos_clbits(2);
    pos_clbits[0] = 1;
    pos_clbits[1] = 3;
    //cout << qf.index_of_matrix("01101000", pos_clbits) << endl << endl;
    qf.apply(hist, 0, 0.1);
    map<string, double> mitigated_hist = qf._mitigated_hist;
    cout << "{";
    for (auto& item: mitigated_hist) {
        cout << "\"" << item.first << "\", " << item.second << endl;
    }
    cout << "}" << endl;
}
