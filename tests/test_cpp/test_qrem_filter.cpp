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

#include "../../cpp/qrem_filter.hpp"
#include "../../cpp/delta_filter.hpp"
#include "../../cpp/least_norm_filter.hpp"
#include "../../cpp/mooney_etal_filter.hpp"
#include "../../cpp/nation_etal_filter.hpp"
#include "./dummy_data.hpp"

using namespace std;
using namespace Eigen;
using namespace libs_qrem;

int main() {
    
    int n = 12;
    vector< vector< vector<double> > > cal_matrices = make_cal_matrices();

    vector< vector<int> > mit_pattern(n, vector<int>(1, 0));
    for (int i = 0; i < n; i++) mit_pattern[i][0] = i;

    vector<int> meas_layout(n, 0);
    for (int i = 0; i < n; i++) meas_layout[i] = i;

    // QREM_Filter* qf = new Least_Norm_Filter(n, cal_matrices, mit_pattern, meas_layout);
    QREM_Filter* qf = new NationEtal_Filter(n, cal_matrices, mit_pattern, meas_layout);

    map<string, int> hist = make_hist();

    Args args;
    args.hist = hist;
    args.d = 0;
    args.method = "lu";
    // args.method = "iterative";

    qf->apply(args);
    
    map<string, double> mitigated_hist = qf->_mitigated_hist;

    /*
    cout << "mitigated hist = {";
    for (auto& item: mitigated_hist) {
        cout << "\"" << item.first << "\", " << item.second << endl;
    }
    cout << "}" << endl;
    */

    /*
    cout << "x_s = {";
    for (auto& item: qf->_x_s) {
        cout << item << endl;
    }
    cout << "}" << endl;
    */
    
    cout << "shots: " << qf->_shots << endl;

    cout << qf->_iterative_one_norm_of_inv_reduced_A << endl;
    qf->iterative_one_norm_of_inv_reduced_A();
    cout << qf->_iterative_one_norm_of_inv_reduced_A << endl;
}
