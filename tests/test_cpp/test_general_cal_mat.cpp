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
#include <fstream>
#include <iomanip>

#include "../../cpp/sgs_algorithm.hpp"
#include "../../cpp/qrem_filter.hpp"
#include "../../cpp/ignis_filter.hpp"
#include "../../cpp/delta_filter.hpp"
#include "../../cpp/least_norm_filter.hpp"
#include "../../cpp/mooney_etal_filter.hpp"
#include "../../cpp/nation_etal_filter.hpp"
#include "./dummy_data12.hpp"
// #include "./dummy_data50.hpp"

using namespace std;
using namespace Eigen;
using namespace libs_qrem;

vector < vector <double> > mat_prod(vector< vector <double> > a, vector< vector <double> > b) {
    vector< vector <double> > ab(a.size(), vector< double >(b[0].size()));
    for (size_t i = 0; i < a.size(); i++) {
        for (size_t j = 0; j < b[0].size(); j++) {
            for (size_t k = 0; j < a[0].size(); k++) {
                ab[i][j] = a[i][k] * b[k][j];
            }
        }
    }
    return ab;
}

vector < vector <double> > mat_tensor(vector< vector <double> > a, vector< vector <double> > b) {
    vector< vector <double> > ab(a.size() * b.size(), vector< double >(a[0].size() * b[0].size()));
    for (size_t i = 0; i < a.size(); i++) {
        for (size_t j = 0; j < a[0].size(); j++) {
            for (size_t ii = 0; ii < b.size(); ii++) {
                for (size_t jj = 0; jj < b[0].size(); jj++) {
                    ab[i * b.size() + ii][j * b[0].size() + jj] = a[i][j] * b[ii][jj];
                }
            }
        }
    }
    return ab;
}


int main() {
    
    int n = 12;
    // int n = 50;
    vector< vector< vector<double> > > local_cal_matrices = make_cal_matrices12();
    // vector< vector< vector<double> > > cal_matrices = make_cal_matrices50();

    vector< vector<int> > mit_pattern;
    mit_pattern.push_back({0,1,2});
    mit_pattern.push_back({3,4});
    mit_pattern.push_back({5});
    mit_pattern.push_back({6});
    mit_pattern.push_back({7,8,9});
    mit_pattern.push_back({10});
    mit_pattern.push_back({11});

    vector< vector< vector<double> > > cal_matrices;
    cal_matrices.push_back(mat_tensor(local_cal_matrices[0], mat_tensor(local_cal_matrices[1], local_cal_matrices[2])));
    cal_matrices.push_back(mat_tensor(local_cal_matrices[3], local_cal_matrices[4]));
    cal_matrices.push_back(local_cal_matrices[5]);
    cal_matrices.push_back(local_cal_matrices[6]);
    cal_matrices.push_back(mat_tensor(local_cal_matrices[7], mat_tensor(local_cal_matrices[8], local_cal_matrices[9])));
    cal_matrices.push_back(local_cal_matrices[10]);
    cal_matrices.push_back(local_cal_matrices[11]);

    cout << "length of cal_matrices: " << cal_matrices.size() << endl;
    for (auto const& mat: cal_matrices) {
        for (auto const& row: mat) {
            for (auto const & elem: row) {
                cout << right << setw(4) << setprecision(2) << elem << " ";
            }
            cout << endl;
        }
        cout << endl;
    }

    vector<int> meas_layout(n, 0);
    for (int i = 0; i < n; i++) meas_layout[i] = i;

    QREM_Filter* qf = new Least_Norm_Filter(n, cal_matrices, mit_pattern, meas_layout);
    // QREM_Filter* qf = new Ignis_Filter(n, cal_matrices, mit_pattern, meas_layout);
    // QREM_Filter* qf = new NationEtal_Filter(n, cal_matrices, mit_pattern, meas_layout);
    // QREM_Filter* qf = new MooneyEtal_Filter(n, cal_matrices, mit_pattern, meas_layout);

    map<string, int> hist = make_hist12();
    
    /*
    map<string, int> hist;
    ifstream ifs("./dummy_hist50.txt");
    string one_line;
    while (getline(ifs, one_line)) {
        hist.insert(make_pair(one_line.substr(0,50), stoi(one_line.substr(50))));
    }
    */

    Args args;
    args.hist = hist;
    args.d = 0;
    args.method = "lu";
    // args.method = "iterative";
    args.threshold = 0.1;
    // args.threshold = 0.01;


    cout << "apply correction" << endl;
    qf->apply(args);
    
    map<string, double> mitigated_hist = qf->_mitigated_hist;

    cout << "mitigated hist = {";
    for (auto& item: mitigated_hist) {
        cout << "\"" << item.first << "\", " << item.second << endl;
    }
    cout << "}" << endl;

    /*
    cout << "x_s = {";
    for (auto& item: qf->_x_s) {
        cout << item << endl;
    }
    cout << "}" << endl;
    */
    
    cout << "shots: " << qf->_shots << endl;

    for (const auto& t: qf->_durations) {
        cout << t.first << ": "  << t.second << " sec " << endl;
    }

    // cout << qf->_iterative_one_norm_of_inv_reduced_A << endl;
    // qf->iterative_one_norm_of_inv_reduced_A();
    cout << qf->_iterative_one_norm_of_inv_reduced_A << endl;

    cout << "sum of x s " << qf->_sum_of_x << endl;
    cout << "sum of x hat " << qf->_sum_of_x_hat << endl;
    cout << "sum of x tilde " << qf->_sum_of_x_tilde << endl;

    double sum_of_x_hat = 0;
    for (size_t i = 0; i < qf->_x_hat.size(); i++) {
        sum_of_x_hat += qf->_x_hat[i];
    }

    cout << "computed sum of x hat " << sum_of_x_hat << endl;

    vector<double> x_tilde = sgs_algorithm(qf->_x_hat, false);
    // cout << x_tilde << endl;
}
