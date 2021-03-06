#ifndef QREM_FILTER_BASE_HPP
#define QREM_FILTER_BASE_HPP

#include <set>
#include <map>
#include <vector>
#include <string>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/SVD>
#include <time.h>

using namespace std;
using namespace Eigen;

namespace libs_qrem {

    // typedef chrono::system_clock::time_point tp_now;

    struct Args {
        map<string, int> hist;
        int d;
        double threshold;
        string method;
        int step;
    };

    /*
    The base class of all qrem filters
    */
    class QREM_Filter {
        public:

            int _num_clbits;
            vector<MatrixXd> _cal_matrices;
            vector< vector<int> > _mit_pattern;
            vector<int> _meas_layout;
            
            vector< JacobiSVD<MatrixXd> > _svd_matrices;
            vector<MatrixXd> _Us;
            vector<MatrixXd> _Sigmas;
            vector<MatrixXd> _Vs;

            vector<MatrixXd> _pinv_matrices;
            vector<MatrixXd> _pinvUs; // unused
            vector<MatrixXd> _pinvSigmas;
            vector<MatrixXd> _pinvVs;

            int _shots;
            size_t _dim;
            vector< vector<double> > _reduced_A;
            vector< vector<double> > _reduced_inv_A; // For proposed methods
            vector< vector<double> > _inv_reduced_A; // For Nation et al.
            double _exact_one_norm_of_reduced_inv_A = -1; // For proposed methods
            double _exact_one_norm_of_inv_reduced_A = -1; // For Nation et al.
            double _iterative_one_norm_of_inv_reduced_A = -1; // For Nation et al.

            vector<int> _qubits_to_clbits;
            vector< vector<int> > _poses_clbits;
            vector< vector<int> > _indices_of_matrices;
            vector<string> _indices_to_keys_vector;

            double _sum_of_x;
            double _sum_of_x_hat;
            double _sum_of_x_tilde;

            map<string, double> _durations;
            vector<double> _x_s;
            vector<double> _x_hat;
            vector<double> _x_tilde;
            map<string, double> _mitigated_hist;

            QREM_Filter(int num_clbits,
                        vector< vector< vector<double> > > cal_matrices,
                        vector< vector<int> > mit_pattern = vector< vector<int> >(0),
                        vector<int> meas_layout = vector<int>(0));

            virtual ~QREM_Filter();

            int index_of_matrix(string state, 
                                vector<int>& pos_clbits);

            int index_of_matrix(int state, 
                                vector<int>& pos_clbits);

            void compute_reduced_A(size_t size);

            void compute_reduced_inv_A(size_t size);

            void exact_one_norm_of_inv_reduced_A();
            void exact_one_norm_of_reduced_inv_A();
            void iterative_one_norm_of_inv_reduced_A(string method = "iterative");

            double mitigate_one_state(int target_index, 
                                      vector<double>& extended_hist, 
                                      vector<string>& indices_to_keys_vector);

            vector<double> col_basis(int col_index, 
                                     vector<MatrixXd>& pinv_mats,
                                     vector<string>& indices_to_keys_vector);

            vector<VectorXd> choose_vecs(string state, 
                                         vector<MatrixXd> matrices);

            double sum_of_tensored_vector(vector<VectorXd> vecs);

            vector<double> preprocess(map<string, int> hist, 
                                      int d);

            void recover_histogram();

            string btos(int target_int, int n);

            virtual void apply(Args args) = 0;

    };
}

#endif