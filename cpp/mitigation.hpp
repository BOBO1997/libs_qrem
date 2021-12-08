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
    class QREM_Filter {
        public:

            int _num_clbits;
            vector<Matrix2d> _cal_matrices;
            // vector< vector< vector <double> > > _cal_matrices;
            vector< vector<int> > _mit_pattern;
            vector<int> _meas_layout;
            
            vector< JacobiSVD< Matrix<double, 2, 2> > > _svd_matrices;
            vector<Matrix2d> _Us;
            vector<Matrix2d> _Sigmas;
            vector<Matrix2d> _Vs;

            vector<Matrix2d> _pinv_matrices;
            vector<Matrix2d> _pinvUs; // unused
            vector<Matrix2d> _pinvSigmas;
            vector<Matrix2d> _pinvVs;

            vector<int> _qubits_to_clbits;

            // MatrixXd A_tilde;
            
            map<string, double> _durations;
            vector< vector<double> > _mitigatd_hists;

            QREM_Filter(int num_clbits,
                        vector< vector< vector<double> > > cal_matrices,
                        vector< vector<int> > mit_pattern,
                        vector<int> meas_layout);

            int index_of_matrix(string state, 
                                vector<int> pos_clbits);
            double mitigate_one_state(string target_state, 
                                      vector<double> extended_hist, 
                                      map<string, int> keys_to_indices);

            vector<double> col_basis(string col_state, 
                                     set<string> labels, 
                                     vector<Matrix2d> pinv_mats,
                                     map<string, int> keys_to_indices);
            // vector<double> col_basis(int col_idx, set<string> labels, vector< vector< vector<double> > > pinv_mats);
            vector<Vector2d> choose_vecs(string state, 
                                         vector<Matrix2d> matrices);
            // vector< vector<double> > choose_vecs(int state_idx, vector< vector< vector<double> > > matrices);
            double sum_of_tensored_vector(vector<Vector2d> vecs);
            // double sum_of_tensored_vector(vector< vector<double> > vecs);

            map<string, double> apply(map<string, int> hist, 
                                      int d, 
                                      double threshold);
    };
}