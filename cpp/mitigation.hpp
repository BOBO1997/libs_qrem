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
            vector< vector<int> > _mit_pattern;
            vector<int> _meas_layout;
            
            vector< JacobiSVD<Matrix2d> > _svd_matrices;
            vector<Matrix2d> _Us;
            vector<Matrix2d> _Sigmas;
            vector<Matrix2d> _Vs;

            vector<Matrix2d> _pinv_matrices;
            vector<Matrix2d> _pinvUs; // unused
            vector<Matrix2d> _pinvSigmas;
            vector<Matrix2d> _pinvVs;

            vector<int> _qubits_to_clbits;
            vector< vector<int> > _poses_clbits;

            vector< vector<int> > _indices_of_matrices;

            map<string, double> _durations;
            map<string, double> _mitigated_hist;

            QREM_Filter(int num_clbits,
                        vector< vector< vector<double> > > cal_matrices,
                        vector< vector<int> > mit_pattern,
                        vector<int> meas_layout);

            int index_of_matrix(string state, 
                                vector<int>& pos_clbits);
            double mitigate_one_state(int target_index, 
                                      vector<double>& extended_hist, 
                                      vector<string>& indices_to_keys_vector);

            vector<double> col_basis(int col_index, 
                                     vector<Matrix2d>& pinv_mats,
                                     vector<string>& indices_to_keys_vector);
            vector<Vector2d> choose_vecs(string state, 
                                         vector<Matrix2d> matrices);
            double sum_of_tensored_vector(vector<Vector2d> vecs);

            map<string, double> apply(map<string, int> hist, 
                                      int d, 
                                      double threshold);
    };
}