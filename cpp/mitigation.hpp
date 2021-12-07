#include <map>
#include <vector>
#include <Eigen/Dense>
#include <time.h>

using namespace std;
using namespace Eigen;

namespace libs_qrem {
    class QREM_Filter {
        public:

            int _num_clbits;
            vector<Matrix2d> _cal_matrices;
            vector<vector<int>> _mit_pattern;
            vector<int> _meas_layout;
            
            // vector<Matrix2d> _svd_matrices;
            // vector<Matrix2d> _Us;
            // vector<Matrix2d> _Sigmas;
            // vector<Matrix2d> _Vs;

            vector<Matrix2d> _pinv_matrices;
            // vector<Matrix2d> _pinvUs;
            // vector<Matrix2d> _pinvSigmas;
            // vector<Matrix2d> _pinvVs;

            vector<int> _qubits_to_clbits;

            QREM_Filter(vector< vector< vector<double> > > cal_matrices);
            void apply(map<string, int> hist);
    };
}