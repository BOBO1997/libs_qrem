#include <map>
#include <vector>
#include <Eigen/Dense>
#include <time.h>

using namespace std;
using Eigen::MatrixXd;

using namespace std;

namespace libs_qrem {
    class QREM_Filter {
        public:

            /*
            int num_clbits;
            vector<MatrixXd> cal_matrices;
            vector<vector<int>> mit_pattern;
            vector<int> meas_layout;
            
            vector<MatrixXd> svd_matrices;
            vector<MatrixXd> Us;
            vector<MatrixXd> Sigmas;
            vector<MatrixXd> Vs;

            vector<MatrixXd> pinv_matrices;
            vector<MatrixXd> pinvUs;
            vector<MatrixXd> pinvSigmas;
            vector<MatrixXd> pinvVs;

            vector<int> qubits_to_clbits;
            */

            // QREM_Filter(vector<MatrixXd> cal_matrices);
            QREM_Filter(vector<double> cal_matrices);


            

            void apply(map<string, double> hist);
    };
}