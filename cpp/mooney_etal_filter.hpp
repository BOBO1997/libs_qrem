#include <set>
#include <map>
#include <vector>
#include <string>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/SVD>
#include <time.h>

#include "qrem_filter.hpp"

using namespace std;
using namespace Eigen;

namespace libs_qrem {

    class MooneyEtal_Filter : public QREM_Filter {
        public:

            MooneyEtal_Filter(int num_clbits,
                        vector< vector< vector<double> > > cal_matrices,
                        vector< vector<int> > mit_pattern,
                        vector<int> meas_layout);
            
            string btos(int target_int, int n);

            int flip_state(int state, int mat_idx, vector<int>& flip_poses);

            virtual void apply(map<string, int> hist, 
                               int d = 0, 
                               double threshold = 0.1);
    };
}