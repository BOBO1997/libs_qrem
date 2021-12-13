#include <set>
#include <map>
#include <vector>
#include <string>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/SVD>
#include <time.h>

#include "qrem_filter_base.hpp"

using namespace std;
using namespace Eigen;

namespace libs_qrem {

    class QREM_Filter_MooneyEtal : public QREM_Filter_Base {
        public:

            QREM_Filter_MooneyEtal(int num_clbits,
                        vector< vector< vector<double> > > cal_matrices,
                        vector< vector<int> > mit_pattern,
                        vector<int> meas_layout);
            
            int flip_state(int state, int mat_idx, vector<int>& flip_poses);

            virtual void apply(map<string, int> hist, 
                               int d, 
                               double threshold);
    };
}