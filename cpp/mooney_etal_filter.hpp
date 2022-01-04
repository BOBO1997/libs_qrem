#include <vector>
#include <string>

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
            
            int flip_state(int state, int mat_idx, vector<int>& flip_poses);
            string flip_state(string state, int mat_idx, vector<int>& flip_poses);
            string flip_binary_string_at_pos(string str, int pos);

            virtual void apply(Args args);
    };
}