#include <vector>

#include "qrem_filter.hpp"

using namespace std;
using namespace Eigen;

namespace libs_qrem {

    class Ignis_Filter : public QREM_Filter {
        public:

            Ignis_Filter(int num_clbits,
                        vector< vector< vector<double> > > cal_matrices,
                        vector< vector<int> > mit_pattern,
                        vector<int> meas_layout);

            int flip_state(int target_idx, int mat_index, vector<int>& flip_poses);

            virtual void apply(Args args);
    };
}