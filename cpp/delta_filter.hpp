#include <vector>

#include "qrem_filter.hpp"

using namespace std;
using namespace Eigen;

namespace libs_qrem {

    class Delta_Filter : public QREM_Filter {
        public:

            Delta_Filter(int num_clbits,
                              vector< vector< vector<double> > > cal_matrices,
                              vector< vector<int> > mit_pattern,
                              vector<int> meas_layout);

            virtual void apply(Args args);
    };
}