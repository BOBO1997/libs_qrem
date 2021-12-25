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

    class NationEtal_Filter : public QREM_Filter {
        public:

            NationEtal_Filter(int num_clbits,
                        vector< vector< vector<double> > > cal_matrices,
                        vector< vector<int> > mit_pattern,
                        vector<int> meas_layout);

            virtual void apply(map<string, int> hist, 
                               int d = 0,
                               double threshold = 0.1);
    };
}