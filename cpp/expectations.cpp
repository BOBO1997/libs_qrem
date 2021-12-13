#include <map>
#include <cstdlib>
#include <bitset>
#include <cmath>

#include "expectations.hpp"

namespace libs_qrem {

    pair<double, double> expval_stddev(map<string, double> hist) {
        int shots = 0;
        int n = 0;
        double expval = 0;
        double sq_expval = 0;
        for (const auto& item: hist) {
            shots += item.second;
            n = (int)item.first.size();
            int sigma_z = 1;
            for (int j = 0; j < n; j++) {
                if (item.first[j] == '1') {
                    sigma_z *= -1;
                }
            }
            expval += (double)sigma_z * item.second;
            sq_expval += item.second * item.second;
        }
        expval /= shots;
        double stddev = sqrt(sq_expval / shots - expval * expval);
        return make_pair(expval, stddev);
    }
}