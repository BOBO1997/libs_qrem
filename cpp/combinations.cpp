#include <functional>
#include <iostream>
#include <vector>

#include "combinations.hpp"

using namespace std;

namespace libs_qrem {

void recursive_comb(vector<int>& indices, int s, int rest, vector< vector<int> >& nCk) {
    if (rest == 0) {
        nCk.push_back(indices);
    } else {
        if (s < 0) return;
        recursive_comb(indices, s - 1, rest, nCk);
        indices[rest - 1] = s;
        recursive_comb(indices, s - 1, rest - 1, nCk);
    }
}

vector< vector<int> > combinations(int n, int k) {
    vector<int> indices = vector<int>(k, 0);
    vector< vector<int> > ans = vector< vector<int> >(0, vector<int>(k, 0));
    recursive_comb(indices, n - 1, k, ans);
    return ans;
}

}