#include <iostream>
#include <map>
#include <set>
#include <vector>
#include <string>

#include "combinations.hpp"
#include "hamming.hpp"

using namespace std;

namespace libs_qrem {

string change_bit_at_poses(string key, vector<int> poses) {

    for (const auto& pos: poses) {
        if (key[pos] == '0') {
            key[pos] = '1';
        }
        else {
            key[pos] = '0';
        }
    }
    return key;
}

set<string> extend_keys(set<string>& original_keys, int max_dist) {
    set<string> extended_key_set;
    for (const auto& key: original_keys) {
        extended_key_set.insert(key);
        int n = key.size();
        for (int d = 0; d < max_dist; d++) {
            vector< vector<int> > combs = combinations(n, d + 1);
            for (const auto& comb: combs) {
                string new_key = change_bit_at_poses(key, comb);
                extended_key_set.insert(new_key);
            }
        }
    }
    return extended_key_set;
}

vector<double> extend_vectors(map<string, double>& y, vector<string>& indices_to_keys_vector) {
    vector<double> extended_y(indices_to_keys_vector.size(), 0);
    for (size_t i = 0; i < indices_to_keys_vector.size(); i++) {
        extended_y[i] = y[ indices_to_keys_vector[i] ];
    }
    return extended_y;
}

}