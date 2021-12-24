#include <iostream>
#include <vector>
#include <queue>
#include <utility>
#include <cassert>

#include "sgs_algorithm.hpp"

using namespace std;

namespace libs_qrem {

vector<double> sgs_algorithm(vector<double>& x, bool make_sum_to_one = false) {
    priority_queue< pair<double, int>, 
                    vector< pair<double, int> >,
                    greater< pair<double, int> > > pq;
    double sum_of_x = 0;
    for (int state_idx = 0; state_idx < (int)x.size(); state_idx++) {
        if (x[state_idx] > 0) {
            pq.push(make_pair(x[state_idx], state_idx));
            sum_of_x += x[state_idx];
        }
    }
    double negative_accumulator = 0;
    if (make_sum_to_one) {
        negative_accumulator = 1 - sum_of_x;
    }
    while (pq.size() > 0) {
        pair<double, int> top = pq.top();
        if (top.first + negative_accumulator / pq.size() < 0) {
            negative_accumulator += top.first;
            pq.pop();
        }
        else { 
            break;
        }
    }
    int denominator = pq.size();
    vector<double> x_tilde = vector<double>(x.size(), 0);
    while (pq.size() > 0) {
        pair<double, int> top = pq.top();
        x_tilde[top.second] = top.first + negative_accumulator / denominator;
        assert(x_tilde[top.second] >= 0);
        pq.pop();
    }
    return x_tilde;
}

}