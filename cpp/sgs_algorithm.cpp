#include <iostream>
#include <vector>
#include <queue>
#include <utility>

#include "sgs_algorithm.hpp"

using namespace std;

namespace libs_qrem {

vector<double> sgs_algorithm(vector<double> x) {
    priority_queue< pair<double, int> > pq;
    double sum_of_x = 0;
    for (int state_idx = 0; state_idx < x.size(); state_idx++) {
        if (x[state_idx] > 0) {
            pq.push(make_pair(x[state_idx], state_idx));
            sum_of_x += x[state_idx];
        }
    }
    double negative_accumulator = 1 - sum_of_x;
    if (negative_accumulator >= 0) {
        cout << "accumulator is positive, we might even ignoring the necessal positive values." << endl;
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
    vector<double> x_tilde = vector<double>(denominator, 0);
    while (pq.size() > 0) {
        pair<double, int> top = pq.top();
        x_tilde[top.second] = top.first + negative_accumulator / denominator;
        pq.pop();
    }
    return x_tilde;
}

}