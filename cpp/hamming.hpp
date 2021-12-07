using namespace std;

namespace libs_qrem
{
    string change_bit_at_poses(string key, vector<int> poses);
    set<string> extend_keys(set<string> original_keys, int max_dist);
    vector<double> extended_vectors(map<string, double> y, map<string, int> keys_to_indices);
    void print_vec1d();
}