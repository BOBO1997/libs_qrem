using namespace std;

namespace libs_qrem
{
    string change_bit_at_poses(string key, vector<int> poses);
    set<string> extend_keys(set<string>& original_keys, int max_dist);
    vector<double> extend_vector(map<string, double>& y, vector<string>& keys_to_indices);
}