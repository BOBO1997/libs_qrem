from dummy_data50 import *

cpp_code = "#include <map>\n"
cpp_code += "#include <string>\n"
cpp_code += "#include <vector>\n"
cpp_code += "#include <iostream>\n"
cpp_code += "\n"
cpp_code += "using namespace std;\n"
cpp_code += "\n"
cpp_code += "map<string, int> hist;\n"
cpp_code += "cout << hist.size() << endl;\n"
# for key, value in hist.items():
#     cpp_code += "hist.insert(make_pair(\"" + key + \
#         "\", " + str(value) + "));\n"

cpp_code += "\n"
cpp_code += "int n = 50;\n"
cpp_code += "vector< vector< vector<double> > > cal_matrices;\n"
for matrix in cal_matrices:
    cpp_code += "cal_matrices.push_back({{" + str(matrix[0, 0]) + "," + str(
        matrix[0, 1]) + "},{" + str(matrix[1, 0]) + "," + str(matrix[1, 1]) + "}});\n"

with open("dummy_data50.cpp", "w") as f:
    print("writing...")
    f.write(cpp_code)
    f.close()
print("finished")