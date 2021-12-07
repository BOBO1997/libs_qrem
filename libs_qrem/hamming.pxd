import cython
cimport cython

from libcpp.map cimport map
from libcpp.string cimport string
from libcpp.vector cimport vector
from libcpp.set cimport set

cdef extern from "../cpp/combinations.hpp" namespace "libs_qrem":
    cdef void recursive_comb(vector[int]&, int s, int rest, vector[ vector[int] ]& nCk)
    cdef vector[ vector[int] ] combinations(int n, int k)

cdef extern from "../cpp/hamming.hpp" namespace "libs_qrem":
    cdef void print_vec1d()
    cdef set[string] extend_keys(set[string] original_keys, int max_dist)
    cdef vector[double] extended_vectors(map[string, double] y, map[string, int] keys_to_indices)