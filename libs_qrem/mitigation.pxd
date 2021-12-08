import cython
cimport cython

from libcpp.vector cimport vector
from libcpp.set cimport set
from libcpp.map cimport map
from libcpp.string cimport string

cdef extern from "../cpp/mitigation.hpp" namespace "libs_qrem":
    cdef cppclass QREM_Filter:
        map[string, double] _durations
        QREM_Filter(int num_clbits,
                    vector[vector[vector[double]]] cal_matrices,
                    vector[vector[int]] mit_pattern,
                    vector[int] meas_layout)
        map[string, double] apply(map[string, int] hist, int d, double threshold)