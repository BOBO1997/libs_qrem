import cython
cimport cython

from libcpp.vector cimport vector
from libcpp.set cimport set
from libcpp.map cimport map
from libcpp.string cimport string

cdef extern from "../cpp/mitigation.hpp" namespace "libs_qrem":
    cdef cppclass QREM_Filter:
        QREM_Filter(vector[double] cal_matrices)
        # QREM_Filter(vector[MatrixXd] cal_matrices)
        void apply(map[string, double] hist)