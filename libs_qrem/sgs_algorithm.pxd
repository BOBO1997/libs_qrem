import cython
cimport cython

from libcpp.vector cimport vector

cdef extern from "../cpp/mitigation.hpp" namespace "libs_qrem":
    cdef vector[double] sgs_algorithm(vector[double])