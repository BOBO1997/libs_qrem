import cython
cimport cython

from libcpp.map cimport map
from libcpp.string cimport string
from libcpp.vector cimport vector

include "qrem_filter.pyx"

def test_qrem_filter_cython():
    cdef vector[vector[vector[double]]] cal_matrices
    cal_matrices.push_back([[0.95,0.05],[0.1,0.9]])
    cal_matrices.push_back([[0.9,0.1],[0.2,0.8]])
    print("cal_matrices:", cal_matrices)
    meas_filter = QREM_Filter_Cython(2, cal_matrices)

    cdef map[string, int] hist
    hist["00".encode()] = 44
    hist["01".encode()] = 8
    hist["10".encode()] = 12
    hist["11".encode()] = 36
    print("hist: ", hist)
    mitigated_hist = meas_filter.apply(hist)
    print("mitigated_hist:", mitigated_hist)

def test_constructor():
    cdef vector[vector[vector[double]]] cal_matrices
    cal_matrices.push_back([[0.95,0.05],[0.1,0.9]])
    cal_matrices.push_back([[0.9,0.1],[0.2,0.8]])
    print("cal_matrices:", cal_matrices)
    meas_filter = QREM_Filter_Cython(2, cal_matrices)