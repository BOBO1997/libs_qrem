import cython
cimport cython

include "qrem_filter.pyx"

def test():
    cdef vector[double] cal_matrices
    # cdef vector[MatrixXd] cal_matrices
    return QREM_Filter_Cython(cal_matrices)