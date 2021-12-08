import cython
cimport cython

from libcpp.map cimport map
from libcpp.string cimport string

cdef class QREM_Filter_Cython:
    cdef QREM_Filter* ptr

    def __cinit__(self, n, cal_matrices, mit_pattern = [], meas_layout = []):
        self.ptr = new QREM_Filter(n, cal_matrices, mit_pattern, meas_layout)
    
    def __deadaloc(self):
        del self.ptr
    
    def apply(self, hist, d = 0, threshold = 0.1):
        self.ptr.apply(hist, d, threshold)