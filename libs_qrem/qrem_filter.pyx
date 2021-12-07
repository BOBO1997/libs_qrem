import cython
cimport cython

from libcpp.map cimport map
from libcpp.string cimport string

cdef class QREM_Filter_Cython:
    cdef QREM_Filter* ptr

    def __cinit__(self, cal_matrices):
        self.ptr = new QREM_Filter(cal_matrices)
    
    def __deadaloc(self):
        del self.ptr
    
    def apply(self):
        cdef map[string, double] hist
        self.ptr.apply(hist)