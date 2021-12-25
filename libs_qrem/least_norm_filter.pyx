import cython
import numpy as np
cimport numpy as np
cimport cython

from scipy.optimize import minimize

from libcpp.map cimport map
from libcpp.string cimport string

# OK
cdef class LeastNormFilter(BaseFilter):
    cdef Least_Norm_Filter* instance_ptr

    def __cinit__(self, n, cal_matrices, mit_pattern = [], meas_layout = []):
        self.instance_ptr = new Least_Norm_Filter(n, cal_matrices, mit_pattern, meas_layout)
        self.ptr = self.instance_ptr
        self.method = "least norm"
    
    def __dealloc__(self):
        del self.ptr
    
    def apply(self, hist, d = 0, threshold = 0.1, silent = True):
        cdef map[string, int] cpp_hist
        for key, value in hist.items():
            cpp_hist[key.encode('utf-8')] = value
            self.shots += <double>value
        self.ptr.apply(cpp_hist, d, threshold)
        if not silent:
            print("mitigation finished")
            for item in self.ptr._durations:
                print("time of", item.first.decode(), "is", item.second, "msec")
        return self.mitigated_hist()