import cython
import numpy as np
cimport numpy as np
cimport cython

from scipy.optimize import minimize

from libcpp.map cimport map
from libcpp.string cimport string

# OK
cdef class NationEtalFilter(BaseFilter):
    cdef NationEtal_Filter* instance_ptr

    def __cinit__(self, n, cal_matrices, mit_pattern = [], meas_layout = []):
        self.instance_ptr = new NationEtal_Filter(n, cal_matrices, mit_pattern, meas_layout)
        self.ptr = self.instance_ptr
    
    def __dealloc__(self):
        del self.ptr
    
    def apply(self, hist, d = 0, threshold = 0.1, method = "iterative", silent = True):
        cdef map[string, int] cpp_hist
        for key, value in hist.items():
            cpp_hist[key.encode('utf-8')] = value
            self.shots += <double>value
        self.ptr.apply(cpp_hist, d, threshold)
        cdef map[string, double] mitigated_hist
        mitigated_hist = self.ptr._mitigated_hist
        self._x_s.vec = self.ptr._x_s
        self._x_hat.vec = self.ptr._x_hat
        self._x_tilde.vec = self.ptr._x_tilde
        if not silent:
            print("mitigation finished")
            for item in self.ptr._durations:
                print("time of", item.first.decode(), "is", item.second, "msec")
        hist_dict = dict()
        for item in mitigated_hist:
            hist_dict[item.first.decode('utf-8')] = item.second
        return hist_dict