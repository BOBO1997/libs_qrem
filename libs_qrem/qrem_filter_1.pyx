import cython
import numpy as np
cimport numpy as np
cimport cython

from scipy.optimize import minimize

from libcpp.map cimport map
from libcpp.string cimport string

# OK
cdef class QREM_Filter_1:
    cdef QREM_Filter* ptr
    cdef double expval, stddev
    cdef VectorDouble _x_s
    cdef VectorDouble _x_hat
    cdef VectorDouble _x_tilde

    def __cinit__(self, n, cal_matrices, mit_pattern = [], meas_layout = []):
        self.ptr = new QREM_Filter(n, cal_matrices, mit_pattern, meas_layout)
    
    def __dealloc__(self):
        del self.ptr
    
    def sum_of_x(self):
        return self.ptr._sum_of_x

    def sum_of_x_hat(self):
        return self.ptr._sum_of_x_hat

    def sum_of_x_tilde(self):
        return self.ptr._sum_of_x_tilde

    def reduced_inv_A(self):
        return matrix_to_ndarray(self.ptr._reduced_inv_A)

    def one_norm(self):
        return self.ptr._one_norm

    def mitigated_hist(self):
        hist_dict = dict()
        for item in self.ptr._mitigated_hist:
            hist_dict[item.first.decode('utf-8')] = item.second
        return hist_dict
    
    def x_s(self):
        return vector_to_list(self.ptr._x_s)

    def x_hat(self):
        return vector_to_list(self.ptr._x_hat)

    def x_tilde(self):
        return vector_to_list(self.ptr._x_tilde)
    
    def indices_to_keys_vector(self):
        return vector_to_list(self.ptr._indices_to_keys_vector)

    def times(self):
        times = dict()
        for item in self.ptr._durations:
            times[item.first.decode('utf-8')] = item.second
        return times

    def expval_stddev(self):
        self.expval, self.stddev = expval_stddev(self.mitigated_hist())
        return self.expval, self.stddev

    def apply(self, hist, d = 0, threshold = 0.1):
        cdef map[string, int] cpp_hist
        for key, value in hist.items():
            cpp_hist[key.encode('utf-8')] = value
        self.ptr.apply(cpp_hist, d, threshold)
        cdef map[string, double] mitigated_hist
        mitigated_hist = self.ptr._mitigated_hist
        print("mitigation finished")
        self._x_s.vec = self.ptr._x_s
        self._x_hat.vec = self.ptr._x_hat
        self._x_tilde.vec = self.ptr._x_tilde
        for item in self.ptr._durations:
            print("time of", item.first.decode(), "is", item.second, "msec")
        hist_dict = dict()
        for item in mitigated_hist:
            hist_dict[item.first.decode('utf-8')] = item.second
        return hist_dict