import cython
import numpy as np
cimport numpy as np
cimport cython

from scipy.optimize import minimize

from libcpp.map cimport map
from libcpp.string cimport string

# OK
cdef class QREM_Filter_Base:
    cdef QREM_Filter* ptr
    cdef double expval, stddev
    cdef VectorDouble _x_s
    cdef VectorDouble _x_hat
    cdef VectorDouble _x_tilde

    def __cinit__(self, n, cal_matrices, mit_pattern = [], meas_layout = []):
        self.ptr = new QREM_Filter(n, cal_matrices, mit_pattern, meas_layout)
    
    def __deadaloc(self):
        del self.ptr
    
    def sum_of_x(self):
        return self.ptr._sum_of_x

    def sum_of_x_hat(self):
        return self.ptr._sum_of_x_hat

    def mitigated_hist(self):
        hist_dict = dict()
        for item in self.ptr._mitigated_hist:
            hist_dict[item.first.decode('utf-8')] = item.second
        return hist_dict
    
    def x_s(self):
        return np.asarray(self._x_s)

    def x_hat(self):
        return np.asarray(self._x_hat)

    def x_tilde(self):
        return np.asarray(self._x_tilde)

    def times(self):
        times = dict()
        for item in self.ptr._durations:
            times[item.first.decode('utf-8')] = item.second
        return times

    def expval_stddev(self):
        self.expval, self.stddev = expval_stddev(self.mitigated_hist())
        return self.expval, self.stddev