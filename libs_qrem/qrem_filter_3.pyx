import cython
import numpy as np
cimport numpy as np
cimport cython

from scipy.optimize import minimize

from libcpp.map cimport map
from libcpp.string cimport string

# OK
cdef class QREM_Filter_3:
    cdef QREM_Filter_MooneyEtal* ptr

    def __cinit__(self, n, cal_matrices, mit_pattern = [], meas_layout = []):
        self.ptr = new QREM_Filter_MooneyEtal(n, cal_matrices, mit_pattern, meas_layout)
    
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

    def times(self):
        times = dict()
        for item in self.ptr._durations:
            times[item.first.decode('utf-8')] = item.second
        return times

    def apply(self, hist, d = 0, threshold = 0.1):
        cdef map[string, int] cpp_hist
        for key, value in hist.items():
            cpp_hist[key.encode('utf-8')] = value
        self.ptr.apply(cpp_hist, d, threshold)
        cdef map[string, double] mitigated_hist
        mitigated_hist = self.ptr._mitigated_hist
        print("mitigation finished")
        for item in self.ptr._durations:
            print("time of", item.first.decode(), "is", item.second, "msec")
        hist_dict = dict()
        for item in mitigated_hist:
            hist_dict[item.first.decode('utf-8')] = item.second
        return hist_dict