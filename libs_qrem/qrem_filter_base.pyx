import cython
import numpy as np
cimport numpy as np
cimport cython

from scipy.optimize import minimize

from libcpp.map cimport map
from libcpp.string cimport string

# OK
cdef class QREMFilter:
    cdef QREM_Filter* ptr
    cdef double expval, stddev
    cdef double shots

    def __cinit__(self):
        pass
    
    def __dealloc__(self):
        del self.ptr
    
    def sum_of_x(self):
        return self.ptr._sum_of_x

    def sum_of_x_hat(self):
        return self.ptr._sum_of_x_hat

    def sum_of_x_tilde(self):
        return self.ptr._sum_of_x_tilde
    
    def reduced_A(self):
        self.ptr.compute_reduced_A(self.ptr._indices_to_keys_vector)
        return matrix_to_ndarray(self.ptr._reduced_A)

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
        return vector_to_list_double(self.ptr._x_s)

    def x_hat(self):
        return vector_to_list_double(self.ptr._x_hat)

    def x_tilde(self):
        return vector_to_list_double(self.ptr._x_tilde)
    
    def indices_to_keys_vector(self):
        return vector_to_list_string(self.ptr._indices_to_keys_vector)

    def times(self):
        times = dict()
        for item in self.ptr._durations:
            times[item.first.decode('utf-8')] = item.second
        return times

    def expval_stddev(self):
        self.expval, self.stddev = expval_stddev(self.mitigated_hist())
        self.stddev = self.one_norm() / np.sqrt(self.shots)
        return self.expval, self.stddev
