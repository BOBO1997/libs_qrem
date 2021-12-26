import cython
import numpy as np
cimport numpy as np
cimport cython

from scipy.optimize import minimize

from libcpp.map cimport map
from libcpp.string cimport string

# OK
cdef class BaseFilter:
    cdef QREM_Filter* ptr
    cdef double shots
    cdef str method

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
        if self.ptr._reduced_A.size() == 0:
            self.ptr.compute_reduced_A(self.ptr._indices_to_keys_vector.size())
        return matrix_to_ndarray(self.ptr._reduced_A)

    def normalized_reduced_A(self):
        if self.ptr._reduced_A.size() == 0:
            self.ptr.compute_reduced_A(self.ptr._indices_to_keys_vector.size())
            normalize_cols(self.ptr._reduced_A)
        return matrix_to_ndarray(self.ptr._reduced_A)

    def reduced_inv_A(self):
        return matrix_to_ndarray(self.ptr._reduced_inv_A)

    def exact_one_norm_of_inv_reduced_A(self):
        if self.ptr._exact_one_norm_of_inv_reduced_A < 0:
            self.ptr.exact_one_norm_of_inv_reduced_A()
        return self.ptr._exact_one_norm_of_inv_reduced_A

    def exact_one_norm_of_reduced_inv_A(self):
        if self.ptr._exact_one_norm_of_reduced_inv_A < 0:
            self.ptr.exact_one_norm_of_reduced_inv_A()
        return self.ptr._exact_one_norm_of_reduced_inv_A

    def iterative_one_norm_of_inv_reduced_A(self, method="iterative"):
        if self.ptr._iterative_one_norm_of_inv_reduced_A < 0:
            self.ptr.iterative_one_norm_of_inv_reduced_A(method.encode('utf-8'))
        return self.ptr._iterative_one_norm_of_inv_reduced_A

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

    def expval(self):
        expval = 0
        cdef str key
        cdef double count, sigma_z
        for key, count in self.mitigated_hist().items():
            sigma_z = 1
            for s in key:
                if s == "1":
                    sigma_z *= -1
            expval += sigma_z * count
        expval /= self.shots
        return expval

    def mitigation_stddev(self, norm_type = "exact"):
        if self.method == "delta" or self.method == "SLSQP" or self.method == "least norm":
            return self.exact_one_norm_of_reduced_inv_A() / np.sqrt(self.shots)
        elif self.method == "Nation et al.":
            if norm_type == "exact" or norm_type == "direct":
                return self.exact_one_norm_of_inv_reduced_A() / np.sqrt(self.shots)
            else:
                return self.iterative_one_norm_of_inv_reduced_A() / np.sqrt(self.shots)
        elif self.method == "Mooney et al.":
            print("Cannot compute the standard deviation of mitigation.")
