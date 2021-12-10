from time import perf_counter
import cython
import numpy as np
cimport numpy as np
cimport cython

from scipy.optimize import minimize

from libcpp.pair cimport pair 
from libcpp.map cimport map
from libcpp.string cimport string

include "vector_wrapper.pyx"

# x
cdef class QREM_Filter_2:
    cdef QREM_Filter_Nlp* ptr
    cdef VectorDouble x_hat_vector
    # cdef np.ndarray[np.float64_t, ndim=1] x_hat_ndarray

    def __cinit__(self, n, cal_matrices, mit_pattern = [], meas_layout = []):
        self.ptr = new QREM_Filter_Nlp(n, cal_matrices, mit_pattern, meas_layout)
        self.x_hat_vector = VectorDouble(1)

    def __dealloc__(self):
        del self.ptr
    
    def sum_of_x(self):
        return self.ptr._sum_of_x

    def times(self):
        times = dict()
        cdef pair[string, double] item
        for item in self.ptr._durations:
            times[item.first.decode('utf-8')] = item.second
        return times

    def apply(self, hist, d = 0, threshold = 0.1):
        cdef str key
        cdef int value
        cdef map[string, int] cpp_hist
        for key, value in hist.items():
            cpp_hist[key.encode('utf-8')] = value

        # apply inverse
        self.x_hat_vector.vec = self.ptr.apply(cpp_hist, d, threshold)
        times = self.ptr._durations
        print("finished")
        for item in times:
            print("time of", item.first.decode(), "is", item.second, "msec")
        
        cdef int vec_size = self.ptr._indices_to_keys_vector.size()

        # apply optimization
        self.x_hat_ndarray = np.asarray(self.x_hat_vector)
        def fun(np.ndarray[np.float64_t, ndim=1] x):
            return sum((x - self.x_hat_ndarray) ** 2)

        cdef np.ndarray[np.float64_t, ndim=1] x0 = np.random.rand(vec_size)
        x0 = x0 / sum(x0)
        cons = ({'type': 'eq', 'fun': lambda x: 1 - sum(x)})

        cdef double t1, t2
        t1 = perf_counter() * 1000
        cdef np.ndarray[np.float64_t, ndim=1] res = minimize(fun, x0, method='SLSQP', constraints=cons, tol=1e-6)
        t2 = perf_counter() * 1000

        self.ptr_durations["slsqp".encode('utf-8')] = t2 - t1

        # apply sgs_algorithm
        cdef vector[double] x_tilde = sgs_algorithm(res.x)

        t1 = perf_counter() * 1000
        hist_dict = dict()
        cdef int i
        cdef string state
        for i, state in enumerate(self.ptr._indices_to_keys_vector):
            hist_dict[state.decode('utf-8')] = x_tilde[i]
        t2 = perf_counter() * 1000

        self.ptr_durations["sgs_algorithm".encode('utf-8')] = t2 - t1

        return hist_dict