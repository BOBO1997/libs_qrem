from time import perf_counter
import cython
import numpy as np
cimport numpy as np
cimport cython

import scipy

from libcpp.pair cimport pair 
from libcpp.map cimport map
from libcpp.string cimport string

include "vector_wrapper.pyx"

# x
cdef class SLSQPFilter(BaseFilter):
    cdef SLSQP_Filter* instance_ptr
    cdef VectorDouble x_hat_vector
    cdef np.float64_t[:] x_hat_ndarray
    cdef np.float64_t[:] _x_hat

    def __cinit__(self, n, cal_matrices, mit_pattern = [], meas_layout = []):
        self.instance_ptr = new SLSQP_Filter(n, cal_matrices, mit_pattern, meas_layout)
        self.ptr = self.instance_ptr
        self.x_hat_vector = VectorDouble(1)
        self.method = "SLSQP"

    def __dealloc__(self):
        del self.ptr
    
    def apply(self, hist, d = 0, threshold = 0.1, silent = True):
        cdef str key
        cdef int value
        cdef map[string, int] cpp_hist
        for key, value in hist.items():
            cpp_hist[key.encode('utf-8')] = value
            self.shots += <double>value

        # apply inverse
        self.ptr.apply(cpp_hist, d, threshold)
        self.x_hat_vector.vec = self.ptr._x_s
        
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
        # cdef scipy.optimize.OptimizeResult res
        res = scipy.optimize.minimize(fun, x0, method='SLSQP', constraints=cons, tol=1e-6)
        cdef np.ndarray[np.float64_t, ndim=1] res_x
        res_x = res.x
        self._x_hat = res.x
        t2 = perf_counter() * 1000

        cdef pair[string, double] duration
        duration.first = "slsqp".encode('utf-8')
        duration.second = t2 - t1
        self.ptr._durations.insert(duration)

        # apply sgs_algorithm
        self._x_tilde.vec = sgs_algorithm(res_x, False)

        t1 = perf_counter() * 1000
        hist_dict = dict()
        cdef int i
        cdef string state
        for i, state in enumerate(self.ptr._indices_to_keys_vector):
            if not self._x_tilde.vec[i] == 0:
                hist_dict[state.decode('utf-8')] = self._x_tilde.vec[i] * self.shots
        t2 = perf_counter() * 1000

        duration.first = "sgs_algorithm".encode('utf-8')
        duration.second = t2 - t1
        self.ptr._durations.insert(duration)

        duration.first = "total".encode('utf-8')
        duration.second = 0

        cdef pair[string, double] item
        for item in self.ptr._durations:
            duration.second += item.second
        self.ptr._durations.insert(duration)

        if not silent:
            print("mitigation finished")
            for item in self.ptr._durations:
                print("time of", item.first.decode(), "is", item.second, "msec")

        return hist_dict