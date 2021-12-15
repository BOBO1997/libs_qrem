import cython
import numpy as np
cimport numpy as np
cimport cython

from libcpp.pair cimport pair
from libcpp.map cimport map
from libcpp.string cimport string

def expval_stddev(hist):
    cdef double shots
    shots = 0
    cdef double expval, sq_expval
    expval, sq_expval = 0, 0
    cdef str key
    cdef double count
    cdef double sigma_z
    for key, count in hist.items():
        shots += count
        sigma_z = 1
        for s in key:
            if s == "1":
                sigma_z *= -1
        expval += sigma_z ** 2 * count
        sq_expval += count
    expval /= shots
    cdef double stddev
    stddev = np.sqrt(sq_expval / shots - expval ** 2)
    return expval, stddev