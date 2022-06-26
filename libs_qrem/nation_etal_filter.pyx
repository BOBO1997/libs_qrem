import sys
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
        self.method = "Nation et al."
    
    def __dealloc__(self):
        del self.instance_ptr

    def iterations(self):
        return self.instance_ptr._iterations
    
    def error(self):
        return self.instance_ptr._error
    
    def apply_specific(self, hist, d=0, silent=True, threshold=0.1, method="bicgstab"):
        cdef map[string, int] cpp_hist
        for key, value in hist.items():
            cpp_hist[key.encode('utf-8')] = value
            self.shots += <double>value
        cdef Args args
        args.hist = cpp_hist
        args.d = d
        if method == "bicgstab" or method == "iterative" or method == "lu" or method == "exact" or method == "direct":
            args.method = method.encode('utf-8')
        else:
            raise Exception("No such method: " + method)
            sys.exit(1)
        self.ptr.apply(args)
        if not silent:
            print("mitigation finished")
            for item in self.ptr._durations:
                print("time of", item.first.decode(), "is", item.second, "msec")
        return self.mitigated_hist()