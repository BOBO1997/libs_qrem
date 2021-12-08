import cython
cimport cython

from libcpp.map cimport map
from libcpp.string cimport string

cdef class QREM_Filter_Cython:
    cdef QREM_Filter* ptr

    def __cinit__(self, n, cal_matrices, mit_pattern = [], meas_layout = []):
        self.ptr = new QREM_Filter(n, cal_matrices, mit_pattern, meas_layout)
    
    def __deadaloc(self):
        del self.ptr
    
    def apply(self, hist, d = 0, threshold = 0.1):
        cdef map[string, int] cpp_hist
        for key, value in hist.items():
            cpp_hist[key.encode()] = value
        mitigated_hist = self.ptr.apply(cpp_hist, d, threshold)
        times = self.ptr._durations
        print("finished")
        for item in times:
            print("time of", item.first.decode(), "is", item.second, "msec")
        hist_dict = dict()
        for item in mitigated_hist:
            hist_dict[item.first] = item.second
        return hist_dict