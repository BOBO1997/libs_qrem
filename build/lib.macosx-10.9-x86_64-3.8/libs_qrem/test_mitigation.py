# from mitigation import qrem_filter
from mitigation import QREM_Filter_Cython

class qrem_filter:
    def __init__(self, n, cal_matrices, mit_pattern=[], meas_layout=[]) -> None:
        # self.num_clbits = n
        # self.cal_matrices = cal_matrices
        # self.mit_pattern = mit_pattern
        # self.meas_layout = meas_layout
        self.qrem_filter_cython = QREM_Filter_Cython(
            n, cal_matrices, mit_pattern, meas_layout)

    def apply(self, hist, d=0, threshold=0.1):
        return self.qrem_filter_cython(hist, d, threshold)

n = 3
cal_matrices = [[[0.9, 0.1], [0.2, 0.8]], [[0.9, 0.1], [0.2, 0.8]], [[0.9, 0.1], [0.2, 0.8]]]

# meas_filter = qrem_filter(3, cal_matrices)
meas_filter = QREM_Filter_Cython(n, cal_matrices)

# hist = {"000": 50, "001": 3, "010": 4, "011": 2, "100": 0, "101": 50, "110":50,"111":50,}
hist = {"000": 48, "010": 4, "101": 1, "111":43}
print(hist)
mitigated_hist = meas_filter.apply(hist)
print(mitigated_hist)
