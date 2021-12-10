# from mitigation import qrem_filter
from mitigation import QREM_Filter_1, QREM_Filter_2
from dummy_data import n, cal_matrices, hist
from pprint import pprint

class qrem_filter:
    def __init__(self, n, cal_matrices, mit_pattern=[], meas_layout=[]) -> None:
        self.qrem_filter_cython = QREM_Filter_1(
            n, cal_matrices, mit_pattern, meas_layout)

    def apply(self, hist, d=0, threshold=0.1):
        return self.qrem_filter_cython(hist, d, threshold)

meas_filter1 = QREM_Filter_1(n, cal_matrices)
meas_filter2 = QREM_Filter_2(n, cal_matrices)

pprint(hist)
mitigated_hist1 = meas_filter1.apply(hist)
pprint(mitigated_hist1)

mitigated_hist2 = meas_filter2.apply(hist)
pprint(mitigated_hist2)
