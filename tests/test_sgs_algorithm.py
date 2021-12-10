import numpy as np
from mitigation import QREM_Filter_1, QREM_Filter_2
from dummy_data import n, cal_matrices, hist

mitigator = QREM_Filter_1(len(cal_matrices), cal_matrices)

mitigated_hist = mitigator.apply(hist)

from pprint import pprint
pprint(mitigated_hist)