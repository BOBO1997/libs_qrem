from typing import List
import numpy as np
import time
from pprint import pprint
import copy

import sgs_algorithm
from sgs_algorithm import sgs_algorithm

from mitigation_tools import MitigationTools


class InvSLM0SGS(MitigationTools):

    # OK
    def __init__(self,
                 num_clbits: int,
                 cal_matrices: List[np.array],
                 mit_pattern: List[List[int]] = None,
                 meas_layout: List[int] = None) -> None:
        """
        Initialize the TensoredMitigation class
        
        Arguments
            num_clbits: number of measured qubits (int)
            cal_matrices: calibration matrices (list of 2 * 2 numpy array)
            meas_layout: the mapping from classical registers to qubits
        """
        super().__init__(num_clbits=num_clbits,
                         cal_matrices=cal_matrices,
                         mit_pattern=mit_pattern,
                         meas_layout=meas_layout)

    # OK
    def apply(self,
              counts: dict,
              shots: int = None,
              sgs: bool = True,
              rescale: bool = True,
              silent: bool = False) -> dict:
        """
        O(s * s * n) time and O(s) space

        Arguments
            counts: raw counts (dict of str to int)
            shots: total number of shot (int)
        Returns
            mitigated_counts: mitigated counts (dict of str to float)
        """

        if shots is None:
            shots = sum(counts.values())

        # make probability vector (dict), O(s) time
        y = {int(state, 2): counts[state] / shots for state in counts}

        if not silent:
            print("Restriction to labels of y + Lagrange Multiplier + SGS algorithm")

        t1 = time.time()
        # preprocess 1: compute sum of x, O(s * s * n) time in total
        sum_of_x = 0
        x_s = {state_idx: 0 for state_idx in y}  # O(s) space # e basis
        for state_idx in y:  # O(s) time
            sum_of_col = self.sum_of_tensored_vector(self.choose_vecs(state_idx, self.pinv_matrices))  # O(n) time
            sum_of_x += sum_of_col * y.get(state_idx, 0)
            x_s[state_idx] = self.mitigate_one_state(state_idx, y)  # O(n * s) time
        if not silent:
            print("sum of mitigated probability vector x_s:", sum(x_s.values()))

        # preprocess 2: compute the denominator of delta only summing up the first term of the denominator, O(n) time in total
        sum_of_vi = self.sum_of_tensored_vector(self.choose_vecs(0, self.pinvVs))  # O(n) time
        lambda_i = self.sum_of_tensored_vector(self.choose_vecs(0, self.pinvSigmas))  # O(n) time
        delta_denom = (sum_of_vi ** 2) / (lambda_i ** 2)  # O(1) time
        delta_coeff = (1 - sum_of_x) / delta_denom  # O(1) time

        # prepare x_hat_s only using the first term Delta_0 # only choose column 0, O(s * n) time in total
        x_hat_s = copy.deepcopy(x_s)
        sum_of_vi = self.sum_of_tensored_vector(self.choose_vecs(0, self.pinvVs))  # O(n) time
        lambda_i = self.sum_of_tensored_vector(self.choose_vecs(0, self.pinvSigmas))  # O(n) time
        delta_col = sum_of_vi / (lambda_i ** 2)  # O(1) time
        v_col = self.v_basis(0, x_hat_s.keys())  # O(s * n) time
        for state_idx in x_hat_s:  # O(s) time
            x_hat_s[state_idx] += delta_coeff * delta_col * v_col.get(state_idx, 0)  # O(1) time
        if not silent:
            print("sum of mitigated probability vector x_hat_s:", sum(x_hat_s.values()))

        # algorithm by Smolin et al. # O(s * log(s)) time
        # print(x_hat_s)
        x_tilde = sgs_algorithm(x_hat_s, silent=silent) if sgs else x_hat_s

        t2 = time.time()
        self.time = t2 - t1

        if not silent:
            print(t2 - t1, "s")

        if not silent:
            print("main process: Done!")
        mitigated_counts = {format(state, "0"+str(self.num_clbits)+"b"): x_tilde[state] * shots for state in x_tilde} if rescale else x_tilde # rescale to counts
        return mitigated_counts
