from typing import List
import numpy as np
from pprint import pprint
import time

from mitigation_tools import MitigationTools

from scipy.sparse.linalg import gmres

from sgs_algorithm import sgs_algorithm

class NationEtal(MitigationTools):

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

    def prepare_nation_A_tilde(self, keys: list):
        """
        prepare the numpy 2d matrix with the indices selected from keys
        """
        self.prepare_raw_A_tilde(keys)

        for j, col_key in enumerate(keys):
            sum_of_col = 0
            for i, row_key in enumerate(keys):
                sum_of_col += self.A_tilde[i][j]
            for i, row_key in enumerate(keys):
                self.A_tilde[i][j] /= sum_of_col

    def apply_inverse_A(self, extended_y):
        return np.linalg.inv(self.A_tilde) @ extended_y

    # OK
    def apply(self,
              counts: dict,
              shots: int = None,
              d: int = 0,
              rescale: bool = True,
              silent: bool = False) -> dict:
        """
        Arguments
            counts: raw counts (dict of str to int)
            shots: total number of shot (int)

        Returns
            mitigated_counts: mitigated counts (dict of str to float)
        """

        if shots is None:
            shots = sum(counts.values())

        # make probability vector (dict)
        # y = {int(state, 2): counts[state] / shots for state in counts}
        y = {state: counts[state] / shots for state in counts}

        if not silent:
            print("Method by Nation, Kang, Sundaresan, and Gambatta")

        t1 = time.time()

        # Prepare small calibration matrix A tilde
        extended_keys = self.extend_keys(set(y.keys()), d)
        sorted_extended_keys = sorted(extended_keys)
        self.prepare_nation_A_tilde(sorted_extended_keys)

        keys_to_indices = dict()
        for i, key in enumerate(sorted_extended_keys):
            keys_to_indices[key] = i
        
        extended_y = self.extend_vectors(y, keys_to_indices)
        x, _ = gmres(self.A_tilde, extended_y)
        # x = np.linalg.inv(self.A_tilde) @ extended_y
        # x = self.apply_inverse_A(extended_y)
        if not silent:
            print("main process: Done!")
            print("sum of mitigated probability vector x:", np.sum(x))
        
        x_hat = dict()
        for key, index in keys_to_indices.items():
            x_hat[key] = x[index]

        if not silent:
            print("start sgs_algorithm")

        x_tilde = sgs_algorithm(x_hat)

        t2 = time.time()
        self.time = t2 - t1

        if not silent:
            print(t2 - t1, "s")

        mitigated_counts = {state: x_tilde[state] * shots for state in x_tilde} if rescale else x_tilde # rescale to counts
        # mitigated_counts = {format(state, "0"+str(self.num_clbits)+"b"): x_tilde[state] * shots for state in x_tilde} if rescale else x_tilde # rescale to counts
        return mitigated_counts
