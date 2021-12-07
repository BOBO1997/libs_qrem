from typing import List, Union
import numpy as np
import time
from pprint import pprint

from numpy.core.records import array

import sgs_algorithm
from sgs_algorithm import sgs_algorithm

import mitigation_tools
from mitigation_tools import MitigationTools

class MooneyEtal(MitigationTools):

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

    def flip_state(self, state: Union[str,int], mat_index: int, flip_poses: List[int]) -> Union[str, int]:
        """
        Flip the state according to the chosen qubit positions
        mat_index indicates where to flip and where not to flip (as a mask of mat_index)
        flip_poses indicates the indices to be flipped
        """
        flip_poses = [pos for i, pos in enumerate(flip_poses) if (mat_index >> i) & 1] # apply mask
        if isinstance(state, str):
            flip_poses = sorted(flip_poses)
            new_state = ""
            pos = 0
            for flip_pos in flip_poses:
                new_state += state[pos:flip_pos]
                new_state += str(int(state[flip_pos], 2) ^ 1)  # flip the state
                pos = flip_pos + 1
            new_state += state[pos:]
            return new_state
        else:
            for flip_pos in flip_poses:
                state = state ^ ((1 << (self.num_clbits - 1)) >> flip_pos)
            return state

    # OK
    def apply(self,
              counts: dict,
              threshold: float = 0.1,
              shots: int = None,
              sgs: bool = True,
              rescale: bool = True,
              silent: bool = False) -> dict:
        """
        The algorithm by Mooney et al.

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
            print("The heuristcs by Mooney et al. + SGS algorithm")

        t1 = time.time()
        
        for pinv_mat, pos_qubits in zip(self.pinv_matrices, self.mit_pattern): # O(n) time
            x = dict()
            pos_clbits = [self.qubits_to_clbits[qubit] for qubit in pos_qubits]
            for state, _ in y.items(): # O(s) time
                first_index = self.index_of_mat(state, pos_clbits)
                sum_of_count = 0
                for i in range(len(pinv_mat)):  # i is index of pinv_mat # O(1) time
                    source_state = self.flip_state(state, i, pos_clbits)
                    second_index = self.index_of_mat(source_state, pos_clbits)
                    sum_of_count += pinv_mat[first_index, second_index] * y.get(source_state, 0)
                if abs(sum_of_count) >= threshold:
                    x[state] = sum_of_count
            y = x

        sum_of_x = sum(x.values())
        if not silent:
            print(x)
            print("sum of counts:", sum_of_x)
        if sum_of_x < 0:
            print("negative counts")

        # algorithm by Smolin et al. # O(s * log(s)) time
        x_tilde = sgs_algorithm(x, silent=silent) if sgs else x

        t2 = time.time()
        self.time = t2 - t1
        
        if not silent:
            print(t2 - t1, "s")
        
        if not silent:
            print("main process: Done!")
        mitigated_counts = {format(state, "0"+str(self.num_clbits)+"b"): x_tilde[state] * shots for state in x_tilde} if rescale else x_tilde  # rescale to counts
        return mitigated_counts
