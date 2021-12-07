from typing import Union, List
import numpy as np
from pprint import pprint

import svds
from svds import SVDs

from itertools import combinations
import copy

from draw_heatmap import draw_heatmap

class MitigationTools(SVDs):

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

    
        self.A_tilde = None
        self.time = -1

    # OK
    def index_of_mat(self, state: Union[str, int], pos_clbits: List[int]) -> int:
        """
        Compute the index of (pseudo inverse) calibration matrix for the input quantum state
        If input state is string, then this process costs O(n) time, O(n) space
        If input state is integer, then this process costs O(len(pos_clbits)) time, O(1) space
        Using tensored mitigation, len(pos_clbits) == 1, therefore the time complexity would be O(1).
        
        Arguments
            state: the quantum states we focusing
            pos_clbits: the classical qubits of which the matrix is in charge (equivalent to indicating the matrix)
        Returns
            sub_state: the position in the matrix
        """
        if isinstance(state, str):  # O(n) time, O(n) space
            sub_state = ""
            for pos in pos_clbits:
                sub_state += state[pos]
            return int(sub_state, 2)
        else:
            sub_state = 0
            for pos in pos_clbits:  # O(len(pos_clbits)) time, O(1) space
                sub_state <<= 1
                sub_state += (state >> (self.num_clbits - pos - 1)) & 1
            return sub_state

    # OK
    def mitigate_one_state(self, target_state: Union[str, int], counts: dict) -> float:
        """
        Mitigate one state using inverse calibration matrices.
        O(n * shots) time, O(shots) space
        
        Arguments
            target_state: quanutum state to be mitigated (str or int)
            counts: raw counts (dict)
        Returns
            new_count: mitigated count of target state (float)
        """
        new_count = 0
        for source_state in counts:  # O(shots)
            tensor_elem = 1.
            # O(n) time
            for pinv_mat, pos_qubits in zip(self.pinv_matrices, self.mit_pattern):
                # if completely tensored, then len(pos_clbits) == 1 -> O(1) time
                pos_clbits = [self.qubits_to_clbits[qubit] for qubit in pos_qubits]
                first_index = self.index_of_mat(target_state, pos_clbits)  # O(1) time
                second_index = self.index_of_mat(source_state, pos_clbits)  # O(1) time
                tensor_elem *= pinv_mat[first_index, second_index]  # O(1) time
            new_count += tensor_elem * counts[source_state]  # O(1) time
        return new_count

    # OK
    def col_basis(self, col_idx: int, labels: List[int], pinv_mats: List[np.array]) -> dict:
        """
        Coordinate transformation from pinv v basis to e basis
        O(s * n) time, O(s * n) space

        Arguments
            col_idx: 
            labels: 
        Returns
            col_i: a (col_idx)-th vector of pinvVs (default type is dict)
        """
        col_i = {label: 0 for label in labels}
        for label in labels:  # O(s) times
            tensor_elem = 1.
            for pinv_mat, pos_qubits in zip(pinv_mats, self.mit_pattern):  # O(n) time
                # if completely tensored, then len(pos_clbits) == 1 -> O(1) time
                pos_clbits = [self.qubits_to_clbits[qubit] for qubit in pos_qubits]
                first_index = self.index_of_mat(label, pos_clbits)  # O(1) time
                second_index = self.index_of_mat(col_idx, pos_clbits)  # O(1) time
                tensor_elem *= pinv_mat[first_index, second_index]  # O(1) time
            col_i[label] = tensor_elem

        return col_i

    # OK
    def v_basis(self, col_idx: int, labels: List[int]) -> dict:
        return self.col_basis(col_idx, labels, self.pinvVs)

    # OK
    def choose_vecs(self, state_idx: int, matrices: List[np.array]) -> List[np.array]:
        """
        O(n) time, O(n) space
        
        Arguments
            state_idx: the focusing state
            matrices: list of matrices to be choosed its row
        Returns
            vecs: the list of vector
        """
        vecs = []
        for mat, pos_qubits in zip(matrices, self.mit_pattern):
            pos_clbits = [self.qubits_to_clbits[qubit] for qubit in pos_qubits]
            vecs.append(mat[self.index_of_mat(state_idx, pos_qubits)])
        return vecs

    # OK
    def sum_of_tensored_vector(self, vecs: list) -> int:
        """
        O(n) time, O(n) space
        
        Arguments
            vecs: the vectors consisted of the tensor product
        Returns
            sum_val: the sum of tensored vector of input vectors
        """
        sum_val = 1
        for vec in reversed(vecs):
            sum_val *= sum(vec)
        return sum_val

    # ====================================== for nation etal ==================================== #

    # OK, O(n)
    def change_bit_at_poses(self, key: str, poses: int) -> str:
        for pos in poses:
            key = key[:pos] + "1" + key[pos +
                                        1:] if key[pos] == "0" else key[:pos] + "0" + key[pos+1:]
        return key

    # OK, O(n * 2^d)
    def extend_keys(self, original_keys: set, max_dist: int) -> set:

        extended_key_set = copy.deepcopy(original_keys)

        for key in original_keys:
            n = len(key)
            for d in range(max_dist):
                combs = combinations(range(n), d + 1)
                for comb in combs:
                    new_key = self.change_bit_at_poses(key, comb)
                    extended_key_set.add(new_key)
        return extended_key_set

    def extend_vectors(self, y: dict, keys_to_indices: dict):
        extended_y = np.zeros(len(list(keys_to_indices)))
        for key, value in y.items():
            extended_y[keys_to_indices[key]] = value
        return extended_y

    # OK, O(n)
    def compute_mat_elem(self, row_key: str, col_key: str) -> float:
        """
        used as a subrouting in self.prepare_A_tilde
        O(num_clbits)
        """
        ret = 1.0
        for mat_idx, (i, j) in enumerate(zip(row_key, col_key)):
            ret *= self.cal_matrices[mat_idx][int(i)][int(j)]
        return ret

    # O(n * s * s), where s is the number of extended keys
    def prepare_raw_A_tilde(self, keys: list):
        """
        prepare the numpy 2d matrix with the indices selected from keys
        """
        self.A_tilde = np.zeros((len(keys), len(keys)))
        for i, row_key in enumerate(keys):
            for j, col_key in enumerate(keys):
                self.A_tilde[i][j] = self.compute_mat_elem(row_key, col_key)

    # ====================================== for debug ========================================== #
    
    # OK
    def listup_sigmas(self) -> List[float]:
        """
        List up all the sigma_i values to see their distribution
        For debug and inspection.
        """
        self.run_inv_svd()
        # O(n) time
        return [self.sum_of_tensored_vector(self.choose_vecs(state_idx, self.pinvSigmas)) for state_idx in range(2 ** self.num_clbits)]
    
    # OK
    def listup_deltas(self, counts: dict, shots: int = None, basis: str = "v") -> List[float]:
        """
        List up all the Delta_i values to see their distribution
        For debug and inspection.
        """
        if shots is None:
            shots = sum(counts.values())
        y = {int(state, 2): counts[state] / shots for state in counts}
        self.run_inv_svd()

        sum_of_x = 0
        delta_denom = 0
        deltas = [0 for _ in range(2 ** self.num_clbits)]  # v basis
        for state_idx in range(2 ** self.num_clbits):  # O(2^n) time
            # O(n * shots) time
            sum_of_x += self.mitigate_one_state(state_idx, y)
            sum_of_vi = self.sum_of_tensored_vector(
                self.choose_vecs(state_idx, self.pinvVs))  # O(n) time
            lambda_i = self.sum_of_tensored_vector(
                self.choose_vecs(state_idx, self.pinvSigmas))  # O(n) time
            delta_denom += (sum_of_vi ** 2) / (lambda_i ** 2)  # O(1) time
            deltas[state_idx] = sum_of_vi / (lambda_i ** 2)
        delta_coeff = (1 - sum_of_x) / delta_denom  # O(1) time
        deltas = [delta_coeff * delta_i for delta_i in deltas]
        print("sum_of_x: ", sum_of_x)
        print("delta_denom: ", delta_denom)
        print("delta_coeff: ", delta_coeff)

        self.inspect_sum_of_x(counts)

        if basis == "v":
            return deltas
        else:
            e_deltas = np.zeros(2 ** self.num_clbits)  # O(2^n) space
            V = self.kron_matrices(self.pinvVs)  # O(4^n) time, O(4^n) space
            _ = draw_heatmap(V, list(range(2 ** self.num_clbits)),
                             list(range(2 ** self.num_clbits)))
            for i, vi in enumerate(V.T):
                e_deltas += deltas[i] * vi
            return e_deltas.tolist()

    # OK
    def listup_sum_of_cols(self, counts: dict, shots: int = None, method="naive") -> tuple:
        """
        List up all the Delta_i values to see their distribution
        For debug and inspection.
        """
        if shots is None:
            shots = sum(counts.values())
        y = {int(state, 2): counts[state] / shots for state in counts}
        self.run_inv_svd()

        sum_of_pinv_vis = [self.sum_of_tensored_vector(self.choose_vecs(
            state_idx, self.pinvVs)) for state_idx in range(2 ** self.num_clbits)]
        pinv_lambdas = [self.sum_of_tensored_vector(self.choose_vecs(
            state_idx, self.pinvSigmas)) for state_idx in range(2 ** self.num_clbits)]

        return sum_of_pinv_vis, pinv_lambdas

    # OK
    def inspect_sum_of_x(self, counts: dict, shots: int = None) -> None:
        """
        For debug and inspection.
        """
        if shots is None:
            shots = sum(counts.values())
        y = {int(state, 2): counts[state] / shots for state in counts}
        self.run_inv_svd()

        sum_of_x = 0
        for state_idx in range(2 ** self.num_clbits):  # O(2^n) time
            # O(n * shots) time
            sum_of_x += self.mitigate_one_state(state_idx, y)

        sum_x = 0
        for state_idx, value in y.items():  # O(shots) time
            sum_of_col = self.sum_of_tensored_vector(
                self.choose_vecs(state_idx, self.pinv_matrices))  # O(n) time
            sum_x += sum_of_col * value

        print("sum_of_x == sum_x")
        print("sum_of_x: ", sum_of_x)
        print("sum_x: ", sum_x)
        print(sum_of_x == sum_x)
