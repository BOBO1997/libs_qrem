from typing import Union, List
import numpy as np
from pprint import pprint

class SVDs():
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
        self.num_clbits = num_clbits
        self.cal_matrices = cal_matrices
        self.mit_pattern = mit_pattern if mit_pattern is not None else [
            [i] for i in range(self.num_clbits)]
        self.pinv_matrices = None
        self.svd_matrices = None
        self.Us = None
        self.Sigmas = None
        self.Vs = None
        self.pinv_svd_matrices = None
        self.pinvUs = None
        self.pinvSigmas = None
        self.pinvVs = None

        if meas_layout is None:
            meas_layout = [i for i in range(self.num_clbits)]
        meas_layout = meas_layout[::-1]  # reverse endian

        self.qubits_to_clbits = [-1 for _ in range(max(meas_layout) + 1)]
        for i, qubit in enumerate(meas_layout):
            self.qubits_to_clbits[qubit] = i

        # compute pseudo inverse and svd of all calibartion matrices
        self.run_inv_svd()

    # OK
    def kron_matrices(self, matrices: list):
        """
        Tensor product of given matrices (in the given order)
        """
        ret = matrices[0]
        for mat in matrices[1:]:
            ret = np.kron(ret, mat)
        return ret

    # OK
    def run_inv(self) -> None:
        """
        Prepare inverse of calibration matrices.
        
        Internal Method.
        O(n) time, O(n) space
        """
        if self.pinv_matrices is None:
            self.pinv_matrices = list(map(np.linalg.pinv, self.cal_matrices))

    # OK
    def run_svd(self) -> None:
        """
        Prepare svd of calibration matrices.
        
        Internal Method.
        O(n) time, O(n) space
        """
        if self.svd_matrices is None:
            self.svd_matrices = list(map(np.linalg.svd, self.cal_matrices))
            self.Us = [U for U, _, _ in self.svd_matrices]
            self.Sigmas = [np.diag(Sigma) for _, Sigma, _ in self.svd_matrices]
            self.Vs = [V.T for _, _, V in self.svd_matrices]
            for i in range(len(self.mit_pattern)):
                negative_H = - np.array([[1, 1], [1, -1]])
                if (self.Us[i] == np.diag([1] * (2 ** len(self.mit_pattern[i])))).all() and (self.Vs[i] == np.diag([1] * (2 ** len(self.mit_pattern[i])))).all():
                    print(
                        "changed identity matrix into negative Hadamard matrix at mit_pattern[", i, "]")
                    self.Us[i] = self.kron_matrices(
                        [negative_H] * len(self.mit_pattern[i])) / np.sqrt(2 ** len(self.mit_pattern[i]))
                    self.Vs[i] = self.kron_matrices(
                        [negative_H] * len(self.mit_pattern[i])) / np.sqrt(2 ** len(self.mit_pattern[i]))

    # OK
    def run_inv_svd(self) -> None:
        """
        Do singular value decomposition of all inverse calibration matrices.
        
        Internal Method.
        O(n) time, O(n) space
        """
        self.run_inv()
        self.run_svd()
        if self.pinv_svd_matrices is None:
            # self.pinv_svd_matrices = list(map(np.linalg.svd, self.pinv_matrices))
            # self.pinvUs = [U for U , _, _ in self.pinv_svd_matrices]
            # self.pinvSigmas = [np.diag(Sigma) for _, Sigma, _ in self.pinv_svd_matrices]
            # self.pinvVs = [V.T for _ , _, V in self.pinv_svd_matrices]
            self.pinvUs = list(map(np.linalg.pinv, self.Us))
            self.pinvSigmas = list(map(np.linalg.pinv, self.Sigmas))
            self.pinvVs = list(map(np.linalg.pinv, self.Vs))
