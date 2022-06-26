# An Efficient Quantum Readout Error Mitigation for Sparse Measurement Outcomes of Near-term Quantum Devices

`libs_qrem` is a python package which executes efficient quantum readout error mitigation (QREM) written in C++/Cython.
This package mitigates the readout errors in 65 qubit measurement result of GHZ state from ibmq_brooklyn in few seconds.
- Time Complexity: $O(ns^2)$
- Space Complexity: $O(s^2)$ ( Can be reduced into $O(ns)$ )

# Installation

## Requirements
- [C++/Eigen](https://eigen.tuxfamily.org/index.php?title=Main_Page)
- [Cython](https://cython.org/)

## Install via pip
```sh
pip install git+https://github.com/BOBO1997/libs_qrem
```

c.f.) Reinstall via pip
```sh
pip install --upgrade --force-reinstall git+https://github.com/BOBO1997/libs_qrem
```

## Install manually
```sh
git clone https://github.com/BOBO1997/libs_qrem.git
cd libs_qrem
python setup.py install --record install_record.txt
```

## Uninstall via pip

```sh
pip uninstall libs_qrem
```

## Inplace build and clean

- build
```sh
python setup.py build_ext --inplace
```

- clean
```sh
rm -rf dist/ build/ libs_qrem.egg-info/ libs_qrem.cpython-38-darwin.so libs_qrem/*.cpp
```

# Usage

## Classes

There are 4 different classes that support 4 different QREM methods.

1. `DeltaFilter`: Apply inverse matrix for the vector elements in subspace + correct the vector by adding a correction vector "delta" which is approximated through the solution of Lagrange multiplier + apply SGS algorithm.
2. `LeastNormFilter`: Apply inverse matrix for the vector elements in subspace + apply the solution of least norm problem to compute the closest vector that meets all elements are summed up to 1 + apply SGS algorithm.
3. `MooneyEtalFilter`: Method by [Mooney, White, Hill, Hollenberg, 2021](https://arxiv.org/abs/2101.08946) + apply SGS algorithm.
4. `NationEtalFilter`: Method by [Nation, Kang, Sundaresan, Gambetta, 2021](https://arxiv.org/abs/2108.12518) + apply SGS algorithm.
5. `IgnisFilter`: Apply full inverse matrix under tensor product noise model (such as `TensoredFilter` in [qiskit.ignis](https://github.com/Qiskit/qiskit-ignis)) + apply SGS algorithm.

where SGS algorithm is the algorithm proposed by [Smolin, Gambetta, Smith, 2012](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.108.070502).

Each class inherits the base class `BaseFilter` which has the following methods in order to get access to the internal information.

- Matrices
    - `reduced_A()`: returns a `list` of `list` with `double` elements
    - `normalized_reduced_A()`: returns a `list` of `list` with `double` elements
    - `reduced_inv_A()`: returns a `list` of `list` with `double` elements
    - `exact_one_norm_of_inv_reduced_A()`: returns a `double` value
    - `iterative_one_norm_of_inv_reduced_A()`: returns a `double` value
    - `exact_one_norm_of_reduced_inv_A()`: returns a `double` value
- Mitigated vectors and vectors under the procedure
    - `mitigated_hist()`: returns a `dict` with `str` keys and `double` values
    - `x_s()`: returns a `list` with `double` elements
    - `x_hat()`: returns a `list` with `double` elements
    - `x_tilde()`: returns a `list` with `double` elements
- Sum of vectors
    - `sum_of_x()`: returns a `double` value
    - `sum_of_x_hat()`: returns a `double` value
    - `sum_of_x_tilde()`: returns a `double` value
- Other information
    - `indices_to_keys_vector()`: returns a `list` with `str` elements
    - `times()`: returns a `dict` with `str` keys and `double` values
    - `expval()`: returns a `double` value
    - `mitigation_stddev(norm_type = "exact")`: returns a `double` value

## Examples

Each QREM filter in `libs_qrem` takes the number of qubits and calibration matrices in the following way.
```py
from libs_qrem import LeastNormFilter
meas_filter = LeastNormFilter(n, meas_fitter.cal_matrices)
```
Giving a dictionary typed noisy probability distribution or noisy histogram to `LeastNormFilter.apply()` method, it returns a mitigated probability distribution.
```py
mitigated_hist = meas_filter.apply(noisy_hist)
```

For running calibration circuit, please see [qiskit tutorial](https://qiskit.org/documentation/tutorials/noise/3_measurement_error_mitigation.html).

An easy example code becomes as follows.

```py
from qiskit import QuantumRegister, Aer
from qiskit.ignis.mitigation.measurement import tensored_meas_cal
# prepare calibration circuit (same as qiskit tutorial)
n = 5
qr = qiskit.QuantumRegister(n)
mit_pattern = [[2], [3, 4]]
meas_calibs, state_labels = tensored_meas_cal(mit_pattern=mit_pattern, qr=qr, circlabel='mcal')

# run calibration circuit (same as qiskit tutorial)
backend = qiskit.Aer.get_backend('ibmq_brooklyn')
job = qiskit.execute(meas_calibs, backend=backend, shots=5000)
cal_results = job.result()

from qiskit.ignis.mitigation.measurement import TensoredMeasFitter
meas_fitter = TensoredMeasFitter(cal_results, mit_pattern=mit_pattern)

# Create mitigator instance (corresponds to meas_fitter.filter in the tutorial code.)
# this is very similar to the usage of qiskit.ignis.mitigation modules 
# meas_filter = meas_fitter.filter
from libs_qrem import LeastNormFilter
meas_filter = LeastNormFilter(n, meas_fitter.cal_matrices)

# apply mitigation
# Let `noisy_hist` be a dict variable representing a noisy histogram obtained from `.get_counts()` method in `qiskit.result.Result` instance.
# e.g. `noisy_hist = {"000": 50, "101": 20, "111": 30}`
# Then you can mitigate it as follows.
mitigated_hist = meas_filter.apply(noisy_hist)
```

For detailed examples, see [here](https://github.com/BOBO1997/master_thesis/tree/main/test_libs_qrem).

# Publications

## Paper

This package is based on the paper: [An Efficient Quantum Readout Error Mitigation for Sparse Measurement Outcomes of Near-term Quantum Devices](https://arxiv.org/abs/2201.11046).

Demonstrations in the paper are also stored [here](https://github.com/BOBO1997/master_thesis/tree/main/test_libs_qrem).

## International Conferences
<!-- 
- [The 3rd Workshop on Quanutm Software, Information Processing Society of Japan](https://www.ipsj.or.jp/kenkyukai/event/qs3.html) 
-->
- Efficient Readout Error Mitigation Heuristic for Measurement Outcomes with Few
States (Bo Yang, Rudy Raymond and Shumpei Uno) [AQIS2021, Poster Session B22](http://aqis-conf.org/2021/)

## Cite This Package

```
@misc{yang2022efficient,
      title={An Efficient Quantum Readout Error Mitigation for Sparse Measurement Outcomes of Near-term Quantum Devices}, 
      author={Bo Yang and Rudy Raymond and Shumpei Uno},
      year={2022},
      eprint={2201.11046},
      archivePrefix={arXiv},
      primaryClass={quant-ph}
}
```

Now in the preparation for journal submission.

# Package Contributors

- [Bo Yang](https://github.com/BOBO1997)
- [Dongjia Zhang](https://github.com/tokatoka)
