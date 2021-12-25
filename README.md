# An Efficient Quantum Readout Error Mitigation Heuristic for Sparse Probability Distribution

`libs_qrem` is a python package which executes efficient quantum readout error mitigation (QREM) written in C++/Cython.
This package mitigates the readout errors in 65 qubit measurement result of GHZ state from ibmq_brooklyn in few seconds.
- Time Complexity: $O(s^2 + ns)$
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
python3 setup.py install --record install_record.txt
```

## Uninstall via pip

```sh
pip uninstall libs_qrem
```

## Inplace build and clean

- build
```sh
python3 setup.py build_ext --inplace
```

- clean
```sh
rm -r dist/ build/ libs_qrem.egg-info/ libs_qrem.cpython-38-darwin.so libs_qrem/*cpp
```

# Usage

## Classes

There are 4 different classes that support 4 different QREM methods.

1. `DeltaFilter`: Apply inverse matrix for the vector elements in subspace + correct the vector by adding a correction vector "delta" which is approximated through the solution of Lagrange multiplier + apply SGS algorithm.
2. `SLSQPFilter`: Apply inverse matrix for the vector elements in subspace + run `scipy.optimize.minimize` to find the closest vector that meets all elements are summed up to 1 + apply SGS algorithm.
3. `LeastNormFilter`: Apply inverse matrix for the vector elements in subspace + apply the solution of least norm problem to compute the closest vector that meets all elements are summed up to 1 + apply SGS algorithm.
4. `MooneyEtalFilter`: Method by [Mooney, White, Hill, Hollenberg, 2021](https://arxiv.org/abs/2101.08946) + apply SGS algorithm.
5. `NationEtalFilter`: Method by [Nation, Kang, Sundaresan, Gambetta, 2021](https://arxiv.org/abs/2108.12518) + apply SGS algorithm.

where SGS algorithm is the algorithm proposed by [Smolin, Gambetta, Smith, 2012](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.108.070502).

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
mit_pattern = [[2], [3], [4]] # only completely local noise setting is supported currently, but we will extend our mitigator soon
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

For a detailed example, see [here](https://github.com/BOBO1997/qip2021_poster549/blob/main/master_thesis/qrem_benchmarkings/ghz_states/brooklyn_main8192_mit8192/mitigation.ipynb).

# Publications

## Paper

## International Conferences
<!-- 
- [The 3rd Workshop on Quanutm Software, Information Processing Society of Japan](https://www.ipsj.or.jp/kenkyukai/event/qs3.html) 
-->
- Efficient Readout Error Mitigation Heuristic for Measurement Outcomes with Few
States (Bo Yang, Rudy Raymond and Shumpei Uno) [AQIS2021, Poster Session B22](http://aqis-conf.org/2021/)

## Cite this package

To be updated