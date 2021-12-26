import cython
cimport cython

from libcpp.pair cimport pair
from libcpp.vector cimport vector
from libcpp.set cimport set
from libcpp.map cimport map
from libcpp.string cimport string
from libcpp cimport bool

cdef extern from "../cpp/qrem_filter.hpp" namespace "libs_qrem":

    cdef cppclass Args:
        map[string, int] hist
        int d
        double threshold
        string method

    cdef cppclass QREM_Filter:
        double _sum_of_x
        double _sum_of_x_hat
        double _sum_of_x_tilde
        map[string, double] _durations
        vector[double] _x_s
        vector[double] _x_hat
        vector[double] _x_tilde

        vector[vector[double]] _reduced_A
        vector[vector[double]] _reduced_inv_A #  For proposed methods
        vector[vector[double]] _inv_reduced_A #  For Nation et al.
        double _exact_one_norm_of_reduced_inv_A #  For proposed methods
        double _exact_one_norm_of_inv_reduced_A #  For Nation et al.
        double _iterative_one_norm_of_inv_reduced_A #  For Nation et al.

        map[string, double] _mitigated_hist
        vector[string] _indices_to_keys_vector

        QREM_Filter(int num_clbits,
                          vector[vector[vector[double]]] cal_matrices,
                          vector[vector[int]] mit_pattern,
                          vector[int] meas_layout)
        void compute_reduced_A(int size)
        void apply(Args args)

cdef extern from "../cpp/delta_filter.hpp" namespace "libs_qrem":
    cdef cppclass Delta_Filter(QREM_Filter):
        Delta_Filter(int num_clbits,
                          vector[vector[vector[double]]] cal_matrices,
                          vector[vector[int]] mit_pattern,
                          vector[int] meas_layout)
        void apply(Args args)

cdef extern from "../cpp/slsqp_filter.hpp" namespace "libs_qrem":
    cdef cppclass SLSQP_Filter(QREM_Filter):
        SLSQP_Filter(int num_clbits,
                    vector[vector[vector[double]]] cal_matrices,
                    vector[vector[int]] mit_pattern,
                    vector[int] meas_layout)
        void apply(Args args)

cdef extern from "../cpp/least_norm_filter.hpp" namespace "libs_qrem":
    cdef cppclass Least_Norm_Filter(QREM_Filter):
        Least_Norm_Filter(int num_clbits,
                    vector[vector[vector[double]]] cal_matrices,
                    vector[vector[int]] mit_pattern,
                    vector[int] meas_layout)
        void apply(Args args)

cdef extern from "../cpp/mooney_etal_filter.hpp" namespace "libs_qrem":
    cdef cppclass MooneyEtal_Filter(QREM_Filter):
        MooneyEtal_Filter(int num_clbits,
                    vector[vector[vector[double]]] cal_matrices,
                    vector[vector[int]] mit_pattern,
                    vector[int] meas_layout)
        void apply(Args args)

cdef extern from "../cpp/nation_etal_filter.hpp" namespace "libs_qrem":
    cdef cppclass NationEtal_Filter(QREM_Filter):
        NationEtal_Filter(int num_clbits,
                    vector[vector[vector[double]]] cal_matrices,
                    vector[vector[int]] mit_pattern,
                    vector[int] meas_layout)
        void apply(Args args)

cdef extern from "../cpp/sgs_algorithm.hpp" namespace "libs_qrem":
    cdef vector[double] sgs_algorithm(vector[double], bool make_sum_to_one)

cdef extern from "../cpp/eigen_utils.hpp" namespace "libs_qrem":
    cdef void normalize_cols(vector[vector[double]]& A)