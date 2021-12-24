import cython
cimport cython

from libcpp.pair cimport pair
from libcpp.vector cimport vector
from libcpp.set cimport set
from libcpp.map cimport map
from libcpp.string cimport string
from libcpp cimport bool

cdef extern from "../cpp/qrem_filter_delta.hpp" namespace "libs_qrem":
    cdef cppclass QREM_Filter_Delta:
        double _sum_of_x
        double _sum_of_x_hat
        double _sum_of_x_tilde
        map[string, double] _durations
        vector[double] _x_s
        vector[double] _x_hat
        vector[double] _x_tilde
        double _one_norm
        vector[vector[double]] _reduced_A
        vector[vector[double]] _reduced_inv_A
        map[string, double] _mitigated_hist
        vector[string] _indices_to_keys_vector

        QREM_Filter_Delta(int num_clbits,
                          vector[vector[vector[double]]] cal_matrices,
                          vector[vector[int]] mit_pattern,
                          vector[int] meas_layout)
        void compute_reduced_A(vector[string]& indices_to_keys_vector)
        void apply(map[string, int] hist, int d)

cdef extern from "../cpp/qrem_filter_nlp.hpp" namespace "libs_qrem":
    cdef cppclass QREM_Filter_Nlp:
        double _sum_of_x
        double _sum_of_x_hat
        double _sum_of_x_tilde
        vector[string] _indices_to_keys_vector
        map[string, double] _durations
        vector[double] _x_s
        vector[double] _x_hat
        vector[double] _x_tilde
        double _one_norm
        vector[vector[double]] _reduced_A
        vector[vector[double]] _reduced_inv_A
        vector[string] _indices_to_keys_vector

        QREM_Filter_Nlp(int num_clbits,
                    vector[vector[vector[double]]] cal_matrices,
                    vector[vector[int]] mit_pattern,
                    vector[int] meas_layout)
        void compute_reduced_A(vector[string]& indices_to_keys_vector)
        void apply(map[string, int] hist, int d)

cdef extern from "../cpp/qrem_filter_mooney_etal.hpp" namespace "libs_qrem":
    cdef cppclass QREM_Filter_MooneyEtal:
        double _sum_of_x
        double _sum_of_x_hat
        double _sum_of_x_tilde
        map[string, double] _durations
        vector[double] _x_s
        vector[double] _x_hat
        vector[double] _x_tilde
        double _one_norm
        vector[vector[double]] _reduced_A
        vector[vector[double]] _reduced_inv_A
        map[string, double] _mitigated_hist
        vector[string] _indices_to_keys_vector

        QREM_Filter_MooneyEtal(int num_clbits,
                    vector[vector[vector[double]]] cal_matrices,
                    vector[vector[int]] mit_pattern,
                    vector[int] meas_layout)
        void compute_reduced_A(vector[string]& indices_to_keys_vector)
        void apply(map[string, int] hist, int d, double threshold)

cdef extern from "../cpp/qrem_filter_lnp.hpp" namespace "libs_qrem":
    cdef cppclass QREM_Filter_Lnp:
        double _sum_of_x
        double _sum_of_x_hat
        double _sum_of_x_tilde
        map[string, double] _durations
        vector[double] _x_s
        vector[double] _x_hat
        vector[double] _x_tilde
        double _one_norm
        vector[vector[double]] _reduced_A
        vector[vector[double]] _reduced_inv_A
        map[string, double] _mitigated_hist
        vector[string] _indices_to_keys_vector

        QREM_Filter_Lnp(int num_clbits,
                    vector[vector[vector[double]]] cal_matrices,
                    vector[vector[int]] mit_pattern,
                    vector[int] meas_layout)
        void compute_reduced_A(vector[string]& indices_to_keys_vector)
        void apply(map[string, int] hist, int d)

cdef extern from "../cpp/sgs_algorithm.hpp" namespace "libs_qrem":
    cdef vector[double] sgs_algorithm(vector[double], bool make_sum_to_one)