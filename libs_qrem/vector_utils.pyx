# reference
# https://ymd_h.gitlab.io/ymd_blog/posts/sharing_memory_between_numpy_and_vector_by_cython/

from libc.stdlib cimport malloc, free
from libcpp.vector cimport vector

cdef vector_to_list_double(vector[double] vec):
    lst = []
    cdef size_t i
    for i in range(vec.size()):
        lst.append(vec[i])
    return lst

cdef vector_to_list_string(vector[string] vec):
    lst = []
    cdef size_t i
    for i in range(vec.size()):
        lst.append(vec[i].decode())
    return lst

cdef vector_to_dict(vector[double] vec, vector[string] indices_to_keys_vector):
    assert vec.size() == indices_to_keys_vector.size()
    dct = dict()
    cdef size_t i
    for i in range(vec.size()):
        dct[indices_to_keys_vector[i]] = vec[i]
    return dct

cdef matrix_to_ndarray(vector[vector[double]] matrix):
    array = []
    cdef size_t i, j
    for i in range(matrix.size()):
        sub_array = []
        for j in range(matrix[0].size()):
            sub_array.append(matrix[i][j])
        array.append(sub_array)
    return array