# reference
# https://ymd_h.gitlab.io/ymd_blog/posts/sharing_memory_between_numpy_and_vector_by_cython/

from libc.stdlib cimport malloc, free
from libcpp.vector cimport vector

cdef class VectorWrapper:
    cdef Py_ssize_t *shape
    cdef Py_ssize_t *strides
    cdef Py_ssize_t itemsize
    cdef int ndim
    cdef int value_dim

    def __cinit__(self):
        self.shape   = <Py_ssize_t*>malloc(sizeof(Py_ssize_t) * 2)
        self.strides = <Py_ssize_t*>malloc(sizeof(Py_ssize_t) * 2)

    def __dealloc__(self):
        free(self.shape)
        free(self.strides)

    cdef void update_size(self):
        self.shape[0] = <Py_ssize_t>(self.vec_size()//self.value_dim)
        self.strides[self.ndim -1] = <Py_ssize_t> self.itemsize

        if self.ndim is 2:
            self.shape[1] = <Py_ssize_t> (self.value_dim)
            self.strides[0] = self.value_dim * <Py_ssize_t> self.itemsize

    cdef void set_buffer(self,Py_buffer *buffer):
        pass

    def __getbuffer__(self, Py_buffer *buffer, int flags):
        # relevant documentation http://cython.readthedocs.io/en/latest/src/userguide/buffer.html#a-matrix-class

        self.update_size()

        self.set_buffer(buffer)
        buffer.len = self.vec_size() * self.itemsize
        buffer.readonly = 0
        buffer.ndim = self.ndim
        buffer.shape = self.shape
        buffer.strides = self.strides
        buffer.suboffsets = NULL
        buffer.itemsize = self.itemsize
        buffer.internal = NULL
        buffer.obj = self

    def __releasebuffer__(self, Py_buffer *buffer):
        pass

cdef class VectorInt(VectorWrapper):
    cdef vector[int] vec

    def __cinit__(self,value_dim=1):
        self.vec = vector[int]()
        self.itemsize = sizeof(int)

        self.ndim = 1 if value_dim is 1 else 2
        self.value_dim = value_dim

    def vec_size(self):
        return self.vec.size()

    cdef void set_buffer(self,Py_buffer* buffer):
        buffer.buf = <void*>(self.vec.data())
        buffer.format = 'i'

cdef class VectorDouble(VectorWrapper):
    cdef vector[double] vec

    def __cinit__(self,value_dim=1):
        self.vec = vector[double]()
        self.itemsize = sizeof(double)

        self.ndim = 1 if value_dim is 1 else 2
        self.value_dim = value_dim

    def vec_size(self):
        return self.vec.size()

    cdef void set_buffer(self,Py_buffer* buffer):
         buffer.buf = <void*>(self.vec.data())
         buffer.format = 'd'

cdef class VectorULong(VectorWrapper):
    cdef vector[size_t] vec

    def __cinit__(self,value_dim=1):
        self.vec = vector[size_t]()
        self.itemsize = sizeof(size_t)

        self.ndim = 1 if value_dim is 1 else 2
        self.value_dim = value_dim

    def vec_size(self):
        return self.vec.size()

    cdef void set_buffer(self,Py_buffer* buffer):
        buffer.buf = <void*>(self.vec.data())
        buffer.format = 'L'

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
    # cdef np.ndarray[double, ndim=2] array
    array = []
    # array = np.zeros((matrix.size(), matrix[0].size()))
    cdef size_t i, j
    for i in range(matrix.size()):
        sub_array = []
        for j in range(matrix[0].size()):
            sub_array.append(matrix[i][j])
            # array[i][j] = matrix[i][j]
        array.append(sub_array)
    return array