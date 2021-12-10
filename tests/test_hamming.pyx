import cython
cimport cython

# from libcpp.map cimport map
from libcpp.string cimport string
from libcpp.vector cimport vector
from libcpp.set cimport set

cdef void print_vec2d(vector[vector[int]] vec2d):
    for vec1d in vec2d:
        for elem in vec1d:
            print(elem)
        print()

cdef void print_set_of_string(set[string] s):

    # for i in range(s.size()):
        # print(i)
        # print(s[i])

    for i, ss in s:
        print(i)
        print(ss)


def test_extend_keys(key, max_dist):
    cdef set[string] s = extend_keys(key, max_dist)
    print(s)
    # print_vec2d(vec2d)
    print_set_of_string(s)
    print_vec1d()

def aaa(keys, max_dist):
    cdef set[string] cppkeys
    for key in keys:
        print(type(key))
        print(type(key.encode()))
        cppkeys.insert(key.encode())
    print(type(cppkeys))
    print(cppkeys)
    cdef set[string] s = extend_keys(cppkeys, max_dist)
    print(s)
    # extended_keys = {}
    # for ss in s:
        # extended_keys.add(ss.decode())
    # print(extended_keys)