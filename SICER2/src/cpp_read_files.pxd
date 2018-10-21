
from libcpp.pair cimport pair
from libcpp.string cimport string
from libcpp.vector cimport vector
from libcpp.map cimport map
from libc.stdint cimport uint32_t
# cdef from namespace "std":

ctypedef pair[string, char] key
ctypedef vector[uint32_t] intvec


cdef extern from "<unordered_map>" namespace "std":
    cdef cppclass unordered_map[T, T]:
        pass

ctypedef map[key, intvec] genome_map


cdef extern from "read_files.hpp":
    genome_map read_bed(const char *)
