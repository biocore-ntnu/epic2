from libcpp.pair cimport pair
from libcpp.string cimport string
from libcpp.vector cimport vector
from libcpp.map cimport map
from libc.stdint cimport uint32_t

ctypedef pair[string, char] key
ctypedef vector[uint32_t] intvec
ctypedef map[key, intvec] genome_map

cdef struct interval:
        uint32_t start
        uint32_t end

ctypedef vector[interval] interval_vector
ctypedef map[key, intvec] genome_tags
ctypedef map[key, interval_vector] genome_intervals

cdef genome_map read_bed(const char *, uint32_t)
cdef genome_map read_bed_gz(const char *, uint32_t)
cdef genome_map read_bedpe(const char *, uint32_t)
cdef genome_map read_bedpe_gz(const char *, uint32_t)

cdef genome_map intervals_to_midpoint(genome_intervals genome, uint32_t drop_duplicates)
cdef genome_map intervals_to_five_end(genome_intervals genome, uint32_t drop_duplicates)
