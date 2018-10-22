
from cython.operator import dereference, postincrement
from libc.stdint cimport uint32_t, uint16_t, int32_t, int64_t

from cython.operator import dereference
from libcpp.algorithm cimport sort as stdsort
from libcpp.map cimport map as cppmap
from libcpp.algorithm cimport unique
from libcpp.vector cimport vector
from libcpp.string cimport string
from libcpp.pair cimport pair
from libcpp.map cimport map

cdef extern from "<utility>" namespace "std" nogil:
    pair[T,U] make_pair[T,U](T&,U&)

ctypedef vector[uint32_t] intvec
ctypedef pair [string, char] key
ctypedef pair [int, char] intkey
ctypedef map[intkey, intvec] genome_map_int
ctypedef map[key, intvec] genome_map


# include "<htslib>"


# cpdef read_bam(filename):

#     cdef:


import pysam


# from pysam.libcalignedsegment cimport AlignedSegment

# ctypedef struct Alignment:
#     int32_t start
#     uint32_t alen
#     uint32_t reference_id
#     int64_t flag

cpdef read_bam(filename):

    cdef:
        uint32_t flag
        int32_t start
        int32_t end
        int32_t length
        uint32_t is_strand
        char forward = "+"
        char reverse = "-"
        char strand
        string chromosome
        int chromosome_id
        intkey chrom_strand
        key chrom_strand_fixed
        uint32_t five_end
        genome_map_int genome
        genome_map genome_fixed
        # Alignment a

    if filename.endswith(".bam"):
        samfile = pysam.AlignmentFile(filename, "rb")
    else:
        samfile = pysam.AlignmentFile(filename, "r")


    for a in samfile:
        flag = a.flag

        # https://broadinstitute.github.io/picard/explain-flags.html
        if flag & 0x4 or flag & 0x200 or flag & 0x400:
            continue

        is_reverse = flag & 0x10

        start = a.reference_start

        end = start + a.alen

        if start < 0 or end < 0:
            continue

        chromosome_id = a.reference_id
        if is_reverse:
            chrom_strand = make_pair(<int>chromosome_id, <char>reverse)
            genome[chrom_strand].push_back(end + 1)
        else:
            chrom_strand = make_pair(<int>chromosome_id, <char>forward)
            genome[chrom_strand].push_back(start + 1)


    it = genome.begin();

    while (it != genome.end()):
        chromosome_id = dereference(it).first.first
        chromosome = samfile.get_reference_name(chromosome_id).encode("utf-8")

        strand = dereference(it).first.second

        if chr(strand) == "+":
            chrom_strand_fixed = make_pair(<string>chromosome, <char>forward)
            genome_fixed[chrom_strand_fixed] = dereference(it).second
        else:
            chrom_strand_fixed = make_pair(<string>chromosome, <char>reverse)
            genome_fixed[chrom_strand_fixed] = dereference(it).second

        postincrement(it)

    return genome_fixed
