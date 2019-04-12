
from cython.operator import dereference, postincrement
from libc.stdint cimport uint32_t, uint16_t, int32_t, int64_t, uint64_t

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

ctypedef struct interval:
    uint32_t start
    uint32_t end

ctypedef vector[uint32_t] intvec
ctypedef vector[interval] interval_vector
ctypedef pair [string, char] key
ctypedef pair [int, char] intkey
ctypedef map[intkey, intvec] genome_map_int
ctypedef map[key, intvec] genome_map
ctypedef map[intkey, interval_vector] genome_intervals_int
ctypedef map[key, intvec] genome_intervals


# include "<htslib>"


# cpdef read_bam(filename):

#     cdef:


cdef uint32_t compare_start_end(interval lhs, interval rhs):
  if (lhs.start < rhs.start):
    return <uint32_t> 1
  elif (rhs.start < lhs.start):
      return <uint32_t> 0
  elif (lhs.end < rhs.end):
      return <uint32_t> 1
  else:
    return <uint32_t> 0


# uint32_t compare_by_start_end(const interval lhs, const interval rhs){
#   if (lhs.start < rhs.start){
#     return 1;
#   } else if (rhs.start < lhs.start){
#     return 0;
#   } else if (lhs.end < rhs.end){
#     return 1;
#   } else {
#     return 0;
#   };

# }

cdef uint32_t start_end_equal(interval lhs, interval rhs):
  if ((lhs.start == rhs.start) and (lhs.end == rhs.end)):
      return <uint32_t> 1
  else:
      return <uint32_t> 0




# uint32_t start_end_equal(const interval lhs, const interval rhs){
#   if ((lhs.start == rhs.start) && (lhs.end == rhs.end)){
#       return 1;
#   } else {
#     return 0;
#   }
# }


import pysam


from pysam.libcalignmentfile cimport AlignmentFile, AlignedSegment


cpdef read_bam(filename, uint32_t drop_duplicates, uint32_t mapq, uint64_t required_flag, uint64_t filter_flag):

    cdef:
        uint32_t flag
        int32_t start
        int32_t end
        int32_t length
        uint32_t is_strand
        char forward = b"+"
        char reverse = b"-"
        char strand
        string chromosome
        int chromosome_id
        intkey chrom_strand
        key chrom_strand_fixed
        uint32_t five_end
        genome_intervals_int genome
        genome_intervals genome_fixed
        interval _interval
        interval_vector intervals
        intvec five_ends
        # AlignmentFile samfile
        # AlignedSegment a
        uint32_t i = 0

    samfile = pysam.AlignmentFile(filename, "rb")


    for a in samfile:
        flag = a.flag

        # https://broadinstitute.github.io/picard/explain-flags.html

        if a.mapping_quality < mapq:
            continue
        if (flag & required_flag) != required_flag:
            continue
        if (flag & filter_flag) != 0:
            continue

        is_reverse = flag & 0x10

        start = a.reference_start

        end = start + a.alen

        if start < 0 or end < 0:
            continue

        _interval = [<uint32_t> (start + 1), <uint32_t> (end + 1)]

        chromosome_id = a.reference_id
        if is_reverse:
            chrom_strand = make_pair(<int>chromosome_id, <char>reverse)
            genome[chrom_strand].push_back(_interval)
        else:
            chrom_strand = make_pair(<int>chromosome_id, <char>forward)
            genome[chrom_strand].push_back(_interval)


    it = genome.begin();

    while (it != genome.end()):
        chromosome_id = dereference(it).first.first
        chromosome = samfile.get_reference_name(chromosome_id).encode("utf-8")

        strand = dereference(it).first.second
        intervals = dereference(it).second
        five_ends = intvec()

        if drop_duplicates:

            stdsort(intervals.begin(), intervals.end(), compare_start_end)
            intervals.erase(unique(intervals.begin(), intervals.end(), start_end_equal), intervals.end())

        if chr(strand) == "+":
            chrom_strand_fixed = make_pair(<string>chromosome, <char>forward)
            for i in range(intervals.size()):
                five_ends.push_back(intervals[i].start)
        else:
            chrom_strand_fixed = make_pair(<string>chromosome, <char>reverse)
            for i in range(intervals.size()):
                five_ends.push_back(intervals[i].end)

        genome_fixed[chrom_strand_fixed] = five_ends

        postincrement(it)

    return genome_fixed
