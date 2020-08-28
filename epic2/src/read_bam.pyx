from pysam.libcalignmentfile cimport AlignmentFile, AlignedSegment
from cython.operator import dereference, postincrement
from libc.stdint cimport uint32_t, int32_t, uint64_t

from libcpp.string cimport string
from libcpp.pair cimport pair
from libcpp.map cimport map
from libcpp.utility cimport move

cimport epic2.src.cpp_read_files as cr

cdef extern from "<utility>" namespace "std" nogil:
    pair[T,U] make_pair[T,U](T&,U&)

ctypedef pair [int, char] intkey
ctypedef map[intkey, cr.interval_vector] genome_intervals_int


cdef cr.genome_intervals map_chromosome_id_to_name(genome_intervals_int genome, AlignmentFile samfile):
    cdef:
        string chromosome
        char forward = b"+"
        char reverse = b"-"
        char strand
        cr.key chrom_strand
        cr.genome_intervals genome_fixed
        map[intkey, cr.interval_vector].iterator it

    it = genome.begin()
    while it != genome.end():
        chromosome_id = dereference(it).first.first
        chromosome = samfile.get_reference_name(chromosome_id).encode("utf-8")

        strand = dereference(it).first.second

        if strand == forward:
            chrom_strand_fixed = make_pair(<string>chromosome, <char>forward)
        else:
            chrom_strand_fixed = make_pair(<string>chromosome, <char>reverse)

        genome_fixed[chrom_strand_fixed] = move(dereference(it).second)
        postincrement(it)
    return genome_fixed


cpdef cr.genome_map read_bampe(filename, uint32_t drop_duplicates, uint32_t mapq, uint64_t required_flag, uint64_t filter_flag):

    cdef:
        uint32_t flag
        int32_t template_length
        char forward = b"+"
        char reverse = b"-"
        int chromosome_id
        intkey chrom_strand
        cr.interval interval
        genome_intervals_int genome
        AlignmentFile samfile
        AlignedSegment segment

    samfile = AlignmentFile(filename, "rb")

    # Consider every sequenced fragment exactly once, process only correctly paired 5` reads
    # https://broadinstitute.github.io/picard/explain-flags.html
    required_flag = required_flag | 3

    for segment in samfile.fetch():
        flag = segment.flag

        if segment.mapping_quality < mapq:
            continue
        if (flag & required_flag) != required_flag:
            continue
        if (flag & filter_flag) != 0:
            continue

        template_length = segment.template_length

        # template_length must be > 0 for 5` read
        if template_length <= 0:
            continue

        interval.start = segment.reference_start
        if interval.start < 0:
            continue
        interval.end = interval.start + template_length - 1

        is_reverse = flag & 0x10
        chromosome_id = segment.reference_id

        if is_reverse:
            chrom_strand = make_pair(<int>chromosome_id, <char>reverse)
            genome[chrom_strand].push_back(interval)
        else:
            chrom_strand = make_pair(<int>chromosome_id, <char>forward)
            genome[chrom_strand].push_back(interval)
    return cr.intervals_to_midpoint(map_chromosome_id_to_name(genome, samfile), drop_duplicates)


cpdef cr.genome_map read_bam(filename, uint32_t drop_duplicates, uint32_t mapq, uint64_t required_flag, uint64_t filter_flag):

    cdef:
        uint32_t flag
        int32_t start
        int32_t end
        char forward = b"+"
        char reverse = b"-"
        int chromosome_id
        intkey chrom_strand
        cr.interval interval
        genome_intervals_int genome
        AlignmentFile samfile
        AlignedSegment segment

    samfile = AlignmentFile(filename, "rb")


    for segment in samfile:
        flag = segment.flag

        # https://broadinstitute.github.io/picard/explain-flags.html

        if segment.mapping_quality < mapq:
            continue
        if (flag & required_flag) != required_flag:
            continue
        if (flag & filter_flag) != 0:
            continue

        is_reverse = flag & 0x10

        start = segment.reference_start
        end = start + segment.alen

        if start < 0 or end < 0:
            continue

        # It is not clear why + 1 is used here for start and end
        # should be only end - 1 to map intervals [) -> [] ?
        interval = [<uint32_t> (start + 1), <uint32_t> (end + 1)]

        chromosome_id = segment.reference_id
        if is_reverse:
            chrom_strand = make_pair(<int>chromosome_id, <char>reverse)
            genome[chrom_strand].push_back(interval)
        else:
            chrom_strand = make_pair(<int>chromosome_id, <char>forward)
            genome[chrom_strand].push_back(interval)

    return cr.intervals_to_five_end(map_chromosome_id_to_name(genome, samfile), drop_duplicates)
