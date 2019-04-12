from __future__ import print_function

import logging
from collections import defaultdict
from cython.operator import dereference, postincrement



cimport epic2.src.cpp_read_files as cr

import sys
import numpy as np

from libc.stdint cimport uint32_t, uint16_t

from epic2.src.read_bam import read_bam
from epic2.src.genome_info import sniff

from cython.operator import dereference
from libcpp.algorithm cimport sort as stdsort
from libcpp.map cimport map as cppmap
from libcpp.algorithm cimport unique
from libcpp.vector cimport vector
from libcpp.string cimport string
from libcpp cimport bool
from libcpp.map cimport map

cdef extern from "<algorithm>" namespace "std" nogil:
    OutputIter merge[InputIter1, InputIter2, OutputIter] (InputIter1 first1, InputIter1 last1,
                                                          InputIter2 first2, InputIter2 last2,
                                                          OutputIter result)

cimport cython

import numpy as np

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.initializedcheck(False)
cpdef remove_out_of_bounds_bins(count_dict, chromsizes, bin_size):

    # original code: get_bed_coords() in make_graph_file.py

    cdef:
        # uint32_t[::1] bins
        # uint16_t[::1] counts
        uint32_t chromsize
        size_t i


    for chromosome, (bins, counts) in count_dict.items():
        if not chromosome in chromsizes:
            print("Chromosome {} found in file, but not in the chromsizes dict for the genome.".format(chromosome), file=sys.stderr)
            continue

        chromsize = chromsizes[chromosome]

        # sys.stderr.write(chromosome + "\n")
        # sys.stderr.write("{} {}\n".format(bins[len(bins) - 1], chromsize))
        i = len(bins) - 1
        while (i >= 0 and (bins[i]) > chromsize):
            # print("For {} bin {} is out of bounds ({})\n".format(chromosome, bins[i], chromsize))
            # print("It has a count of ({})\n".format(counts[i]))
            # 52280 - 51862 = 418
            i -= 1

        if i != len(bins) - 1:
            # print("-----")
            # print(chromosome)
            # print("length before:", len(counts))
            count_dict[chromosome] = bins[:i+1], counts[:i+1]
            # print("length after:", len(counts[:i + 1]))



@cython.boundscheck(False)
@cython.wraparound(False)
@cython.initializedcheck(False)
cpdef count_reads_per_bin(tags):

    cdef:
        uint32_t[::1] bins
        uint16_t[::1] counts
        Vector32 v
        uint32_t i
        uint32_t last = 1
        uint32_t current
        uint16_t count = 0
        uint32_t nparsed
        uint32_t last_two_equal

    bins_counts = dict()
    for k, v in tags.items():
        nparsed = 0
        count = 0

        bin_arr = np.ones(len(v), dtype=np.uint32)
        bins = bin_arr
        count_arr = np.ones(len(v), dtype=np.uint16)
        counts = count_arr

        if len(v) >= 1:
            last = v.wrapped_vector[0]

        for i in range(0, len(v)):
            current = v.wrapped_vector[i]

            if current != last:
                bins[nparsed] = last
                counts[nparsed] = count
                last = current
                count = 1
                nparsed += 1
            else:
                count += 1

        last_two_equal = v.wrapped_vector[len(v) - 2] == v.wrapped_vector[len(v) - 1]
        if last_two_equal:
            bins[nparsed] = v.wrapped_vector[len(v) - 1]
            counts[nparsed] = count
            nparsed += 1
        else:
            bins[nparsed] = v.wrapped_vector[len(v) - 1]
            counts[nparsed] = 1
            nparsed += 1


        #chrM    0       16399   52280   87916   1       0.419037470139  1
        # sys.stderr.write("Found {} in count_reads_per_bin\n".format(sum(count_arr[:nparsed])))
        bins_counts[k] = (bin_arr[:nparsed], count_arr[:nparsed])

    return bins_counts


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.initializedcheck(False)
cdef class Vector32:

    cdef vector[uint32_t] wrapped_vector

    cdef push_back(self, uint32_t num):
        self.wrapped_vector.push_back(num)

    def sort(self):
        stdsort(self.wrapped_vector.begin(), self.wrapped_vector.end())

    def unique(self):
        self.wrapped_vector.erase(unique(self.wrapped_vector.begin(), self.wrapped_vector.end()), self.wrapped_vector.end())

    def __str__(self):
        return "[" + ", ".join([str(i) for i in self.wrapped_vector]) + "]"

    def __repr__(self):
        return str(self)

    def __len__(self):
        return self.wrapped_vector.size()

    def __iter__(self):
        # slow, only implemented to ease testing
        return (v for v in self.wrapped_vector)

    def merge(self, Vector32 other):

        cdef vector[uint32_t] o = vector[uint32_t](len(self) + len(other))
        merge(self.wrapped_vector.begin(), self.wrapped_vector.end(),
              other.wrapped_vector.begin(), other.wrapped_vector.end(),
              o.begin())

        cdef Vector32 output = Vector32()
        output.wrapped_vector = o

        return output




@cython.cdivision(True)
@cython.boundscheck(False)
@cython.wraparound(False)
@cython.initializedcheck(False)
cpdef files_to_bin_counts(files, args, datatype):

    cdef:
        uint32_t bin_size = args["bin_size"]
        uint32_t half_fragment_size = args["fragment_size"] / 2
        uint32_t drop_duplicates = args["drop_duplicates"]
        Vector32 v
        Vector32 v2
        long[::1] bin_arr
        bytes py_bytes
        char* c_string
        cr.genome_map cpp_tags
        map[cr.key, cr.intvec].iterator it
        bool paired_end = False
        str file_format


    sum_tags = defaultdict(list)
    logging.info("Parsing {} file(s):".format(datatype))
    sys.stderr.flush()
    tags = dict()

    # uint32_t mapq, uint64_t required_flag, uint64_t filter_flag

    for f in files:

        paired_end = False

        logging.info("  " + f)
        sys.stderr.flush()

        py_bytes = f.encode()
        c_string = py_bytes

        file_format = sniff(f)

        if file_format == "bed":
            cpp_tags = cr.read_bed(c_string, drop_duplicates)
        elif file_format == "bedpe":
            cpp_tags = cr.read_bedpe(c_string, drop_duplicates)
            paired_end = True
        elif file_format == "bam": # sam also okay here
            cpp_tags = read_bam(f, drop_duplicates, args["mapq"], args["required_flag"], args["filter_flag"])
        elif file_format == "bed.gz":
            cpp_tags = cr.read_bed_gz(c_string, drop_duplicates)
        elif file_format == "bedpe.gz":
            cpp_tags = cr.read_bedpe_gz(c_string, drop_duplicates)
            paired_end = True

        it = cpp_tags.begin();

        while (it != cpp_tags.end()):
            chromosome = dereference(it).first.first.decode()

            if chromosome not in args["chromsizes_"].keys():
                postincrement(it)
                continue

            strand = chr(dereference(it).first.second)

            v = Vector32()
            v.wrapped_vector = dereference(it).second
            tags[chromosome, strand] = v

            postincrement(it)

        if not paired_end:
            for (chromosome, strand), v in tags.items():

                v.sort() # needs to be done again, since extracting the 5' end might make tags wrong order

                if strand == "+":
                    for i in range(len(v)):
                        v.wrapped_vector[i] = v.wrapped_vector[i] + half_fragment_size
                        v.wrapped_vector[i] = v.wrapped_vector[i] - (v.wrapped_vector[i] % bin_size)
                else:
                    i = 0
                    while i < len(v) and v.wrapped_vector[i] < half_fragment_size:
                        v.wrapped_vector[i] = 0
                        i += 1

                    for i in range(i, len(v)):
                        v.wrapped_vector[i] = v.wrapped_vector[i] - half_fragment_size
                        v.wrapped_vector[i] = v.wrapped_vector[i] - (v.wrapped_vector[i] % bin_size)

                # for i in range(len(v)):
                #     v.wrapped_vector[i] = v.wrapped_vector[i] - (v.wrapped_vector[i] % bin_size)

                # sys.stderr.write("Found {} for {} {}\n".format(i, chromosome, strand))
                if chromosome not in sum_tags:
                    sum_tags[chromosome] = v
                else:
                    v2 = sum_tags[chromosome]
                    sum_tags[chromosome] = v.merge(v2)

        else:

            for (chromosome, strand), v in tags.items():

                v.sort() # needs to be done again, since extracting the 5' end might make tags wrong order

                for i in range(len(v)):
                    v.wrapped_vector[i] = v.wrapped_vector[i] - (v.wrapped_vector[i] % bin_size)

                if chromosome not in sum_tags:
                    sum_tags[chromosome] = v
                else:
                    v2 = sum_tags[chromosome]
                    sum_tags[chromosome] = v.merge(v2)


    sys.stderr.flush()

    bins_counts = count_reads_per_bin(sum_tags)

    # import pandas as pd
    # bins, counts = [pd.Series(s) for s in bins_counts["chrY"]]
    # print(bins.tail(), file=sys.stderr)
    # print(counts.tail(), file=sys.stderr)

    count = sum([sum(counts) for _, counts in bins_counts.values()])
    remove_out_of_bounds_bins(bins_counts, args["chromsizes_"], bin_size)

    # bins, counts = [pd.Series(s) for s in bins_counts["chrY"]]
    # print(bins.tail())
    # print(counts.tail())

    # print(bins.tail(), file=sys.stderr)
    # print(counts.tail(), file=sys.stderr)

    return bins_counts, count



cpdef add_reads_to_dict(f, chromosomes):

    genome = dict()
    cdef Vector32 v

    for line in open(f):
        chromosome, left, right, _, _, strand = line.split()

        if chromosome not in chromosomes:
            continue

        if strand == "+":
            five_end = <uint32_t> int(left)
        else:
            five_end = <uint32_t> int(right)

        if (chromosome, strand) in genome:
            v = genome[chromosome, strand]
            v.wrapped_vector.push_back(five_end)
        else:
            v = Vector32()
            v.wrapped_vector.push_back(five_end)
            genome[chromosome, strand] = v

    return genome



######## this was actually slower!

# cpdef add_reads_to_dict(f):

#     genome = dict()
#     cdef:
#         Vector v
#         FILE *f_handle
#         char chromosome [10]
#         char strand [2]
#         uint32_t left
#         uint32_t right

#     fp = fopen(f.encode(), "r")

#     while (
#             fscanf(fp, "%s\t%d\t%d\t%*s\t%*d\t%s\n", chromosome, &left, &right, strand) != EOF
#     ):


#         # pruint32_t("----")
#         # pruint32_t("chromosome is", chromosome)
#         # pruint32_t("strand is", strand)
#         # pruint32_t("left is", left)
#         # pruint32_t("right is", right)
#         if strand == b"+":
#             five_end = left
#         else:
#             five_end = right

#         # pruint32_t("five end is ", five_end)

#         if (chromosome, strand) in genome:
#             v = genome[chromosome, strand]
#             v.wrapped_vector.push_back(five_end)
#         else:
#             v = Vector()
#             v.wrapped_vector.push_back(five_end)
#             genome[chromosome, strand] = v

#     return genome
