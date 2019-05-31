from __future__ import print_function

from scipy.stats import poisson, rankdata
import sys
import logging
import cython
from libcpp cimport bool

from libcpp.algorithm cimport sort as stdsort
from libcpp.string cimport string
from libc.math cimport log2
import numpy as np

from libc.stdint cimport uint32_t, uint16_t, int32_t
from epic2.src.statistics import compute_window_score
from libcpp.vector cimport vector
# from libc.math cimport min as cmin

from natsort import natsorted


cdef struct island:
    string chromosome
    uint32_t start
    uint32_t end
    uint32_t chip_count
    uint32_t input_count
    double score
    double p_value
    double fold_change


# cdef bool compare_p_value(const island &lhs, const island &rhs):

#     return lhs.p_value < rhs.p_value



cdef class IslandVector():

    cdef vector[island] wrapped_vector

    cdef push_back(self, island i):
        self.wrapped_vector.push_back(i)

    def __str__(self):
        return "[" + ", ".join([str(i) for i in self.wrapped_vector]) + "]"

    def __repr__(self):
        return str(self)

    def __len__(self):

        return int(self.wrapped_vector.size())

    def __iter__(self):
        # slow, only implemented to ease testing
        return (v for v in self.wrapped_vector)

    # def p_value_sort(self):
    #     stdsort(self.wrapped_vector.begin(), self.wrapped_vector.end(), &compare_p_value)



@cython.boundscheck(False)
@cython.wraparound(True)
@cython.initializedcheck(False)
def find_islands(bins_counts, int gaps_allowed, int bin_size, float score_threshold, uint32_t min_tags_in_window, float average_window_readcount):

    _poisson = poisson(average_window_readcount)

    cdef:
        double slightly_less = score_threshold - 0.0000000001
        double score
        int i
        int length
        uint32_t count
        uint32_t _bin
        uint32_t[::1] bins = np.ones(1, dtype=np.uint32)
        uint16_t[::1] counts = np.ones(1, dtype=np.uint16)
        int32_t dist # distance can be negative on the very first island in each chromosome, hence not uint
        island current_island
        IslandVector v
        int32_t distance_allowed = (gaps_allowed * bin_size) + 2

    chromosomes = natsorted(set(bins_counts.keys()))

    island_dict = dict()

    for chromosome in chromosomes:
        v = IslandVector()
        # print("-------")
        # print(chromosome)
        i = 0

        bins, counts = bins_counts[chromosome]

        length = len(bins)

        # print("length", length)
        while i < length and (counts[i] < min_tags_in_window):
            i += 1

        # print("first i, counts", i, counts[i])

        count = 0
        _bin = bins[i]
        score = compute_window_score(count, _poisson)

        chromosome = chromosome.encode("utf-8")
        current_island = [chromosome, _bin, _bin + bin_size - 1, 0, 0, score, 0, 0]
        while i < length:

            if (counts[i] < min_tags_in_window):
                i += 1
                continue

            _bin = bins[i]


            dist = _bin - current_island.end

            # # _print = chromosome == b"chrM"# and _bin >= 107040800
            # if _print:
            #     print("_bin", _bin, "i", i)
            #     print("_count", counts[i])
            #     print("dist", dist)
            #     print("current_island", current_island)
            #     print("dist allowed", distance_allowed)

            if dist <= distance_allowed:

                # if _print:
                #     print("we are in dist")

                current_island.end = bins[i] + bin_size - 1
                # current_island.chip_count += counts[i]
                current_island.score += compute_window_score(counts[i], _poisson)
            else:
                if current_island.score > slightly_less: # and current_island.chip_count > min_tags_in_windw:
                    # if _print:
                    #     print("pushed back", current_island)
                    v.push_back(current_island)

                # here need to iterate until we get to a
                # bin with at least min_tags_in_window number of reads
                while i < length and (counts[i] < min_tags_in_window):
                    i += 1
                # print("second i, counts", i, counts[i])

                # if _print:
                #     print("current_island_before", current_island)

                score = compute_window_score(counts[i], _poisson)
                current_island = [chromosome, _bin, bins[i] + bin_size - 1, 0, 0, score, 0, 0]

                # if _print:
                #     print("current_island_after", current_island)
                # else: # then length == i


            i += 1

        # takes care of last or potential single bin
        if current_island.score > slightly_less: # and current_island.chip_count > min_tags_in_windw:
            # print("current_island last", current_island)
            v.push_back(current_island)

        island_dict[chromosome.decode()] = v
        # print(chromosome, v.wrapped_vector[0])

    return island_dict


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.initializedcheck(False)
@cython.cdivision(False)
def add_chip_count_to_islands(islands, c_bins_counts, args):

    cdef:
        int i
        int j = 0
        island current_island
        IslandVector _islands
        IslandVector updated_islands
        uint32_t[::1] bins = np.ones(1, dtype=np.uint32)
        uint16_t[::1] counts = np.ones(1, dtype=np.uint16)


    new_islands = {}
    chromosomes = natsorted(set(islands.keys()))

    for chromosome in chromosomes:

        updated_islands = IslandVector()
        j = 0
        _islands = islands[chromosome]
        # print("chr", chromosome, len(_islands))
        # print(chromosome, c_bins_counts)

        if chromosome not in c_bins_counts:
            # print("chromosome", chromosome, "not in", c_bins_counts)
            continue

        bins, counts = c_bins_counts[chromosome]
        number_bins = len(bins)

        for i in range(len(_islands)):
            # _island = &_islands.wrapped_vector[i]
            _island = _islands.wrapped_vector[i]

            # not overlapping
            while j < number_bins and (bins[j] < _island.start):
                # print("not overlapping:")
                # print("bins[j]", bins[j], "counts[j]", counts[j], "_island.start", _island)
                j += 1

            # overlapping
            while j < number_bins and (bins[j] < _island.end and bins[j] >= _island.start):
                # print("overlapping:")
                # print("bins[j]", bins[j], "counts[j]", counts[j], "_island.start", _island)
                _island.chip_count += counts[j]
                j += 1

            updated_islands.wrapped_vector.push_back(_island)

        new_islands[chromosome] = updated_islands

        if not args.get("df"):
            del c_bins_counts[chromosome]

    return new_islands

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.initializedcheck(False)
@cython.cdivision(True)
def compute_fdr(islands, b_bins_counts, int chip_library_size, int control_library_size, double effective_genome_fraction, double fdr_cutoff, output):

    cdef:
        int i
        int j
        int number_bins
        island _island
        IslandVector _islands
        IslandVector all_islands = IslandVector()
        uint16_t count
        uint32_t _bin
        uint32_t[::1] bins = np.ones(1, dtype=np.uint32)
        uint16_t[::1] counts = np.ones(1, dtype=np.uint16)
        double scaling_factor = float(chip_library_size) / float(control_library_size)
        double zero_scaler = (control_library_size / effective_genome_fraction)
        double fdr
        double average
        int num_islands
        int counter

    sf = poisson.sf
    chromosomes = natsorted(set(islands.keys()))

    num_islands_per_chrom = {}
    for c in natsorted(islands):
        _islands = islands[c]
        num_islands_per_chrom[c] = len(_islands)

    for chromosome in chromosomes:

        j = 0
        _islands = islands[chromosome]
        # print("chr", chromosome, len(_islands), _islands)

        if chromosome not in b_bins_counts:
            continue

        bins, counts = b_bins_counts[chromosome]
        number_bins = len(bins)

        for i in range(len(_islands)):

            # not overlapping
            while j < number_bins and (bins[j] < _islands.wrapped_vector[i].start):
                # print("bins[j]", bins[j], "_islands.wrapped_vector[i].start", _islands.wrapped_vector[i].start)
                j += 1

            # overlapping
            while j < number_bins and (bins[j] < _islands.wrapped_vector[i].end and bins[j] >= _islands.wrapped_vector[i].start):
                _islands.wrapped_vector[i].input_count += counts[j]
                j += 1

            # if _islands.wrapped_vector[i].start == 117299800:
            #     _print = True

            # else:
            #     _print = False

            ##### COMPUTE P ########
            if _islands.wrapped_vector[i].input_count > 0:
                average = _islands.wrapped_vector[i].input_count * scaling_factor
                # if _print:
                #     print("input_count", _islands.wrapped_vector[i].input_count)
                #     print("scaling factor", scaling_factor)
                #     print("average", average)
            else:
                average = (_islands.wrapped_vector[i].end - _islands.wrapped_vector[i].start + 1) * zero_scaler
                average = min(0.25, average) * scaling_factor

            _islands.wrapped_vector[i].fold_change = log2(_islands.wrapped_vector[i].chip_count / average)

            if _islands.wrapped_vector[i].chip_count > average:
                _islands.wrapped_vector[i].p_value = sf(_islands.wrapped_vector[i].chip_count, average)
                # if _print:
                #     print("p_val:", sf(_islands.wrapped_vector[i].chip_count, average))
                #     print("chip_count:", _islands.wrapped_vector[i].chip_count)
            else:
                _islands.wrapped_vector[i].p_value = 1


            all_islands.push_back(_islands.wrapped_vector[i])

        del b_bins_counts[chromosome]

    p_values = np.zeros(len(all_islands))

    for i in range(len(all_islands)):
        p_values[i] = all_islands.wrapped_vector[i].p_value

    ranks = rankdata(p_values)

    counter = 0
    j = 0
    num_islands = len(all_islands)

    if not output:
        handle = sys.stdout
    else:
        handle = open(output, "w+")

    print("\t".join(["#Chromosome", "Start", "End", "PValue", "Score", "Strand", "ChIPCount", "InputCount", "FDR", "log2FoldChange"]), file=handle)
    for i in range(num_islands):
            # _island = all_islands.wrapped_vector[i]
            fdr = all_islands.wrapped_vector[i].p_value * num_islands / ranks[i]
            if fdr > 1:
                fdr = 1

            if fdr <= fdr_cutoff:
                chromosome = all_islands.wrapped_vector[i].chromosome.decode()
                print("\t".join(str(e) for e in [chromosome, all_islands.wrapped_vector[i].start, all_islands.wrapped_vector[i].end, all_islands.wrapped_vector[i].p_value,
                                                 min(1000, all_islands.wrapped_vector[i].fold_change * 100), ".", all_islands.wrapped_vector[i].chip_count, all_islands.wrapped_vector[i].input_count, fdr, all_islands.wrapped_vector[i].fold_change]), file=handle)

    if output:
        handle.close()


def write_islands(islands, float average_window_readcount, float fdr_cutoff, int e_value, output):

    cdef:
        IslandVector _islands
        IslandVector all_islands = IslandVector()
        int i
        int num_islands = 0
        int counter


    chromosomes = natsorted(set(islands.keys()))

    num_islands_per_chrom = {}
    for c in natsorted(islands):
        _islands = islands[c]
        num_islands_per_chrom[c] = len(_islands)


    _poisson = poisson(average_window_readcount)

    if not output:
        handle = sys.stdout
    else:
        handle = open(output, "w+")

    print("\t".join(["#Chromosome", "Start", "End", "ChIPCount", "Score", "Strand"]), file=handle)
    for chromosome in chromosomes:
        _islands = islands[chromosome]
        num_islands += len(_islands)

        for i in range(len(_islands)):
            print("\t".join(str(e) for e in [chromosome, _islands.wrapped_vector[i].start, _islands.wrapped_vector[i].end, _islands.wrapped_vector[i].chip_count, float(_islands.wrapped_vector[i].score), "."]),
                  file=handle)

    if output:
        handle.close()

    logging.info("Empirical estimate of FDR is: {} ({}/{})".format(min(1, float(e_value)/num_islands), e_value, num_islands))
    # ranks = rankdata(np.array([_islands.wrapped_vector[i].p_value for _islands.wrapped_vector[i] in all_islands.wrapped_vector[i]s], dtype=np.int))

    # p_values = np.zeros(len(all_islands))

    # for i in range(len(all_islands)):
    #     p_values[i] = all_islands.wrapped_vector[i].p_value

    # ranks = rankdata(p_values)

    # counter = 0
    # num_islands = len(all_islands)

    # for chromosome, chromosome_size in natsorted(num_islands_per_chrom.items()):

    #     i = 0
    #     for i in range(int(chromosome_size)):

    #         #### compute fdr #####
    #         # _islands.wrapped_vector[i] = all_islands.wrapped_vector[i + counter]

    #         fdr = all_islands.wrapped_vector[i + counter].p_value * num_islands / ranks[i + counter]
    #         # print("fdr cutoff", fdr_cutoff)
    #         if fdr > 1:
    #             fdr = 1
    #         if fdr <= fdr_cutoff:
    #             print("\t".join(str(e) for e in [chromosome, all_islands.wrapped_vector[i].start, all_islands.wrapped_vector[i].end, all_islands.wrapped_vector[i].p_value, all_islands.wrapped_vector[i].chip_count, ".", fdr]))


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.initializedcheck(False)
@cython.cdivision(False)
def differential_count_reads_on_islands(dfs, c_bins_counts):

    # want count per island
    # must iterate over islands, per chromosome
    # add counts per islands

    # need vector to hold counts

    cdef:
        int i, number_bins, number_islands
        int j = 0
        uint32_t[::1] bins = np.ones(1, dtype=np.uint32)
        uint16_t[::1] counts = np.ones(1, dtype=np.uint16)
        uint32_t[::1] starts
        uint32_t[::1] ends

    new_islands = {}
    chromosomes = natsorted(dfs)

    for chromosome in chromosomes:
        df = dfs[chromosome]

        j = 0
        i = 0

        number_islands = len(df)
        outcounts = np.zeros(number_islands, dtype=np.uint32)
        starts, ends = df.Start.astype(np.uint32).values, df.End.astype(np.uint32).values

        if chromosome not in c_bins_counts:
            new_islands[chromosome] = outcounts
            continue

        bins, counts = c_bins_counts[chromosome]
        number_bins = len(bins)

        for i in range(number_islands):
            # _island = &_islands.wrapped_vector[i]

            # not overlapping
            while j < number_bins and (bins[j] < starts[i]):
                # print("not overlapping:")
                # print("bins[j]", bins[j], "counts[j]", counts[j], "_island.start", _island)
                j += 1

            # overlapping
            while j < number_bins and (bins[j] < ends[i] and bins[j] >= starts[i]):
                # print("overlapping:")
                # print("bins[j]", bins[j], "counts[j]", counts[j], "_island.start", _island)
                outcounts[i] += counts[j]
                j += 1


        new_islands[chromosome] = outcounts

        # if not args.get("df"):
        #     del c_bins_counts[chromosome]

    return new_islands
