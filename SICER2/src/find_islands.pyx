from scipy.stats import poisson, rankdata
import cython
from libcpp cimport bool

from libcpp.algorithm cimport sort as stdsort
from libc.math cimport log2
import numpy as np

from libc.stdint cimport uint32_t, uint16_t
from SICER2.src.statistics import compute_window_score
from libcpp.vector cimport vector
# from libc.math cimport min as cmin

from natsort import natsorted


cdef struct island:
    uint32_t start
    uint32_t end
    uint16_t chip_count
    uint16_t input_count
    float score
    float p_value
    float fold_change


# cdef bool compare_p_value(const island &lhs, const island &rhs):

#     return lhs.p_value < rhs.p_value



cdef class IslandVector:

    cdef vector[island] wrapped_vector

    cdef push_back(self, island i):
        self.wrapped_vector.push_back(i)

    def __str__(self):
        return "[" + ", ".join([str(i) for i in self.wrapped_vector]) + "]"

    def __repr__(self):
        return str(self)

    def __len__(self):
        return self.wrapped_vector.size()

    def __iter__(self):
        # slow, only implemented to ease testing
        return (v for v in self.wrapped_vector)

    # def p_value_sort(self):
    #     stdsort(self.wrapped_vector.begin(), self.wrapped_vector.end(), &compare_p_value)



@cython.boundscheck(False)
@cython.wraparound(False)
@cython.initializedcheck(False)
def find_islands(bins_counts, int gaps_allowed, int bin_size, float score_threshold, uint32_t island_enriched_threshold, float average_window_readcount):

    _poisson = poisson(average_window_readcount)

    cdef:
        pass
        float slightly_less = score_threshold - 0.0000000001
        float score
        int i
        int j
        uint16_t count
        uint32_t _bin
        uint32_t[::1] bins = np.ones(1, dtype=np.uint32)
        uint16_t[::1] counts = np.ones(1, dtype=np.uint16)
        uint32_t dist
        island _island
        island current_island
        IslandVector v
        uint32_t distance_allowed = (gaps_allowed * bin_size) + 2
        # uint32_t island_enriched_threshold = int(_island_enriched_threshold)


    # bins_arr = np.ones(1, dtype=np.uint32)
    # counts_arr = np.ones(1, dtype=np.uint32)
    chromosomes = set(bins_counts.keys())

    island_dict = dict()

    for chromosome in chromosomes:
        v = IslandVector()
        i = 0
        j = 4294967295

        # TODO: if chromo not in chromsizes, remove
        bins, counts = bins_counts[chromosome]

        count = counts[0]
        _bin = bins[0]
        score = compute_window_score(count, _poisson)

        current_island = [_bin, _bin + bin_size - 1, count, 0, score, 0, 0]
        for i in range(1, len(bins)):

            count = counts[i]

            _bin = bins[i]

            dist = _bin - current_island.end

            if dist <= distance_allowed:
                current_island.end = bins[i] + bin_size - 1
                current_island.chip_count += counts[i]
                current_island.score += compute_window_score(counts[i], _poisson)
            else:
                if current_island.score > slightly_less and current_island.chip_count > island_enriched_threshold:
                    v.push_back(current_island)

                score = compute_window_score(counts[i], _poisson)
                current_island = [_bin, bins[i] + bin_size - 1, count, 0, score, 0, 0]

        # takes care of last or potential single bin
        if current_island.score > slightly_less and current_island.chip_count > island_enriched_threshold:
            v.push_back(current_island)

        island_dict[chromosome] = v

        del bins_counts[chromosome]

    return island_dict

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.initializedcheck(False)
@cython.cdivision(True)
def compute_fdr(islands, b_bins_counts, int chip_library_size, int control_library_size, float effective_genome_fraction):

    cdef:
        int i
        int j
        int number_bins
        # island *_island
        island _island
        IslandVector _islands
        IslandVector all_islands = IslandVector()
        uint16_t count
        uint32_t _bin
        uint32_t[::1] bins = np.ones(1, dtype=np.uint32)
        uint16_t[::1] counts = np.ones(1, dtype=np.uint16)
        float scaling_factor = float(chip_library_size) / float(control_library_size)
        float zero_scaler = (control_library_size / effective_genome_fraction)
        float fdr
        float average
        int num_islands
        int counter

    # print("chip_library_size", chip_library_size)
    # print("control_library_size", control_library_size)
    # print("float(chip_library_size) / float(control_library_size)", float(chip_library_size) / float(control_library_size))
    sf = poisson.sf
    chromosomes = natsorted(set(islands.keys()))
    num_islands_per_chrom = [len(v) for _, v in natsorted(islands.items())]
    # print("scaling_factor", scaling_factor)
    # print("zero_scaler", zero_scaler)
    for chromosome in chromosomes:
        j = 0
        _islands = islands[chromosome]
        # print("chr", chromosome, len(_islands), _islands)

        if chromosome not in b_bins_counts:
            continue

        bins, counts = b_bins_counts[chromosome]
        number_bins = len(bins)

        for i in range(len(_islands)):
            # _island = &_islands.wrapped_vector[i]
            _island = _islands.wrapped_vector[i]

            # not overlapping
            while (bins[j] < _island.start and j < number_bins):
                # print("bins[j]", bins[j], "_island.start", _island.start)
                j += 1

            # overlapping
            while (bins[j] < _island.end and bins[j] >= _island.start and j < number_bins):
                _island.input_count += counts[j]
                j += 1

            ##### COMPUTE P ########
            if _island.input_count > 0:
                average = _island.input_count * scaling_factor
            else:
                average = (_island.end - _island.start + 1) * zero_scaler
                average = min(0.25, average) * scaling_factor

            _island.fold_change = log2(_island.chip_count / average)

            if _island.chip_count > average:
                _island.p_value = sf(_island.chip_count, average)
            else:
                _island.p_value = 1

            all_islands.push_back(_island)

        del b_bins_counts[chromosome]


    ranks = rankdata(np.array([_island.p_value for _island in all_islands], dtype=np.int))

    counter = 0
    num_islands = len(all_islands)

    print("\t".join(["Chromosome", "Start", "End", "PValue", "Score", "Strand", "ChIPCount", "InputCount", "FDR", "log2(FoldChange)"]))
    for chromosome, chromosome_size in zip(chromosomes, num_islands_per_chrom):

        i = 0
        for i in range(chromosome_size):

            #### compute fdr #####
            _island = all_islands.wrapped_vector[i + counter]
            fdr = _island.p_value * num_islands / ranks[i + counter]
            if fdr > 1:
                fdr = 1

            print("\t".join(str(e) for e in [chromosome, _island.start, _island.end, _island.p_value,
                                             min(1000, _island.fold_change * 100), ".", _island.chip_count, _island.input_count, fdr, _island.fold_change]))

        counter += chromosome_size


def write_islands(islands, average_window_readcount):

    cdef:
        island _island
        IslandVector _islands
        IslandVector all_islands = IslandVector()
        int i
        int num_islands
        int counter


    chromosomes = natsorted(set(islands.keys()))
    num_islands_per_chrom = [len(v) for _, v in natsorted(islands.items())]
    _poisson = poisson(average_window_readcount)
    print("\t".join(["Chromosome", "Start", "End", "PValue", "ChIPCount", "Strand", "FDR"]))
    for chromosome in chromosomes:
        _islands = islands[chromosome]
        for _island in _islands:
            _island.p_value = _poisson.pmf(_island.chip_count)
            all_islands.wrapped_vector.push_back(_island)

    ranks = rankdata(np.array([_island.p_value for _island in all_islands], dtype=np.int))
    counter = 0
    num_islands = len(all_islands)
    for chromosome, chromosome_size in zip(chromosomes, num_islands_per_chrom):

        i = 0
        for i in range(chromosome_size):

            #### compute fdr #####
            _island = all_islands.wrapped_vector[i + counter]

            fdr = _island.p_value * num_islands / ranks[i + counter]
            if fdr > 1:
                fdr = 1
            print("\t".join(str(e) for e in [chromosome, _island.start, _island.end, _island.p_value, _island.chip_count, ".", fdr]))
