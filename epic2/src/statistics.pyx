import sys
from libcpp.vector cimport vector
import numpy as np
from scipy.stats import poisson
from numpy import log
from itertools import count

from libc.math cimport log as clog
from libc.math cimport round as cround

try:
    from functools import lru_cache
except ImportError:
    from functools32 import lru_cache


@lru_cache()
def compute_window_score(int i, _poisson):

    cdef float window_score
    cdef float p_value

    if i < _poisson.mean():
        # sys.stderr.write(" ".join([str(e) for e in [ "i of", i, "gives poisson score of", 0, "\n" ]]))
        return 0

    p_value = _poisson.pmf(i)

    if p_value > 0:
        # sys.stderr.write(" ".join([str(e) for e in ["i of", i, "gives poisson score of", -log(p_value), "and p_value of", p_value, "\n"] ]))
        window_score = -log(p_value)
    else:
        # log of zero not defined
        # sys.stderr.write(" ".join([str(e) for e in ["i of", i, "gives poisson score of", -log(p_value), "and p_value of", p_value, "\n"] ]))
        window_score = 1000

    # sys.stderr.write("Returning window score of " + str(window_score) + "\n")
    return window_score




def generate_cumulative_distribution(vector[float] island_expectations):

    cdef:
        int i
        int l = island_expectations.size()
        double [::1] cumulative
        float partial_sum = 0


    cumulative_arr = np.zeros(l, dtype=float)
    cumulative = cumulative_arr

    for i in range(1, l + 1):

        compliment = l - i
        partial_sum += island_expectations[compliment]
        cumulative[compliment] = partial_sum

    partial_sum += island_expectations[l]
    cumulative[1] = partial_sum

    return cumulative


def update_island_expectations(vector[float] island_expectations, int scaled_score, int bin_size, float average_window_readcount, int island_enriched_threshold, float gap_contribution):

    cdef:
        float WINDOW_P_VALUE = 0.20
        float BIN_SIZE = 0.001
        int E_VALUE = 1000
        float E_VALUE_THRESHOLD = E_VALUE * .0000001
        int i = island_enriched_threshold
        float temp
        int index

    _poisson = poisson(average_window_readcount)
    current_max_scaled_score = island_expectations.size() - 1
    if scaled_score > current_max_scaled_score:
        #index is the scaled_score
        for index in range(current_max_scaled_score + 1, scaled_score + 1):
            temp=0.0


            current_island = int(cround(index - compute_window_score(i, _poisson) / BIN_SIZE))

            while (current_island >= 0):
                island_expectation = island_expectations[current_island]

                temp += _poisson.pmf(i) * island_expectation
                i += 1
                if i == 500:
                    break

                current_island = int(cround(index - compute_window_score(i, _poisson) / BIN_SIZE))

                temp *= gap_contribution
            island_expectations.push_back(temp)


    return island_expectations


def compute_score_threshold(average_window_readcount,
                            island_enriched_threshold,
                            gap_contribution, boundary_contribution,
                            genome_length_in_bins, bin_size):
    # type: (float, int, float, float, float) -> float
    """
    What does island_expectations do?
    """

    cdef:
        int current_scaled_score
        int l
        int interval
        float partial_cumu
        float e

        float WINDOW_P_VALUE = 0.20
        float BIN_SIZE = 0.001
        int E_VALUE = 1000
        float E_VALUE_THRESHOLD = E_VALUE * .0000001


    required_p_value = poisson.pmf(island_enriched_threshold,
                                   average_window_readcount)
    # print("required_p_value", required_p_value)
    prob = boundary_contribution * required_p_value
    # print("prob", prob)

    score = -log(required_p_value)
    # print("score", score)

    current_scaled_score = int(round(score / BIN_SIZE))
    # print("current_scaled_score", current_scaled_score)

    cdef vector[float] island_expectations = vector[float](current_scaled_score + 1)

    island_expectations[0] = boundary_contribution * genome_length_in_bins / gap_contribution

    island_expectations[current_scaled_score] = prob * genome_length_in_bins
    # print("len(island_expectations)", island_expectations.size())
    # print(island_expectations)

    current_max_scaled_score = current_scaled_score

    interval = 1000 # 1 / BIN_SIZE)
    partial_cumu = 0.0

    while (partial_cumu > E_VALUE_THRESHOLD or partial_cumu < 1e-100):

        current_scaled_score += interval
        island_expectations = update_island_expectations(island_expectations, current_scaled_score, bin_size, average_window_readcount, island_enriched_threshold, gap_contribution)
        current_expectation = island_expectations[current_scaled_score]

        l = island_expectations.size()
        if l > interval:
            partial_cumu = sum(island_expectations[i] for i in range(l - interval + 1, l))
        else:
            partial_cumu = sum(island_expectations)

    cumulative = generate_cumulative_distribution(island_expectations)
    # print(list(cumulative))

    score_threshold = 0
    for (i, e) in enumerate(cumulative):
        if e < E_VALUE:
            score_threshold = (i - 1) * BIN_SIZE
            break

    return score_threshold



def compute_enriched_threshold(average_window_readcount):
    # type: (float) -> int
    """
    Computes the minimum number of tags required in window for an island to be enriched.
    """

    cdef:

        float WINDOW_P_VALUE = 0.20
        float BIN_SIZE = 0.001
        int E_VALUE = 1000
        float E_VALUE_THRESHOLD = E_VALUE * .0000001

    current_threshold, survival_function = 0, 1
    while True:
        survival_function -= poisson.pmf(current_threshold,
                                         average_window_readcount)
        current_threshold += 1
        if survival_function <= WINDOW_P_VALUE:
            break

    island_enriched_threshold = current_threshold

    return island_enriched_threshold


def compute_gap_factor(island_enriched_threshold,
                       gap_intervals_allowed, poisson_distribution_parameter):
    # type: (int, int, float) -> float

    max_gap_score = 1.0
    gap_factor = single_gap_factor(island_enriched_threshold,
                                   poisson_distribution_parameter)
    max_gap_score += sum([pow(gap_factor, i)
                          for i in range(1, gap_intervals_allowed + 1)])
    return max_gap_score


def single_gap_factor(island_enriched_threshold,
                      poisson_distribution_parameter):
    # type: (int, float) -> float

    poisson_scores = [poisson.pmf(i, poisson_distribution_parameter)
                      for i in range(island_enriched_threshold)]
    return sum(poisson_scores)


def compute_boundary(island_enriched_threshold, gap_intervals_allowed,
                     average):
    # type: (int, int, float) -> float

    single_gap = single_gap_factor(island_enriched_threshold, average)
    single_boundary_score = pow(single_gap, gap_intervals_allowed + 1)
    start_and_end_score = single_boundary_score * single_boundary_score

    return start_and_end_score


def compute_background_probabilities(total_chip_count, bin_size, effective_genome_fraction, gaps_allowed):
    # type: (int, Namespace) -> Tuple[float, int, float]

    print("total_chip_count", total_chip_count)
    print("bin_size", bin_size)
    print("effective_genome_fraction", effective_genome_fraction)

    average_window_readcount = total_chip_count * (bin_size / float(effective_genome_fraction))
    print("average_window_readcount", average_window_readcount)

    island_enriched_threshold = compute_enriched_threshold(average_window_readcount)
    print("island_enriched_threshold", island_enriched_threshold)

    gap_contribution = compute_gap_factor(island_enriched_threshold, gaps_allowed, average_window_readcount)
    print("gap_contribution", gap_contribution)

    boundary_contribution = compute_boundary(island_enriched_threshold, gaps_allowed, average_window_readcount)
    print("boundary_contribution", boundary_contribution)

    genome_length_in_bins = effective_genome_fraction / bin_size
    print("genome_length_in_bins", genome_length_in_bins)

    score_threshold = compute_score_threshold(average_window_readcount, island_enriched_threshold, gap_contribution, boundary_contribution, genome_length_in_bins, bin_size)
    print("score_threshold", score_threshold)

    return score_threshold, island_enriched_threshold, average_window_readcount

# ('total_chip_count', 9999)
# ('bin_size', 200)
# ('effective_genome_fraction', 2476555186.4)
# ('average_window_readcount', 0.0008074926054472355)
# ('island_enriched_threshold', 1)
# ('gap_contribution', 3.9951596055200813)
# ('boundary_contribution', 0.9935608797169523)
# ('genome_length_in_bins', 12382775.932)
# ('required_p_value', 0.0008068408243290163)
# ('prob', 0.0008016454792118883)
# ('score', 7.122384152855723)
# current_scaled_score
# ('len(island_expectations)', 7124)
# (7122, 9926.5966796875)
# (7123, 3079487.0)
