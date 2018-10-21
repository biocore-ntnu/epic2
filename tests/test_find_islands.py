import pytest
import numpy as np

from SICER2.src.find_islands import find_islands
from SICER2.src.find_islands import compute_fdr



def test_find_islands():

    bins_counts = {"chr7": (np.array([0, 600, 10000, 10200, 11000], dtype=np.uint32), np.array([1, 2, 3, 4, 5], dtype=np.uint16))}

    islands = find_islands(bins_counts, 3, 200, 3.14, 1, 0.0005)

    print("islands", islands)

    expected = {'chr7': [{'start': 0, 'end': 799, 'chip_count': 3, 'input_count': 0, 'score': 23.49685287475586, 'p_value': 0.0, 'fold_change': 0.0}, {'start': 10000, 'end': 11199, 'chip_count': 12, 'input_count': 0, 'score': 100.96963500976562, 'p_value': 0.0, 'fold_change': 0.0}]}["chr7"]

    for i, island in enumerate(expected):
        result = list(islands["chr7"])[i]
        assert expected[i] == result


# def test_fdr():

#     result = compute_fdr({"chr7": },
#                          {"chr7": ([0, 200, 1200], [1, 1, 1])}, 3, 1, 1000)
#     print(result)
#     assert 0
