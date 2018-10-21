
import pytest

from SICER2.src.reads_to_bins import add_reads_to_dict, files_to_bin_counts


args = {"bin_size": 200, "fragment_size": 150, "drop_duplicates": True}


def test_files_to_bins():

    d = files_to_bin_counts(["tests/test.bed"], args)
    print(d)

    print(d["chr7"][0])
    assert list(d["chr7"][0]) == [20246600, 39036600]
    assert list(d["chr7"][1]) == [2, 2]



def test_files_to_bins2():

    d = files_to_bin_counts(["tests/test2.bed"], args)
    print(d)

    print(d["chr7"][0])
    assert list(d["chr7"][0]) == [20246600, 39036600, 39046600]
    assert list(d["chr7"][1]) == [1, 1, 1]

def test_files_to_bins3():

    d = files_to_bin_counts(["tests/test3.bed"], args)
    print(d)

    print(d["chr7"][0])
    assert list(d["chr7"][0]) == [20246600, 39036600, 39046600]
    assert list(d["chr7"][1]) == [2, 1, 1]


def test_files_to_bins123():

    d = files_to_bin_counts(["tests/test.bed", "tests/test2.bed", "tests/test3.bed"], args)

    print(d)

    print(d["chr7"][0])
    assert list(d["chr7"][0]) == [20246600, 39036600, 39046600]
    assert list(d["chr7"][1]) == [5, 4, 2]
# def test_reads_to_dict():

#     print(open("tests/test.bed").readlines())

#     d = add_reads_to_dict("tests/test.bed")
#     print(d)

#     assert list(d["chr7", "+"]) == [20246668, 20246668]
