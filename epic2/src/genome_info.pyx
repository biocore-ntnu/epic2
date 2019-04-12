import logging
import gzip
from os.path import basename
from natsort import natsorted
from collections import OrderedDict

from libc.stdint cimport uint32_t, uint16_t
import numpy as np
import sys
from re import search, IGNORECASE
import pkg_resources



def sniff(f):

    try:
        import pysam
        pysam.AlignmentFile(f, "rb")
        return "bam"
    except:
        pass

    try:
        import gzip
        l = gzip.open(f).readline().decode().split()
        if len(l) == 6:
            return "bed.gz"
        elif len(l) == 10:
            return "bedpe.gz"
        else:
            raise ValueError("Wrong number columns in {}: {}.\nLine looks like: {}".format(f, len(l), l))
    except:
        pass

    try:
        l = open(f).readline().split()
        if len(l) == 6:
            return "bed"
        elif len(l) == 10:
            return "bedpe"
        else:
            raise ValueError("Wrong number columns in {}: {}.\nLine looks like: {}".format(f, len(l), l))

    except:
        raise Exception("Could not find format of file {}.".format(f))

def get_closest_readlength(estimated_readlength):
    # type: (int) -> int
    """Find the predefined readlength closest to the estimated readlength.

    In the case of a tie, choose the shortest readlength."""

    readlengths = [36, 50, 75, 100]
    differences = [abs(r - estimated_readlength) for r in readlengths]
    min_difference = min(differences)
    index_of_min_difference = [i
                               for i, d in enumerate(differences)
                               if d == min_difference][0]

    return readlengths[index_of_min_difference]


def find_readlength(args):
    # type: (Namespace) -> int
    """Estimate length of reads based on 1000 first."""

    # from subprocess import check_output
    # import pandas as pd
    # from io import BytesIO
    _file = args["treatment"][0]
    file_format = sniff(_file)

    arr = np.zeros(100, dtype=np.uint32)

    cdef:
        int i = 0
        uint32_t start
        uint32_t end

    if file_format == "bed":
        for line in open(_file):
            _start, _end = line.split()[1:3]
            arr[i] = int(_end) - int(_start)
            i += 1
            if i == 100:
                break


    elif file_format == "bed.gz":
        for line in gzip.open(_file):
            line = line.decode()
            _start, _end = line.split()[1:3]
            arr[i] = int(_end) - int(_start)
            i += 1
            if i == 100:
                break


    elif file_format == "bedpe":
        for line in open(_file):
            ls = line.split()
            _start = ls[1]
            _end = ls[5]
            arr[i] = int(_end) - int(_start)
            i += 1
            if i == 100:
                break

    elif file_format == "bedpe.gz":
        for line in gzip.open(_file):
            line = line.decode()
            ls = line.split()
            _start = ls[1]
            _end = ls[5]
            arr[i] = int(_end) - int(_start)
            i += 1
            if i == 100:
                break


    elif file_format == "bam":

        import pysam

        if _file.endswith(".bam"):
            samfile = pysam.AlignmentFile(_file, "rb")
        elif _file.endswith(".sam"):
            samfile = pysam.AlignmentFile(_file, "r")

        for a in samfile:
            if a.alen == None: continue

            arr[i] = a.alen
            i += 1
            if i == 100:
                break

    else:
        raise IOError("Cannot recognize file extension of: " + _file + ". Must be bed, bedpe, bam, sam, bed.gz or bedpe.gz")

        # print(arr)
    median = np.median(arr)
    logging.info("Found a median readlength of {}\n".format(median))


    return get_closest_readlength(median)



def get_genome_size_file(genome):

    genome_names = pkg_resources.resource_listdir("epic2", "chromsizes")
    name_dict = {n.lower().replace(".chromsizes", ""): n for n in genome_names}

    # # No try/except here, because get_egs would already have failed if genome
    # # did not exist
    genome_exact = name_dict[genome.lower()]

    return pkg_resources.resource_filename(
        "epic2", "chromsizes/{}".format(genome_exact))


def create_genome_size_dict(genome):
    # type: (str) -> Dict[str,int]
    """Creates genome size dict from string containing data."""

    size_file = get_genome_size_file(genome)
    size_lines = open(size_file).readlines()

    size_dict = {}
    for line in size_lines:
        genome, length = line.split()
        size_dict[genome] = int(length)

    return size_dict


def create_genome_size_dict_custom_genome(chromsizes):

    chromosome_lengths = [l.split() for l in open(chromsizes).readlines()]

    od = OrderedDict()          # type: OrderedDict[str, int]

    for c, l in natsorted(chromosome_lengths):
        od[c] = int(l)

    return od


def get_effective_genome_fraction(genome, read_length):

    genome_names = pkg_resources.resource_listdir("epic2",
                                                  "effective_sizes")

    name_dict = {n.lower(): n for n in genome_names}

    # # No try/except here, because get_egs would already have failed if genome
    # # did not exist
    genome_exact = name_dict[genome.lower() + "_" + str(read_length) + ".txt"]


    try:
        egf = pkg_resources.resource_string( # type: ignore
            "epic2", "effective_sizes/{}".format(
                genome_exact)).split()[-1].decode()
    except KeyError:

        genome_list = "\n".join([basename(g) for g in genome_names])

        raise Exception("Genome " + genome +
            " not found.\n These are the available genomes: " + genome_list +
            "\nIf yours is not there, please request it at github.com/endrebak/epic2 .")

    return float(egf)



def get_chroms_from_bam(bamfile, remove_pattern):
    """Read bam file and return dict of {chrom_name: length}
    """
    import pysam
    inbam = pysam.AlignmentFile(bamfile)
    lengths = inbam.lengths
    references = inbam.references
    chroms = {}

    removed = []

    if not remove_pattern:
        for x, l in zip(references, lengths):
            chroms[x] = l
    else:
        import re
        pattern = re.compile(remove_pattern)
        for x, l in zip(references, lengths):
            if not pattern.search(x):
                chroms[x] = l
            else:
                removed.append(x)

    logging.info("Chromosomes included in analysis: {}".format(", ".join(natsorted(chroms))))
    if removed:
        logging.info("Chromosomes removed from analysis: {} (and {} more.)\n".format(", ".join(natsorted(removed[:5])), max(0, len(removed) - 5)))
    else:
        logging.info("No chromosomes matching removal pattern: {}\n".format(remove_pattern))


    return chroms

def egl_and_chromsizes(args):

    read_length = find_readlength(args)

    if args["autodetect_chroms"]:
        chromsizes = get_chroms_from_bam([f for f in args["treatment"] if f.endswith(".bam")][0], args["discard_chromosomes_pattern"])
        egf = get_effective_genome_fraction(args["genome"], read_length)
        genome_length = sum(chromsizes.values())
        egl = egf * genome_length

    else:

        have_chromsizes = args["chromsizes"] != None
        have_effective_genome_fraction = args["effective_genome_fraction"] != None

        if have_chromsizes and have_effective_genome_fraction:
            # assert have_chromsizes
            # assert have_effective_genome_fraction

            chromsizes = create_genome_size_dict_custom_genome(args["chromsizes"])
            egf = args["effective_genome_fraction"]
            egl = egf * sum(chromsizes.values())
            genome_length = sum(chromsizes.values())
        elif have_chromsizes:
            chromsizes = create_genome_size_dict_custom_genome(args["chromsizes"])
            egf = get_effective_genome_fraction(args["genome"], read_length)
            egl = egf * sum(chromsizes.values())
            genome_length = sum(chromsizes.values())
        elif have_effective_genome_fraction:
            chromsizes = create_genome_size_dict(args["genome"])
            egf = args["effective_genome_fraction"]
            genome_length = sum(chromsizes.values())
            egl = egf * genome_length

        else:
            chromsizes = create_genome_size_dict(args["genome"])
            egf = get_effective_genome_fraction(args["genome"], read_length)
            genome_length = sum(chromsizes.values())
            egl = egf * genome_length

    logging.info("Using genome {}.\n".format(args["genome"]))
    logging.info("Using an effective genome length of ~{} * 1e6\n".format(int(egl/1e6)))

    return egl, chromsizes
