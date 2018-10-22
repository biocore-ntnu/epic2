from os.path import basename
from natsort import natsorted
from collections import OrderedDict

from libc.stdint cimport uint32_t, uint16_t
import numpy as np
import sys
from re import search, IGNORECASE
import pkg_resources


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


    arr = np.zeros(100, dtype=np.uint32)

    cdef:
        int i = 0
        uint32_t start
        uint32_t end

    if _file.endswith(".bed"):
        for line in open(_file):
            _start, _end = line.split()[1:3]
            arr[i] = int(_end) - int(_start)
            i += 1
            if i == 100:
                break

    elif _file.endswith(".bam") or _file.endswith(".sam"):

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

        print(arr)


    return get_closest_readlength(np.median(arr))


    # filereader = "cat "
    # if bed_file.endswith(".gz"):
    #     filereader = "gzcat "
    # elif bed_file.endswith(".bz2"):
    #     filereader = "bzgrep "
    # elif bed_file.endswith(".bam"):
    #     filereader = "samtools view "

    # command = filereader + "{} | head -10000".format(bed_file)
    # output = check_output(command, shell=True)

    # df = pd.read_table(
    #     BytesIO(output),
    #     header=None,
    #     usecols=[1, 2],
    #     sep="\t",
    #     names=["Start", "End"])

    # readlengths = df.End - df.Start
    # mean_readlength = readlengths.mean()
    # median_readlength = readlengths.median()
    # max_readlength = readlengths.max()
    # min_readlength = readlengths.min()

    # logging.info((
    #     "Used first 10000 reads of {} to estimate a median read length of {}\n"
    #     "Mean readlength: {}, max readlength: {}, min readlength: {}.").format(
    #         bed_file, median_readlength, mean_readlength, max_readlength,
    #         min_readlength))




def get_genome_size_file(genome):

    # genome_names = pkg_resources.resource_listdir("epic", "scripts/chromsizes")
    # name_dict = {n.lower().replace(".chromsizes", ""): n for n in genome_names}

    # # No try/except here, because get_egs would already have failed if genome
    # # did not exist
    # genome_exact = name_dict[genome.lower()]

    return pkg_resources.resource_filename(
        "SICER2", "chromsizes/{}.chromsizes".format(genome.lower()))


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

    genome_names = pkg_resources.resource_listdir("SICER2",
                                                  "effective_sizes")

    try:
        egf = pkg_resources.resource_string( # type: ignore
            "SICER2", "effective_sizes/{}_{}.txt".format(
                genome.lower(), read_length)).split()[-1].decode()
    except KeyError:

        genome_list = "\n".join([basename(g) for g in genome_names])

        raise Exception("Genome " + genome +
            " not found.\n These are the available genomes: " + genome_list +
            "\nIf yours is not there, please request it at github.com/endrebak/SICER2 .")

    return float(egf)


def egl_and_chromsizes(args):

    have_chromsizes = args["chromsizes"] != None
    have_effective_genome_fraction = args["effective_genome_fraction"] != None

    if have_chromsizes or have_effective_genome_fraction:
        assert have_chromsizes
        assert have_effective_genome_fraction

        chromsizes = create_genome_size_dict_custom_genome(args["chromsizes"])
        egf = args["effective_genome_fraction"]
        egl = egf * sum(chromsizes.values())
        genome_length = sum(chromsizes.values())
    else:
        read_length = find_readlength(args)
        chromsizes = create_genome_size_dict(args["genome"])
        egf = get_effective_genome_fraction(args["genome"], read_length)
        genome_length = sum(chromsizes.values())
        egl = egf * genome_length

    sys.stderr.write("Using an effective genome length of ~{} * 1e6 \n\n".format(int(egl/1e6)))

    return egl, chromsizes
