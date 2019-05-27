import pandas as pd
import pyranges as pr
from epic2.src.genome_info import sniff


def file_to_unbinned_ranges(f, args)():

    chromosome_data = {}

    file_format = sniff(f)

    names = "Chromosome Start End Strand".split()
    if file_format == "bed" or "bed.gz":
        df = pd.read_csv(f, sep="\t", header=None, usecols=[0, 1, 2, 5], names=names)
        gr = pr.PyRanges(df)
    elif file_format == "bedpe":
        df = pd.read_csv(f, sep="\t", header=None, usecols=[0, 1, 5, 9], names=names)
        gr = pr.PyRanges(df)
    else:
        gr = pr.read_bam(f)

    return gr



def file_to_binned_ranges(f, args):

    chromosome_data = {}
    bin_counts, _ = files_to_bin_counts([f], args, datatype)

    for c, v in bin_counts.items():
        starts, scores = bin_counts
        ends = starts + args["bin_size"]
        df = pd.DataFrame({"Start": starts, "End": end, "Score": scores})
        df.insert(0, "Chromosome", chromosome)
        chromosome_data[c] = df

    return pr.PyRanges(chromosome_data)





def files_to_coverage(files):

    # have same writing-mechanism for both
    # just need to make bins into pyranges before
    # bigwig
    if args["bin"]:
        file_to_ranges = file_to_binned_ranges
    else:
        file_to_ranges = file_to_unbinned_ranges


    # read data and make ranges of coverage
    file_ranges = {}
    for f in file:
        cr = file_to_coverage_ranges(f)
        file_ranges[f] = file_to_ranges(f, args)


    # create coverage
    file_coverage = {}
    for n, gr in file_ranges.items():
        tmp[n] = gr

    return file_coverage


def main():

    args = {}

    treatment_ranges = files_to_coverage(args["treatment"], args)

    if "control" in args:
        control_ranges = files_to_coverage(args["treatment"], args)
        control_ranges = pr.concat(control_ranges)

    treatment_sum = pr.concat(treatment_ranges)
    control_sum = pr.concat(control_ranges)
