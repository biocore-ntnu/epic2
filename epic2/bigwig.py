import pandas as pd
import pyranges as pr
from epic2.src.genome_info import sniff
from epic2.src.reads_to_bins import files_to_bin_counts

# note that


def file_to_unbinned_ranges(f, args, _):

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

    if args["drop_duplicates"]:
        gr = gr.drop_duplicate_positions()

    # print("args", args)# ["chromsizes"])
    return gr #.remove_out_of_genome_bounds_intervals()



def file_to_binned_ranges(f, args, datatype):

    chromosome_data = {}
    # counts = 0
    bin_counts, _ = files_to_bin_counts([f], args, datatype)

    for c, v in bin_counts.items():
        starts, scores = v
        ends = starts + args["bin_size"]
        df = pd.DataFrame({"Start": starts, "End": ends, "Score": scores})
        df.insert(0, "Chromosome", c)
        chromosome_data[c] = df
        # counts += len(df)

    return pr.PyRanges(chromosome_data)





def files_to_coverage(files, args):

    # have same writing-mechanism for both
    # just need to make bins into pyranges before
    # bigwig
    if args["raw"]:
        file_to_ranges = file_to_unbinned_ranges
    else:
        file_to_ranges = file_to_binned_ranges


    # read data and make ranges of coverage
    file_ranges = {}
    for f in files:
        ranges = file_to_ranges(f, args, "ChIP")
        file_ranges[f] = ranges


    # create coverage
    file_coverage = {}
    for n, gr in file_ranges.items():
        file_coverage[n] = gr

    return file_coverage


def main(args):

    # TODO: need to create coverage of file if raw
    # else cluster Scores of binned

    treatment_ranges = files_to_coverage(args["treatment"], args)

    if args.get("control"):
        control_ranges = files_to_coverage(args["control"], args)
        control_sum = pr.concat(control_ranges.values())

    treatment_sum = pr.concat(treatment_ranges.values())

    print("treatment_sum", treatment_sum)
