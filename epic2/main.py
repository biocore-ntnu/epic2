from sys import argv
import sys
import os

from epic2.src.reads_to_bins import files_to_bin_counts
from epic2.src.find_islands import find_islands, compute_fdr, write_islands, add_chip_count_to_islands
from epic2.src.genome_info import egl_and_chromsizes

from collections import OrderedDict


def _main(args):

    args["drop_duplicates"] = int(not args["keep_duplicates"])

    import logging

    if args["quiet"]:
        level = logging.CRITICAL
    else:
        level = logging.INFO

    logging.basicConfig(
        level=level,
        format='%(message)s ',
        datefmt='%a, %d %b %Y %H:%M:%S',
        stream=sys.stderr)

    output = args["output"]
    if output:
        outfile_directory = os.path.dirname(output)
        if outfile_directory and not os.path.exists(outfile_directory):
            os.makedirs(outfile_directory)

    effective_genome_length, chromsizes = egl_and_chromsizes(args)
    args["chromsizes_"] = chromsizes
    args["effective_genome_size"] = effective_genome_length

    c_bins_counts, chip_count_before = files_to_bin_counts(
        args["treatment"], args, "ChIP")
    chip_count = sum(sum(counts) for _, counts in c_bins_counts.values())

    logging.info(
        "\nValid ChIP reads: {} ({} before out of bounds removal)\n".format(
            chip_count, chip_count_before))

    if args["original_statistics"]:
        from epic2.src.SICER_stats import compute_score_threshold
    else:
        from epic2.src.SICER_stats2 import compute_score_threshold

    score_threshold, island_enriched_threshold, average_window_readcount = compute_score_threshold(
        chip_count, args["bin_size"], effective_genome_length,
        args["gaps_allowed"] * args["bin_size"], args["e_value"])

    logging.info(
        "Number of tags in a window: {}\n".format(island_enriched_threshold))

    islands = find_islands(c_bins_counts, args["gaps_allowed"],
                           args["bin_size"], score_threshold,
                           island_enriched_threshold, average_window_readcount)
    new_islands = add_chip_count_to_islands(islands, c_bins_counts, args)

    logging.info("Number of islands found: {}\n".format(
        sum(len(i) for i in islands.values())))

    if args["control"]:

        b_bins_counts, background_count_before = files_to_bin_counts(
            args["control"], args, "Input")
        background_count = sum(
            sum(counts) for _, counts in b_bins_counts.values())
        logging.info(
            "\nValid Background reads: {} ({} before out of bounds removal)\n".
            format(background_count, background_count_before))

        if args["original_algorithm"]:
            cc, bc = chip_count_before, background_count_before
        else:
            cc, bc = chip_count, background_count

        compute_fdr(new_islands, b_bins_counts, cc, bc,
                    effective_genome_length,
                    args["false_discovery_rate_cutoff"], args["output"])
    else:
        write_islands(new_islands, average_window_readcount,
                      args["false_discovery_rate_cutoff"], args["e_value"],
                      args["output"])

        if args["original_algorithm"]:
            cc = chip_count_before
        else:
            cc = chip_count

    return c_bins_counts
