
def remove_out_of_bounds_bins(count_dict, chromsizes, bin_size):


    for chromosome, (bins, counts) in count_dict.items():
        if not chromosome in chromsizes:
            print("Chromosome {} found in file, but not in the chromsizes dict for the genome.".format(chromosome), file=sys.stderr)
            continue

        chromsize = chromsizes[chromosome]

        try:
            i = len(bins) - 1
            while (i >= 0 and (bins[i]) > chromsize):
                # print("For {} bin {} is out of bounds ({})\n".format(chromosome, bins[i], chromsize))
                # print("It has a count of ({})\n".format(counts[i]))
                # 52280 - 51862 = 418
                i -= 1
        except OverflowError as e:
            print("\nAdditional info:\n")
            print("Chromosome:", chromosome)
            print("Chromsize:", chromsize)
            raise e
        # sys.stderr.write(chromosome + "\n")
        # sys.stderr.write("{} {}\n".format(bins[len(bins) - 1], chromsize))

        if i != len(bins) - 1:
            # print("-----")
            # print(chromosome)
            # print("length before:", len(counts))
            count_dict[chromosome] = bins[:i+1], counts[:i+1]
            # print("length after:", len(counts[:i + 1]))
