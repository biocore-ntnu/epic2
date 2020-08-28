from libc.stdint cimport uint32_t


cdef extern from "epic2/src/read_files.cpp":
    cdef struct interval:
        uint32_t start
        uint32_t end

    genome_map read_bed(const char *, uint32_t)
    genome_map read_bed_gz(const char *, uint32_t)
    genome_map read_bedpe(const char *, uint32_t)
    genome_map read_bedpe_gz(const char *, uint32_t)

    genome_map intervals_to_midpoint(genome_intervals genome, uint32_t drop_duplicates)
    genome_map intervals_to_five_end(genome_intervals genome, uint32_t drop_duplicates)
# cdef extern from "../../lib/htslib/htslib/hts.h":
#     hts_open(const char *fn, const char *mode)



# cdef extern from "lib/htslib/htslib/sam.h":
#     ctypedef struct bam_hdr_t:
#         pass

#     ctypedef struct samFile:
#         pass

    pass
    # bam_hdr_t *sam_hdr_read(samFile *fp)
#     bam_init1()

# cdef extern from "SICER2/src/"



# typedef struct {
#     int32_t n_targets, ignore_sam_err;
#     uint32_t l_text;
#     uint32_t *target_len;
#     int8_t *cigar_tab;
#     char **target_name;
#     char *text;
#     void *sdict;
# } bam_hdr_t;
