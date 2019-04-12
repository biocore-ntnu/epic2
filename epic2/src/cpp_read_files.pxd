
from libcpp.pair cimport pair
from libcpp.string cimport string
from libcpp.vector cimport vector
from libcpp.map cimport map
from libc.stdint cimport uint32_t
# cdef from namespace "std":

ctypedef pair[string, char] key
ctypedef vector[uint32_t] intvec

ctypedef map[key, intvec] genome_map


cdef extern from "epic2/src/read_files.cpp":
    genome_map read_bed(const char *, uint32_t)
    genome_map read_bed_gz(const char *, uint32_t)
    genome_map read_bedpe(const char *, uint32_t)
    genome_map read_bedpe_gz(const char *, uint32_t)
    genome_map read_bam(const char *, uint32_t)




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
