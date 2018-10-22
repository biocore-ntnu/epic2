#include <iostream>
#include <algorithm>
#include <fstream>
#include <map>
#include <vector>
#include <string>
#include <vector>
#include <iterator>
#include <unordered_map>

// struct pair_hash {
//   template <class T1, class T2>
//   std::size_t operator () (const std::pair<T1,T2> &p) const {
//     auto h1 = std::hash<T1>{}(p.first);
//     auto h2 = std::hash<T2>{}(p.second);

//     // Mainly for demonstration purposes, i.e. works but is overly simple
//     // In the real world, use sth. like boost.hash_combine
//     return h1 ^ h2;
//   }
// };


typedef std::vector<uint32_t> intvec;
typedef std::pair <std::string, char> key;
typedef std::map<key, intvec> genome_map;
// typedef std::unordered_map<std::string, intvec> bin_map;



genome_map read_bed(char const* fileName)
{
  std::ifstream file(fileName);

  std::string   chromosome;
  std::string   junk;
  std::uint32_t start;
  std::uint32_t end;
  char strand;

  key chrom_strand;
  std::uint32_t five_end;
  genome_map genome;

  while(file >> chromosome >> start >> end >> junk >> junk >> strand)
    {

      chrom_strand = std::make_pair(chromosome, strand);

      if (strand == '+'){
        five_end = start;
      } else {
        five_end = end;
      }

      genome[chrom_strand].push_back(five_end);

    }

  return genome;
}

// #include <htslib/sam.h>
// #include <htslib/hts.h>
// #include <samtools/samtools.h>

// #include <htslib/hts.h>
// #include <seqan/bam_io.h>


// genome_map read_bam(char const* fileName)
// {
//   samFile *fp_in = hts_open(fileName, "r");
// 	bam_hdr_t *bamHdr = sam_hdr_read(fp_in); //read header
// 	bam1_t *aln = bam_init1(); //initialize an alignment
//   std::uint32_t pos;
//   std::uint32_t len;
//   std::string chr;
//   char forward = '+';
//   char reverse = '-';
//   key chrom_strand;
//   genome_map genome;

//   while(sam_read1(fp_in,bamHdr,aln) > 0){

// 		pos = aln->core.pos + 1; //left most position of alignment in zero based coordianate (+1)
// 		chr = std::string(bamHdr->target_name[aln->core.tid]); //contig name (chromosome)
// 		len = aln->core.l_qseq; //length of the read.

//     if (aln->core.flag & 0x10){ // reverse
//       chrom_strand = std::make_pair(chr, reverse);
//       genome[chrom_strand].push_back(pos + len);
//     } else {
//       chrom_strand = std::make_pair(chr, forward);
//       genome[chrom_strand].push_back(pos);
//     }

//   }

//   return genome;
// }
