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


genome_map read_bedpe(char const* fileName)
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

  while(file >> chromosome >> start >> junk >> junk >> junk >> end >> junk >> junk >> strand >> junk)
    {

      chrom_strand = std::make_pair(chromosome, strand);

      five_end = start + ((end - start) / 2);

      genome[chrom_strand].push_back(five_end);

    }

  return genome;
}

// #include <api/BamReader.h>
// using namespace BamTools;
// #include <htslib/sam.h>

// // #include <samtools/samtools.h>

// genome_map read_bam(char const* fileName){

//   genome_map genome;
//   samFile *fp_in = hts_open(fileName,"r"); //open bam file
//   bam_hdr_t *bamHdr = sam_hdr_read(fp_in); //read header
//   bam1_t *aln = bam_init1(); //initialize an alignment

//   char *chrom = "chr2";
//   int locus = 5;
//   int comp ;

//   printf("%s\t%d\n", chrom, locus);

//   //header parse
//   //uint32_t *tar = bamHdr->text ;
//   //uint32_t *tarlen = bamHdr->target_len ;

//   //printf("%d\n",tar);

//   while(sam_read1(fp_in,bamHdr,aln) > 0){

//     int32_t pos = aln->core.pos +1; //left most position of alignment in zero based coordianate (+1)
//     char *chr = bamHdr->target_name[aln->core.tid] ; //contig name (chromosome)
//     uint32_t len = aln->core.l_qseq; //length of the read.

//     uint8_t *q = bam_get_seq(aln); //quality string
//     uint32_t q2 = aln->core.qual ; //mapping quality


//     char *qseq = (char *)malloc(len);

//     // for(int i=0; i< len ; i++){
//     //   qseq[i] = seq_nt16_str[bam_seqi(q,i)]; //gets nucleotide id and converts them into IUPAC id.
//     // }

//     //printf("%s\t%d\t%d\t%s\t%s\t%d\n",chr,pos,len,qseq,q,q2);

//     if(strcmp(chrom, chr) == 0){

//       if(locus > pos+len){
//         // printf("%s\t%d\t%d\t%s\t%s\t%d\n",chr,pos,len,qseq,q,q2);
//         printf("%s\t%d\t%d\t%s\t%d\n",chr,pos,len,q,q2);
//       }
//     }
//   }

//   bam_destroy1(aln);
//   sam_close(fp_in);

//   return genome;
// }

//   // samFile *fp_in = hts_open(fileName,"r"); //open bam file
//   // bam_hdr_t *bamHdr = sam_hdr_read(fp_in); //read header
// 	// bam1_t *aln = bam_init1(); //initialize an alignment

// genome_map read_bam(char const* fileName){

//   genome_map genome;
//   BamTools::BamReader reader = BamTools::BamReader();

//   if ( !reader.Open(fileName) ) {
//     std::cerr << "Could not open input BAM files" << " " << fileName << std::endl;
//     return genome;
//   }

//   BamTools::BamAlignment al;
//   const BamTools::SamHeader header = reader.GetHeader();
//   const BamTools::RefVector references = reader.GetReferenceData();

//   while ( reader.GetNextAlignmentCore(al) ) {
//     // if ( al.MapQuality >= 90 )
//     std::cout << al.MapQuality << std::endl;

//       // writer.SaveAlignment(al);
//   }

//   reader.Close();

//   return genome;
// }

// gcc test.c /mnt/work/me/software/anaconda/pkgs/htslib-1.9-hc238db4_4/lib/libhts.a  -I /mnt/work/me/software/anaconda/pkgs/htslib-1.9-hc238db4_4/include/
