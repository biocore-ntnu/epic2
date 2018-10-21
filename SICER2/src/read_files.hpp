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
