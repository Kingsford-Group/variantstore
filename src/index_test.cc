#include <sdsl/bit_vectors.hpp>
#include <sdsl/int_vector.hpp>
#include <sdsl/util.hpp>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <time.h>
#include <map>
#include <iterator>

int main()
{
    std::cout << "Construct bit_vector ..." << std::endl;
    const uint64_t node_sz = 306009792;
    const uint64_t REF_GENOME_LEN = 3099706404;
    const uint16_t BLOCK_SIZE = 127;
    sdsl::bit_vector b(REF_GENOME_LEN, 0);
    sdsl::int_vector<>node_list(node_sz, 0, 64);

    std::cout << "Construct index ..." << std::endl;
    // Create random position
    std::vector<uint64_t> pos;
    for (uint64_t i = 0; i < REF_GENOME_LEN; i ++)
      pos.push_back(i);

    std::srand (unsigned(std::time(0)));
    std::random_shuffle(pos.begin(), pos.end());

    std::vector<uint64_t> node_pos;
    for (uint64_t i = 0; i < node_sz; i ++)
      node_pos.push_back(pos[i]);

    std::sort(node_pos.begin(), node_pos.end());

    // Create random node_list
    std::vector<uint64_t> double_nodes;
    for (uint64_t i = 0; i < node_sz * 2; i ++)
      double_nodes.push_back(i);

    std::random_shuffle(double_nodes.begin(), double_nodes.end());
    std::vector<uint64_t> node_ids;
    for (uint64_t i = 0; i < node_sz; i ++)
      node_ids.push_back(double_nodes[i]);

    std::sort(node_ids.begin(), node_ids.end());

    std::cout << "Construct index map ..." << std::endl;

    std::map<uint64_t, uint64_t> idx_map;

    for (uint64_t i = 0; i < node_sz; i ++)
    {
      b[node_pos[i]] = 1;
      node_list[i] = node_ids[i];

      idx_map[node_pos[i]] = node_ids[i];
    }


    sdsl::rrr_vector<BLOCK_SIZE> rrrb;
    sdsl::util::assign(rrrb, sdsl::rrr_vector<BLOCK_SIZE>(b));
    sdsl::rrr_vector<BLOCK_SIZE>::rank_1_type rank_rrrb(&rrrb);

    std::cout<< "size of rrrb in MB: " << size_in_mega_bytes(rrrb)<< std::endl;
    std::cout<< "size of rank_rrrb in MB: " << size_in_mega_bytes(rank_rrrb)<< std::endl;

    sdsl::util::bit_compress(node_list);
    std::cout<< "size of compressed int_vec in MB: " << size_in_mega_bytes(node_list) << std::endl;

    std::cout << "Correctness test ... ";

    for (uint64_t i = 0; i < REF_GENOME_LEN; i ++)
    {
      uint64_t rank = rank_rrrb(i);
      auto it = idx_map.lower_bound(i);


      if (rank == 0 || it == idx_map.begin())
      {
        if (rank != 0 || it != idx_map.begin())
        {
          std::cout << "failed at pos " << i  << std::endl;
          std::cout << "rank in bit vector is: " << rank << std::endl;
          return 1;
        }
        continue;
      }

      uint64_t id1 = node_list[rank-1];
      it--;
      uint64_t id2 = it->second;

      if (id1 != id2)
      {
        std::cout << "failed at pos " << i << " with id " << id1 << "," << id2 << std::endl;
        return 1;
      }
    }

    std:: cout << "successfull" << std::endl;

    std:: cout << "Serialization..." << std::endl;
    std::ofstream out("index.sdsl");
    serialize(rrrb, out);
    serialize(node_list, out);
    return 0;

}
