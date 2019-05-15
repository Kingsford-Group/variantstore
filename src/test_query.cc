#include <iostream>
#include <fstream>
#include <vector>
#include <time.h>
#include <string>

#include "variant_graph.h"
#include "dot_graph.h"
#include "query.h"

using namespace variantdb;

int
main ( int argc, char *argv[] )
{
	if (argc < 2) {
		fprintf(stderr, "Please specify the reference fasta file and vcf file.\n");
		exit(1);
	}

	//Construct vg & index
	PRINT("Creating variant graph");
	std::string ref_file(argv[1]);
	std::string vcf_file(argv[2]);
	std::vector<std::string> vcfs = {vcf_file};
	VariantGraph vg(ref_file, vcfs);

	PRINT("Creating Index");
	Index idx(&vg);

/*
  PRINT("TEST get_prev_vertex_with_sample");
  int pos = 1;
  std::string sample_id;
  PRINT("Please specify position: ");
  std::cin >> pos;
  PRINT("Please specify sample id: ");
  std::cin.ignore();
  std::getline(std::cin, sample_id);
  while (pos != -1)
  {
    uint64_t ref_pos = 0;
    uint64_t sample_pos = 0;
    PRINT("Query...");
    get_prev_vertex_with_sample(&vg, &idx, pos, sample_id, ref_pos, sample_pos);
    PRINT("Please specify position: ");
    std::cin >> pos;
    PRINT("Please specify sample id: ");
    std::cin.ignore();
    std::getline(std::cin, sample_id);
  }
*/
  PRINT("TEST query_sample_from_ref");
  int pos_x = 1;
  int pos_y = 1;

  while (pos_x != -1)
  {
    std::string sample_id;
    PRINT("Please specify position_x: ");
    std::cin >> pos_x;
    PRINT("Please specify position_y: ");
    std::cin >> pos_y;
    PRINT("Please specify sample id: ");
    std::cin.ignore();
    std::getline(std::cin, sample_id);
    PRINT("Query...");
    std::string seq = query_sample_from_ref(&vg, &idx, pos_x, pos_y, sample_id);
    PRINT("Sample " << sample_id << " has sequence " << seq
          << " in reference coordinate.");
  }

}



// test_get_prev_vertex_with_sample
