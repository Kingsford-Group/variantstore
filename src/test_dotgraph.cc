/*
 * ============================================================================
 *
 *       Filename:  test_graphcontainer.cc
 *
 *         Author:  Prashant Pandey (), ppandey2@cs.cmu.edu
 										Yinjie Gao <yinjieg@andrew.cmu.edu>
 *   Organization:  Carnegie Mellon University
 *
 * ============================================================================
 */


#include "dot_graph.h"
// Create vg & output dot format graph
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

#include "index.h"
#include "variant_graph.h"

using namespace variantstore;
std::shared_ptr<spdlog::logger> console;

void print_index_info(const Index &idx, uint64_t len)
{
	PRINT("Sequence length is");
	PRINT("Pos\tNode_id");

	for (uint64_t i=0; i<len; i++)
	{
		uint64_t node_id = idx.find(i);
		PRINT(i << ":" << node_id);
	}
	return;
}

void print_vg_info(VariantGraph& vg, std::string& vcf_file) {
	PRINT("Graph stats:");
	PRINT("Chromosome: " << vg.get_chr() << " #Vertices: " << vg.get_num_vertices()
				<< " Seq length: " << vg.get_seq_length());

	PRINT("Variant Graph nodes:");
	for (int i = 0; i < 8; i++) {
		auto bfs = vg.find(0, i);
		std::cout << "radius: " << i << " -- ";
		while (!bfs.done()) {
			std::cout << (*bfs)->vertex_id() << " ";
			++bfs;
		}
		PRINT("");
	}

	auto itr = vg.find("ref");
	PRINT("Ref nodes:");
	while (!itr.done()) {
		//vg.print_vertex_info(**itr);
		std::cout << vg.get_sequence(**itr);
		++itr;
	}

	PRINT("");
	PRINT("Samples:");
	vcflib::VariantCallFile variantFile;
	variantFile.open(vcf_file);
	vcflib::Variant var(variantFile);

	// get all samples
	for (auto sample : variantFile.sampleNames) {
		PRINT("Sample: " << sample);
		auto itr = vg.find(sample);
		while (!itr.done()) {
			//vg.print_vertex_info(**itr);
			std::cout << vg.get_sequence(**itr);
			++itr;
		}
		PRINT("");
	}
	PRINT("");
}


	int
main ( int argc, char *argv[] )
{
	if (argc < 2) {
		fprintf(stderr, "Please specify the reference fasta file and vcf file.\n");
		exit(1);
	}

	console = spdlog::default_logger();
#ifdef DEBUG_MODE
	console->set_level(spdlog::level::debug);
#endif
	//std::string filename = argv[1];

	PRINT("Creating variant graph");
	std::string ref_file(argv[1]);
	std::string vcf_file(argv[2]);
	std::vector<std::string> vcfs = {vcf_file};
	VariantGraph vg(ref_file, vcfs, READ_COMPLETE_GRAPH);

	print_vg_info(vg, vcf_file);
	createDotGraph(&vg, "graph.dot");
	return 0;
}
