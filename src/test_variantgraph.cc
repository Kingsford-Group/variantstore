/*
 * ============================================================================
 *
 *       Filename:  test_variantgraph.cc
 *
 *         Author:  Prashant Pandey (), ppandey2@cs.cmu.edu
 *   Organization:  Carnegie Mellon University
 *
 * ============================================================================
 */

#include <stdlib.h>
#include <assert.h>

#include <iostream>
#include <string>
#include <openssl/rand.h>

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>

#include "variant_graph.h"

using namespace variantdb;

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
	}
	PRINT("");
}

#include	<stdlib.h>

/* 
 * ===  FUNCTION  =============================================================
 *         Name:  main
 *  Description:  
 * ============================================================================
 */
	int
main ( int argc, char *argv[] )
{
	if (argc < 2) {
		fprintf(stderr, "Please specify the log of the number of slots in the CQF.\n");
		exit(1);
	}

	//std::string filename = argv[1];

	PRINT("Creating variant graph");
	std::string ref_file(argv[1]);
	std::string vcf_file(argv[2]);
	std::vector<std::string> vcfs = {vcf_file};
	VariantGraph vg(ref_file, vcfs);

	print_vg_info(vg, vcf_file);

	PRINT("Serialiing variant graph to disk");
	vg.serialize("./ser");

	PRINT("Loading variant graph from disk");
	VariantGraph file_vg("./ser");

	PRINT("Printing variant graph info from file vg");
	print_vg_info(file_vg, vcf_file);
	
	return EXIT_SUCCESS;
}				/* ----------  end of function main  ---------- */

