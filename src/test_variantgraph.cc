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

#include "spdlog/spdlog.h"

#include "dot_graph.h"
#include "variant_graph.h"

using namespace variantdb;
std::shared_ptr<spdlog::logger> console;

void print_vg_info(VariantGraph& vg, std::string& vcf_file,
									 std::vector<std::string> sampleNames) {
	PRINT("Graph stats:");
	PRINT("Chromosome: " << vg.get_chr() << " #Vertices: " << vg.get_num_vertices()
				<< " Seq length: " << vg.get_seq_length() << " Ref length: " <<
				vg.get_ref_length());

	PRINT("Variant Graph nodes:");
	//for (int i = 0; i < 8; i++) {
	auto bfs = vg.find(0);
	//std::cout << "radius: " << i << " -- ";
	while (!bfs.done()) {
		std::cout << (*bfs)->vertex_id() << " ";
		++bfs;
	}
	PRINT("");
	//}

	auto itr = vg.find("ref");
	PRINT("Ref nodes:");
	while (!itr.done()) {
		//vg.print_vertex_info(**itr);
		std::cout << vg.get_sequence(**itr);
		++itr;
	}

	PRINT("");
	PRINT("Samples:");

	// get all samples
	for (auto sample : sampleNames) {
		PRINT("Sample: " << sample);
		auto itr = vg.find(sample);
		while (!itr.done()) {
			//vg.print_vertex_info(**itr);
			std::cout << vg.get_sequence(**itr);
			++itr;
		}
		PRINT("");
	}
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

	console = spdlog::default_logger();
#ifdef DEBUG_MODE
	console->set_level(spdlog::level::debug);
#endif

	PRINT("Creating variant graph");
	std::string ref_file(argv[1]);
	std::string vcf_file(argv[2]);
	VariantGraph vg(ref_file, vcf_file, "./ser");

	vcflib::VariantCallFile variantFile;
	variantFile.open(vcf_file);
	vcflib::Variant var(variantFile);

	print_vg_info(vg, vcf_file, variantFile.sampleNames);

	PRINT("Serialiing variant graph to disk");
	vg.serialize();

	PRINT("Loading variant graph from disk");
	VariantGraph file_vg("./ser", READ_COMPLETE_GRAPH);

	PRINT("Printing variant graph info from file vg");
	print_vg_info(file_vg, vcf_file, variantFile.sampleNames);

	createDotGraph(&file_vg, "./ser");

	return EXIT_SUCCESS;
}				/* ----------  end of function main  ---------- */

