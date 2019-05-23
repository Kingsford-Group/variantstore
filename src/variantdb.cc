/*
 * ============================================================================
 *
 *       Filename:  variantdb.cc
 *
 *         Author:  Prashant Pandey (), ppandey2@cs.cmu.edu
 *   Organization:  Carnegie Mellon University
 *
 * ============================================================================
 */

#include <stdlib.h>

#include "spdlog/spdlog.h"

#include "dot_graph.h"
#include "index.h"
#include "variant_graph.h"

using namespace variantdb;

std::shared_ptr<spdlog::logger> console;

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

	console = spdlog::default_logger();
#ifdef DEBUG_MODE
	console->set_level(spdlog::level::debug);
#endif

	std::string ref_file(argv[1]);
	std::string vcf_file(argv[2]);
	
	console->info("Creating variant graph");
	std::vector<std::string> vcfs = {vcf_file};
	VariantGraph vg(ref_file, vcfs);

	console->info("Graph stats:");
	console->info("Chromosome: {} #Vertices: {} #Edges: {} Seq length: {}",
								vg.get_chr(), vg.get_num_vertices() , vg.get_num_edges(),
								vg.get_seq_length());
	console->info("Serializing variant graph to disk");
	vg.serialize(argv[3]);

	console->info("Creating Index");
	Index idx(&vg);
	console->info("Serializing index to disk");
	idx.serialize(argv[3]);

	return EXIT_SUCCESS;
}				/* ----------  end of function main  ---------- */
