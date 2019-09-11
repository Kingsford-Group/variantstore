/*
 * ============================================================================
 *
 *       Filename:  construct.cc
 *
 *         Author:  Prashant Pandey (), ppandey2@cs.cmu.edu
 *   Organization:  Carnegie Mellon University
 *
 * ============================================================================
 */

#include <stdlib.h>

#include "spdlog/spdlog.h"

#include "index.h"
#include "variant_graph.h"
#include "progopts.h"


using namespace variantdb;

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  main
 *  Description:  
 * =====================================================================================
 */
	int
construct_main (ConstructOpts &opts)
{
	console->info("Creating variant graph");
	VariantGraph vg(opts.ref, opts.vcf, opts.prefix);

	console->info("Graph stats:");
	console->info("Chromosome: {} #Vertices: {} #Edges: {} Seq length: {}",
								vg.get_chr(), vg.get_num_vertices() , vg.get_num_edges(),
								vg.get_seq_length());
	console->info("Serializing variant graph to disk");
	vg.serialize();

	console->info("Creating Index");
	Index idx(&vg);
	console->info("Serializing index to disk");
	idx.serialize(opts.prefix);

	console->info("Loading variant graph");
	VariantGraph file_vg(opts.prefix, READ_COMPLETE_GRAPH);
	console->info("Graph stats:");
	console->info("Chromosome: {} #Vertices: {} #Edges: {} Seq length: {}",
								file_vg.get_chr(), file_vg.get_num_vertices() , file_vg.get_num_edges(),
								file_vg.get_seq_length());

	return EXIT_SUCCESS;
}				/* ----------  end of function main  ---------- */
