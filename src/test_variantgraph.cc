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

	std::cout << "Creating VG" << '\n';
	std::string ref_file(argv[1]);
	std::string vcf_file(argv[2]);
	std::vector<std::string> vcfs = {vcf_file};
	VariantGraph vg(ref_file, vcfs);
	std::cout << vg.get_chr() << ' ' << vg.get_num_vertices() << " " <<
		vg.get_seq_length() << '\n';

	return EXIT_SUCCESS;
}				/* ----------  end of function main  ---------- */

