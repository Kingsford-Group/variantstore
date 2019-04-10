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

//#include "variant_graph.h"

#include "vcflib/Variant.h"

//using namespace variantdb;


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

	vcflib::VariantCallFile variantFile;
	std::string filename = argv[1];
	variantFile.open(filename);
	vcflib::Variant var(variantFile);

	int sampleSize = variantFile.sampleNames.size();
	for (auto sample : variantFile.sampleNames)
		std::cout << sample << " ";
	std::cout << "\n";

	long int count = 0;
	while (variantFile.getNextVariant(var)) {
		count+= 1;
		std::cout << var << "\n";
	}

	return EXIT_SUCCESS;
}				/* ----------  end of function main  ---------- */

