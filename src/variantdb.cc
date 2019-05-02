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

#include "dot_graph.h"
#include "index.h"
#include "variant_graph.h"

using namespace variantdb;

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  main
 *  Description:  
 * =====================================================================================
 */
	int
main ( int argc, char *argv[] )
{
	if (argc < 2) {
		fprintf(stderr, "Please specify the log of the number of slots in the CQF.\n");
		exit(1);
	}

	std::string ref_file(argv[1]);
	std::string vcf_file(argv[2]);
	
	PRINT("Creating variant graph");
	std::vector<std::string> vcfs = {vcf_file};
	VariantGraph vg(ref_file, vcfs);

	PRINT("Creating Index");
	Index idx(&vg);

	PRINT("Serialing variant graph to disk");
	vg.serialize("./ser");

	PRINT("Serialing index to disk");
	idx.serialize("./ser");

	return EXIT_SUCCESS;
}				/* ----------  end of function main  ---------- */
