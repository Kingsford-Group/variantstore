/*
 * ============================================================================
 *
 *       Filename:  main.cc
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

#include "gqf_cpp.h"

#include "vcflib/Variant.h"

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
	uint64_t qbits = atoi(argv[1]);
	uint64_t keybits = 32;
	uint64_t nslots = (1ULL << qbits);
	uint64_t nvals = 500*nslots/1000;

	CQF<KeyObject> graph(nslots, keybits, keybits, QF_HASH_INVERTIBLE);

	PRINT(graph.total_slots());	

	/* Generate random values */
	uint32_t *vals = (uint32_t*)malloc(nvals*sizeof(vals[0]));
	RAND_bytes((unsigned char *)vals, sizeof(*vals) * nvals);
	for (uint32_t i = 0; i < nvals; i++) {
		vals[i] = (1 * vals[i]) % graph.range();
	}

	srand(time(NULL));
	for (uint32_t i = 0; i < nvals;) {
		uint32_t key = vals[i];
		uint32_t nedges = rand() % 4;		// up to 3 edges
		
		for (uint32_t j = 0; j < nedges; j++) {
			uint32_t tonode = vals[rand() % nvals];
			graph.insert(KeyObject(key, tonode), QF_NO_LOCK);
			i++;
		}
	}

	CQF<KeyObject>::Iterator it = graph.begin();	

	while (!it.done()) {
		KeyObject edge = *it;
		PRINT(edge.key << " - " << edge.value);
		++it;
	}

#if 0	
	vcflib::VariantCallFile variantFile;
	std::string filename = argv[1];
	variantFile.open(filename);
	vcflib::Variant var(variantFile);

	long int count = 0;
	while (variantFile.getNextVariant(var)) {
		count+= 1;
		std::cout << var << "\n";
	}
#endif
	return EXIT_SUCCESS;
}				/* ----------  end of function main  ---------- */
