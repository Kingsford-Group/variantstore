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
#include <vector>
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
	uint64_t nvals = 250*nslots/1000;

	CQF<KeyObject> graph(nslots, keybits, 1, QF_HASH_INVERTIBLE);

	PRINT(graph.total_slots());	

	/* Generate random values */
	uint32_t *vals = (uint32_t*)malloc(nvals*sizeof(vals[0]));
	RAND_bytes((unsigned char *)vals, sizeof(*vals) * nvals);
	for (uint32_t i = 0; i < nvals; i++) {
		vals[i] = (1 * vals[i]) % graph.range();
	}

	srand(time(NULL));
	std::vector<std::vector<uint32_t>> edgelist;
	for (uint32_t i = 0; i < nvals; i++) {
		uint32_t key = vals[i];
		uint32_t nedges = rand() % 4;		// up to 3 outgoing edges from a node
		
		if (nedges > 1)
			PRINT("Higher degree node: " << key);
		else if (nedges == 1)
			PRINT("Single degree node: " << key);

		for (uint32_t j = 0; j < nedges; j++) {
			uint32_t tonode = vals[rand() % nvals];
			uint64_t is_inplace = 0;
			uint32_t existing_node = graph.query_key(KeyObject(key, 0, 0),
																							 &is_inplace, QF_NO_LOCK);
			if (existing_node == 0) 	// new node
				graph.insert(KeyObject(key, 1, tonode), QF_NO_LOCK);
			else {		// existing node
				if (is_inplace == 1) { // there's only one outgoing edge
					// create a new vector and add to-nodes
					std::vector<uint32_t> tonodes;
					tonodes.emplace_back(existing_node);
					tonodes.emplace_back(tonode);
					edgelist.emplace_back(tonodes);
					uint32_t pointer = edgelist.size(); 
					// the pointer is always increamented by 1
					graph.replace_key(KeyObject(key, 1, existing_node), KeyObject(key, 0,
																																		pointer),
												QF_NO_LOCK);
				} else {
					edgelist[existing_node - 1].emplace_back(tonode);
				}
			}
		}
	}

	CQF<KeyObject>::Iterator it = graph.begin();	

	while (!it.done()) {
		KeyObject edge = *it;
		if (edge.value == 1)
			PRINT(edge.key << " - " << edge.count);
		else { 
			std::cout << edge.key << " - ";
			for (const auto tonode : edgelist[edge.count-1])
				std::cout << tonode << " ";
			std::cout << '\n';
		}
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
