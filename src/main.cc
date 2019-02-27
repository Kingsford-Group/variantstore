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
#include "graph.h"

#include "vcflib/Variant.h"

using namespace variantdb;

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
	uint64_t nslots = (1ULL << qbits);
	uint64_t nvals = 10*nslots/100;

	Graph graph(nslots);

	/* Generate random values */
	uint32_t *vals = (uint32_t*)malloc(nvals*sizeof(vals[0]));
	RAND_bytes((unsigned char *)vals, sizeof(*vals) * nvals);
	for (uint32_t i = 0; i < nvals; i++) {
		vals[i] = (1 * vals[i]) % UINT32_MAX;
	}

	srand(time(NULL));
	// to check the correctness of our graph implementation.
	std::unordered_map<uint32_t, std::unordered_set<uint32_t>> adj_list;
	for (uint32_t i = 0; i < nvals; i++) {
		uint32_t key = vals[i];
		uint32_t nedges = rand() % 4;		// up to 3 outgoing edges from a node
		
		//if (nedges > 1)
			//PRINT("Higher degree node: " << key << " " << nedges);
		//else if (nedges == 1)
			//PRINT("Single degree node: " << key);
		
		std::unordered_set<uint32_t> vec;
		for (uint32_t j = 0; j < nedges; j++) {
			uint32_t tonode = vals[rand() % nvals];
			graph.add_edge(key, tonode);
			vec.insert(tonode);
		}
		adj_list[key].merge(vec);
	}

	PRINT("Num vertices: " << graph.get_num_vertices());
	PRINT("Num edges: " << graph.get_num_edges());

	for (auto it = adj_list.begin(); it != adj_list.end(); ++it) {
		Graph::vertex_set neighbors = graph.out_neighbors(it->first);

		if (it->second != neighbors) {
			std::cout << it->first << " - ";
			std::cout << '\n';
			for (auto const n : neighbors) {
				std::cout << n << " ";
			}

			std::cout << '\n';
			for (auto const n : it->second) {
				std::cout << n << " ";
			}
			std::cout << '\n';
			ERROR("correctness test failed!");
			return EXIT_FAILURE;
		}
	}

	Graph::VertexIterator itr = graph.begin_vertex();

	while (!itr.done()) {
		Graph::vertex v = *itr;
		Graph::vertex_set neighbors = graph.out_neighbors(v);

		if (adj_list[v] != neighbors) {
			std::cout << v << " - ";
			std::cout << '\n';
			for (auto const n : neighbors) {
				std::cout << n << " ";
			}

			std::cout << '\n';
			for (auto const n : adj_list[v]) {
				std::cout << n << " ";
			}
			std::cout << '\n';
			ERROR("correctness test failed!");
			return EXIT_FAILURE;
		}
		//assert(adj_list[v] == neighbors);
		++itr;
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
