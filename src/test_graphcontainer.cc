/*
 * ============================================================================
 *
 *       Filename:  test_graphcontainer.cc
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

#include "graph.h"

using namespace variantstore;
std::shared_ptr<spdlog::logger> console;

int test_correctness(const Graph& graph, std::unordered_map<uint32_t,
											std::unordered_set<uint32_t>>& adj_list,
											std::set<std::pair<uint32_t, uint32_t>>& edge_list) {
	// check graph correctness by iterating over @adj_list
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
	PRINT("Correctness check passed: by iterating over @adj_list!");

	// check graph correctness by iterating over Graph::VertexIterator
	Graph::VertexIterator vitr = graph.begin_vertices();
	while (!vitr.done()) {
		Graph::vertex v = *vitr;
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
		++vitr;
	}
	PRINT("Correctness check passed: by iterating over Graph::VertexIterator!");

	// check graph correctness by iterating over @edge_list
	for (auto it = edge_list.begin(); it != edge_list.end(); ++it) {
		if (!graph.is_edge(*it)) {
			PRINT("Graph edge: " << (*it).first << " - " << (*it).second);
			ERROR("correctness test failed!");
			return EXIT_FAILURE;
		}
	}
	PRINT("Correctness check passed: by iterating over @edge_list!");

	// check graph correctness by iterating over Graph::EdgeIterator
	Graph::EdgeIterator eitr = graph.begin_edges();
	while(!eitr.done()) {
		Graph::edge e = *eitr;
		if (edge_list.find(e) == edge_list.end()) {
			PRINT("Graph edge: " << e.first << " - " << e.second);
			ERROR("correctness test failed!");
			return EXIT_FAILURE;
		}
		++eitr;
	}
	PRINT("Correctness check passed: by iterating over Graph::EdgeIterator!");

	return EXIT_SUCCESS;
}

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
	uint64_t nvals = 95*nslots/100;
	//nslots *= 16;

	// create a typedef for the Graph type
	//typedef boost::adjacency_list<boost::setS, boost::vecS, boost::directedS>
		//BoostGraph;

	console = spdlog::default_logger();
	Graph graph(nslots);
	//BoostGraph bg(nslots);

	/* Generate random values */
	uint32_t *vals = (uint32_t*)malloc(nvals*sizeof(vals[0]));
	RAND_bytes((unsigned char *)vals, sizeof(*vals) * nvals);
	for (uint32_t i = 0; i < nvals; i++) {
		vals[i] = (1 * vals[i]) % nvals + 1;
	}

	srand(time(NULL));
	// to check the correctness of our graph implementation.
	std::unordered_map<uint32_t, std::unordered_set<uint32_t>> adj_list;
	std::set<std::pair<uint32_t, uint32_t>> edge_list;
	for (uint32_t i = 0; i < nvals; i++) {
		uint32_t key = vals[i];
		uint32_t nedges = rand() % 4 + 1;		// up to 4 outgoing edges from a node

		//if (nedges > 1)
			//PRINT("Higher degree node: " << key << " " << nedges);
		//else if (nedges == 1)
			//PRINT("Single degree node: " << key);
		
		std::unordered_set<uint32_t> vec;
		for (uint32_t j = 0; j < nedges; j++) {
			uint32_t tonode = vals[rand() % nvals];
			graph.add_edge(key, tonode);
			//boost::add_edge(key, tonode, bg);
			vec.insert(tonode);
			edge_list.insert(std::make_pair(key, tonode));
		}
		adj_list[key].merge(vec);
	}

	PRINT("Num vertices: " << graph.get_num_vertices());
	PRINT("Num edges: " << graph.get_num_edges());

	if (test_correctness(graph, adj_list, edge_list) == EXIT_FAILURE)
		return EXIT_FAILURE;

	PRINT("Serializing graph");
	graph.serialize("./ser");

	PRINT("Loading graph from disk");
	Graph file_graph("./ser");

	PRINT("Testing correctness for file graph");
	if (test_correctness(file_graph, adj_list, edge_list) == EXIT_FAILURE)
		return EXIT_FAILURE;

	return EXIT_SUCCESS;
}				/* ----------  end of function main  ---------- */
