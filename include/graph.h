/*
 * ============================================================================
 *
 *       Filename:  graph.h
 *
 *         Author:  Prashant Pandey (), ppandey2@cs.cmu.edu
 *   Organization:  Carnegie Mellon University
 *
 * ============================================================================
 */

#ifndef __GRAPH_H__
#define __GRAPH_H__

#include <stdlib.h>

#include <vector>

#include "gqf_cpp.h"

namespace variantdb {

	// An efficient way to store graph topology for sparse (and skewed) DAGs. 
	// Vertexes are defined as 32-bit integers.
	// It uses a variable-length encoding from the QF (counter encoding) data
	// structure to endcode the neighbors of a vertex. 
	// If there's only one neighbor then it is stored next to the source vertex
	// in the QF. If there are more than one neighbors then they are stored in a
	// vector and a pointer is stored next to the source vertex in the QF.
	template <class key_obj>
		class Graph {
			public:
				typedef uint32_t vertex;
				typedef std::pair<uint32_t, uint32_t> edge;
				typedef vector<vertex> vertex_set;
				typedef vector<edge> edge_set;

				Graph();
				Graph(uint64_t size);

				int add_edge(vertex s, vertex d);

				vertex_set get_neighbors(vertex s);

				class Iterator {

				};


			private:
				CQF<key_obj> adj_list;
		};
}

#endif	// __GRAPH_H_
