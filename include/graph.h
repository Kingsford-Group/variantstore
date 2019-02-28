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

#include <unordered_set>
#include <vector>

#include "gqf_cpp.h"

namespace variantdb {

#define DEFAULT_SIZE (1ULL << 16)
#define KEYBITS 32

	// An efficient way to store graph topology for sparse (and skewed) DAGs. 
	// Vertices are defined as 32-bit integers.
	// There can be only one edge between two vertices, i.e., simple graphs.
	// It uses a variable-length encoding from the QF (counter encoding) data
	// structure to encode the neighbors of a vertex. 
	// If there's only one neighbor then it is stored next to the source vertex
	// in the QF. If there are more than one neighbors then they are stored in a
	// vector and a pointer is stored next to the source vertex in the QF.
	class Graph {
		public:
			typedef uint32_t vertex;
			typedef std::pair<uint32_t, uint32_t> edge;
			typedef std::unordered_set<vertex> vertex_set;
			typedef std::unordered_set<edge> edge_set;

			Graph();
			Graph(uint32_t size);

			int add_edge(vertex s, vertex d);
			int remove_edge(vertex s, vertex d);

			vertex_set out_neighbors(vertex s) const;
			uint32_t out_degree(vertex v) const;

			vertex_set in_neighbors(vertex s) const;
			uint32_t in_degree(vertex v) const;

			uint32_t get_num_vertices() const;
			uint32_t get_num_edges() const;

			// Iterate over all nodes in the graph
			class VertexIterator {
				public:
					VertexIterator(CQF<KeyObject>::Iterator it);
					vertex operator*(void) const;
					void operator++(void);
					bool done(void) const;

				private:
					CQF<KeyObject>::Iterator adj_list_itr;
			};

			VertexIterator begin_vertex(void) const;
			VertexIterator end_vertex(void) const;

			// Iterate over all edges in the graph
			class EdgeIterator {
				public:
					EdgeIterator(CQF<KeyObject>::Iterator it);
					edge operator*(void) const;
					void operator++(void);
					bool done(void) const;

				private:
					CQF<KeyObject>::Iterator adj_list_itr;
			};

			EdgeIterator begin_edge(void) const;
			EdgeIterator end_edge(void) const;

		private:
			CQF<KeyObject> adj_list;
			std::vector<std::unordered_set<vertex>> aux_vertex_list;
			uint32_t num_edges{0};
	};

	Graph::Graph() {
		adj_list = CQF<KeyObject>(DEFAULT_SIZE, KEYBITS, 1, QF_HASH_INVERTIBLE);
		adj_list.set_auto_resize();
	}

	Graph::Graph(uint32_t size) {
		adj_list = CQF<KeyObject>(size, KEYBITS, 1, QF_HASH_INVERTIBLE);
		adj_list.set_auto_resize();
	}

	int Graph::add_edge(vertex s, vertex d) {
		uint64_t is_inplace {0};
		vertex val = adj_list.query_key(KeyObject(s, 0, 0), &is_inplace,
																		QF_NO_LOCK);
		if (d == 0)
			return 0;
		if (is_inplace == 1 && val == d)	// edge already exist
			return 0;	
		else if (val == 0) 	{	// new vertex
			num_edges++;
			return adj_list.insert(KeyObject(s, 1, d), QF_NO_LOCK);
		}
		else {		// existing node
			if (is_inplace == 1) { // there's only one outgoing edge
				// create a new vector and add outgoing vertices 
				std::unordered_set<vertex> neighbors;
				neighbors.insert(val);
				neighbors.insert(d);
				aux_vertex_list.emplace_back(neighbors);
				uint32_t pointer = aux_vertex_list.size();	// it is an offset in the vector. 
				num_edges++;
				// the pointer is always increamented by 1
				return adj_list.replace_key(KeyObject(s, 1, val), KeyObject(s, 0, pointer),
																		QF_NO_LOCK);
			} else {
				if (aux_vertex_list[val - 1].insert(d).second)
					num_edges++;
			}
		}
		return 0;
	}

	int Graph::remove_edge(vertex s, vertex d) {
		uint64_t is_inplace {0};
		vertex val = adj_list.query_key(KeyObject(s, 0, 0), &is_inplace,
																		QF_NO_LOCK);
		if (val == 0)
			return 0;
		else {
			if (is_inplace == 1) {
				return adj_list.delete_key(KeyObject(s, 1, d), QF_NO_LOCK);
			} else {
				//uint32_t idx = 0;
				for (auto const vertex : aux_vertex_list[val - 1]) {
					if (vertex == d)
						aux_vertex_list[val - 1].erase(aux_vertex_list[val - 1].begin());
				}
			}
		}
		return 0;
	}

	Graph::vertex_set Graph::out_neighbors(vertex s) const {
		Graph::vertex_set neighbor_set;
		uint64_t is_inplace {0};
		vertex val = adj_list.query_key(KeyObject(s, 0, 0), &is_inplace,
																		QF_NO_LOCK);
		if (val == 0)
			return neighbor_set;
		else {
			if (is_inplace == 1) {
				neighbor_set.insert(val);
			} else {
				neighbor_set = aux_vertex_list[val - 1];
			}
			return neighbor_set;
		}
	}

	uint32_t Graph::out_degree(vertex v) const {
		uint64_t is_inplace {0};
		vertex val = adj_list.query_key(KeyObject(v, 0, 0), &is_inplace,
																		QF_NO_LOCK);
		if (val == 0)
			return 0;
		else {
			if (is_inplace == 1) {
				return 1;	
			} else {
				return aux_vertex_list[val - 1].size();
			}
		}
		return 0;
	}

	uint32_t Graph::get_num_vertices() const {
		return adj_list.dist_elts();
	}

	uint32_t Graph::get_num_edges() const {
		return num_edges;
	}

	Graph::VertexIterator Graph::begin_vertex(void) const {
		return VertexIterator(adj_list.begin());
	}

	Graph::VertexIterator Graph::end_vertex(void) const {
		return VertexIterator(adj_list.end());
	}

	Graph::VertexIterator::VertexIterator(CQF<KeyObject>::Iterator it) :
		adj_list_itr(it) {};

	Graph::vertex Graph::VertexIterator::operator*(void) const {
		return (*adj_list_itr).key;
	}

	void Graph::VertexIterator::operator++(void) {
		++adj_list_itr;
	}

	bool Graph::VertexIterator::done(void) const {
		return adj_list_itr.done();
	}

}

#endif	// __GRAPH_H_
