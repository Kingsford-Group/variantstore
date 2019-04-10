/*
 * ============================================================================
 *
 *       Filename:  variant_graph.h
 *
 *         Author:  Prashant Pandey (), ppandey2@cs.cmu.edu
 *   Organization:  Carnegie Mellon University
 *
 * ============================================================================
 */

#ifndef __VARIANT_GRAPH_H__
#define __VARIANT_GRAPH_H__

#include <iostream>
#include <string>
#include <vector>
#include <map>

namespace variantdb {

	// Construction:
	// Create a variant graph based on a reference genome.
	// Add vcf files during construction or iteratively later on.
	// The constructor takes a reference genome and zero or more vcf files.
	// All vcf files should be based on the reference genome specified during
	// the construction. Otherwise, it will throw an error.
	// The graph can be serialized to location on disk and deserialized from a
	// location on disk.
	//
	// Comments on the construction process:
	// During constructions we need to support predecessor queries to map a
	// given index in the reference genome to the node_id. This structure is
	// also a dynamic position-based index.
	//
	// Querying:
	// Given a node Id it supports adjacency queries to traverse the graph.
	// Given a specific tuple <node_id, sample_id>, returns an iterator
	// which traverses the path taken by <sample_id> in the variant graph.
	//
	// Indexing:
	// Given a variant graph object a position-based is built which is a static
	// data structure.
	//
	// Schema:
	// Node list schema: <Vertex_ID, Offset, Length, Index, Sample_ID>
	// Sequence buffer
	// Graph topology structure

	class VariantGraph {
		public:
			// construct variant graph using a reference genome and zero or more vcf
			// files.
			VariantGraph(const std::string ref, const std::vector<std::string> vcfs);
			// read variant graph from disk
			VariantGraph(const std::string prefix);

			~VariantGraph();

			// add one or more vcf files to the variant graph
			void add_vcfs(const std::vector<std::string> vcfs);

			// persist variant graph to disk
			void serialize(const std::string prefix);

			// structure of vertex in variant graph
			// currently this is a naive structure.
			// TODO potential scope for space optimization.
			typedef struct VariantGraphVertex {
				uint64_t vertex_id;
				uint64_t offset;
				uint64_t length;
				uint64_t index;
				uint64_t sample_id;
			} VariantGraphVertex;

			// iterator traversing a specific path in the variant graph
			class VariantGraphPathIterator {
				public:
					VariantGraphIterator(VertexIterator vitr);
					vertex operator*(void) const;
					void operator++(void);
					bool done(void) const;

				private:
					VertexIterator vitr;
			};

			std::vector<VariantGraphVertex> out_neighbors(uint64_t vertex_id);
			std::vector<VariantGraphVertex> in_neighbors(uint64_t vertex_id);

			VariantGraphIterator find(uint64_t vertex_id, uint64_t sample_id);
			// iterator will be positioned at the start of the path.
			VariantGraphIterator find(uint64_t sample_id);

			
		private:
			void add_mutation(const std::string org, const std::string mut, uint64_t
												pos);
			void split_vertex(uint64_t pos);
			void split_vertex(uint64_t pos1, uint64_t pos2);
			void add_vertex(const std::string seq, uint64_t index, uint64_t sample_id);

			std::map<uint64_t, uint64_t> idx_vertex_id;

			// structures to persist when serializing variant graph.
			std::vector<uint64_t> sample_id_list;
			Graph topology;
			std::map<uint64_t, VariantGraphVertex> vertex_list;
			sdsl::int_vector<> seq_buffer;
	};
}

#endif //__VARIANT_GRAPH_H__
