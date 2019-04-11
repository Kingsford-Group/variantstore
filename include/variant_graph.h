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

#include <sdsl/bit_vectors.hpp>
#include <sdsl/int_vector.hpp>
#include <sdsl/util.hpp>

#include "vcflib/Variant.h"

#include "variantgraphvertex.pb.h"
#include "graph.h"

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
			// can't create a VariantGraph without a reference genome and zero or
			// more vcf files.
			VariantGraph() = delete; 
			// construct variant graph using a reference genome and zero or more vcf
			// files.
			VariantGraph(const std::string ref_file, const std::vector<std::string>&
									 vcfs = std::vector<std::string>());
			// read variant graph from disk
			//VariantGraph(const std::string prefix);

			~VariantGraph();

			// add one or more vcf files to the variant graph
			void add_vcfs(const std::vector<std::string>& vcfs);

			// persist variant graph to disk
			void serialize(const std::string prefix);

			uint64_t get_num_vertices() const;
			uint64_t get_seq_length() const;

			// iterator traversing a specific path in the variant graph
			class VariantGraphPathIterator {
				public:
					VariantGraphPathIterator(Graph::vertex v, const std::string sample_id);
					VariantGraphVertex operator*(void) const;
					void operator++(void);
					bool done(void) const;
		
				private:
					VariantGraphVertex cur;
			};

			std::vector<VariantGraphVertex> out_neighbors(uint64_t vertex_id);
			std::vector<VariantGraphVertex> in_neighbors(uint64_t vertex_id);

			VariantGraphPathIterator find(uint64_t vertex_id, uint64_t sample_id);
			// iterator will be positioned at the start of the path.
			VariantGraphPathIterator find(uint64_t sample_id);


		private:
			void add_mutation(const std::string org, const std::string mut, uint64_t
												pos);
			void split_vertex(uint64_t pos);
			void split_vertex(uint64_t pos1, uint64_t pos2);
			// returns the vertex_id of the new vertex
			uint64_t add_vertex(const std::string& seq, uint64_t index, const
													std::string& sample_id);

			uint64_t seq_length{0};
			uint64_t num_vertices{0};
			std::string chr;
			std::map<uint64_t, uint64_t> idx_vertex_id;

			// structures to persist when serializing variant graph.
			Graph topology;
			VariantGraphVertexList vertex_list;
			sdsl::int_vector<> seq_buffer;
	};

	VariantGraph::VariantGraph(const std::string ref_file, const
														 std::vector<std::string>& vcfs) {
		std::cout << "In VG constructor" << "\n";
		std::string ref;
		read_fasta(ref_file, chr, ref);
		// initialize the seq buffer
		sdsl::util::assign(seq_buffer, sdsl::int_vector<>(ref.size(), 0, 3));

		// add ref node
		add_vertex(ref, 0, "ref");

		// Add vcf files
		add_vcfs(vcfs);
	}

	VariantGraph::~VariantGraph() {
	}

	void VariantGraph::add_vcfs(const std::vector<std::string>& vcfs) {
		for (auto vcf : vcfs) {
			vcflib::VariantCallFile variantFile;
			variantFile.open(vcf);
			vcflib::Variant var(variantFile);

			for (auto sample : variantFile.sampleNames)
				std::cout << sample << " ";
			std::cout << "\n";

			long int count = 0;
			while (variantFile.getNextVariant(var)) {
				count+= 1;
				std::cout << var << "\n";
			}
		}
	}

	uint64_t VariantGraph::add_vertex(const std::string& seq, uint64_t index,
																		const std::string& sample_id) {
		// resize the seq_buffer.
		seq_buffer.resize(seq_buffer.size() + seq.size());
		uint64_t start_offset = seq_length;
		// Add seq to seq_buffer
		for (const auto c : seq) {
			seq_buffer[seq_length] = c;
			seq_length++;
		}
		// create vertex object and add to vertex_list
		//VariantGraphVertex v = {	num_vertices,
			//start_offset,
			//seq_length - start_offset,
			//index,
			//sample_id};
		//vertex_list[v.vertex_id] = v;
		//// add to idx-vertex map
		//idx_vertex_id[index] = v.vertex_id;

		// increment vertex count
		num_vertices++;

		//return v.vertex_id;
		return 0;
	}

	uint64_t VariantGraph::get_num_vertices(void) const {
		return num_vertices;
	}

	uint64_t VariantGraph::get_seq_length(void) const {
		return seq_length;
	}

}

#endif //__VARIANT_GRAPH_H__
