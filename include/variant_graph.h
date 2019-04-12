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

			uint64_t get_num_vertices(void) const;
			uint64_t get_seq_length(void) const;
			std::string get_chr(void) const;

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
			enum MUTATION_TYPE {
				INSERTION,
				DELETION,
				SUBSTITUTION
			};

			void update_idx_vertex_id_map(const VariantGraphVertex& v);
			VariantGraphVertex* create_vertex(uint64_t id, uint64_t offset, uint64_t
																				length,
																				const std::vector<VariantGraphVertex::sample_info>&
																				samples);
			VariantGraphVertex::sample_info* create_sample_info(uint64_t index,
																													const std::string
																													sample_id, bool gt1,
																													bool gt2); 
			void add_mutation(const std::string ref, const std::string alt, uint64_t
												pos, std::string sample_id, bool gt1, bool gt2);
			// we only split vertices from the ref.
			// splits the vertex into two. Connects the cur vertex and
			// the new one. sets the vertex_id of the new vertex.
			void split_vertex(uint64_t vertex_id, uint64_t pos, uint64_t*
												new_vertex);
			// we only split vertices from the ref.
			// splits the vertex into three vertices. Connects the cur vertex and
			// the new ones. sets the vertex_id of the new vertices.
			void split_vertex(uint64_t vertex_id, uint64_t pos1, uint64_t pos2,
												uint64_t* new_vertex_1,  uint64_t* new_vertex_2);
			// returns the vertex_id of the new vertex
			VariantGraphVertex* add_vertex(const std::string& seq, uint64_t index,
																		 const std::string& sample_id, bool gt1,
																		 bool gt2);

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
		// Verify that the version of the library that we linked against is
		// compatible with the version of the headers we compiled against.
		GOOGLE_PROTOBUF_VERIFY_VERSION;

		std::string ref;
		read_fasta(ref_file, chr, ref);
		// initialize the seq buffer
		sdsl::util::assign(seq_buffer, sdsl::int_vector<>(ref.size(), 0, 3));

		// add ref node
		VariantGraphVertex *v = add_vertex(ref, 0, "ref", 0, 0);
		
		// update idx->vertex_id map
		update_idx_vertex_id_map(*v);

		// Add vcf files
		add_vcfs(vcfs);
	}

	VariantGraph::~VariantGraph() {
		// Optional:  Delete all global objects allocated by libprotobuf.
		google::protobuf::ShutdownProtobufLibrary();
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

	void VariantGraph::update_idx_vertex_id_map(const VariantGraphVertex& v) {
		for (int i = 0; i < v.s_info_size(); i++) {
			const VariantGraphVertex::sample_info& s = v.s_info(i);
			idx_vertex_id[s.index()] = v.vertex_id();
		}
	}

	VariantGraphVertex* VariantGraph::create_vertex(uint64_t id, uint64_t
																									offset, uint64_t length,
																									const std::vector<VariantGraphVertex::sample_info>&
																									samples) {
		// create vertex object and add to vertex_list
		VariantGraphVertex* v = vertex_list.add_vertex();
		v->set_vertex_id(id);
		v->set_offset(offset);
		v->set_length(length);
		for (const auto sample : samples) {
			VariantGraphVertex::sample_info* s = v->add_s_info();
			s->set_index(sample.index());
			s->set_sample_id(sample.sample_id());
			s->set_gt_1(sample.gt_1());
			s->set_gt_2(sample.gt_2());
		}

		return v;
	}

	void VariantGraph::add_mutation(const std::string ref, const std::string
																	alt, uint64_t pos, std::string sample_id,
																	bool gt1, bool gt2) {
		// find the type mutatuin
		enum MUTATION_TYPE mutation;
		if (ref.size() == alt.size())
			mutation = SUBSTITUTION;
		else if (ref.size() > alt.size()) 
			mutation = DELETION;
		else
			mutation = INSERTION;

		uint64_t ref_vertex_id;
		auto fit = idx_vertex_id.lower_bound(pos);
		if (fit != idx_vertex_id.end()) {
			// get to the vertex which contain @pos
			while (fit->first > pos - 1 && fit != idx_vertex_id.begin())	
				--fit;
			ref_vertex_id = fit->second;
		} else {
			auto rit = idx_vertex_id.rbegin();
			while (rit->first > pos - 1 && rit != idx_vertex_id.rend())
				++rit;
			ref_vertex_id = rit->second;
		}

		VariantGraphVertex v = vertex_list.vertex(ref_vertex_id); 

		// check if need to split
		if (v.s_info(0).index() < pos - 1 && v.s_info(0).index() + v.length() >
				pos - 1) { // need splitting
			// split the ref node
			// for insertion split the node once. 
			// for substitution and deletion split twice.
			uint64_t v1, v2;
			if (mutation == INSERTION)
				split_vertex(ref_vertex_id, pos - v.s_info(0).index() - 1, &v1);
			else
				split_vertex(ref_vertex_id, pos - v.s_info(0).index() - 1, pos +
										 alt.size() - v.s_info(0).index() - 1, &v1, &v2);
		}

		//VariantGraphVertex *v = add_vertex(alt, [>index<], sample_id, gt1, gt2);


	}

	void VariantGraph::split_vertex(uint64_t vertex_id, uint64_t pos, uint64_t*
																	new_vertex) {
		VariantGraphVertex cur_vertex = vertex_list.vertex(vertex_id);

		// create vertex object and add to vertex_list
		uint64_t offset = cur_vertex.offset() + pos - 1;
		uint64_t length = cur_vertex.length() - (pos - cur_vertex.offset()) + 1;
		VariantGraphVertex::sample_info s;
		s.set_index(cur_vertex.s_info(0).index() + pos);
		s.set_sample_id(cur_vertex.s_info(0).sample_id());
		s.set_gt_1(0);
		s.set_gt_2(0);
		std::vector<VariantGraphVertex::sample_info> samples = {s};
		VariantGraphVertex *v = create_vertex(num_vertices, offset, length, samples);

		*new_vertex = v->vertex_id();
		update_idx_vertex_id_map(*v);

		// update length of the seq in the cur_vertex
		cur_vertex.set_length(pos - cur_vertex.offset() - 1);

		// add the edge
		topology.add_edge(vertex_id, *new_vertex);

		// increment vertex count
		num_vertices++;
	}

	void VariantGraph::split_vertex(uint64_t vertex_id, uint64_t pos1, uint64_t
																	pos2, uint64_t* new_vertex_1,  uint64_t*
																	new_vertex_2) {
		uint64_t v1, v2;
		split_vertex(vertex_id, pos1, &v1);
		split_vertex(v1, pos2-pos1, &v2);
	}

	VariantGraphVertex* VariantGraph::add_vertex(const std::string& seq,
																							 uint64_t index, const
																							 std::string& sample_id, bool
																							 gt1, bool gt2) {
		// resize the seq_buffer.
		seq_buffer.resize(seq_buffer.size() + seq.size());
		uint64_t start_offset = seq_length;
		// Add seq to seq_buffer
		for (const auto c : seq) {
			seq_buffer[seq_length] = c;
			seq_length++;
		}
		// create vertex object and add to vertex_list
		VariantGraphVertex::sample_info s;
		s.set_index(index);
		s.set_sample_id(sample_id);
		s.set_gt_1(gt1);
		s.set_gt_2(gt2);
		std::vector<VariantGraphVertex::sample_info> samples = {s};
		VariantGraphVertex *v = create_vertex(num_vertices, start_offset,
																					seq_length, samples);

		// increment vertex count
		num_vertices++;

		return v;
	}

	uint64_t VariantGraph::get_num_vertices(void) const {
		return num_vertices;
	}

	uint64_t VariantGraph::get_seq_length(void) const {
		return seq_length;
	}

	std::string VariantGraph::get_chr(void) const {
		return chr;
	}

}

#endif //__VARIANT_GRAPH_H__
