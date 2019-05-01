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

#include <assert.h>

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
			VariantGraph(const std::string& ref_file, const std::vector<std::string>&
									 vcfs);

			// read variant graph from disk
			VariantGraph(const std::string& prefix);

			~VariantGraph();

			// add one or more vcf files to the variant graph
			void add_vcfs(const std::vector<std::string>& vcfs);

			// persist variant graph to disk
			void serialize(const std::string& prefix);

			uint64_t get_num_vertices(void) const;
			uint64_t get_seq_length(void) const;
			const std::string get_chr(void) const;
			uint64_t get_ref_length(void) const;

			void print_vertex_info(const VariantGraphVertex& v) const;
			const std::string get_sequence(const VariantGraphVertex& v) const;

			bool get_sample_from_vertex_if_exists(Graph::vertex v, const std::string
																						sample_id,
																						VariantGraphVertex::sample_info&
																						sample) const;

			// iterator for a breadth-first traversal in the variant graph
			class VariantGraphIterator {
				public:
					VariantGraphIterator(const VariantGraph* g, Graph::vertex v,
															 uint64_t r);
					const VariantGraphVertex* operator*(void) const;
					void operator++(void);
					bool done(void) const;

				private:
					const VariantGraph *vg;
					Graph::GraphIterator itr;
			};

			VariantGraphIterator find(Graph::vertex v = 0, uint64_t radius =
																UINT64_MAX);

			// iterator traversing a specific path in the variant graph
			class VariantGraphPathIterator {
				public:
					VariantGraphPathIterator(const VariantGraph* g, Graph::vertex v,
																	 const std::string sample_id);
					const VariantGraphVertex* operator*(void) const;
					void operator++(void);
					bool done(void) const;

				private:
					const VariantGraph* vg;
					VariantGraphVertex cur;
					std::string s_id;
					bool is_done;
			};

			//std::vector<VariantGraphVertex> out_neighbors(uint64_t vertex_id);
			//std::vector<VariantGraphVertex> in_neighbors(uint64_t vertex_id);

			VariantGraph::VariantGraphPathIterator find(uint64_t vertex_id, const
																									std::string sample_id)
				const;
			// iterator will be positioned at the start of the path.
			VariantGraph::VariantGraphPathIterator find(const std::string sample_id)
				const;

		private:
			enum MUTATION_TYPE {
				INSERTION,
				DELETION,
				SUBSTITUTION
			};

			const std::string mutation_string(MUTATION_TYPE m) const;

			void update_idx_vertex_id_map(const VariantGraphVertex& v);
			uint64_t find_sample_index(Graph::vertex ref_v_id, std::string
																 sample_id) const;
			bool get_neighbor_vertex(Graph::vertex id, const std::string sample_id,
															 Graph::vertex* v) const;
			void add_sample_to_vertex(Graph::vertex id, uint64_t sample_idx, const
																std::string sample_id, bool gt1, bool gt2);
			bool check_if_mutation_exists(Graph::vertex prev, Graph::vertex next,
																		const std::string alt, Graph::vertex* v)
				const;
			bool check_if_mutation_exists(Graph::vertex prev, uint64_t offset,
																		uint64_t length, Graph::vertex* v) const;
			VariantGraphVertex* create_vertex(uint64_t id, uint64_t offset, uint64_t
																				length,
																				const std::vector<VariantGraphVertex::sample_info>&
																				samples);
			VariantGraphVertex::sample_info* create_sample_info(uint64_t index,
																													const std::string
																													sample_id, bool gt1,
																													bool gt2); 
			void add_mutation(std::string ref, std::string alt, uint64_t pos,
												std::string sample_id, bool gt1, bool gt2);
			// we only split vertices from the ref.
			// splits the vertex into two. Connects the cur vertex and
			// the new one. sets the vertex_id of the new vertex.
			void split_vertex(uint64_t vertex_id, uint64_t pos, Graph::vertex*
												new_vertex);
			// we only split vertices from the ref.
			// splits the vertex into three vertices. Connects the cur vertex and
			// the new ones. sets the vertex_id of the new vertices.
			void split_vertex(uint64_t vertex_id, uint64_t pos1, uint64_t pos2,
												Graph::vertex* new_vertex_1,  Graph::vertex*
												new_vertex_2);
			// returns the vertex_id of the new vertex
			VariantGraphVertex* add_vertex(uint64_t offset, uint64_t length,
																		 uint64_t index, const std::string&
																		 sample_id, bool gt1, bool gt2);
			// returns the vertex_id of the new vertex
			VariantGraphVertex* add_vertex(const std::string& seq, uint64_t index,
																		 const std::string& sample_id, bool gt1,
																		 bool gt2);

			uint64_t seq_length{0};
			uint64_t num_vertices{0};
			std::string chr;
			uint64_t ref_length;
			std::map<uint64_t, uint64_t> idx_vertex_id;

			// structures to persist when serializing variant graph.
			VariantGraphVertexList vertex_list;
			sdsl::int_vector<> seq_buffer;
			Graph topology;
	};

	VariantGraph::VariantGraph(const std::string& ref_file, const
														 std::vector<std::string>& vcfs) {
		// Verify that the version of the library that we linked against is
		// compatible with the version of the headers we compiled against.
		GOOGLE_PROTOBUF_VERIFY_VERSION;

		std::string ref;
		read_fasta(ref_file, chr, ref);
		ref_length = ref.size();
		// initialize the seq buffer
		sdsl::util::assign(seq_buffer, sdsl::int_vector<>(ref.size(), 0, 3));

		// add ref node
		// we set the index to 1.
		VariantGraphVertex *v = add_vertex(ref, 1, "ref", 0, 0);

		// update idx->vertex_id map
		update_idx_vertex_id_map(*v);

		// Add vcf files
		add_vcfs(vcfs);
	}

	VariantGraph::VariantGraph(const std::string& prefix) : topology(prefix) {
		// load vertex list
		std::string vertex_list_name = prefix + "/vertex_list.proto";
		fstream input(vertex_list_name, ios::in | ios::binary);
		if (!vertex_list.ParseFromIstream(&input)) {
			ERROR("Failed to parse vertex list.");
			abort();
		}
		// load seq buffer
		std::string seq_buffer_name = prefix + "/seq_buffer.sdsl";
		sdsl::load_from_file(seq_buffer, seq_buffer_name);
		// load topology
		//topology = Graph(prefix);	
		num_vertices = topology.get_num_vertices() + 1;
		seq_length = seq_buffer.size();
	}

	VariantGraph::~VariantGraph() {
		// Optional:  Delete all global objects allocated by libprotobuf.
		google::protobuf::ShutdownProtobufLibrary();
	}

	void VariantGraph::serialize(const std::string& prefix) {
		// serialize vertex list
		std::string vertex_list_name = prefix + "/vertex_list.proto";
		std::fstream output(vertex_list_name, ios::out | ios::trunc | ios::binary);
		if (!vertex_list.SerializeToOstream(&output)) {
			ERROR("Failed to write vertex list.");
			abort();
		}	
		// serialize seq buffer
		std::string seq_buffer_name = prefix + "/seq_buffer.sdsl";
		seq_buffer.resize(seq_buffer.size());
		sdsl::store_to_file(seq_buffer, seq_buffer_name);
		// serialize topology
		topology.serialize(prefix);
	}

	void VariantGraph::add_vcfs(const std::vector<std::string>& vcfs) {
		for (auto vcf : vcfs) {
			vcflib::VariantCallFile variantFile;
			variantFile.open(vcf);
			vcflib::Variant var(variantFile);

			// substitutions are always only one base at a time.
			// Insertions/deletions can be multiple bases.
			long int num_mutations = 0;
			while (variantFile.getNextVariant(var)) {
				num_mutations += 1;
				for (const auto alt : var.alt) {
					for (const auto sample : var.samples) {
						auto gt = sample.second.find("GT");
						auto gt_vec = *gt;
						std::string gt_info = gt_vec.second[0];
						//assert(gt_info.size() == 3);

						bool gt1, gt2;
						bool add = false;
						if (gt_info.size() == 3) {
							// extract gt info
							const char *str = gt_info.c_str();
							int first = str[0] - '0';
							char phase = str[1];
							int second = str[2] - '0';
							if (phase == '|') {
								if (first > 0 || second > 0) {
									if (first > 0)
										gt1 = true;
									else 
										gt1 = false;
									if (second > 0)
										gt2 = true;
									else 
										gt2 = false;

									add = true;
								}	
							} else if (phase == '/') {
								if (first > 0 && second > 0) {
									gt1 = true; gt2 = true;
									add = true;
								} else if (first == '1' || second == '1') {
									gt1 = false; gt2 = false;
									add = true;
								}
							}
						} else if (gt_info.size() == 1) {
							int present = stoi(gt_info);
							if (present) {
								add = true;
								gt1 = true;
								gt2 = false;
							}
						}
						if (add)
							add_mutation(var.ref, alt, var.position, sample.first, gt1, gt2);
					}
				}
			}
			PRINT("Num mutations: " << num_mutations);
		}
	}

	VariantGraphVertex* VariantGraph::add_vertex(uint64_t offset, uint64_t length,
																							 uint64_t index, const
																							 std::string& sample_id, bool
																							 gt1, bool gt2) {
		// create vertex object and add to vertex_list
		VariantGraphVertex::sample_info s;
		s.set_index(index);
		s.set_sample_id(sample_id);
		s.set_gt_1(gt1);
		s.set_gt_2(gt2);
		std::vector<VariantGraphVertex::sample_info> samples = {s};
		VariantGraphVertex *v = create_vertex(num_vertices, offset, length,
																					samples);

		// increment vertex count
		num_vertices++;

		return v;
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
			seq_buffer[seq_length] = map_base(c);
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
																					seq.size(), samples);

		// increment vertex count
		num_vertices++;

		return v;
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

	void VariantGraph::split_vertex(uint64_t vertex_id, uint64_t pos,
																	Graph::vertex* new_vertex) {
		DEBUG("Spliting vertex: " << vertex_id << " " << pos);
		const VariantGraphVertex cur_vertex = vertex_list.vertex(vertex_id);

		// create vertex object and add to vertex_list
		uint64_t offset = cur_vertex.offset() + pos - 1;
		uint64_t length = cur_vertex.length() - pos  + 1;
		VariantGraphVertex::sample_info s;
		s.set_index(cur_vertex.s_info(0).index() + pos - 1);
		s.set_sample_id(cur_vertex.s_info(0).sample_id());
		s.set_gt_1(cur_vertex.s_info(0).gt_1());
		s.set_gt_2(cur_vertex.s_info(0).gt_2());
		std::vector<VariantGraphVertex::sample_info> samples = {s};
		VariantGraphVertex *v = create_vertex(num_vertices, offset, length, samples);

		*new_vertex = v->vertex_id();
		update_idx_vertex_id_map(*v);

		// update length of the seq in the cur_vertex
		vertex_list.mutable_vertex(vertex_id)->set_length(cur_vertex.length() -
																											length);

		// move outgoing connections from old_node to the new_node
		for (const auto v : topology.out_neighbors(vertex_id)) {
			topology.add_edge(*new_vertex, v);
			topology.remove_edge(vertex_id, v);
		}

		// add the edge
		topology.add_edge(vertex_id, *new_vertex);

		// increment vertex count
		num_vertices++;
	}

	void VariantGraph::split_vertex(uint64_t vertex_id, uint64_t pos1, uint64_t
																	pos2, Graph::vertex* new_vertex_1,
																	Graph::vertex* new_vertex_2) {
		split_vertex(vertex_id, pos1, new_vertex_1);
		split_vertex(*new_vertex_1, pos2 - pos1 + 1, new_vertex_2);
	}

	uint64_t VariantGraph::get_num_vertices(void) const {
		return num_vertices;
	}

	uint64_t VariantGraph::get_seq_length(void) const {
		return seq_length;
	}

	const std::string VariantGraph::get_chr(void) const {
		return chr;
	}

	uint64_t VariantGraph::get_ref_length(void) const {
		return ref_length;
	}

	const std::string VariantGraph::mutation_string(MUTATION_TYPE m) const {
		if (m == INSERTION)
			return std::string("INSERTION");
		else if (m == DELETION)
			return std::string("DELETION");
		else
			return std::string("SUBSTITUTION");
	}

	void VariantGraph::print_vertex_info(const VariantGraphVertex& v) const {
		std::string samples;
		for (int i = 0; i < v.s_info_size(); i++) {
			const VariantGraphVertex::sample_info& s = v.s_info(i);
			samples.append("Sample id: " + s.sample_id());
			samples.append(" Index: " +  std::to_string((int)s.index()) + " ");
		}
		PRINT("ID: " << v.vertex_id() << " Offset: " << v.offset() <<
					" length: " << v.length() << " Samples: " << samples); 
	}

	const std::string VariantGraph::get_sequence(const VariantGraphVertex& v)
		const {
			std::string seq;
			for (uint64_t i = v.offset(); i < v.offset() + v.length(); i++) {
				seq += map_int((uint8_t)seq_buffer[i]);
			}
			return seq;
		}

	// the idx map should only contain information about the ref vertices.
	void VariantGraph::update_idx_vertex_id_map(const VariantGraphVertex& v) {
		for (int i = 0; i < v.s_info_size(); i++) {
			const VariantGraphVertex::sample_info& s = v.s_info(i);
			if (s.sample_id() == "ref")
				idx_vertex_id[s.index()] = v.vertex_id();
		}
	}

	bool VariantGraph::get_sample_from_vertex_if_exists(Graph::vertex v, const
																											std::string sample_id,
																											VariantGraphVertex::sample_info&
																											sample) const {
		VariantGraphVertex cur_vertex = vertex_list.vertex(v);
		for (int i = 0; i < cur_vertex.s_info_size(); i++) {
			const VariantGraphVertex::sample_info& s = cur_vertex.s_info(i);
			if (s.sample_id() == sample_id) {
				sample = s;
				return true;
			}
		}
		return false;
	}

	uint64_t VariantGraph::find_sample_index(Graph::vertex ref_v_id,
																					 const std::string sample_id) const {
		Graph::vertex cur_vertex_id = ref_v_id;
		VariantGraphVertex cur_vertex = vertex_list.vertex(cur_vertex_id);
		// check if the cur vertex contains the sample
		VariantGraphVertex::sample_info sample;
		if (get_sample_from_vertex_if_exists(cur_vertex_id, sample_id, sample))
			return sample.index() + cur_vertex.length();

		// else traverse back in the graph to find a vertex with @sample_id
		get_sample_from_vertex_if_exists(ref_v_id, "ref", sample);
		uint64_t cur_index = sample.index();
		auto temp_itr = idx_vertex_id.lower_bound(cur_index);
		--temp_itr;
		while (cur_index > 1) {
			// move to prev vertex.
			cur_index = temp_itr->first; 
			cur_vertex_id = temp_itr->second;
			cur_vertex = vertex_list.vertex(cur_vertex_id);

			if (get_sample_from_vertex_if_exists(cur_vertex_id, sample_id, sample))
				break;

			Graph::vertex neighbor_vertex_id;
			if (get_neighbor_vertex(cur_vertex_id, sample_id, &neighbor_vertex_id)) {
				cur_vertex_id = neighbor_vertex_id;
				cur_vertex = vertex_list.vertex(neighbor_vertex_id);
				break;
			}
			--temp_itr;
		}

		// init cur distance
		uint64_t cur_distance = 0;
		if (get_sample_from_vertex_if_exists(cur_vertex_id, sample_id, sample)) {
			cur_distance += sample.index();
		} else if (get_sample_from_vertex_if_exists(cur_vertex_id, "ref", sample)) {
			cur_distance += sample.index();
		}
		// traverse the graph forward till you find the @ref_v_id
		auto it = this->find(cur_vertex_id, sample_id);
		while (!it.done() && (*it)->vertex_id() != ref_v_id) {
			cur_distance += (*it)->length();
			++it;
		}

		return cur_distance += vertex_list.vertex(ref_v_id).length();
	}

	// if there is not neighbor with @sample_id then set @v to "ref"
	bool VariantGraph::get_neighbor_vertex(Graph::vertex id, const std::string
																				 sample_id, Graph::vertex* v)
		const {
			std::vector<Graph::vertex> neighbors;
			for (const auto v_id : topology.out_neighbors(id)) {
				VariantGraphVertex vertex = vertex_list.vertex(v_id);
				for (int i = 0; i < vertex.s_info_size(); i++) {
					const VariantGraphVertex::sample_info& s = vertex.s_info(i);
					if (sample_id == "ref" && s.sample_id() == sample_id) { // handle special case for "ref" path.
						neighbors.emplace_back(v_id);
					} else {
						if (s.sample_id() == sample_id) {
							*v = v_id; 
							return true;
						} else if (s.sample_id() == "ref") {	// if sample_id is not found follow "ref" path
							*v = v_id;
						}
					}
				}
			}
			if (sample_id == "ref" && neighbors.size() > 0) {
				//if (neighbors.size() == 0) {
				//ERROR("Can't find the vertex with sample id: " <<  sample_id <<
				//" at vertex id: " << id);
				//abort();
				//}
				if (neighbors.size() > 2) {
					ERROR("More than two neighbors in ref path at vertex id: " << id);
					abort();
				}
				if (neighbors.size() == 1) {
					*v = neighbors[0];
					return true;
				}
				if (topology.is_edge(std::make_pair(neighbors[0], neighbors[1]))) {
					*v = neighbors[0];
					return true;
				} else if (topology.is_edge(std::make_pair(neighbors[1],
																									 neighbors[0]))) {
					*v = neighbors[1];
					return true;
				} else {
					ERROR("Two neighbors are not connected in ref path at vertex id: "
								<< id);
					abort();
				}
			}

			return false;
		}

	void VariantGraph::add_sample_to_vertex(Graph::vertex id, uint64_t
																					sample_idx, const std::string
																					sample_id, bool gt1, bool gt2) {
		VariantGraphVertex::sample_info* s =
			vertex_list.mutable_vertex(id)->add_s_info();
		s->set_index(sample_idx);
		s->set_sample_id(sample_id);
		s->set_gt_1(gt1);
		s->set_gt_2(gt2);
	}

	bool VariantGraph::check_if_mutation_exists(Graph::vertex prev,
																							Graph::vertex next, const
																							std::string alt,
																							Graph::vertex* v) const {
		for (const auto v_id : topology.out_neighbors(prev)) {
			if (topology.is_edge(std::make_pair(v_id, next))) {
				const VariantGraphVertex neighbor = vertex_list.vertex(v_id);
				const std::string seq = get_sequence(neighbor);
				if (seq == alt) {
					*v = v_id;
					return true;
				}
			}
		}
		return false;
	}

	bool VariantGraph::check_if_mutation_exists(Graph::vertex prev,
																							uint64_t offset, uint64_t length,
																							Graph::vertex* v) const {
		for (const auto v_id : topology.out_neighbors(prev)) {
			const VariantGraphVertex neighbor = vertex_list.vertex(v_id);
			if (neighbor.offset() == offset && neighbor.length() == length) {
				*v = v_id;
				return true;
			}
		}
		return false;
	}

	void VariantGraph::add_mutation(std::string ref, std::string alt, uint64_t
																	pos, std::string sample_id,
																	bool gt1, bool gt2) {
		// find the type mutatuin
		enum MUTATION_TYPE mutation;
		if (ref.size() == alt.size())
			mutation = SUBSTITUTION;
		else if (ref.size() > alt.size()) 
			mutation = DELETION;
		else
			mutation = INSERTION;

		DEBUG("Adding mutation: " << mutation_string(mutation) << " " << ref <<
					" " << alt << " " << pos << " " << sample_id);
		// update pos and alt/ref if it's an insertion/deletion.
		if (mutation == INSERTION) {
			pos = pos + ref.size();
			alt = alt.substr(ref.size());
		} else if (mutation == DELETION) {
			pos = pos + alt.size();
			ref = ref.substr(alt.size());
		}

		// find the vertex in the ref corresponding to @pos
		auto ref_idx_itr = idx_vertex_id.lower_bound(pos);
		if (ref_idx_itr->first != pos) {	// the iterator is positioned at an index greater than @pos
			--ref_idx_itr;	// move to the prev node
			if (ref_idx_itr->first >= pos) {
				std::cerr << "The prev ref vertex has an index greater than pos." <<
					"\n";
				abort();
			}
		}
		uint64_t ref_vertex_idx = ref_idx_itr->first;
		Graph::vertex ref_vertex_id = ref_idx_itr->second;
		VariantGraphVertex ref_vertex = vertex_list.vertex(ref_vertex_id);
		// check if the mapping is consistant with the vertex structure.
		//assert(ref_vertex_idx == get_sample_from_vertex(ref_vertex_id,
																										//"ref").index());

		if (mutation == SUBSTITUTION) {
			Graph::vertex prev_ref_vertex_id = 0, next_ref_vertex_id = 0;
			// Either we got the vertex where there's already a substition 
			// or the vertex contains the seq with @pos
			if (ref_vertex_idx == pos && ref_vertex.length() == ref.size()) { // splitting not needed
				// find the prev ref vertex
				auto temp_itr = idx_vertex_id.lower_bound(ref_vertex_idx);
				if (temp_itr->first != ref_vertex_idx) {
					ERROR("Vertex id not found in the map");
					abort();
				}
				--temp_itr;
				prev_ref_vertex_id = temp_itr->second;

				// find the next ref vertex
				get_neighbor_vertex(ref_vertex_id, "ref", &next_ref_vertex_id);		
			} else if (ref_vertex_idx == pos && ref_vertex.length() > ref.size()) { // vertex needs to be split into two
				// find the prev ref vertex
				auto temp_itr = idx_vertex_id.lower_bound(ref_vertex_idx);
				if (temp_itr->first != ref_vertex_idx) {
					ERROR("Vertex id not found in the map");
					abort();
				}
				--temp_itr;
				prev_ref_vertex_id = temp_itr->second;

				// split the vertex
				split_vertex(ref_vertex_id, ref.size() + 1, &next_ref_vertex_id);
			} else if (ref_vertex_idx < pos && ref_vertex.length() > ref.size()) { // vertex needs to be split into three
				// split the vertex
				uint64_t split_pos = pos - ref_vertex.offset();
				split_vertex(ref_vertex_id, split_pos, split_pos + ref.size(),
										 &prev_ref_vertex_id, &next_ref_vertex_id);
				std::swap(ref_vertex_id, prev_ref_vertex_id);
			}
			// find the index of the sample at this position in the graph.
			uint64_t sample_idx = find_sample_index(prev_ref_vertex_id, sample_id);
			Graph::vertex exist_vertex = 0;
			if (check_if_mutation_exists(prev_ref_vertex_id, next_ref_vertex_id, alt,
																	 &exist_vertex)) {
				add_sample_to_vertex(exist_vertex, sample_idx, sample_id, gt1, gt2);
			} else {
				VariantGraphVertex* sample_vertex = add_vertex(alt, sample_idx,
																											 sample_id, gt1, gt2);

				// make connections for the new vertex in the graph
				topology.add_edge(prev_ref_vertex_id, sample_vertex->vertex_id());
				topology.add_edge(sample_vertex->vertex_id(), next_ref_vertex_id);
			}
		} else if (mutation == INSERTION) {
			Graph::vertex prev_ref_vertex_id = 0, next_ref_vertex_id = 0;
			// Either we got the vertex where there's already an insertion
			// or the vertex contains the seq with @pos
			if (ref_vertex_idx == pos && ref_vertex.length() == 1) { // splitting not needed
				// find the prev ref vertex
				auto temp_itr = idx_vertex_id.lower_bound(ref_vertex_idx);
				if (temp_itr->first != ref_vertex_idx) {
					ERROR("Vertex id not found in the map");
					abort();
				}
				--temp_itr;
				prev_ref_vertex_id = temp_itr->second;

				// find the next ref vertex
				get_neighbor_vertex(ref_vertex_id, "ref", &next_ref_vertex_id);		
			} else if (ref_vertex_idx <= pos) {	// vertex needs to be split into two
				// split the vertex
				split_vertex(ref_vertex_id, pos - ref_vertex_idx + 1,
										 &next_ref_vertex_id);
				prev_ref_vertex_id = ref_vertex_id;
			}
			// find the index of the sample at this position in the graph.
			uint64_t sample_idx = find_sample_index(prev_ref_vertex_id, sample_id);
			Graph::vertex exist_vertex = 0;
			if (check_if_mutation_exists(prev_ref_vertex_id, next_ref_vertex_id, alt,
																	 &exist_vertex)) {
				add_sample_to_vertex(exist_vertex, sample_idx, sample_id, gt1, gt2);
			} else {
				VariantGraphVertex* sample_vertex = add_vertex(alt, sample_idx,
																											 sample_id, gt1, gt2);

				// make connections for the new vertex in the graph
				topology.add_edge(prev_ref_vertex_id, sample_vertex->vertex_id());
				topology.add_edge(sample_vertex->vertex_id(), next_ref_vertex_id);
			}
		} else { // it's deletion
			Graph::vertex prev_ref_vertex_id = 0, next_ref_vertex_id = 0;
			// Either we got the vertex where there's already a deletion
			// or the vertex contains the seq with @pos
			if (ref_vertex_idx == pos && ref_vertex.length() == ref.size()) { // splitting not needed
				// find the prev ref vertex
				auto temp_itr = idx_vertex_id.lower_bound(ref_vertex_idx);
				if (temp_itr->first != ref_vertex_idx) {
					ERROR("Vertex id not found in the map");
					abort();
				}
				--temp_itr;
				prev_ref_vertex_id = temp_itr->second;

				// find the next ref vertex
				get_neighbor_vertex(ref_vertex_id, "ref", &next_ref_vertex_id);	
			} else if (ref_vertex_idx == pos && ref_vertex.length() > ref.size()) { // vertex needs to be split into two
				// find the prev ref vertex
				auto temp_itr = idx_vertex_id.lower_bound(ref_vertex_idx);
				if (temp_itr->first != ref_vertex_idx) {
					ERROR("Vertex id not found in the map");
					abort();
				}
				--temp_itr;
				prev_ref_vertex_id = temp_itr->second;

				// split the vertex
				split_vertex(ref_vertex_id, ref.size() + 1, &next_ref_vertex_id);
			} else if (ref_vertex_idx < pos && ref_vertex.length() > ref.size()) { // vertex needs to be split into three
				// split the vertex
				uint64_t split_pos = pos - ref_vertex.offset();
				split_vertex(ref_vertex_id, split_pos, split_pos + ref.size(),
										 &prev_ref_vertex_id, &next_ref_vertex_id);
				std::swap(ref_vertex_id, prev_ref_vertex_id);
			}
			// find the index of the sample at this position in the graph.
			uint64_t sample_idx = find_sample_index(prev_ref_vertex_id, sample_id);
			add_sample_to_vertex(next_ref_vertex_id, sample_idx, sample_id, gt1,
													 gt2);
			// make connections for the new vertex in the graph
			topology.add_edge(prev_ref_vertex_id, next_ref_vertex_id);
		}
	}

	VariantGraph::VariantGraphPathIterator::VariantGraphPathIterator(const
																																	 VariantGraph*
																																	 g,
																																	 Graph::vertex
																																	 v, const
																																	 std::string
																																	 sample_id)
	{
		vg = g; 
		cur = vg->vertex_list.vertex(v);
		s_id = sample_id;
		is_done = false;
	}

	const VariantGraphVertex*
		VariantGraph::VariantGraphPathIterator::operator*(void) const {
			return &cur;
		}

	void VariantGraph::VariantGraphPathIterator::operator++(void) {
		Graph::vertex next_vertex = 0; // there shoould be no incoming edge to 0
		if (!vg->get_neighbor_vertex(cur.vertex_id(), s_id, &next_vertex) &&
				next_vertex == 0) {
			is_done = true;
		}
		cur = vg->vertex_list.vertex(next_vertex);
	}

	bool VariantGraph::VariantGraphPathIterator::done(void) const {
		return is_done;
	}

	VariantGraph::VariantGraphPathIterator VariantGraph::find(uint64_t
																														vertex_id, const
																														std::string
																														sample_id) const {
		return VariantGraphPathIterator(this, vertex_id, sample_id);	
	}

	// iterator will be positioned at the start of the path.
	VariantGraph::VariantGraphPathIterator VariantGraph::find(const std::string
																														sample_id) const {
		return VariantGraphPathIterator(this, 0, sample_id);	
	}

	VariantGraph::VariantGraphIterator::VariantGraphIterator(const VariantGraph*
																													 g, Graph::vertex v,
																													 uint64_t r) :
		vg(g), itr(&g->topology, v, r) {};

	const VariantGraphVertex*
		VariantGraph::VariantGraphIterator::operator*(void) const {
			return &vg->vertex_list.vertex(*itr);
	}

	void VariantGraph::VariantGraphIterator::operator++(void) {
		++itr;
	}

	bool VariantGraph::VariantGraphIterator::done(void) const {
		return itr.done();
	}
	
	VariantGraph::VariantGraphIterator VariantGraph::find(Graph::vertex v,
																												uint64_t radius) {
		return VariantGraphIterator(this, v, radius);
	}

}

#endif //__VARIANT_GRAPH_H__
