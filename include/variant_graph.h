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
#include <unordered_map>
#include <regex>

#include <assert.h>

#include <sdsl/bit_vectors.hpp>
#include <sdsl/int_vector.hpp>
#include <sdsl/util.hpp>

#include "vcflib/Variant.h"
#include "lru/lru.hpp"
#include "stream.hpp"

#include "variantgraphvertex.pb.h"
#include "graph.h"

namespace variantdb {

#define CACHE_SIZE 256

	// to map a sample --> vertex ids.
	using Cache = LRU::Cache<uint32_t, Graph::vertex>;

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

	typedef struct {
		uint32_t sample_id;
		bool gt1;
		bool gt2;
	} sample_struct;

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
			uint64_t get_num_edges(void) const;
			uint64_t get_seq_length(void) const;
			const std::string get_chr(void) const;
			uint64_t get_ref_length(void) const;
			std::string get_sample_name(uint32_t id) const;
			double get_cache_hit_rate(void) const;

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
					const VariantGraph *vg;
					VariantGraphVertex cur;
					uint32_t s_id;
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
			bool get_sample_from_vertex_if_exists(Graph::vertex v, uint32_t sample_id,
																						VariantGraphVertex::sample_info&
																						sample) const;
			const std::string mutation_string(MUTATION_TYPE m) const;

			void update_idx_vertex_id_map(const VariantGraphVertex& v);
			uint64_t find_sample_index(Graph::vertex ref_v_id, uint32_t sample_id)
				const;
			bool get_neighbor_vertex(Graph::vertex id, uint32_t sample_id,
															 Graph::vertex* v) const;
			void add_sample_to_vertex(Graph::vertex id, uint64_t sample_idx, uint32_t
																sample_id, bool gt1, bool gt2);
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
												std::vector<sample_struct>& sample_list);
			void fix_sample_indexes(void);

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
																		 uint64_t index, uint32_t sample_id, bool
																		 gt1, bool gt2);
			// returns the vertex_id of the new vertex
			VariantGraphVertex* add_vertex(const std::string& seq, uint64_t index,
																		 uint32_t sample_id, bool gt1, bool gt2);

			// iterator for a breadth-first traversal in the variant graph
			class VariantGraphMutableIterator {
				public:
					VariantGraphMutableIterator(VariantGraph* g, Graph::vertex v,
																			uint64_t r);
					VariantGraphVertex* operator*(void);
					void operator++(void);
					bool done(void) const;

				private:
					VariantGraph *vg;
					Graph::GraphIterator itr;
			};

			VariantGraphMutableIterator mutable_find(Graph::vertex v = 0, uint64_t
																							 radius = UINT64_MAX);

			class VariantGraphPathMutableIterator {
				public:
					VariantGraphPathMutableIterator(VariantGraph* g, Graph::vertex v,
																					const std::string sample_id);
					VariantGraphVertex* operator*(void);
					void operator++(void);
					bool done(void) const;

				private:
					VariantGraph *vg;
					VariantGraphVertex *cur;
					uint32_t s_id;
					bool is_done;
			};
			
			// iterator will be positioned at the start of the path.
			VariantGraph::VariantGraphPathMutableIterator mutable_find(const
																																 std::string
																																 sample_id);

			uint64_t seq_length{0};
			uint64_t num_vertices{0};
			std::string chr;
			uint64_t ref_length;
			std::map<uint64_t, uint64_t> idx_vertex_id;
			std::unordered_map<uint32_t, std::string> idsample_map;
			Cache cache;

			// structures to persist when serializing variant graph.
			VariantGraphVertexList vertex_list;
			sdsl::int_vector<> seq_buffer;
			Graph topology;
			std::unordered_map<std::string, uint32_t> sampleid_map;
	};

	VariantGraph::VariantGraph(const std::string& ref_file, const
														 std::vector<std::string>& vcfs) :
		cache(CACHE_SIZE) {
			// Verify that the version of the library that we linked against is
			// compatible with the version of the headers we compiled against.
			GOOGLE_PROTOBUF_VERIFY_VERSION;

#ifdef DEBUG_MODE
			cache.monitor();
#endif

			std::string ref;
			read_fasta(ref_file, chr, ref);
			ref_length = ref.size();
			// initialize the seq buffer
			sdsl::util::assign(seq_buffer, sdsl::int_vector<>(0, 0, 3));

			// add ref node
			// we set the index to 1.
			sampleid_map.insert(std::make_pair("ref", sampleid_map.size()));
			idsample_map.insert(std::make_pair(0, "ref"));
			// ref id is 0
			VariantGraphVertex *v = add_vertex(ref, 1, 0, 0, 0);

			// update idx->vertex_id map
			update_idx_vertex_id_map(*v);

			// Add vcf files
			add_vcfs(vcfs);

			// fix sample indexes
			fix_sample_indexes();
		}

	VariantGraph::VariantGraph(const std::string& prefix) : topology(prefix) {
		// load vertex list
		std::string vertex_list_name = prefix + "/vertex_list.proto";
		fstream input(vertex_list_name, ios::in | ios::binary);

		function<void(VariantGraphVertex&)> lambda = [this](VariantGraphVertex& v)
		{
			VariantGraphVertex* vertex = vertex_list.add_vertex();
			*vertex = v;
    };

		if (!stream::for_each(input, lambda)) {
			ERROR("Failed to parse vertex list.");
			abort();
		}

		// load seq buffer
		std::string seq_buffer_name = prefix + "/seq_buffer.sdsl";
		sdsl::load_from_file(seq_buffer, seq_buffer_name);
		num_vertices = topology.get_num_vertices() + 1;
		seq_length = seq_buffer.size();

		//load sampleid map
		std::string sampleid_map_name = prefix + "/sampleid_map.lst";
		std::ifstream sampleid_file(sampleid_map_name);
		std::string sample;
		uint32_t id;
		while (sampleid_file >> sample >> id)  {
			sampleid_map.insert(std::make_pair(sample, id));
			idsample_map.insert(std::make_pair(id, sample));
		}

		sampleid_file.close();
	}

	VariantGraph::~VariantGraph() {
		// Optional:  Delete all global objects allocated by libprotobuf.
		google::protobuf::ShutdownProtobufLibrary();
	}

	void VariantGraph::serialize(const std::string& prefix) {
		// serialize vertex list
		std::string vertex_list_name = prefix + "/vertex_list.proto";
		std::fstream output(vertex_list_name, ios::out | ios::trunc | ios::binary);

		std::function<VariantGraphVertex(uint64_t)> lambda =
			[this](uint64_t n) {
			return vertex_list.vertex(n);
		};
	
		if (!stream::write(output, vertex_list.vertex_size(), lambda)) {
			ERROR("Failed to write vertex list.");
			abort();
		}	

		// serialize seq buffer
		std::string seq_buffer_name = prefix + "/seq_buffer.sdsl";
		seq_buffer.resize(seq_buffer.size());
		sdsl::util::bit_compress(seq_buffer);
		sdsl::store_to_file(seq_buffer, seq_buffer_name);
		// serialize topology
		topology.serialize(prefix);

		// serialize sampleid_map
		std::string sampleid_map_name = prefix + "/sampleid_map.lst";
		std::ofstream sampleid_file(sampleid_map_name);
		for (const auto sample : sampleid_map)
			sampleid_file << sample.first << " " << sample.second << "\n";
		sampleid_file.close();
	}

	void VariantGraph::add_vcfs(const std::vector<std::string>& vcfs) {
		for (auto vcf : vcfs) {
			vcflib::VariantCallFile variantFile;
			variantFile.open(vcf);
			vcflib::Variant var(variantFile);

			PRINT("Adding mutations from: " << vcf << " #Samples:" <<
						variantFile.sampleNames.size());
			// insert samples into sampleid_map
			for (const auto sample : variantFile.sampleNames) {
				sampleid_map.insert(std::make_pair(sample, sampleid_map.size()));
				idsample_map.insert(std::make_pair(sampleid_map.size() - 1, sample));
			}

			long int num_mutations = 0;
			uint32_t num_samples_in_mutation = 0;
			while (variantFile.getNextVariant(var)) {
				num_mutations += 1;
				if (num_mutations % 10000 == 0) {
					PRINT("Mutations added: " << num_mutations << " #Vertices: " <<
								get_num_vertices() << " #Edges: " << get_num_edges());
					PRINT("Average num samples in mutations: " <<
								num_samples_in_mutation / (double)10000);
					num_samples_in_mutation = 0;
				}
				for (const auto alt : var.alt) {
					std::vector<sample_struct> sample_list;
					if (std::regex_match(alt, std::regex("^[ACTG]+$"))) {
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
								int present = 0;
								try {
									present = stoi(gt_info);
								} catch (const std::invalid_argument& ia) {}
								if (present) {
									add = true;
									gt1 = true;
									gt2 = false;
								}
							}
							if (add) {
								auto it = sampleid_map.find(sample.first);
								if (it == sampleid_map.end()) {
									ERROR("Unkown sample: " << sample.first);
									abort();
								} else {
									sample_struct s = {it->second, gt1, gt2};
									sample_list.emplace_back(s);
									num_samples_in_mutation++;
								}
							}
						}
					} else {
						//ERROR("Unsupported variant allele: " << alt);
					}
					if (sample_list.size() > 0)
						add_mutation(var.ref, alt, var.position, sample_list);
				}
			}
			PRINT("Num mutations: " << num_mutations);
		}
	}

	VariantGraphVertex* VariantGraph::add_vertex(uint64_t offset, uint64_t length,
																							 uint64_t index, uint32_t
																							 sample_id, bool
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
																							 uint64_t index, uint32_t
																							 sample_id, bool
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
		//DEBUG("Splitting vertex: " << vertex_id << " " << pos);
		const VariantGraphVertex cur_vertex = vertex_list.vertex(vertex_id);

		// create vertex object and add to vertex_list
		uint64_t offset = cur_vertex.offset() + pos - 1;
		uint64_t length = cur_vertex.length() - pos  + 1;
		if (length == 0) {
			ERROR("Vertex length is 0");
			abort();
		}
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
		return topology.get_num_vertices();
	}

	uint64_t VariantGraph::get_num_edges(void) const {
		return topology.get_num_edges();
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

	double VariantGraph::get_cache_hit_rate(void) const {
#ifdef DEBUG_MODE
		return cache.stats().hit_rate();
#endif
		return 0;
	}

	std::string VariantGraph::get_sample_name(uint32_t id) const {
		auto it = idsample_map.find(id);
		if (it == idsample_map.end()) {
			ERROR("Unknown sample id: " << id);
		}
		return it->second;
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
			samples.append("Sample id: " + std::to_string((int)s.sample_id()));
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
			// ref id is 0
			if (s.sample_id() == 0)
				idx_vertex_id[s.index()] = v.vertex_id();
		}
	}

	bool VariantGraph::get_sample_from_vertex_if_exists(Graph::vertex v,
																											uint32_t sample_id,
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

	bool VariantGraph::get_sample_from_vertex_if_exists(Graph::vertex v,
																											const std::string
																											sample_id,
																											VariantGraphVertex::sample_info&
																											sample) const {
		auto map_it = sampleid_map.find(sample_id);
		if (map_it == sampleid_map.end()) {
			ERROR("Sample not found: " << sample_id);
		}

		return get_sample_from_vertex_if_exists(v, map_it->second, sample);	
	}

	uint64_t VariantGraph::find_sample_index(Graph::vertex ref_v_id,
																					 uint32_t sample_id) const {
		Graph::vertex cur_vertex_id = ref_v_id;
		VariantGraphVertex cur_vertex = vertex_list.vertex(cur_vertex_id);
		// check if the cur vertex contains the sample
		VariantGraphVertex::sample_info sample;
		if (get_sample_from_vertex_if_exists(cur_vertex_id, sample_id, sample))
			return sample.index() + cur_vertex.length();

		// check in the cache
		if (cache.contains(sample_id))
			cur_vertex_id = cache.lookup(sample_id);
		else {	// else traverse back in the graph to find a vertex with @sample_id
			// ref id is 0
			get_sample_from_vertex_if_exists(ref_v_id, 0, sample);
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
		}

		// init cur distance
		uint64_t cur_distance = 0;
		if (get_sample_from_vertex_if_exists(cur_vertex_id, sample_id, sample)) {
			cur_distance += sample.index();
			// ref id is 0
		} else if (get_sample_from_vertex_if_exists(cur_vertex_id, 0, sample)) {
			cur_distance += sample.index();
		}
		// traverse the graph forward till you find the @ref_v_id
		auto map_it = idsample_map.find(sample_id);
		if (map_it == idsample_map.end()) {
			ERROR("Sample id not found: " << sample_id);
		}
		auto it = this->find(cur_vertex_id, map_it->second);
		while (!it.done() && (*it)->vertex_id() != ref_v_id) {
			cur_distance += (*it)->length();
			++it;
		}

		return cur_distance += vertex_list.vertex(ref_v_id).length();
	}

	// if there is not neighbor with @sample_id then set @v to "ref"
	bool VariantGraph::get_neighbor_vertex(Graph::vertex id, uint32_t sample_id,
																				 Graph::vertex* v) const {
			uint32_t min_idx = UINT32_MAX;
			for (const auto v_id : topology.out_neighbors(id)) {
				VariantGraphVertex vertex = vertex_list.vertex(v_id);
				for (int i = 0; i < vertex.s_info_size(); i++) {
					const VariantGraphVertex::sample_info& s = vertex.s_info(i);
					if (s.sample_id() != 0 && s.sample_id() == sample_id) {
						*v = v_id; 
						return true;
					} else if (s.sample_id() == 0) {	// if sample_id is not found follow "ref" path
						// if there are multiple outgoing ref vertexes then we return the
						// one with the smallest index.
						if (min_idx > s.index()) {
							*v = v_id;
							min_idx = s.index();
						}
					}
				}
			}
			if (*v != 0)
				return true;
			return false;
		}

	void VariantGraph::add_sample_to_vertex(Graph::vertex id, uint64_t
																					sample_idx, uint32_t  sample_id,
																					bool gt1, bool gt2) {
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
																	pos, std::vector<sample_struct>& sample_list)
	{
		// find the type mutation
		enum MUTATION_TYPE mutation;
		if (ref.size() == alt.size())
			mutation = SUBSTITUTION;
		else if (ref.size() > alt.size()) 
			mutation = DELETION;
		else
			mutation = INSERTION;

		//DEBUG("Adding mutation: " << mutation_string(mutation) << " " << ref <<
		//" " << alt << " " << pos << " " << sample_list.size());
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
				get_neighbor_vertex(ref_vertex_id, 0, &next_ref_vertex_id);		
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
			} else if (ref_vertex_idx == pos && ref_vertex.length() < ref.size()) { // substitution will skip 1 or more vertexes.
				// find the next ref vertex
				VariantGraphVertex next_ref_vertex;
				auto temp_itr = idx_vertex_id.lower_bound(ref_vertex_idx);
				do {
					++temp_itr;
					next_ref_vertex = vertex_list.vertex(temp_itr->second);
				} while (temp_itr->first < pos && temp_itr->first +
								 next_ref_vertex.length() < pos + ref.size());
				if (temp_itr->first + next_ref_vertex.length() == pos +
						ref.size()) {
					get_neighbor_vertex(temp_itr->second, 0, &next_ref_vertex_id);
				} else {
					// split the vertex
					split_vertex(temp_itr->second, pos + ref.size() - temp_itr->first + 1,
											 &next_ref_vertex_id);
				}
				temp_itr = idx_vertex_id.lower_bound(ref_vertex_idx);
				if (temp_itr->first != ref_vertex_idx) {
					ERROR("Vertex id not found in the map");
					abort();
				}
				--temp_itr;
				prev_ref_vertex_id = temp_itr->second;
			} else if (ref_vertex_idx < pos && ref_vertex_idx + ref_vertex.length()
								 > pos + ref.size()) { // vertex needs to be split into three
				// split the vertex
				uint64_t split_pos = pos - ref_vertex.offset();
				split_vertex(ref_vertex_id, split_pos, split_pos + ref.size(),
										 &prev_ref_vertex_id, &next_ref_vertex_id);
				std::swap(ref_vertex_id, prev_ref_vertex_id);
			} else if (ref_vertex_idx < pos && ref_vertex_idx + ref_vertex.length()
								 < pos + ref.size()) { // split the current vertex to get prev_vertex
				prev_ref_vertex_id = ref_vertex_id;
				split_vertex(prev_ref_vertex_id, pos - ref_vertex_idx + 1,
										 &ref_vertex_id);	
				// find the next ref vertex
				VariantGraphVertex next_ref_vertex;
				auto temp_itr = idx_vertex_id.lower_bound(ref_vertex_idx);
				do {
					++temp_itr;
					next_ref_vertex = vertex_list.vertex(temp_itr->second);
				} while (temp_itr->first < pos && temp_itr->first +
								 next_ref_vertex.length() < pos + ref.size());
				if (temp_itr->first + next_ref_vertex.length() == pos +
						ref.size()) {
					get_neighbor_vertex(temp_itr->second, 0, &next_ref_vertex_id);
				} else {
					// split the vertex
					split_vertex(temp_itr->second, pos + ref.size() - temp_itr->first + 1,
											 &next_ref_vertex_id);
				}
			}
			// create a vertex for the mutation using the first sample from the
			// list.
			uint32_t sample_id = sample_list[0].sample_id;
			bool gt1 = sample_list[0].gt1;
			bool gt2 = sample_list[0].gt2;
			//uint64_t sample_idx = find_sample_index(prev_ref_vertex_id, sample_id);
			VariantGraphVertex* sample_vertex = add_vertex(alt, 0,
																										 sample_id, gt1, gt2);
			// make connections for the new vertex in the graph
			topology.add_edge(prev_ref_vertex_id, sample_vertex->vertex_id());
			topology.add_edge(sample_vertex->vertex_id(), next_ref_vertex_id);

			// remove the first sample from the list
			sample_list.erase(sample_list.begin());

			// add rest of the samples to the vertex.
			for (auto sample : sample_list) {
				//sample_idx = find_sample_index(prev_ref_vertex_id,
																								//sample.sample_id);
				add_sample_to_vertex(sample_vertex->vertex_id(), 0,
														 sample.sample_id, sample.gt1, sample.gt2);
			}
		} else if (mutation == INSERTION) {
			Graph::vertex prev_ref_vertex_id = 0, next_ref_vertex_id = 0;
			// Either we got the vertex where there's already an insertion
			// or the vertex contains the seq with @pos
			if (ref_vertex_idx == pos) { // splitting not needed
				// find the prev ref vertex
				auto temp_itr = idx_vertex_id.lower_bound(ref_vertex_idx);
				if (temp_itr->first != ref_vertex_idx) {
					ERROR("Vertex id not found in the map");
					abort();
				}
				--temp_itr;
				prev_ref_vertex_id = temp_itr->second;

				// find the next ref vertex
				next_ref_vertex_id = ref_vertex_id;
			} else if (ref_vertex_idx < pos && ref_vertex_idx + ref_vertex.length()
								 > pos) {	// vertex needs to be split into two
				// split the vertex
				split_vertex(ref_vertex_id, pos - ref_vertex_idx + 1,
										 &next_ref_vertex_id);
				prev_ref_vertex_id = ref_vertex_id;
			} else { // to handle insertions after ref seq length. 
				prev_ref_vertex_id = ref_vertex_id;
			}
			// create a vertex for the mutation using the first sample from the
			// list.
			uint32_t sample_id = sample_list[0].sample_id;
			bool gt1 = sample_list[0].gt1;
			bool gt2 = sample_list[0].gt2;
			//uint64_t sample_idx = find_sample_index(prev_ref_vertex_id, sample_id);
			VariantGraphVertex* sample_vertex = add_vertex(alt, 0,
																										 sample_id, gt1, gt2);
			// make connections for the new vertex in the graph
			topology.add_edge(prev_ref_vertex_id, sample_vertex->vertex_id());
			if (next_ref_vertex_id != 0)
				topology.add_edge(sample_vertex->vertex_id(), next_ref_vertex_id);

			// remove the first sample from the list
			sample_list.erase(sample_list.begin());

			// add rest of the samples to the vertex.
			for (auto sample : sample_list) {
				//sample_idx = find_sample_index(prev_ref_vertex_id,
																								//sample.sample_id);
				add_sample_to_vertex(sample_vertex->vertex_id(), 0,
														 sample.sample_id, sample.gt1, sample.gt2);
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
				get_neighbor_vertex(ref_vertex_id, 0, &next_ref_vertex_id);	
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
			} else if (ref_vertex_idx == pos && ref_vertex.length() < ref.size()) { // deletion will skip 1 or more vertexes
				// find the next ref vertex
				VariantGraphVertex next_ref_vertex;
				auto temp_itr = idx_vertex_id.lower_bound(ref_vertex_idx);
				do {
					++temp_itr;
					next_ref_vertex = vertex_list.vertex(temp_itr->second);
				} while (temp_itr->first < pos && temp_itr->first +
								 next_ref_vertex.length() < pos + ref.size());
				if (temp_itr->first + next_ref_vertex.length() == pos +
						ref.size()) {
					get_neighbor_vertex(temp_itr->second, 0, &next_ref_vertex_id);
				} else {
					// split the vertex
					split_vertex(temp_itr->second, pos + ref.size() - temp_itr->first + 1,
											 &next_ref_vertex_id);
				}
				temp_itr = idx_vertex_id.lower_bound(ref_vertex_idx);
				if (temp_itr->first != ref_vertex_idx) {
					ERROR("Vertex id not found in the map");
					abort();
				}
				--temp_itr;
				prev_ref_vertex_id = temp_itr->second;
			} else if (ref_vertex_idx < pos && ref_vertex_idx + ref_vertex.length()
								 > pos + ref.size()) { // vertex needs to be split into three
				// split the vertex
				uint64_t split_pos = pos - ref_vertex.offset();
				split_vertex(ref_vertex_id, split_pos, split_pos + ref.size(),
										 &prev_ref_vertex_id, &next_ref_vertex_id);
				std::swap(ref_vertex_id, prev_ref_vertex_id);
			} else if (ref_vertex_idx < pos && ref_vertex_idx + ref_vertex.length()
								 < pos + ref.size()) { // split the current vertex to get prev_vertex
				prev_ref_vertex_id = ref_vertex_id;
				split_vertex(prev_ref_vertex_id, pos - ref_vertex_idx + 1,
										 &ref_vertex_id);	
				// find the next ref vertex
				VariantGraphVertex next_ref_vertex;
				auto temp_itr = idx_vertex_id.lower_bound(ref_vertex_idx);
				do {
					++temp_itr;
					next_ref_vertex = vertex_list.vertex(temp_itr->second);
				} while (temp_itr->first < pos && temp_itr->first +
								 next_ref_vertex.length() < pos + ref.size());
				if (temp_itr->first + next_ref_vertex.length() == pos +
						ref.size()) {
					get_neighbor_vertex(temp_itr->second, 0, &next_ref_vertex_id);
				} else {
					// split the vertex
					split_vertex(temp_itr->second, pos + ref.size() - temp_itr->first + 1,
											 &next_ref_vertex_id);
				}
			}
			// add samples to the vertex.
			for (auto sample : sample_list) {
				//uint64_t sample_idx = find_sample_index(prev_ref_vertex_id,
																								//sample.sample_id);
				add_sample_to_vertex(next_ref_vertex_id, 0, sample.sample_id,
														 sample.gt1, sample.gt2);
			}
			// make connections for the new vertex in the graph
			topology.add_edge(prev_ref_vertex_id, next_ref_vertex_id);
		}
	}

	void VariantGraph::fix_sample_indexes(void) {
#if 0
		for (const auto sample : sampleid_map) {
			if (sample.first == "ref")
				continue;
			uint32_t cur_index = 1;
			auto path = mutable_find(sample.first);
			while (!path.done()) {
				auto vertex = *path;
				for (int i = 0; i < vertex->s_info_size(); ++i) {
					VariantGraphVertex::sample_info* s = vertex->mutable_s_info(i);
					if (s->sample_id() == sample.second) {
						s->set_index(cur_index);
					}
				}
				cur_index += vertex->length();
				++path;
			}
		}
#else	
		// Optimized solution.
		// Perform a BFS on the graph. From each vertex, explore the neighbors and
		// is a sample_info is stored at the vertex then update the sample index
		// as vertex_index + vertex_len + delta
		// We update the delta if we travel from a vertex with sample info to a
		// ref vertex. Or from a ref vertex to a sample vertex where sample index
		// is already been set.
		PRINT("Fixing indexes");
		// map to keep track of delta for samples.
		std::unordered_map<uint32_t, int32_t> sampleid_delta;

		// init the index map with all deltas as 0
		for (const auto sample : idsample_map)
			sampleid_delta.insert(std::make_pair(sample.first, 0));

		// init BFS iterator
		auto it = mutable_find();
		while(!it.done()) {
			auto cur_vertex = *it;
			// get neighbors of the vertex
			VariantGraphVertex::sample_info ref_sample; 
			if (get_sample_from_vertex_if_exists(cur_vertex->vertex_id(), 0,
																					 ref_sample)) { // this is a ref vertex
				uint32_t ref_index = ref_sample.index();
				for (auto neighbor_id :
						 topology.out_neighbors(cur_vertex->vertex_id())) {
					auto cur_neighbor = vertex_list.mutable_vertex(neighbor_id);
					// for each sample other than "ref" in the vertex update the index
					for (int i = 0; i < cur_neighbor->s_info_size(); ++i) {
						VariantGraphVertex::sample_info* s = cur_neighbor->mutable_s_info(i);
						if (s->sample_id() != 0 && s->index() == 0) {	// 0 is the ref sample
							auto map_it = sampleid_delta.find(s->sample_id());
							if (map_it == sampleid_delta.end()) {
								ERROR("Unknown sample id: " << s->sample_id());
								abort();
							}
							int32_t delta = map_it->second;
							int32_t sample_index = ref_index + cur_vertex->length() + delta;
							if (sample_index < 0) {
								ERROR("Sample index is less the 0: " <<
											get_sample_name(s->sample_id()) << " " <<
											cur_vertex->vertex_id());
								abort();
							}
							s->set_index(sample_index);
						} else if (s->sample_id() != 0 && s->index() != 0) {	// reset the sample index
							sampleid_delta[s->sample_id()] =  s->index() - (ref_index +
																															cur_vertex->length());
						}
					}
				}
			} else { // this is not a ref vertex
				// only sample vertexes should only have one outgoing edge.
				if (topology.out_neighbors(cur_vertex->vertex_id()).size() > 1) {
					ERROR("Sample vertex has more than 1 neighbor: " <<
								cur_vertex->vertex_id());
					abort();
				}
				for (auto neighbor_id :
						 topology.out_neighbors(cur_vertex->vertex_id())) {
					auto cur_neighbor = vertex_list.mutable_vertex(neighbor_id);
					// for each sample other than "ref" in the vertex update the index
					uint32_t cur_length = cur_vertex->length();
					// for each sample update the delta
					for (int i = 0; i < cur_vertex->s_info_size(); ++i) {
						const VariantGraphVertex::sample_info s = cur_vertex->s_info(i);
						uint32_t cur_index = s.index();
						VariantGraphVertex::sample_info ref_sample; 
						if (!get_sample_from_vertex_if_exists(cur_neighbor->vertex_id(), 0,
																								 ref_sample)) { // this is a ref vertex
							ERROR("Ref vertex not found as a neighbor from sample vertex. "
										<< cur_neighbor->vertex_id());
							abort();
						}
						sampleid_delta[s.sample_id()] = cur_index + cur_length -
							ref_sample.index();
					}
				}
			}
			++it;
		}
#endif
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
		auto it = vg->sampleid_map.find(sample_id);
		if (it == vg->sampleid_map.end()) {
			ERROR("Sample not found");
			abort();
		} else {
			s_id = it->second;
		}

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

	VariantGraph::VariantGraphPathMutableIterator::VariantGraphPathMutableIterator(
									VariantGraph* g, Graph::vertex v, const std::string
									sample_id) {
		vg = g; 
		cur = vg->vertex_list.mutable_vertex(v);
		auto it = vg->sampleid_map.find(sample_id);
		if (it == vg->sampleid_map.end()) {
			ERROR("Sample not found");
			abort();
		} else {
			s_id = it->second;
		}

		is_done = false;
	}

	VariantGraphVertex*
		VariantGraph::VariantGraphPathMutableIterator::operator*(void) {
			return cur;
		}

	void VariantGraph::VariantGraphPathMutableIterator::operator++(void) {
		Graph::vertex next_vertex = 0; // there shoould be no incoming edge to 0
		if (!vg->get_neighbor_vertex(cur->vertex_id(), s_id, &next_vertex) &&
				next_vertex == 0) {
			is_done = true;
		}
		cur = vg->vertex_list.mutable_vertex(next_vertex);
	}

	bool VariantGraph::VariantGraphPathMutableIterator::done(void) const {
		return is_done;
	}

	// iterator will be positioned at the start of the path.
	VariantGraph::VariantGraphPathMutableIterator VariantGraph::mutable_find(const
																																	 std::string
																																	 sample_id) {
		return VariantGraphPathMutableIterator(this, 0, sample_id);	
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

	VariantGraph::VariantGraphMutableIterator::VariantGraphMutableIterator(VariantGraph*
																																				 g,
																																				 Graph::vertex
																																				 v,
																																				 uint64_t
																																				 r) :
		vg(g), itr(&g->topology, v, r) {};

	VariantGraphVertex*
		VariantGraph::VariantGraphMutableIterator::operator*(void) {
			return vg->vertex_list.mutable_vertex(*itr);
		}

	void VariantGraph::VariantGraphMutableIterator::operator++(void) {
		++itr;
	}

	bool VariantGraph::VariantGraphMutableIterator::done(void) const {
		return itr.done();
	}

	VariantGraph::VariantGraphMutableIterator
		VariantGraph::mutable_find(Graph::vertex v, uint64_t radius) {
		return VariantGraphMutableIterator(this, v, radius);
	}

}

#endif //__VARIANT_GRAPH_H__
