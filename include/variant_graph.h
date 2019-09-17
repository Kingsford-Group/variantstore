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

#include "gqf/hashutil.h"
#include "util.h"
#include "variantgraphvertex.pb.h"
#include "graph.h"
#include "variantdb_fs.h"

namespace variantdb {

#define CACHE_SIZE 256
//#define NUM_VERTEXES_IN_BLOCK 1000000
#define NUM_VERTEXES_IN_BLOCK 10000

	// to map a sample --> vertex ids.
	using Cache = LRU::Cache<uint32_t, Graph::vertex>;

	static inline uint64_t word_rank(uint64_t val) {
		asm("popcnt %[val], %[val]"
				: [val] "+r" (val)
				:
				: "cc");
		return val;
	}

	// Returns the position of the rank'th 1.  (rank = 0 returns the 1st 1)
	// Returns 64 if there are fewer than rank+1 1s.
	static inline uint64_t word_select(uint64_t val, int rank) {
		uint64_t i = 1ULL << rank;
		asm("pdep %[val], %[mask], %[val]"
				: [val] "+r" (val)
				: [mask] "r" (i));
		asm("tzcnt %[bit], %[index]"
				: [index] "=r" (i)
				: [bit] "g" (val)
				: "cc");
		return i;
	}

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

	enum READ_TYPE {
		READ_COMPLETE_GRAPH,
		READ_INDEX_ONLY
	};

	enum MODE {
		READ_WRITE_MODE,
		READ_ONLY_MODE
	};

	class VariantGraph {
		public:
			// can't create a VariantGraph without a reference genome and zero or
			// more vcf files.
			VariantGraph() = delete; 
			// construct variant graph using a reference genome and zero or more vcf
			// files.
			VariantGraph(const std::string& ref_file, const std::string& vcf_file,
									 const std::string& prefix);

			// read variant graph from disk
			VariantGraph(const std::string& prefix, enum READ_TYPE type);

			~VariantGraph();

			// add one or more vcf files to the variant graph
			void add_vcfs(const std::string& vcf_file);

			// persist variant graph to disk
			void serialize();

			uint64_t get_num_vertices(void) const;
			uint64_t get_num_edges(void) const;
			uint64_t get_seq_length(void) const;
			const std::string get_chr(void) const;
			uint64_t get_ref_length(void) const;
			uint64_t get_num_sample_classes(void) const;
			std::string get_sample_name(uint32_t id) const;
			//double get_cache_hit_rate(void) const;

			std::string print_vertex_info(const VariantGraphVertex& v) const;
			const std::string get_sequence(const VariantGraphVertex& v) const;

			bool get_sample_from_vertex_if_exists(Graph::vertex v, const std::string
																						sample_id,
																						VariantGraphVertex::sample_info&
																						sample) const;
			uint32_t get_sample_id(uint32_t sampleclass_id, uint32_t index) const;
			//uint32_t get_sample_id_fast(uint32_t sampleclass_id, uint32_t index)
				//const;

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
																UINT64_MAX) const;

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

			const VariantGraphVertex& get_vertex(Graph::vertex id) const;
			VariantGraphVertex* get_mutable_vertex(Graph::vertex id);
			void serialize_vertex_list(uint64_t index) const;
			void load_vertex_list(uint64_t index) const;

			void flush_in_memory_vertex_lists(void) const;
			void add_or_replace_in_memory_vertex_list(uint32_t index,
																								VariantGraphVertexList* v) const;

			bool get_sample_from_vertex_if_exists(Graph::vertex v, uint32_t sample_id,
																						VariantGraphVertex::sample_info&
																						sample) const;
			const std::string mutation_string(MUTATION_TYPE m) const;

			void update_idx_vertex_id_map(const VariantGraphVertex& v);
			uint64_t find_sample_index(Graph::vertex ref_v_id, uint32_t sample_id)
				const;
			bool get_neighbor_vertex(Graph::vertex id, uint32_t sample_id,
															 Graph::vertex* v) const;
			const std::string get_sequence(uint64_t start, uint32_t length) const;
			void add_sample_to_vertex(Graph::vertex id, uint64_t sample_idx, bool
																gt1, bool gt2);
			bool check_if_mutation_exists(Graph::vertex prev, Graph::vertex next,
																		const std::string alt, Graph::vertex* v)
				const;
			bool check_if_mutation_exists(Graph::vertex prev, uint64_t offset,
																		uint64_t length, Graph::vertex* v) const;
			void validate_ref_path_edge(Graph::vertex src, Graph::vertex dest) const;
			VariantGraphVertex* create_vertex(uint64_t id, uint64_t offset, uint64_t
																				length, uint32_t sampleclass_id,
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
																		 uint64_t index, uint32_t sampleclass_id,
																		 bool gt1, bool gt2);
			// returns the vertex_id of the new vertex
			VariantGraphVertex* add_vertex(const std::string& seq, uint64_t index,
																		 uint32_t sampleclass_id, bool gt1, bool
																		 gt2);
			void add_sample_vector(const sdsl::bit_vector& vector, uint64_t
														 sampleclass_id);
			uint32_t find_sample_vector_or_add(const std::vector<sample_struct>&
																				 sample_list);
			bool update_vertex_sample_class(Graph::vertex vertex_id, const
																			std::vector<sample_struct>& sample_list);
			uint32_t get_popcnt(uint32_t sampleclass_id) const;
			sdsl::bit_vector get_bit_vector(uint32_t sampleclass_id) const;

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
			std::string prefix;
			uint64_t ref_length{0};
			uint64_t num_samples{0};
			enum MODE mode;
			enum READ_TYPE load_type;
			std::map<uint64_t, uint64_t> idx_vertex_id;
			std::unordered_map<uint32_t, std::string> idsample_map;
			std::unordered_map<uint64_t, uint32_t> sampleclass_map;
			sdsl::rrr_vector<SDSL_BITVECTOR_BLOCK_SIZE> rrr_sample_vector;
			mutable std::map<uint32_t, VariantGraphVertexList*> in_memory_vertex_lists;

			/* structures to persist when serializing variant graph. */
			// a vector of VariantGraphVertexLists. Each list object is a block
			// containing vertexes.
			std::vector<VariantGraphVertexList> vertex_block_list;
			sdsl::int_vector<> seq_buffer;
			sdsl::bit_vector sample_vector;
			Graph topology;
			std::unordered_map<std::string, uint32_t> sampleid_map;
	};

	VariantGraph::VariantGraph(const std::string& ref_file, const
														 std::string& vcf_file,const std::string& prefix) :
		prefix(prefix), mode(READ_WRITE_MODE), load_type(READ_INDEX_ONLY) {
			// Verify that the version of the library that we linked against is
			// compatible with the version of the headers we compiled against.
			GOOGLE_PROTOBUF_VERIFY_VERSION;

			std::string ref;
			read_fasta(ref_file, chr, ref);
			ref_length = ref.size();
			// initialize the seq buffer
			sdsl::util::assign(seq_buffer, sdsl::int_vector<>(0, 0, 3));
			// initialize the sample vector
			sdsl::util::assign(sample_vector, sdsl::bit_vector(0,0));

			// add ref node
			// we start the index to 1.
			sampleid_map.insert(std::make_pair("ref", sampleid_map.size()));
			idsample_map.insert(std::make_pair(0, "ref"));
			// ref id is 0
			VariantGraphVertex *v = add_vertex(ref, 1, 0, 0, 0);

			// update idx->vertex_id map
			update_idx_vertex_id_map(*v);

			// Add vcf files
			add_vcfs(vcf_file);

			// flush remaining vertex lists.
			flush_in_memory_vertex_lists();

			// fix sample indexes
			fix_sample_indexes();
		}

	VariantGraph::VariantGraph(const std::string& prefix, enum READ_TYPE type) :
		prefix(prefix), mode(READ_ONLY_MODE), load_type(type), topology(prefix) {

		// Read all proto files and sort them based on ids.
		std::vector<std::string> proto_files = fs::GetFilesExt(prefix.c_str(), ".proto");
		std::map<int, std::string> sorted_proto_files;
		for (auto file : proto_files) {
			int id = std::stoi(file.substr(file.find_last_of('_') + 1,
																		 file.find_last_of('.')));
			sorted_proto_files[id] = file;
		}

		if (load_type == READ_COMPLETE_GRAPH) {
			function<void(VariantGraphVertexList&, uint32_t)> lambda =
				[this](VariantGraphVertexList& v, uint32_t index)
				{
					vertex_block_list.emplace_back(v);
					//VariantGraphVertex* vertex = vertex_list.add_vertex();
					//*vertex = v;
				};

			for (auto file_obj : sorted_proto_files) {
				// load vertex list
				std::string vertex_list_name = prefix + "/" + file_obj.second;
				fstream input(vertex_list_name, ios::in | ios::binary);
				if (!stream::for_each(input, 0, lambda)) {	// passing 0. Index is not used here.
					console->error("Failed to parse vertex list {}.", vertex_list_name);
					abort();
				}
			}
		}

		// load seq buffer
		std::string seq_buffer_name = prefix + "/seq_buffer.sdsl";
		if (!sdsl::load_from_file(seq_buffer, seq_buffer_name)) {
			console->error("Failed to load seq buffer {}.", seq_buffer_name);
			abort();
		}
		num_vertices = topology.get_num_vertices() + 1;
		seq_length = seq_buffer.size();

		// load sample bit vector
		std::string sample_vector_name = prefix + "/sample_vector.sdsl";
		if (!sdsl::load_from_file(rrr_sample_vector, sample_vector_name)) {
			console->error("Failed to load sample vector {}.", sample_vector_name);
			abort();
		}

		//load sampleid map
		std::string sampleid_map_name = prefix + "/sampleid_map.lst";
		std::ifstream sampleid_file(sampleid_map_name);
		if (!sampleid_file.good()) {
			console->error("Failed to open sampleid map file {}.", sampleid_map_name);
			abort();
		}
		std::string sample;
		uint32_t id;
		// read chr name
		sampleid_file >> chr >> ref_length;
		// read num samples
		sampleid_file >> sample >> num_samples;
		if (num_samples <= 0) {
			console->error("Num samples is less or equal to 0.");
			abort();
		}
		while (sampleid_file >> sample >> id)  {
			sampleid_map.insert(std::make_pair(sample, id));
			idsample_map.insert(std::make_pair(id, sample));
		}
		if (num_samples != sampleid_map.size()) {
			console->error("Num samples is not equal to num entries in samples file.");
			abort();
		}

		sampleid_file.close();
	}

	VariantGraph::~VariantGraph() {
		// Optional:  Delete all global objects allocated by libprotobuf.
		google::protobuf::ShutdownProtobufLibrary();
	}

	void VariantGraph::serialize_vertex_list(uint64_t index) const {
		std::function<VariantGraphVertexList(uint64_t)> lambda =
			[this](uint64_t n) {
				if (load_type == READ_COMPLETE_GRAPH) {
					return vertex_block_list[n];
				} else {
					auto itr = in_memory_vertex_lists.find(n);
					if (itr == in_memory_vertex_lists.end()) {
						console->error("Can't find the vertex list for id: {}", n);
						abort();
					}
					return *(itr->second);
				}
				//return vertex_list.vertex(n);
			};

		// serialize vertex list
		std::string vertex_list_name = prefix + "/vertex_list_" +
			std::to_string(index) + ".proto";
		std::fstream output(vertex_list_name, ios::out | ios::trunc | ios::binary);
		if (!stream::write(output, index, lambda)) {
			console->error("Failed to write vertex list.", vertex_list_name);
			abort();
		}
	}

	void VariantGraph::load_vertex_list(uint64_t index) const {
		function<void(VariantGraphVertexList&, uint32_t)> lambda =
			[this](VariantGraphVertexList& v, uint32_t index)
			{
				VariantGraphVertexList *list = new VariantGraphVertexList();
				*list = v;
				add_or_replace_in_memory_vertex_list(index, list);
				//VariantGraphVertex* vertex = vertex_list.add_vertex();
				//*vertex = v;
			};

		// load vertex list
		std::string vertex_list_name = prefix + "/vertex_list_" +
			std::to_string(index) + ".proto";
		fstream input(vertex_list_name, ios::in | ios::binary);
		if (!stream::for_each(input, index, lambda)) {
			console->error("Failed to parse vertex list {}.", vertex_list_name);
			abort();
		}
	}

	void VariantGraph::serialize() {
		if (mode == READ_ONLY_MODE) {
			console->error("Serialization not allowed. VariantDB is loaded in READ ONLY mode");
			abort();
		}

		if (load_type == READ_COMPLETE_GRAPH) {
			for (uint64_t i = 0; i < vertex_block_list.size(); i++)
				serialize_vertex_list(i);
		} else {
			for (auto itr : in_memory_vertex_lists)
				serialize_vertex_list(itr.first);
		}
		// clear in-memory vertex lists
		in_memory_vertex_lists.clear();

		// serialize seq buffer
		std::string seq_buffer_name = prefix + "/seq_buffer.sdsl";
		seq_buffer.resize(seq_buffer.size());
		sdsl::util::bit_compress(seq_buffer);
		if (!sdsl::store_to_file(seq_buffer, seq_buffer_name)) {
			console->error("Failed to serialize seq buffer {}.", seq_buffer_name);
			abort();
		}

		// serialize topology
		topology.serialize(prefix);

		console->info("Number of sample vector classes: {}", sampleclass_map.size());
		// serialize sample vector
		sdsl::util::assign(rrr_sample_vector,
											 sdsl::rrr_vector<SDSL_BITVECTOR_BLOCK_SIZE>(sample_vector));
		std::string sample_vector_name = prefix + "/sample_vector.sdsl";
		if (!sdsl::store_to_file(rrr_sample_vector, sample_vector_name)) {
			console->error("Failed to serialize compressed sample vector {}.",
										 sample_vector_name);
			abort();
		}

		// serialize sampleid_map
		std::string sampleid_map_name = prefix + "/sampleid_map.lst";
		std::ofstream sampleid_file(sampleid_map_name);
		if (!sampleid_file.good()) {
			console->error("Failed to open sampleid file {}.", sampleid_map_name);
			abort();
		}

		// write the chromosome name and length
		sampleid_file << chr << " " << std::to_string(ref_length) << "\n";
		// write num samples
		sampleid_file << chr << " " << std::to_string(num_samples) << "\n";

		for (const auto sample : sampleid_map)
			sampleid_file << sample.first << " " << sample.second << "\n";
		sampleid_file.close();
	}

	void VariantGraph::add_vcfs(const std::string& vcf_file) {
		std::string m_vcf_file = vcf_file;
		vcflib::VariantCallFile variantFile;
		variantFile.open(m_vcf_file);
		vcflib::Variant var(variantFile);

		console->info("Adding mutations from: {} #Samples: {}", vcf_file,
									variantFile.sampleNames.size());
		// insert samples into sampleid_map
		for (const auto sample : variantFile.sampleNames) {
			sampleid_map.insert(std::make_pair(sample, sampleid_map.size()));
			idsample_map.insert(std::make_pair(sampleid_map.size() - 1, sample));
		}
		num_samples += sampleid_map.size();

		long int num_mutations = 0;
		long int num_mutations_samples = 0;
		uint32_t num_samples_in_mutation = 0;
		while (variantFile.getNextVariant(var)) {
			// verify mutation.
			if (var.sequenceName != chr || var.position < 1 ||
					(uint64_t)var.position > ref_length || var.ref !=
					get_sequence(var.position - 1, var.ref.size())) {
				console->error("Unsupported mutation: {} {} {}", var.sequenceName,
											 var.position, var.ref);
				continue;
			}
			num_mutations += 1;
			if (num_mutations % 100000 == 0) {
				console->debug("Mutations added: {} #Vertices: {} #Edges: ",
											 num_mutations, get_num_vertices(), get_num_edges());
				console->debug("Average num samples in mutations: {}",
											 num_samples_in_mutation / (double)100000);
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
								console->error("Unkown sample: {}", sample.first);
								abort();
							} else {
								sample_struct s = {it->second, gt1, gt2};
								sample_list.emplace_back(s);
								num_samples_in_mutation++;
							}
						}
					}
				} else {
					//console->error("Unsupported variant allele: {} {}", var.position,
					//alt);
					continue;
				}
				if (sample_list.size() > 0) {
					add_mutation(var.ref, alt, var.position, sample_list);
					num_mutations_samples += sample_list.size();
				}
			}
		}
		console->info("Num mutations: {}", num_mutations_samples);
	}

	VariantGraphVertex* VariantGraph::add_vertex(uint64_t offset, uint64_t length,
																							 uint64_t index, uint32_t
																							 sampleclass_id, bool
																							 gt1, bool gt2) {
		// create vertex object and add to vertex_list
		VariantGraphVertex::sample_info s;
		s.set_index(index);
		//s.set_sample_id(sample_id);
		s.set_gt_1(gt1);
		s.set_gt_2(gt2);
		std::vector<VariantGraphVertex::sample_info> samples = {s};
		VariantGraphVertex *v = create_vertex(num_vertices, offset, length,
																					sampleclass_id, samples);

		// increment vertex count
		num_vertices++;

		return v;
	}

	VariantGraphVertex* VariantGraph::add_vertex(const std::string& seq,
																							 uint64_t index, uint32_t
																							 sampleclass_id, bool
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
		//s.set_sample_id(sample_id);
		s.set_gt_1(gt1);
		s.set_gt_2(gt2);
		std::vector<VariantGraphVertex::sample_info> samples = {s};
		VariantGraphVertex *v = create_vertex(num_vertices, start_offset,
																					seq.size(), sampleclass_id, samples);

		// increment vertex count
		num_vertices++;

		return v;
	}

	void VariantGraph::add_sample_vector(const sdsl::bit_vector& vector,
																			 uint64_t sampleclass_id) {
		if (sampleclass_id < 1) {
			console->error("Sample class is smaller than 1.");
			abort();
		}
		sample_vector.resize(sample_vector.size() + num_samples);
		uint64_t start_idx = (sampleclass_id - 1) * num_samples;
		for (uint32_t i = 0; i < num_samples/64*64; i+=64)
			sample_vector.set_int(start_idx+i, vector.get_int(i, 64), 64);
		if (num_samples%64)
			sample_vector.set_int(start_idx+num_samples/64*64,
														vector.get_int(num_samples/64*64, num_samples%64),
														num_samples%64);
	}

	uint32_t VariantGraph::find_sample_vector_or_add(const
																									 std::vector<sample_struct>&
																									 sample_list) {
		sdsl::bit_vector vector(num_samples, 0);
		for (const auto sample : sample_list) {
			vector[sample.sample_id] = 1;
		}

		uint64_t vec_hash = MurmurHash64A((void*)vector.data(),
																			vector.capacity()/8, 2038074743);
		auto it = sampleclass_map.find(vec_hash);
		if (it == sampleclass_map.end()) {
			uint32_t sampleclass_id = sampleclass_map.size() + 1; // 0 is reserved for ref sample
			sampleclass_map.insert(std::make_pair(vec_hash, sampleclass_id));
			add_sample_vector(vector, sampleclass_id);	// bit vectors are stored 0-based.
			// validate: query each sample id from the bit vector.
			//int i = 0;
			//for (const auto sample : sample_list) {
			//uint32_t id = get_sample_id(sampleclass_id, i++);
			//if (id != sample.sample_id) {
			//console->error("Bit vector construction/query failed. Sample class id: {}",
			//sampleclass_id);
			//abort();
			//}
			//}
			return sampleclass_id;
		} else {
			return it->second;
		}
	}

	bool VariantGraph::update_vertex_sample_class(Graph::vertex vertex_id, const
																								std::vector<sample_struct>&
																								sample_list) {
		std::map<uint32_t, sample_struct> sample_indexes;
		std::vector<sample_struct> list;  
		for (const auto sample : sample_list)
			sample_indexes.insert(std::make_pair(sample.sample_id, sample));

		VariantGraphVertex *v = get_mutable_vertex(vertex_id);
		uint64_t ref_index = 0;
		for (int i = v->s_info_size() - 1; i >= 0; i--) {
			sample_struct s = {get_sample_id(v->sampleclass_id(), i), 0, 0};
			// save ref index
			if (s.sample_id == 0)
				ref_index = v->s_info(i).index();
			sample_indexes.insert(std::make_pair(s.sample_id, s));
		}
		for (const auto sample : sample_indexes)
			list.emplace_back(sample.second);

		uint32_t class_id = find_sample_vector_or_add(list);
		v->set_sampleclass_id(class_id);
		//validate list size and popcnt of the vector.
		if (list.size() != get_popcnt(class_id)) {
			console->error("Sample class update failed. list size: {} popcnt: {}",
										 list.size(), get_popcnt(class_id));
			return false;
		}

		// clear s_info and reinsert sample indexes in the order of sample_ids.
		v->clear_s_info();
		for (const auto sample : list) {
			if (sample.sample_id == 0)
				add_sample_to_vertex(vertex_id, ref_index, sample.gt1, sample.gt2);
			else
				add_sample_to_vertex(vertex_id, 0, sample.gt1, sample.gt2);
		}

		return true;
	}

	uint32_t VariantGraph::get_sample_id(uint32_t sampleclass_id, uint32_t
																			 index) const {
		if (sampleclass_id == 0) { // only ref sample
			return 0;
		} else {
			uint64_t start_idx = (sampleclass_id - 1) * num_samples;
			uint64_t vector_offset = start_idx;
			uint32_t rank = index + 1;
			for (uint32_t i = 0; i < num_samples/64*64; i+=64) {
				uint64_t word;
				if (mode == READ_WRITE_MODE)
					word = sample_vector.get_int(start_idx+i, 64);
				else
					word = rrr_sample_vector.get_int(start_idx+i, 64);

				if (word_rank(word) >= rank) {
					return word_select(word, rank - 1) + start_idx + i - vector_offset;
				} else {
					rank -= word_rank(word);
				}
			}

			if (num_samples%64) {
				uint64_t word;
				if (mode == READ_WRITE_MODE)
					word = sample_vector.get_int(start_idx+num_samples/64*64,
																			 num_samples%64);
				else
					word = rrr_sample_vector.get_int(start_idx+num_samples/64*64,
																					 num_samples%64);

				if (word_rank(word) >= rank) {
					return word_select(word, rank - 1) + num_samples/64*64;
				} else {
					console->error("Index {} passed is outside the bounds for sample class {} popcnt: {}.",
												 index, sampleclass_id, get_popcnt(sampleclass_id));
					abort();
				}
			}
		}
		return UINT32_MAX;
	}

	sdsl::bit_vector VariantGraph::get_bit_vector(uint32_t sampleclass_id) const
	{
		uint64_t start_idx = (sampleclass_id - 1) * num_samples;
		sdsl::bit_vector sampleclass_vector(num_samples, 0);
		for (uint32_t i = 0; i < num_samples/64*64; i+=64)
			sampleclass_vector.set_int(i, sample_vector.get_int(start_idx+i, 64),
																 64);
		if (num_samples%64)
			sampleclass_vector.set_int(num_samples/64*64,
																 sample_vector.get_int(start_idx+num_samples/64*64,
																											 num_samples%64),
																 num_samples%64);
		return sampleclass_vector;
	}

	//uint32_t VariantGraph::get_sample_id(uint32_t sampleclass_id, uint32_t
																			 //index) const {
		//if (sampleclass_id == 0) { // only ref sample
			//return 0;
		//} else {
			//// extract sample class vector for sampleclass_id
			//sdsl::bit_vector sampleclass_vector = get_bit_vector(sampleclass_id);

			//// select based on index.
			////sdsl::rrr_vector<SDSL_BITVECTOR_BLOCK_SIZE> rrr_vec(sampleclass_vector);
			////sdsl::rrr_vector<SDSL_BITVECTOR_BLOCK_SIZE>::select_1_type select_vec(&rrr_vec);
			//sdsl::bit_vector::select_1_type select_vec(&sampleclass_vector);


			//// index 0 means the first 1 in the select vector.
			//uint32_t sample_id = select_vec(index + 1);
			//uint32_t sample_id_fast = get_sample_id_fast(sampleclass_id, index);
			//if (sample_id != sample_id_fast) {
				//std::cout << sampleclass_vector << '\n';
				////uint32_t sample_id_fast = get_sample_id_fast(sampleclass_id, index);
				//console->error("Sample ids don't match old: {} new: {} index: {}",
											 //sample_id, sample_id_fast, index);
				//abort();
			//}
			//if (sample_id >= num_samples) {
				//console->error("Index {} passed is outside the bounds for sample class {} popcnt: {}.",
											 //index, sampleclass_id, get_popcnt(sampleclass_id));
				//abort();
			//}
			//return sample_id;
		//}
	//}

	uint32_t VariantGraph::get_popcnt(uint32_t sampleclass_id) const {
		if (sampleclass_id == 0) { // only ref sample
			return 1;
		} else {
			uint64_t start_idx = (sampleclass_id - 1) * num_samples;
			uint64_t popcnt = 0;
			for (uint32_t i = 0; i < num_samples/64*64; i+=64) {
				uint64_t word = sample_vector.get_int(start_idx+i, 64);
				popcnt += word_rank(word);
			}

			if (num_samples%64) {
				uint64_t word = sample_vector.get_int(start_idx+num_samples/64*64,
																							num_samples%64);
				popcnt += word_rank(word);
			} 
			return popcnt;
		}
	}

	void VariantGraph::add_or_replace_in_memory_vertex_list(uint32_t index,
																													VariantGraphVertexList*
																													v) const {
		if (in_memory_vertex_lists.size() < 2) {
			in_memory_vertex_lists[index] = v;
		} else {
			auto itr = in_memory_vertex_lists.begin();
			uint32_t id1 = itr->first;
			++itr;
			uint32_t id2 = itr->first;
			uint32_t remove_id = id1 < id2 ? id1 : id2;
			serialize_vertex_list(remove_id);
			delete in_memory_vertex_lists.find(remove_id)->second;
			in_memory_vertex_lists.erase(in_memory_vertex_lists.find(remove_id));
			in_memory_vertex_lists[index] = v;
		}
	}

	void VariantGraph::flush_in_memory_vertex_lists(void) const {
		for (auto it : in_memory_vertex_lists)
			serialize_vertex_list(it.first);
		in_memory_vertex_lists.clear();
	}

	VariantGraphVertex* VariantGraph::create_vertex(uint64_t id, uint64_t
																									offset, uint64_t length,
																									uint32_t sampleclass_id,
																									const std::vector<VariantGraphVertex::sample_info>&
																									samples) {
		//if (samples.size() != get_popcnt(sampleclass_id)) {
			//console->error("Num of samples: {} is not equal to num of 1s {} in the sample class.",
										 //samples.size(), get_popcnt(sampleclass_id));
			//abort();
		//}
		// check if we need to create a new partition
		if (id == 0 || id % NUM_VERTEXES_IN_BLOCK == 0) {	// id is the number of vertices including the new vertex.
			VariantGraphVertexList *list = new VariantGraphVertexList();
			if (load_type == READ_COMPLETE_GRAPH) {
				vertex_block_list.emplace_back(*list);
			} else { // remove the vertex list with smaller id and add the new vertex list.
				uint32_t block_index = id / NUM_VERTEXES_IN_BLOCK;
				add_or_replace_in_memory_vertex_list(block_index, list);
			}
		}
		uint32_t block_index = id / NUM_VERTEXES_IN_BLOCK;
		// create vertex object and add to vertex_list
		VariantGraphVertex* v;
		if (load_type == READ_COMPLETE_GRAPH) {
			v = vertex_block_list[block_index].add_vertex();
		} else {
			v = in_memory_vertex_lists[block_index]->add_vertex();
		}
		v->set_vertex_id(id);
		v->set_offset(offset);
		v->set_length(length);
		v->set_sampleclass_id(sampleclass_id);
		for (const auto sample : samples) {
			VariantGraphVertex::sample_info* s = v->add_s_info();
			s->set_index(sample.index());
			//s->set_sample_id(sample.sample_id());
			s->set_gt_1(sample.gt_1());
			s->set_gt_2(sample.gt_2());
		}

		return v;
	}

	void VariantGraph::split_vertex(uint64_t vertex_id, uint64_t pos,
																	Graph::vertex* new_vertex) {
		//console->debug("Splitting vertex: {} {}", vertex_id, pos);
		const VariantGraphVertex cur_vertex = get_vertex(vertex_id);
		if (pos > cur_vertex.length()) {
			console->error("Split position is greater than vertex length. {} {} {}",
										 vertex_id, cur_vertex.length(), pos);
			abort();
		}

		// create vertex object and add to vertex_list
		uint64_t offset = cur_vertex.offset() + pos - 1;
		uint64_t length = cur_vertex.length() - pos  + 1;
		VariantGraphVertex::sample_info s;
		s.set_index(cur_vertex.s_info(0).index() + pos - 1);
		//s.set_sample_id(cur_vertex.s_info(0).sample_id());
		s.set_gt_1(cur_vertex.s_info(0).gt_1());
		s.set_gt_2(cur_vertex.s_info(0).gt_2());
		std::vector<VariantGraphVertex::sample_info> samples = {s};
		VariantGraphVertex *v = create_vertex(num_vertices, offset, length,
																					0,
																					samples);

		*new_vertex = v->vertex_id();
		update_idx_vertex_id_map(*v);

		// update length of the seq in the cur_vertex
		get_mutable_vertex(vertex_id)->set_length(cur_vertex.length() - length);

		// move outgoing connections from old_node to the new_node
		for (const auto v : topology.out_neighbors(vertex_id)) {
			validate_ref_path_edge(*new_vertex, v);
			topology.add_edge(*new_vertex, v);
			topology.remove_edge(vertex_id, v);
		}

		// add the edge
		validate_ref_path_edge(vertex_id, *new_vertex);
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

	const VariantGraphVertex& VariantGraph::get_vertex(Graph::vertex id) const {
		uint64_t block_index = id / NUM_VERTEXES_IN_BLOCK;
		uint64_t block_offset = id % NUM_VERTEXES_IN_BLOCK;
		if (load_type == READ_COMPLETE_GRAPH) {
			return vertex_block_list[block_index].vertex(block_offset);
		} else {
			auto itr = in_memory_vertex_lists.find(block_index);
			if (itr == in_memory_vertex_lists.end()) {	// load the vertex list
				load_vertex_list(block_index);
				itr = in_memory_vertex_lists.find(block_index);
			}
			return itr->second->vertex(block_offset);
		}
	}

	VariantGraphVertex* VariantGraph::get_mutable_vertex(Graph::vertex id) {
		uint64_t block_index = id / NUM_VERTEXES_IN_BLOCK;
		uint64_t block_offset = id % NUM_VERTEXES_IN_BLOCK;
		if (load_type == READ_COMPLETE_GRAPH) {
			return vertex_block_list[block_index].mutable_vertex(block_offset);
		} else {
			auto itr = in_memory_vertex_lists.find(block_index);
			if (itr == in_memory_vertex_lists.end()) {	// load the vertex list
				load_vertex_list(block_index);
				itr = in_memory_vertex_lists.find(block_index);
			}
			return itr->second->mutable_vertex(block_offset);
		}
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

	uint64_t VariantGraph::get_num_sample_classes(void) const {
		return sampleclass_map.size();
	}

	//double VariantGraph::get_cache_hit_rate(void) const {
//#ifdef DEBUG_MODE
		//return cache.stats().hit_rate();
//#endif
		//return 0;
	//}

	std::string VariantGraph::get_sample_name(uint32_t id) const {
		auto it = idsample_map.find(id);
		if (it == idsample_map.end()) {
			console->error("Unknown sample id: {}", id);
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

	std::string VariantGraph::print_vertex_info(const VariantGraphVertex& v) const {
		std::string samples;
		for (int i = 0; i < v.s_info_size(); i++) {
			const VariantGraphVertex::sample_info& s = v.s_info(i);
			samples.append("Sample id: " +
										 get_sample_name(get_sample_id(v.sampleclass_id(),
																									 i)));
			samples.append(" Index: " +  std::to_string((int)s.index()) + " ");
		}
		return "ID: " + std::to_string(v.vertex_id()) + " Offset: " +
			std::to_string(v.offset()) + " length: " + std::to_string(v.length()) +
			" sample_class: " + std::to_string(v.sampleclass_id()) + " Samples: " +
			samples; 
	}

	const std::string VariantGraph::get_sequence(const VariantGraphVertex& v)
		const {
			std::string seq;
			for (uint64_t i = v.offset(); i < v.offset() + v.length(); i++) {
				seq += map_int((uint8_t)seq_buffer[i]);
			}
			return seq;
		}

	const std::string VariantGraph::get_sequence(uint64_t start, uint32_t length)
		const {
			std::string seq;
			for (uint64_t i = start; i < start + length; i++) {
				seq += map_int((uint8_t)seq_buffer[i]);
			}
			return seq;
		}

	// the idx map should only contain information about the ref vertices.
	void VariantGraph::update_idx_vertex_id_map(const VariantGraphVertex& v) {
		for (int i = 0; i < v.s_info_size(); i++) {
			const VariantGraphVertex::sample_info& s = v.s_info(i);
			// ref id is 0
			if (get_sample_id(v.sampleclass_id(), i) == 0)
				idx_vertex_id[s.index()] = v.vertex_id();
		}
	}

	bool VariantGraph::get_sample_from_vertex_if_exists(Graph::vertex v,
																											uint32_t sample_id,
																											VariantGraphVertex::sample_info&
																											sample) const {
		const VariantGraphVertex cur_vertex = get_vertex(v);
		for (int i = 0; i < cur_vertex.s_info_size(); i++) {
			const VariantGraphVertex::sample_info& s = cur_vertex.s_info(i);
			if (get_sample_id(cur_vertex.sampleclass_id(), i) == sample_id) {
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
			console->error("Sample not found: {}", sample_id);
		}

		return get_sample_from_vertex_if_exists(v, map_it->second, sample);	
	}

	uint64_t VariantGraph::find_sample_index(Graph::vertex ref_v_id,
																					 uint32_t sample_id) const {
		Graph::vertex cur_vertex_id = ref_v_id;
		const VariantGraphVertex cur_vertex = get_vertex(cur_vertex_id);
		// check if the cur vertex contains the sample
		VariantGraphVertex::sample_info sample;
		if (get_sample_from_vertex_if_exists(cur_vertex_id, sample_id, sample))
			return sample.index() + cur_vertex.length();

		// check in the cache
		//if (cache.contains(sample_id))
			//cur_vertex_id = cache.lookup(sample_id);
		//else {
			// else traverse back in the graph to find a vertex with @sample_id
			// ref id is 0
			get_sample_from_vertex_if_exists(ref_v_id, 0, sample);
			uint64_t cur_index = sample.index();
			auto temp_itr = idx_vertex_id.lower_bound(cur_index);
			--temp_itr;
			while (cur_index > 1) {
				// move to prev vertex.
				cur_index = temp_itr->first; 
				cur_vertex_id = temp_itr->second;
				//cur_vertex = get_vertex(cur_vertex_id);

				if (get_sample_from_vertex_if_exists(cur_vertex_id, sample_id, sample))
					break;

				Graph::vertex neighbor_vertex_id;
				if (get_neighbor_vertex(cur_vertex_id, sample_id, &neighbor_vertex_id)) {
					cur_vertex_id = neighbor_vertex_id;
					//cur_vertex = get_vertex(neighbor_vertex_id);
					break;
				}
				--temp_itr;
			}
		//}

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
			console->error("Sample id not found: {}", sample_id);
		}
		auto it = this->find(cur_vertex_id, map_it->second);
		while (!it.done() && (*it)->vertex_id() != ref_v_id) {
			cur_distance += (*it)->length();
			++it;
		}

		return cur_distance += get_vertex(ref_v_id).length();
	}

	// if there is not neighbor with @sample_id then set @v to "ref"
	bool VariantGraph::get_neighbor_vertex(Graph::vertex id, uint32_t sample_id,
																				 Graph::vertex* v) const {
		uint32_t min_idx = UINT32_MAX;
		for (const auto v_id : topology.out_neighbors(id)) {
			const VariantGraphVertex vertex = get_vertex(v_id);
			for (int i = 0; i < vertex.s_info_size(); i++) {
				const VariantGraphVertex::sample_info& s = vertex.s_info(i);
				uint32_t s_id = get_sample_id(vertex.sampleclass_id(), i);
				if (s_id != 0 && s_id == sample_id) {
					*v = v_id; 
					return true;
				} else if (s_id == 0) {	// if sample_id is not found follow "ref" path
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
																					sample_idx, bool gt1, bool gt2) {
		VariantGraphVertex::sample_info* s =
			get_mutable_vertex(id)->add_s_info();
		s->set_index(sample_idx);
		//s->set_sample_id(sample_id);
		s->set_gt_1(gt1);
		s->set_gt_2(gt2);
	}

	bool VariantGraph::check_if_mutation_exists(Graph::vertex prev,
																							Graph::vertex next, const
																							std::string alt,
																							Graph::vertex* v) const {
		for (const auto v_id : topology.out_neighbors(prev)) {
			if (topology.is_edge(std::make_pair(v_id, next))) {
				const VariantGraphVertex neighbor = get_vertex(v_id);
				const std::string seq = get_sequence(neighbor);
				if (seq == alt) {
					*v = v_id;
					return true;
				}
			}
		}
		return false;
	}

	void VariantGraph::validate_ref_path_edge(Graph::vertex src, Graph::vertex
																						dest) const {
		VariantGraphVertex::sample_info src_sample, dest_sample;
		if (!get_sample_from_vertex_if_exists(src, 0, src_sample) &&
				!get_sample_from_vertex_if_exists(dest, 0, dest_sample)) {
			if (src_sample.index() >= dest_sample.index()) {
				console->error("Source ref index is not smaller than dest ref index {}:{} {}:{}",
											 src, src_sample.index(), dest, dest_sample.index());
				abort();
			}
		}
	}

	bool VariantGraph::check_if_mutation_exists(Graph::vertex prev,
																							uint64_t offset, uint64_t length,
																							Graph::vertex* v) const {
		for (const auto v_id : topology.out_neighbors(prev)) {
			const VariantGraphVertex neighbor = get_vertex(v_id);
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

		//console->debug("Adding mutation: {} {} {} {} {}",
									 //mutation_string(mutation), ref, alt, pos,
									 //sample_list.size());
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
		VariantGraphVertex ref_vertex = get_vertex(ref_vertex_id);

		if (mutation == SUBSTITUTION) {
			Graph::vertex prev_ref_vertex_id = 0, next_ref_vertex_id = 0;
			// Either we got the vertex where there's already a substition 
			// or the vertex contains the seq with @pos
			if (ref_vertex_idx == pos && ref_vertex.length() == ref.size()) {
				// splitting not needed. mutation spans the whole vertex.
				// find the prev ref vertex
				auto temp_itr = idx_vertex_id.lower_bound(ref_vertex_idx);
				if (temp_itr->first != ref_vertex_idx) {
					console->error("Vertex id not found in the map");
					abort();
				}
				--temp_itr;
				prev_ref_vertex_id = temp_itr->second;

				// find the next ref vertex
				get_neighbor_vertex(ref_vertex_id, 0, &next_ref_vertex_id);		
			} else if (ref_vertex_idx == pos && ref_vertex.length() > ref.size()) {
				// vertex needs to be split into two. mutation is contained in the
				// vertex.
				// find the prev ref vertex
				auto temp_itr = idx_vertex_id.lower_bound(ref_vertex_idx);
				if (temp_itr->first != ref_vertex_idx) {
					console->error("Vertex id not found in the map");
					abort();
				}
				--temp_itr;
				prev_ref_vertex_id = temp_itr->second;

				// split the vertex
				split_vertex(ref_vertex_id, ref.size() + 1, &next_ref_vertex_id);
			} else if (ref_vertex_idx == pos && ref_vertex.length() < ref.size()) {
				// splitting needed. mutation spans one or more vertexes.
				// find the next ref vertex which contains the pos + ref.size()
				VariantGraphVertex next_ref_vertex;
				auto temp_itr = idx_vertex_id.lower_bound(ref_vertex_idx);
				do {
					++temp_itr;
					next_ref_vertex = get_vertex(temp_itr->second);
				} while (temp_itr->first + next_ref_vertex.length() < pos + ref.size());
				if (temp_itr->first + next_ref_vertex.length() == pos +
						ref.size()) {
					get_neighbor_vertex(temp_itr->second, 0, &next_ref_vertex_id);
				} else {
					// split the vertex
					split_vertex(temp_itr->second, pos + ref.size() - temp_itr->first + 1,
											 &next_ref_vertex_id);
				}

				// find the prev vertex.
				temp_itr = idx_vertex_id.lower_bound(ref_vertex_idx);
				if (temp_itr->first != ref_vertex_idx) {
					console->error("Vertex id not found in the map");
					abort();
				}
				--temp_itr;
				prev_ref_vertex_id = temp_itr->second;
			} else if (ref_vertex_idx < pos && ref_vertex_idx + ref_vertex.length()
								 > pos + ref.size()) {
				// vertex needs to be split into three. mutation contained in the
				// vertex.
				// split the vertex
				uint64_t split_pos = pos - ref_vertex.offset();
				split_vertex(ref_vertex_id, split_pos, split_pos + ref.size(),
										 &prev_ref_vertex_id, &next_ref_vertex_id);
				std::swap(ref_vertex_id, prev_ref_vertex_id);
			} else if (ref_vertex_idx < pos && ref_vertex_idx + ref_vertex.length()
								 < pos + ref.size()) {
				// split the current vertex to get prev_vertex. mutation spans one or
				// more vertexes.
				prev_ref_vertex_id = ref_vertex_id;
				split_vertex(prev_ref_vertex_id, pos - ref_vertex_idx + 1,
										 &ref_vertex_id);

				// find the next ref vertex
				VariantGraphVertex next_ref_vertex;
				auto temp_itr = idx_vertex_id.lower_bound(ref_vertex_idx);
				do {
					++temp_itr;
					next_ref_vertex = get_vertex(temp_itr->second);
				} while (temp_itr->first + next_ref_vertex.length() < pos + ref.size());
				if (temp_itr->first + next_ref_vertex.length() == pos +
						ref.size()) {
					get_neighbor_vertex(temp_itr->second, 0, &next_ref_vertex_id);
				} else {
					// split the vertex
					split_vertex(temp_itr->second, pos + ref.size() - temp_itr->first + 1,
											 &next_ref_vertex_id);
				}
			} else if (ref_vertex_idx < pos && ref_vertex_idx + ref_vertex.length()
								 == pos + ref.size()) { // split the current vertex to get prev_vertex
				prev_ref_vertex_id = ref_vertex_id;
				split_vertex(prev_ref_vertex_id, pos - ref_vertex_idx + 1,
										 &ref_vertex_id);	

				// find the next ref vertex
				get_neighbor_vertex(ref_vertex_id, 0, &next_ref_vertex_id);		
			}
			// create a vertex for the mutation using the first sample from the
			// list.
			//uint32_t sample_id = sample_list[0].sample_id;
			bool gt1 = sample_list[0].gt1;
			bool gt2 = sample_list[0].gt2;
			//uint64_t sample_idx = find_sample_index(prev_ref_vertex_id, sample_id);
			uint32_t sampleclass_id = find_sample_vector_or_add(sample_list);
			VariantGraphVertex* sample_vertex = add_vertex(alt, 0,
																										 sampleclass_id, gt1, gt2);
			// make connections for the new vertex in the graph
			topology.add_edge(prev_ref_vertex_id, sample_vertex->vertex_id());
			topology.add_edge(sample_vertex->vertex_id(), next_ref_vertex_id);

			// remove the first sample from the list
			sample_list.erase(sample_list.begin());

			// add rest of the samples to the vertex.
			for (auto sample : sample_list) {
				//sample_idx = find_sample_index(prev_ref_vertex_id,
				//sample.sample_id);
				add_sample_to_vertex(sample_vertex->vertex_id(), 0, sample.gt1,
														 sample.gt2);
			}
			// validate popcnt and s_info size.
			if ((uint32_t)sample_vertex->s_info_size() !=
					get_popcnt(sample_vertex->sampleclass_id())) {
				console->error("Num of samples: {} is not equal to num of 1s: {} in the sample class.",
											 sample_vertex->s_info_size(),
											 get_popcnt(sampleclass_id));
				abort();
			}
		} else if (mutation == INSERTION) {
			Graph::vertex prev_ref_vertex_id = 0, next_ref_vertex_id = 0;
			// Either we got the vertex where there's already an insertion
			// or the vertex contains the seq with @pos
			if (ref_vertex_idx == pos) { // splitting not needed
				// find the prev ref vertex
				auto temp_itr = idx_vertex_id.lower_bound(ref_vertex_idx);
				if (temp_itr->first != ref_vertex_idx) {
					console->error("Vertex id not found in the map");
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
			} else if (ref_vertex_idx + ref_vertex.length() == pos) { // splitting not needed
				prev_ref_vertex_id = ref_vertex_id;

				// find the next ref vertex
				get_neighbor_vertex(ref_vertex_id, 0, &next_ref_vertex_id);		
			} else { // to handle insertions after ref seq length. 
				prev_ref_vertex_id = ref_vertex_id;
			}
			// create a vertex for the mutation using the first sample from the
			// list.
			//uint32_t sample_id = sample_list[0].sample_id;
			bool gt1 = sample_list[0].gt1;
			bool gt2 = sample_list[0].gt2;
			//uint64_t sample_idx = find_sample_index(prev_ref_vertex_id, sample_id);
			uint32_t sampleclass_id = find_sample_vector_or_add(sample_list);
			VariantGraphVertex* sample_vertex = add_vertex(alt, 0,
																										 sampleclass_id, gt1, gt2);
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
				add_sample_to_vertex(sample_vertex->vertex_id(), 0, sample.gt1,
														 sample.gt2);
			}
			// validate popcnt and s_info size.
			if ((uint32_t)sample_vertex->s_info_size() !=
					get_popcnt(sample_vertex->sampleclass_id())) {
				console->error("Num of samples: {} is not equal to num of 1s: {} in the sample class.",
											 sample_vertex->s_info_size(),
											 get_popcnt(sampleclass_id));
				abort();
			}
		} else { // it's deletion
			Graph::vertex prev_ref_vertex_id = 0, next_ref_vertex_id = 0;
			// Either we got the vertex where there's already a deletion
			// or the vertex contains the seq with @pos
			if (ref_vertex_idx == pos && ref_vertex.length() == ref.size()) {
				// splitting not needed
				// find the prev ref vertex
				auto temp_itr = idx_vertex_id.lower_bound(ref_vertex_idx);
				if (temp_itr->first != ref_vertex_idx) {
					console->error("Vertex id not found in the map");
					abort();
				}
				--temp_itr;
				prev_ref_vertex_id = temp_itr->second;

				// find the next ref vertex
				get_neighbor_vertex(ref_vertex_id, 0, &next_ref_vertex_id);	
			} else if (ref_vertex_idx == pos && ref_vertex.length() > ref.size()) { 
				// vertex needs to be split into two
				// find the prev ref vertex
				auto temp_itr = idx_vertex_id.lower_bound(ref_vertex_idx);
				if (temp_itr->first != ref_vertex_idx) {
					console->error("Vertex id not found in the map");
					abort();
				}
				--temp_itr;
				prev_ref_vertex_id = temp_itr->second;

				// split the vertex
				split_vertex(ref_vertex_id, ref.size() + 1, &next_ref_vertex_id);
			} else if (ref_vertex_idx == pos && ref_vertex.length() < ref.size()) {
				// deletion will skip 1 or more vertexes
				// find the next ref vertex
				VariantGraphVertex next_ref_vertex;
				auto temp_itr = idx_vertex_id.lower_bound(ref_vertex_idx);
				do {
					++temp_itr;
					next_ref_vertex = get_vertex(temp_itr->second);
				} while (temp_itr->first + next_ref_vertex.length() < pos + ref.size());
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
					console->error("Vertex id not found in the map");
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
					next_ref_vertex = get_vertex(temp_itr->second);
				} while (temp_itr->first + next_ref_vertex.length() < pos + ref.size());
				if (temp_itr->first + next_ref_vertex.length() == pos + ref.size()) {
					get_neighbor_vertex(temp_itr->second, 0, &next_ref_vertex_id);
				} else {
					// split the vertex
					split_vertex(temp_itr->second, pos + ref.size() - temp_itr->first + 1,
											 &next_ref_vertex_id);
				}
			} else if (ref_vertex_idx < pos && ref_vertex_idx + ref_vertex.length()
								 == pos + ref.size()) { // split the current vertex to get prev_vertex
				prev_ref_vertex_id = ref_vertex_id;
				split_vertex(prev_ref_vertex_id, pos - ref_vertex_idx + 1,
										 &ref_vertex_id);	

				// find the next ref vertex
				get_neighbor_vertex(ref_vertex_id, 0, &next_ref_vertex_id);	
			}
			// update sample class.
			// if it's a new ref vertex update sample class.
			if (!update_vertex_sample_class(next_ref_vertex_id, sample_list)) {
				console->error("Unsupported mutation: {} {} {} {} {}",
											 mutation_string(mutation), ref, alt, pos,
											 sample_list.size());
			}
			// validate popcnt and s_info size.
			const VariantGraphVertex next_ref_vertex =
				get_vertex(next_ref_vertex_id);
			if ((uint32_t)next_ref_vertex.s_info_size() !=
					get_popcnt(next_ref_vertex.sampleclass_id())) {
				console->error("Num of samples: {} is not equal to num of 1s: {} in the sample class.",
											 next_ref_vertex.s_info_size(),
											 get_popcnt(next_ref_vertex.sampleclass_id()));
				abort();
			}
			// make connections for the new vertex in the graph
			validate_ref_path_edge(prev_ref_vertex_id, next_ref_vertex_id);
			topology.add_edge(prev_ref_vertex_id, next_ref_vertex_id);
		}
	}

	void VariantGraph::fix_sample_indexes(void) {
#if 0
		// Naive solution.
		// Traverse sample paths for each sample in the database.
		// Update sample indexes during path traversals.
		for (const auto sample : sampleid_map) {
			if (sample.first == "ref")
				continue;
			uint32_t cur_index = 1;
			auto path = mutable_find(sample.first);
			while (!path.done()) {
				auto vertex = *path;
				for (int i = 0; i < vertex->s_info_size(); ++i) {
					VariantGraphVertex::sample_info* s = vertex->mutable_s_info(i);
					if (get_sample_id(vertex->sampleclass_id(), i) == sample.second) {
						s->set_index(cur_index);
					}
				}
				cur_index += vertex->length();
				++path;
			}
		}
#else
		// Optimized solution.
		// Initialize delta for each sample to 0.
		// Perform a BFT on the graph.
		// For each vertex:
		// If it is a ref vertex, for each neighbor
		// 		If the sample index is 0
		// 				update the sample index
		// 				sample_index = vertex_index + vertex_len + delta.
		//		Else if sample index is not 0
		// 				update delta
		// 				sample_delta = sample_index - (ref_index + vertex_len)
		// If it is a sample vertex
		// 		For each sample in the vertex
		// 			sample_delta = (sample_index + vertex_len) - neighbour_vertex_ref_index
		//
		console->info("Fixing sample indexes in the graph");
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
					auto cur_neighbor = get_mutable_vertex(neighbor_id);
					// for each sample other than "ref" in the vertex update the index
					for (int i = 0; i < cur_neighbor->s_info_size(); ++i) {
						VariantGraphVertex::sample_info* s = cur_neighbor->mutable_s_info(i);
						uint32_t s_id = get_sample_id(cur_neighbor->sampleclass_id(), i);
						if (s_id != 0 && s->index() == 0) {	// sample index is 0. Update index
							auto map_it = sampleid_delta.find(s_id);
							if (map_it == sampleid_delta.end()) {
								console->error("Unknown sample id: {}", s_id);
								abort();
							}
							int32_t delta = map_it->second;
							int32_t sample_index = ref_index + cur_vertex->length() + delta;
							if (sample_index < 0) {
								console->error("Sample index is less than 0: {} {}",
															 get_sample_name(s_id),
															 cur_vertex->vertex_id());
								abort();
							}
							s->set_index(sample_index);
						} else if (s_id != 0 && s->index() != 0) {	// Sample index is not 0. Update delta
							sampleid_delta[s_id] =  s->index() - (ref_index +
																										cur_vertex->length());
						}
					}
				}
			} else { // this is not a ref vertex
				// Sample vertexes should only have one outgoing edge.
				if (topology.out_neighbors(cur_vertex->vertex_id()).size() > 1) {
					console->error("Sample vertex has more than 1 neighbor: {}",
												 cur_vertex->vertex_id());
					abort();
				}
				for (auto neighbor_id :
						 topology.out_neighbors(cur_vertex->vertex_id())) {	// should be a single iteration loop.
					auto cur_neighbor = get_vertex(neighbor_id);
					VariantGraphVertex::sample_info ref_sample;
					if (!get_sample_from_vertex_if_exists(cur_neighbor.vertex_id(), 0,
																								ref_sample)) { // this is a ref vertex
						console->error("Ref vertex not found as a neighbor from sample vertex. {}",
													 cur_neighbor.vertex_id());
						abort();
					}
					// for each sample other than "ref" in the vertex update the index
					uint32_t cur_length = cur_vertex->length();
					// for each sample update the delta
					for (int i = 0; i < cur_vertex->s_info_size(); ++i) {
						const VariantGraphVertex::sample_info s = cur_vertex->s_info(i);
						uint32_t cur_index = s.index();
						uint32_t s_id = get_sample_id(cur_vertex->sampleclass_id(), i);
						sampleid_delta[s_id] = cur_index + cur_length - ref_sample.index();
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
		cur = vg->get_vertex(v);
		auto it = vg->sampleid_map.find(sample_id);
		if (it == vg->sampleid_map.end()) {
			console->error("Sample not found");
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
		cur = vg->get_vertex(next_vertex);
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
		cur = vg->get_mutable_vertex(v);
		auto it = vg->sampleid_map.find(sample_id);
		if (it == vg->sampleid_map.end()) {
			console->error("Sample not found");
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
		cur = vg->get_mutable_vertex(next_vertex);
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
			return &vg->get_vertex(*itr);
		}

	void VariantGraph::VariantGraphIterator::operator++(void) {
		++itr;
	}

	bool VariantGraph::VariantGraphIterator::done(void) const {
		return itr.done();
	}

	VariantGraph::VariantGraphIterator VariantGraph::find(Graph::vertex v,
																												uint64_t radius) const
	{
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
			return vg->get_mutable_vertex(*itr);
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
