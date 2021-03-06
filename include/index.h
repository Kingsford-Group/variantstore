/*
 * ============================================================================
 *
 *       Filename:  index.h
 *
 *         Author:  Prashant Pandey <ppandey2@cs.cmu.edu>
 *									Yinjie Gao <yinjieg@andrew.cmu.edu>
 *   Organization:  Carnegie Mellon University
 *
 * ============================================================================
 */

#ifndef __INDEX_H__
#define __INDEX_H__

#include <sdsl/bit_vectors.hpp>
#include <sdsl/int_vector.hpp>
#include <sdsl/util.hpp>
#include <vector>

#include "dot_graph.h"
#include "util.h"
#include "variant_graph.h"

namespace variantstore {
	class Index {
		// A static index is build based on given variant graph
		// An Index supports query node(in reference) at a given position
		// Adding nodes in reference change the index content

	public:
		// Construct Index given Variant Graph
		Index(const VariantGraph *vg);
		// Read Index from disk
		Index(const std::string prefix);
		~Index();
		// Query node_id given position
		Graph::vertex find(const uint64_t pos) const;
		Graph::vertex find(const uint64_t pos, uint64_t &ref_node_rank) const;
		bool is_empty(const uint64_t pos_x, const uint64_t pos_y) const;
		Graph::vertex previous(const uint64_t ref_node_rank) const;
		void serialize(const std::string prefix);

	private:
		//bit_vector bound_vec;
		sdsl::rrr_vector<SDSL_BITVECTOR_BLOCK_SIZE> rrrb;
		sdsl::rrr_vector<SDSL_BITVECTOR_BLOCK_SIZE>::rank_1_type rank_rrrb;
		sdsl::rrr_vector<SDSL_BITVECTOR_BLOCK_SIZE>::select_1_type select_rrrb;
		sdsl::int_vector<> node_list;
	};


	Index::Index(const VariantGraph *vg)
	{
		// Construct a bit vector of size = genome length
		sdsl::bit_vector b(vg->get_ref_length(), 0);
		// Construct a int vector of node_list_sz
		uint64_t init_sz = 1;
		sdsl::util::assign(node_list, sdsl::int_vector<>(init_sz, 0,
																										 sizeof(Graph::vertex)*8));
		uint64_t node_list_sz = 0;
		// Iterate nodes folloing path in REF
		// modify bit vector & node list
		VariantGraph::VariantGraphPathIterator it = vg->find("ref"); //  Iterate through reference
		while(!it.done())
		{
			uint64_t node_id = (*it)->vertex_id();
			Graph::vertex v = node_id;
			VariantGraphVertex::sample_info sample;
			// ref sample id is 0
			if (!vg->get_sample_from_vertex_if_exists(v, "ref", sample)) {
				console->error("Ref sample not found in the vertex: {}", v);
				abort();
			}
			uint64_t idx = 	sample.index();

			//console->debug("At index {} has node {}", idx, node_id);

			// extra check to determine loops in the ref path through the graph.
			//if (b[idx-1] == 1) {
				//console->error("Found a loop in the ref path at node: {} and index: {}"
											 //, node_id, idx - 1);
				//createDotGraph(vg, "./", node_id-5, 10);
				//abort();
			//}
				if (b[idx - 1] != 1) { // only update if not already set.
				node_list_sz++;
				if (node_list_sz > node_list.size()) {node_list.resize(node_list_sz);}
				node_list[node_list_sz - 1] = node_id;
			}
			b[idx-1] = 1;
			++it;
		}
		
		//std::cout << b << std::endl;
		//std::cout << node_list << std::endl;

		// Compress it & Construct rank support vector from vector
		sdsl::util::assign(rrrb, sdsl::rrr_vector<SDSL_BITVECTOR_BLOCK_SIZE>(b));
		sdsl::util::assign(rank_rrrb,
			sdsl::rrr_vector<SDSL_BITVECTOR_BLOCK_SIZE>::rank_1_type(&rrrb));
		sdsl::util::assign(select_rrrb,
			sdsl::rrr_vector<SDSL_BITVECTOR_BLOCK_SIZE>::select_1_type(&rrrb));
		sdsl::util::bit_compress(node_list);
		return;
	} // Index(const VariantGraph vg)

	Index::Index(const std::string prefix)
	{
		sdsl::load_from_file(rrrb, prefix + "/index.sdsl");
		sdsl::load_from_file(node_list, prefix + "/ref_node_id.sdsl");
		sdsl::util::assign(rank_rrrb,
			sdsl::rrr_vector<SDSL_BITVECTOR_BLOCK_SIZE>::rank_1_type(&rrrb));
		sdsl::util::assign(select_rrrb,
			sdsl::rrr_vector<SDSL_BITVECTOR_BLOCK_SIZE>::select_1_type(&rrrb));
		return;
	} // Index(const std::string filename)

	Graph::vertex Index::find(const uint64_t pos) const
	{
		if (pos < 1) {
			console->error("Can't find node corresponding to pos {}", pos);
			abort();
		}
		if ( pos >= rank_rrrb.size())
			return node_list[node_list.size()-1];

		uint64_t node_idx = rank_rrrb(pos);
		if ( node_idx == 0 )
			return node_list[0];

		return node_list[node_idx - 1];
	} // find(uint64_t pos)

	Graph::vertex Index::find(const uint64_t pos, uint64_t &ref_node_rank) const
	{
		if ( pos >= rank_rrrb.size())
		{
			ref_node_rank = node_list.size()-1;
			return node_list[node_list.size()-1];
		}
		uint64_t node_idx = rank_rrrb(pos);
		if ( node_idx == 0 )
			return node_list[0];

		ref_node_rank = node_idx - 1;
		return node_list[node_idx-1];
	} // find(uint64_t pos)

	bool Index::is_empty(const uint64_t pos_x, const uint64_t pos_y) const {
		if (pos_x < 1) {
			console->error("Can't find node corresponding to pos {}", pos_x);
			abort();
		}
		if (pos_x > rank_rrrb.size())
			return true;

		uint64_t pos_x_rank = rank_rrrb(pos_x);
		uint64_t index_x = select_rrrb(pos_x_rank);
		uint64_t index_y = select_rrrb(pos_x_rank + 1);

		if (index_x <= pos_x && index_y <= pos_y)
			return false;

		return true;
	}

	Graph::vertex Index::previous(const uint64_t ref_node_rank) const
	{
		if (ref_node_rank == 0) { return node_list[0]; }
		return node_list[ref_node_rank-1];
	}

	void Index::serialize(const std::string prefix)
	{
		sdsl::store_to_file(rrrb, prefix + "/index.sdsl");
		sdsl::store_to_file(node_list, prefix + "/ref_node_id.sdsl");
		return;
	} // serialize(const std::string prefix)

	Index::~Index(){}

}

#endif // __INDEX_H__
