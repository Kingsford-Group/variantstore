/*
 * ============================================================================
 *
 *       Filename:  index.h
 *
 *         Author:  Prashant Pandey (), ppandey2@cs.cmu.edu
 										Yinjie Gao <yinjieg@andrew.cmu.edu>
 *   Organization:  Carnegie Mellon University
 *
 * ============================================================================
 */


#ifndef __INDEX_H__
#define __INDEX_H_

#include "variant_graph.h"
#include <sdsl/bit_vectors.hpp>
#include <set>

uint64_t REF_GENOME_LEN = 3099706404;
uint16_t BLOCK_SIZE = 64

namespace variantdb {
	class Index {
		// A static index is build based on given variant graph
		// An Index supports query node(in reference) at a given position
		// Adding nodes in reference change the index content

	private:
		//bit_vector bound_vec;
		sdsl::rrr_vector<BLOCK_SIZE> rrrb;
		sdsl::rrr_vector<BLOCK_SIZE>::rank_1_type rank_rrrb;
		std::vector<uint64_t> node_list;

	public:
		// Construct Index given Variant Graph
		Index(const VariantGraph vg);
		~Index();
		// Query node_id given position
		uint64_t find(uint64_t pos);

		void serialize(const std::string prefix);
	};


	Index::Index(const VariantGraph vg)
	{
		// Construct a bit vector of size = genome length
		sdsl::bit_vector b(REF_GENOME_LEN, 0);

		uint64_t node_id; // ?Need node id here?
		uint64_t idx;

		// Iterate nodes folloing path in REF and modify bit vector & node list
		// doneend() for iterator
		// ?Reference takes the sample id 0?
		for (VariantGraphIterator node_it = vg.find(0); node_it != vg.done(0); node_it++)
		{
			node_id; // ??
			idx = node_it->index;
			b[idx] = 1;
			node_list.push_back(node_id);
		}

		// Compress it & Construct rank support vector from vector
		sdsl::util::assign(rrrb, sdsl::rrr_vector<BLOCK_SIZE>(b));
		sdsl::util::assign(rank_rrrb, sdsl::rrr_vector<BLOCK_SIZE>::rank_1_type(&rrrb));
		return;
	}

	uint64_t Index::find(uint64_t pos)
	{
		uint64_t node_idx = rank_rrrb(pos);
		return node_list[node_idx];
	}

}



#endif // __INDEX_H__
