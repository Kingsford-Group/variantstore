/*
 * ============================================================================
 *
 *       Filename:  index.h
 *
 *         Author:  Prashant Pandey <ppandey@cs.stonybrook.edu>
 										Yinjie Gao <yinjieg@andrew.cmu.edu>
 *   Organization:  Carnegie Mellon University
 *
 * ============================================================================
 */

#ifndef __INDEX_H__
#define __INDEX_H_

#include "variant_graph.h"
#include "variantgraphvertex.pb.h"

#include <sdsl/bit_vectors.hpp>
#include <sdsl/int_vector.hpp>
#include <sdsl/util.hpp>
#include <vector>

const uint16_t BLOCK_SIZE = 127;

namespace variantdb {
	class Index {
		// A static index is build based on given variant graph
		// An Index supports query node(in reference) at a given position
		// Adding nodes in reference change the index content

	private:
		//bit_vector bound_vec;
		sdsl::rrr_vector<BLOCK_SIZE> rrrb;
		sdsl::rrr_vector<BLOCK_SIZE>::rank_1_type rank_rrrb;
		sdsl::int_vector<> node_list;

	public:
		// Construct Index given Variant Graph
		Index(const VariantGraph *vg);
		// Read Index from disk
		Index(const std::string filename);
		~Index();
		// Query node_id given position
		uint64_t find(uint64_t pos);

		void serialize(const std::string prefix);
	};


	Index::Index(const VariantGraph *vg)
	{
		// Construct a bit vector of size = genome length
		sdsl::bit_vector b(vg->get_ref_length(), 0);
		// Construct a int vector of) node_list_sz
		sdsl::util::assign(node_list, sdsl::int_vector<>(0, 0, 64));
		uint64_t node_list_sz = 0;
		// Iterate nodes folloing path in REF
		// modify bit vector & node list
		VariantGraph::VariantGraphPathIterator it = vg->find(0); //  Iterate through reference
		while(!it.done())
		{
			node_list_sz++;
			if (node_list_sz > node_list.size()) {node_list.resize(node_list_sz);}
			uint64_t node_id = (*it)->vertex_id();
			uint64_t idx = (*it)->offset();
			//uint64_t idx = (*node_it).s_info(0).index();
			b[idx] = 1;
			node_list[node_list_sz-1] = node_id;
			++it;
		}

		// Compress it & Construct rank support vector from vector
		sdsl::util::assign(rrrb, sdsl::rrr_vector<BLOCK_SIZE>(b));
		sdsl::util::assign(rank_rrrb,
			sdsl::rrr_vector<BLOCK_SIZE>::rank_1_type(&rrrb));
		sdsl::util::bit_compress(node_list);
		return;
	} // Index(const VariantGraph vg)

	Index::Index(const std::string filename)
	{
		std::ifstream in(filename);
		rrrb.load(in);
		sdsl::util::assign(rank_rrrb,
			sdsl::rrr_vector<BLOCK_SIZE>::rank_1_type(&rrrb));
    node_list.load(in);
		return;
	} // Index(const std::string filename)

	uint64_t Index::find(uint64_t pos)
	{
		uint64_t node_idx = rank_rrrb(pos);
		return node_list[node_idx];
	} // find(uint64_t pos)

	void Index::serialize(const std::string prefix)
	{
		std::ofstream out(prefix + ".sdsl");
		sdsl::serialize(rrrb, out);
		sdsl::serialize(node_list, out);
		return;
	} // serialize(const std::string prefix)

	Index::~Index(){}

}

#endif // __INDEX_H__
