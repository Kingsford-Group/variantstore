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
#include "graph.h"
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

	public:
		// Construct Index given Variant Graph
		Index(const VariantGraph *vg);
		// Read Index from disk
		Index(const std::string prefix);
		~Index();
		// Query node_id given position
		uint64_t find(uint64_t pos) const;
		void serialize(const std::string prefix);

	private:
		//bit_vector bound_vec;
		sdsl::rrr_vector<BLOCK_SIZE> rrrb;
		sdsl::rrr_vector<BLOCK_SIZE>::rank_1_type rank_rrrb;
		sdsl::int_vector<> node_list;
	};


	Index::Index(const VariantGraph *vg)
	{
		// Construct a bit vector of size = genome length
		sdsl::bit_vector b(vg->get_ref_length(), 0);
		// Construct a int vector of node_list_sz
		uint64_t init_sz = 1;
		sdsl::util::assign(node_list, sdsl::int_vector<>(init_sz, 0,
																										 sizeof(Graph::vertex)));
		uint64_t node_list_sz = 0;
		// Iterate nodes folloing path in REF
		// modify bit vector & node list
		VariantGraph::VariantGraphPathIterator it = vg->find("ref"); //  Iterate through reference
		while(!it.done())
		{
			node_list_sz++;
			if (node_list_sz > node_list.size()) {node_list.resize(node_list_sz);}
			uint64_t node_id = (*it)->vertex_id();
			Graph::vertex v = node_id;
			VariantGraphVertex::sample_info sample;
			if (!vg->get_sample_from_vertex_if_exists(v, "ref", sample)) {
				ERROR("Ref sample not found in the vertex: " << v);
				abort();
			}
			uint64_t idx = 	sample.index();

			DEBUG("At index " << idx
								<< " has node " << node_id);

			b[idx] = 1;
			node_list[node_list_sz-1] = node_id;
			++it;
		}

		std::cout << std::endl;
		// Compress it & Construct rank support vector from vector
		sdsl::util::assign(rrrb, sdsl::rrr_vector<BLOCK_SIZE>(b));
		sdsl::util::assign(rank_rrrb,
			sdsl::rrr_vector<BLOCK_SIZE>::rank_1_type(&rrrb));
		sdsl::util::bit_compress(node_list);
		return;
	} // Index(const VariantGraph vg)

	Index::Index(const std::string prefix)
	{
		sdsl::load_from_file(rrrb, prefix + "/index.sdsl");
		sdsl::load_from_file(node_list, prefix + "/ref_node_id.sdsl");
		sdsl::util::assign(rank_rrrb,
			sdsl::rrr_vector<BLOCK_SIZE>::rank_1_type(&rrrb));
		return;
	} // Index(const std::string filename)

	uint64_t Index::find(uint64_t pos) const
	{
		if ( pos >= rank_rrrb.size())
			return node_list[node_list.size()-1];

		uint64_t node_idx = rank_rrrb(pos);
		if ( node_idx == 0 )
			return node_list[0];

		return node_list[node_idx-1];
	} // find(uint64_t pos)

	void Index::serialize(const std::string prefix)
	{
		sdsl::store_to_file(rrrb, prefix + "/index.sdsl");
		sdsl::store_to_file(node_list, prefix + "/ref_node_id.sdsl");
		return;
	} // serialize(const std::string prefix)

	Index::~Index(){}

}

#endif // __INDEX_H__
