/*
 * ============================================================================
 *
 *       Filename:  query.h
 *
 *         Author:  Prashant Pandey <ppandey@cs.stonybrook.edu>
 										Yinjie Gao <yinjieg@andrew.cmu.edu>
 *   Organization:  Carnegie Mellon University
 *
 * ============================================================================
 */
#include <string>
#include <cassert>
#include "index.h"
#include "variant_graph.h"
#include "vcflib/Variant.h"
#include "variantgraphvertex.pb.h"

std::string REF = "ref";

namespace variantdb {
// TODO: check validity of inputs

// traverse back until a node containg sample_id is found
// return such node, record the pos in ref and sample
Graph::vertex get_prev_vertex_with_sample ( VariantGraph *vg, Index *idx,
                                            const uint64_t pos,
                                            const string sample_id,
                                            uint64_t &ref_pos,
                                            uint64_t &sample_pos) {
  uint64_t cur_pos = pos;
  Graph::vertex v = idx->find(cur_pos); // the node right before cur_pos
  Graph::vertex v_find = v;
  VariantGraphVertex::sample_info sample;
  VariantGraphVertex::sample_info sample_find;
  bool sample_found = false;

  while (true) { // do until find sample or to the front of graph
    Graph::vertex cur_v = v;
    Graph::vertex new_v = cur_v;

    // find the node before cur_v in the reference
    while (new_v == cur_v && cur_pos > 1) {
      cur_pos = cur_pos - 1;
      new_v = idx->find(cur_pos);
    }
    v = new_v;

    // if reach the first node in the graph, return the first node
    if (cur_pos <= 1) {
      ref_pos = 1;
      v_find = v;
      assert(v_find == 0);
      vg->get_sample_from_vertex_if_exists(v_find, REF, sample_find);
      sample_pos = sample_find.index();
      break;
    }

    // check all outgoing edges of v by BFS, report of a node containing sample
    // is found. ref_pos is updated whenever a node with ref is encountered
    VariantGraph::VariantGraphIterator it = vg->find(v, 1);
    ++it;
    while(!it.done()) {
      v = (*it)->vertex_id();

      // TODO: consider consecutive mutation
      if (vg->get_sample_from_vertex_if_exists(v, REF, sample))
      { ref_pos = sample.index(); }
			if (vg->get_sample_from_vertex_if_exists(v, sample_id, sample_find))
      {
        v_find = v;
        sample_found = true;
        sample_pos = sample_find.index();
      }
      ++it;
    }

    if (sample_found == true) {break;}
  }

  DEBUG("The closest node before " << pos << " contain sample "
      << sample_id << " is node " << v_find << ". Ref Index: " << cur_pos);

  return v_find;
} // get_prev_vertex_with_sample()


// return sequence of sample_id in ref position [pos_x, pos_y)
std::string query_sample_from_ref ( VariantGraph *vg, Index *idx,
                                    const uint64_t pos_x, const uint64_t pos_y,
                                    const string sample_id) {
  std::string seq = "";
  uint64_t ref_pos;
  uint64_t sample_pos;
  // get the closest node containing sample_id upstream pos_x
  Graph::vertex closest_v = get_prev_vertex_with_sample(vg, idx, pos_x,
                                                        sample_id,
                                                        ref_pos, sample_pos);

  // Do DFS on the path of the sample record seq from pos_x to pos_y in ref
  VariantGraph::VariantGraphPathIterator it = vg->find(closest_v, sample_id);
  bool record_seq = false;
  std::string temp;

  while (!it.done()) {
    temp.assign(vg->get_sequence(*(*it))); // store seq of current node
    DEBUG("Check node: " << (*it)->vertex_id() << ", reference Index: "
          << ref_pos << ", seq: " << temp);
    uint64_t l = (*it)->length();

    // if there no outgoing node containing ref
    // TODO: consider consecutive mutation
    uint64_t next_ref_pos = ref_pos + l;

    // Check outgoing nodes, update ref_pos to the index of node containing ref
    VariantGraph::VariantGraphIterator bfs_it = vg->find((*it)->vertex_id(), 1);
    ++bfs_it;
    while (!bfs_it.done()) {
      Graph::vertex v = (*bfs_it)->vertex_id();
      VariantGraphVertex::sample_info sample;
      if (vg->get_sample_from_vertex_if_exists(v, REF, sample)){
        next_ref_pos = sample.index();
        break;
      }
      ++bfs_it;
    }

    DEBUG("Next reference index: " << next_ref_pos);


    if (record_seq == true && next_ref_pos < pos_y) {
      // already start recoding and no reach pos_y at this node
      seq += temp;
    } else if (record_seq == true && next_ref_pos >= pos_y) {
      // already start recoding and no reach pos_y at this node
      seq += temp.substr(0, pos_y - ref_pos);
      break;
    } else if (next_ref_pos >= pos_x && next_ref_pos < pos_y) {
      // reach pos_x at this node and not reach pos_y
      record_seq = true;
      seq += temp.substr(pos_x - ref_pos);
    } else if (next_ref_pos >= pos_x && next_ref_pos >= pos_y) {
      // This node covers [pos_x, pos_y)
      seq = temp.substr(pos_x - ref_pos, pos_y - pos_x);
      break;
    }

    ++it;
    ref_pos = next_ref_pos; // update ref_pos
  }

  return seq;
} // query_sample_from_ref()


// return sequence of sample_id in sample's position [pos_x, pos_y)
std::string query_sample_from_sample ( VariantGraph *vg, Index *idx,
                                       const uint64_t pos_x,
                                       const uint64_t pos_y,
                                       const string sample_id ){
  std::string seq = "";
  uint64_t ref_pos;
  uint64_t sample_pos;

  Graph::vertex closest_v = get_prev_vertex_with_sample(vg, idx, pos_x,
                                                        sample_id,
                                                        ref_pos, sample_pos);
  // look for prev node containing sample until such node's sample pos < pos_x
  while (sample_pos >= pos_x && closest_v > 0)
  {
    uint64_t pos = ref_pos;
    closest_v = get_prev_vertex_with_sample(vg, idx, pos, sample_id,
                                            ref_pos, sample_pos);
  }

  DEBUG("Closest notes containing sample with sample_pos < "
        << pos_x << ": " << closest_v);

  // DFS from such node and record seq btw [pos_x, pos_y) in sample's coordinate
  VariantGraph::VariantGraphPathIterator it = vg->find(closest_v, sample_id);
  bool record_seq = false;
  std::string temp;

  while (!it.done()) {
    temp.assign(vg->get_sequence(*(*it)));
    DEBUG("Check node: " << (*it)->vertex_id() << ", sample Index: "
          << sample_pos << ", seq " << temp);
    uint64_t l = (*it)->length();
    uint64_t next_sample_pos = sample_pos + l;

    if (record_seq == true && next_sample_pos < pos_y) {
      seq += temp;
    } else if (record_seq == true && next_sample_pos >= pos_y) {
      seq += temp.substr(0, pos_y - sample_pos);
      break;
    } else if (next_sample_pos >= pos_x && next_sample_pos < pos_y) {
      record_seq = true;
      seq += temp.substr(pos_x - sample_pos);
    } else if (next_sample_pos >= pos_x && next_sample_pos >= pos_y) {
      seq = temp.substr(pos_x - sample_pos, pos_y - pos_x);
      break;
    }

    ++it;
    sample_pos = next_sample_pos;
  }
  return seq;
} // query_sample_from_sample()


vcflib::Variant closest_mutation (uint64_t pos){
  vcflib::Variant var;
  return var;
}
}
