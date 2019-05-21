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
#include <vector>
#include <map>
#include <cassert>
#include "index.h"
#include "variant_graph.h"
#include "vcflib/Variant.h"
#include "variantgraphvertex.pb.h"

std::string REF = "ref";

namespace variantdb {
// TODO: check validity of inputs

/* ----------------------------------------------------------------------------
Support Func: traverse back until a node containg sample_id is found
return such node, record the pos in ref and sample
*/
Graph::vertex get_prev_vertex_with_sample ( VariantGraph *vg, Index *idx,
                                            const uint64_t pos,
                                            const string sample_id,
                                            uint64_t &ref_pos,
                                            uint64_t &sample_pos) {
  // the node right before cur_pos
  uint64_t cur_pos = pos;
  uint64_t cur_ref_node_idx;
  Graph::vertex v = idx->find(cur_pos, cur_ref_node_idx);
  Graph::vertex v_find = v;
  VariantGraphVertex::sample_info sample;
  VariantGraphVertex::sample_info sample_find;
  bool sample_found = false;

  while (true) { // do until find sample or to the front of graph
    v = idx->previous(cur_ref_node_idx);

    DEBUG("cur_ref_node_idx: " << cur_ref_node_idx);
    // if reach the first node in the graph, return the first node
    if (cur_ref_node_idx <= 1) {
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
      cur_ref_node_idx--;
    }

    if (sample_found == true) {break;}
  }

  DEBUG("The closest node before " << pos << " contain sample "
      << sample_id << " is node " << v_find << ". Ref Index: " << ref_pos
      << ". Sample index: " << sample_pos);

  return v_find;
} // get_prev_vertex_with_sample()


/* ----------------------------------------------------------------------------
Give positional range and sample id
return sequence of sample_id in ref position [pos_x, pos_y)
*/
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

/* ----------------------------------------------------------------------------
Give positional range and sample id
return sequence of sample_id in sample's position [pos_x, pos_y)
*/
std::string query_sample_from_sample ( VariantGraph *vg, Index *idx,
                                       const uint64_t pos_x,
                                       const uint64_t pos_y,
                                       const string sample_id ) {
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


/* ----------------------------------------------------------------------------
Support Func: get samples in a nodes
Return true if the node contains samples other than ref
*/
bool get_samples ( const VariantGraphVertex* v, VariantGraph *vg,
                   std::vector <std::string> &sample_ids)
{
  bool is_var = false;
  sample_ids = {};
  for (int i = 0; i < v->s_info_size(); ++i) {
    const VariantGraphVertex::sample_info& s = v->s_info(i);
    std::string sample_id = vg->get_sample_name(s.sample_id());
    sample_ids.push_back(sample_id);
    if (sample_id != REF) {is_var = true;}
  }
  return is_var;
}


typedef std::map <std::string, std::vector <std::string>> AltSamplesMap;

struct Variant {
  uint64_t var_pos;
  std::string ref;
  std::vector <std::string> alts;
  AltSamplesMap alt_sample_map;
};


/* ----------------------------------------------------------------------------
Support Func: Return true if found the closest variants at or after the given pos. Return false is no variant is found reaching the end of the graph
*/
bool next_variant_in_ref ( VariantGraph *vg, Index *idx, const uint64_t pos,
                           Variant &var) {
  bool found_var = false;
  Graph::vertex v = idx->find(pos); // the node before pos
  VariantGraph::VariantGraphPathIterator it = vg->find(v, REF);
  while (!it.done()) {
    PRINT("Check vertes: " << (*it)->vertex_id());
    VariantGraph::VariantGraphIterator bfs_it=vg->find((*it)->vertex_id(), 1);
    ++bfs_it;
    uint32_t num_out_nodes = 0;

    // check all outgoing edges and record alternatives
    while (!bfs_it.done())
    {
      PRINT("Check bfs vertes: " << (*bfs_it)->vertex_id());
      num_out_nodes++;
      Graph::vertex out_v = (*bfs_it)->vertex_id();
      std::vector <std::string> sample_ids;
      // got variant nodes' samples_ids
      if (get_samples((*bfs_it), vg, sample_ids)) {
        std::string alt = vg->get_sequence(*(*bfs_it));
        var.alts.push_back(alt);
        found_var = true;
        var.alt_sample_map[alt] = sample_ids;
      }
      ++bfs_it;
    }
    DEBUG("Num outgoing edges for node " << v << " is: " << num_out_nodes);

    if (found_var == true) {
      ++it;
      var.ref.assign(vg->get_sequence(*(*it)));
      VariantGraphVertex::sample_info sample;
      if (vg->get_sample_from_vertex_if_exists((*it)->vertex_id(),REF,sample))
        var.var_pos = sample.index();
      else
      {
        ERROR((*it)->vertex_id() << "is not in the reference path");
        abort();
      }
      break;
    }

    ++it;
  } // end DFS in ref
  return found_var;
}


/* ----------------------------------------------------------------------------
Find a variant near given position
*/
bool closest_var ( VariantGraph *vg, Index *idx, const uint64_t pos,
                   Variant &var) {
  bool found_var = false;
  // find next_var_in_ref
  // find prev var in ref ???

  return found_var;
}


/* ----------------------------------------------------------------------------
return all variants in sample_id occuring in sample's coordinate
[pos_x, pos_y)
*/
std::vector <Variant> get_sample_var_in_sample ( VariantGraph *vg, Index *idx,
                                                 const uint64_t pos_x,
                                                 const uint64_t pos_y,
                                                 const string sample_id) {
  // traverse in sample's coordinate from pos_x to pos_y
  // report if the sample's node is a variant node
  std::vector <Variant> vars;
  uint64_t ref_pos;
  uint64_t sample_pos;
  Graph::vertex closest_v = get_prev_vertex_with_sample(vg, idx, pos_x,
                                                        sample_id,
                                                        ref_pos, sample_pos);

  // look for prev node containing sample until such node's sample pos < pos_x
  while (sample_pos >= pos_x && closest_v > 0)
  {
    uint64_t pos = ref_pos;
    closest_v = get_prev_vertex_with_sample(vg, idx, pos, sample_id, ref_pos,
                                            sample_pos);
  }

  // start from the ref node right before the node contains sample
  closest_v = idx->find(ref_pos);
  VariantGraphVertex::sample_info sample;
  uint64_t seq_len;

  if (vg->get_sample_from_vertex_if_exists(closest_v, REF, sample)){
     seq_len = ref_pos - sample.index();
     ref_pos = sample.index();
  } else {
    ERROR("reference node is expected to be found!");
  }

  sample_pos = sample_pos - seq_len;

  DEBUG("Closest note before node containing sample with sample_pos < "
        << pos_x << ": " << closest_v);

  // DFS from such node and record seq btw [pos_x, pos_y) in sample's coordinate
  VariantGraph::VariantGraphPathIterator it = vg->find(closest_v, sample_id);
  std::string cur_ref;

  while (!it.done()) {
    Graph::vertex cur_v = (*it)->vertex_id();
    std::vector <std::string> samples_in_node;
    // get next ref_pos & sample_pos
    uint64_t l = (*it)->length();
    // if there no outgoing node containing ref
    // TODO: consider consecutive mutation
    uint64_t next_ref_pos = ref_pos + l;
    uint64_t next_sample_pos = sample_pos + l;
    std::string next_ref;
    // Check outgoing nodes, update ref_pos to the index of node containing ref
    VariantGraph::VariantGraphIterator bfs_it = vg->find((*it)->vertex_id(), 1);
    ++bfs_it;
    while (!bfs_it.done()) {
      Graph::vertex v = (*bfs_it)->vertex_id();
      VariantGraphVertex::sample_info sample;
      if (vg->get_sample_from_vertex_if_exists(v, REF, sample)) {
        next_ref_pos = sample.index();
        next_ref = vg->get_sequence(*(*bfs_it));
      }

      ++bfs_it;
    }

    if (sample_pos >= pos_x && sample_pos < pos_y &&
        vg->get_sample_from_vertex_if_exists(cur_v, sample_id, sample)){
      Variant var;
      std::string alt = vg->get_sequence(*(*it));
      var.alts.push_back(alt);
      var.ref = cur_ref;
      get_samples((*it), vg, var.alt_sample_map[alt]);
      vars.push_back(var);
    } // Not variant or variant node does not contain sample_id

    cur_ref = next_ref;
    ref_pos = next_ref_pos;
    sample_pos = next_sample_pos;
    ++it;
  }
} // get_sample_var_in_sample()


/* ---------------------------------------------------------------------------- Given the position range and sample_id
Return all variants in sample_id occuring in ref coordinate [pos_x, pos_y)
*/
std::vector <Variant> get_sample_var_in_ref ( VariantGraph *vg, Index *idx,
                                              const uint64_t pos_x,
                                              const uint64_t pos_y,
                                              const string sample_id) {
  // traverse in sample's coordinate from pos_x to pos_y
  // report if the sample's node is a variant node
} // mut_in_ref()


// Given a variant (ref & alt) and the position
// Return samples has such variant in the reference coordinate
std::vector <std::string> samples_has_var ( VariantGraph *vg, Index *idx,
                                       const uint64_t pos, const string ref,
                                       const string alt) {

} // samples_has_var()

}
