#include "dot_graph.h"
#include "index.h"
#include "variant_graph.h"
#include <string>
#include "vcflib/Variant.h"
#include "variantgraphvertex.pb.h"

namespace variantdb {
// TODO: check validity of inputs

// Do BFS of radius 1 to prev node in ref until sample is found
Graph::vertex get_prev_vertex_with_sample (VariantGraph *vg, Index *idx,
                                          uint64_t pos, std::string sample_id,
                                          uint64_t &ref_pos, uint64_t &sample_pos) {
  uint64_t cur_pos = pos;
  Graph::vertex v = idx->find(cur_pos);
  Graph::vertex v_find = v;
  VariantGraphVertex::sample_info sample;
  VariantGraphVertex::sample_info sample_find;

  bool find_sample = false;

  while (true) {
    Graph::vertex prev_v = v;
    Graph::vertex new_v = prev_v;
    while (new_v == prev_v && cur_pos > 1) {
      cur_pos = cur_pos - 1;
      new_v = idx->find(cur_pos);
    }

    v = new_v;

    if (cur_pos <= 1) {
      ref_pos = 1;
      v_find = v;
      vg->get_sample_from_vertex_if_exists(v, "ref", sample_find);
      break;
    }

    VariantGraph::VariantGraphIterator it = vg->find(v, 1);
    while(!it.done())
		{
      v = (*it)->vertex_id();
      //DEBUG("check vertex: " << v);
      if (vg->get_sample_from_vertex_if_exists(v, "ref", sample))
      { ref_pos = sample.index(); }
			if (vg->get_sample_from_vertex_if_exists(v, sample_id, sample_find))
      {
        v_find = v;
        find_sample = true;
      }
      ++it;
    }

    if (find_sample == true) {break;}
  }

  sample_pos = sample_find.index();
  DEBUG("The closest node before " << pos << " contain sample "
      << sample_id << " is node " << v_find << ". Ref Index: " << cur_pos);
  return v_find;
}

// return sequence of sample_id from ref position [x, y)
std::string query_sample_from_ref (VariantGraph *vg, Index *idx,
  uint64_t pos_x, uint64_t pos_y, std::string sample_id) {
  std::string seq = "";
  uint64_t ref_pos;
  uint64_t sample_pos;
  Graph::vertex closest_v = get_prev_vertex_with_sample(vg, idx, pos_x,
    sample_id, ref_pos, sample_pos);
  // Do DFS on the path of the sample, record ref idx along the way
  VariantGraph::VariantGraphPathIterator it = vg->find(closest_v, sample_id);
  bool record_seq = false;
  std::string temp;
  while (!it.done()) {
    temp.assign(vg->get_sequence(*(*it)));
    DEBUG("Check node: " << (*it)->vertex_id() << ", Reference Index: "
          << ref_pos << ", Seq " << temp);
    uint64_t l = (*it)->length();
    //BFS to find the next ref pos
    uint64_t next_ref_pos = ref_pos + l; //TODO: consider other cases

    VariantGraph::VariantGraphIterator temp_it = vg->find((*it)->vertex_id(), 1);
    ++temp_it;
    while (!temp_it.done())
    {
      Graph::vertex v = (*temp_it)->vertex_id();
      VariantGraphVertex::sample_info sample;
      if (vg->get_sample_from_vertex_if_exists(v, "ref", sample))
      { next_ref_pos = sample.index(); }
      ++temp_it;
    }

    DEBUG("Next reference index: " << next_ref_pos);

    if (record_seq == true && next_ref_pos < pos_y) {
      seq += temp;
    } else if (record_seq == true && next_ref_pos >= pos_y) {
      seq += temp.substr(0, pos_y - ref_pos);
      break;
    } else if (next_ref_pos >= pos_x && next_ref_pos < pos_y) {
      record_seq = true;
      seq += temp.substr(pos_x - ref_pos);
    } else if (next_ref_pos >= pos_x && next_ref_pos >= pos_y) {
      DEBUG(pos_x - ref_pos << " " << pos_y - pos_x);
      seq = temp.substr(pos_x - ref_pos, pos_y - pos_x);
      break;
    }

    ++it;
    ref_pos = next_ref_pos;
  }
  // start append to ref from pos_x to pos_y
  return seq;
}

std::string query_sample_from_sample (uint64_t pos_x,
  uint64_t pos_y, std::string sample_id){
  std::string seq;
  return seq;
}

vcflib::Variant closest_mutation (uint64_t pos){
  vcflib::Variant var;
  return var;
}
}
