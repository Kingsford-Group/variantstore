#include "index.h"
#include "variant_graph.h"
#include <string>
#include "vcflib/Variant.h"
#include "variantgraphvertex.pb.h"

namespace variantdb {
// TODO: check validity of inputs

// Do BFS of radius 1 to prev node in ref until sample is found
Graph::vertex get_prev_vertex_with_sample ( const VariantGraph *vg,
                                            const Index *idx,
                                            const uint64_t pos,
                                            const uint32_t sample_id_idx,
                                            uint64_t &ref_pos,
                                            uint64_t &sample_pos) {
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
      vg->get_sample_from_vertex_if_exists(v, 0, sample_find);
      break;
    }

    VariantGraph::VariantGraphIterator it = vg->find(v, 1);
    while(!it.done())
		{
      v = (*it)->vertex_id();
      //DEBUG("check vertex: " << v);
      if (vg->get_sample_from_vertex_if_exists(v, 0, sample))
      { ref_pos = sample.index(); }
			if (vg->get_sample_from_vertex_if_exists(v, sample_id_idx, sample_find))
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
std::string query_sample_from_ref ( const VariantGraph *vg, const Index *idx,
                                    const uint64_t pos_x, const uint64_t pos_y,
                                    const string sample_id) {
  std::string seq = "";
  uint64_t ref_pos;
  uint64_t sample_pos;
  uint32_t sample_id_idx; //TODO
  Graph::vertex closest_v = get_prev_vertex_with_sample(vg, idx, pos_x,
                                                        sample_id_idx,
                                                        ref_pos, sample_pos);
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
      if (vg->get_sample_from_vertex_if_exists(v, 0, sample))
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
      seq = temp.substr(pos_x - ref_pos, pos_y - pos_x);
      break;
    }

    ++it;
    ref_pos = next_ref_pos;
  }
  // start append to ref from pos_x to pos_y
  return seq;
}

std::string query_sample_from_sample ( const VariantGraph *vg, const Index *idx,
                                       const uint64_t pos_x,
                                       const uint64_t pos_y,
                                       const string sample_id ){
  std::string seq = "";
  uint64_t ref_pos;
  uint64_t sample_pos;
  uint32_t sample_id_idx; //TODO
  Graph::vertex closest_v = get_prev_vertex_with_sample(vg, idx, pos_x,
                                                        sample_id_idx,
                                                        ref_pos, sample_pos);
  while (sample_pos >= pos_x)
  {
    uint64_t pos = ref_pos;
    closest_v = get_prev_vertex_with_sample(vg, idx, pos, sample_id_idx,
                                            ref_pos, sample_pos);
  }

  VariantGraph::VariantGraphPathIterator it = vg->find(closest_v, sample_id);
  bool record_seq = false;
  std::string temp;

  while (!it.done()) {
    temp.assign(vg->get_sequence(*(*it)));
    DEBUG("Check node: " << (*it)->vertex_id() << ", Reference Index: "
          << sample_pos << ", Seq " << temp);
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
}

vcflib::Variant closest_mutation (uint64_t pos){
  vcflib::Variant var;
  return var;
}
}
