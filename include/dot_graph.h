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

#ifndef __DOT_GRAPH_H__
#define __DOT_GRAPH_H__

#include "variant_graph.h"
#include "variantgraphvertex.pb.h"
#include "graph.h"
#include "util.h"

#include <string>
#include <fstream>
#include <iostream>

/*
* Sample Output: granph.dot
* digraph {
*     subgraph cluster_0 {
*        label="reference";
*        0 -> 1;
*        1 -> 2;
*        4 -> 5;
*    }
*
*     1 -> 3;
*     3 -> 5;
*
* }
*/
using namespace variantdb;

namespace variantdb {
  bool is_ref_node(const VariantGraphVertex* v) {
    for (int i = 0; i < v->s_info_size(); i++) {
			const VariantGraphVertex::sample_info& s = v->s_info(i);
    //  DEBUG("Node " << v->vertex_id() << " contains sample " << s.sample_id());
			if (s.sample_id() == "ref") {
        return true;
      }
		}
    return false;
  }

	std::string get_samples(const VariantGraphVertex* v) {
		std::string samples;
		samples += "[ label=\"" + std::to_string(v->vertex_id()) + " l:" +
			std::to_string(static_cast<int>(v->length())) + "\n(";
    for (int i = 0; i < v->s_info_size(); ++i) {
			const VariantGraphVertex::sample_info& s = v->s_info(i);
			samples += s.sample_id() + " i:" +
				std::to_string(static_cast<int>(s.index()));
			if (i < v->s_info_size() - 1)
				samples += "\n";
		}
		samples += ")\"]";
    return samples;
  }

  bool createDotGraph(VariantGraph *vg, std::string prefix, Graph::vertex v = 0,
                      uint64_t radius = UINT64_MAX)
  {
    std::ofstream dotfile;
    dotfile.open (prefix + "/graph.dot");
		std::string labels = "";
    std::string ref = "";
    std::string sample = "";

    if (!dotfile.is_open())
    {
      ERROR ("Can't open the file graph output. ");
      return false;
    }


    VariantGraph::VariantGraphIterator it = vg->find(v, radius);
    while(!it.done()) {
			labels += std::to_string((*it)->vertex_id()) + get_samples(*it) + "\n";
      ++it;
		}

    dotfile << "digraph {\n";
		dotfile << labels;

    it = vg->find(v, radius);
    ref += "\tsubgraph cluster_0 {\n";
    ref += "\t\tlabel=\"reference\";\n";

    while(!it.done())
    {
      uint64_t node_id = (*it)->vertex_id();

      // Find neibors of node
      VariantGraph::VariantGraphIterator node = vg->find(node_id, 1);
      ++node;

      while (!node.done()) {
        uint64_t neibor_id = (*node)->vertex_id();
        if (is_ref_node(*node) && is_ref_node(*it)) {
          //DEBUG("Node " << node_id << " is on the reference path");
          ref += "\t\t" + std::to_string(node_id) + " -> " ;
          ref += std::to_string(neibor_id) + "\n";
        }
        else {
          sample += "\t" + std::to_string(node_id) + " -> " ;
          sample += std::to_string(neibor_id) + "\n";
        }
        ++node;
      }

      ++it;
    }

    ref += "\t}\n";
    dotfile << ref;
    dotfile << sample;
    dotfile << "}";
    dotfile.close();
    return true;
  } // DotGraph

}

#endif
