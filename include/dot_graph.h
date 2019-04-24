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
  bool createDotGraph(VariantGraph *vg, Graph::vertex v = 0,
          uint64_t radius = UINT64_MAX, std::string prefix="graph")
  {
    std::ofstream dotfile;
    dotfile.open (prefix + ".dot");

    if (!dotfile.is_open())
    {
      ERROR ("Can't open the file graph output. ");
      return false;
    }

    VariantGraph::VariantGraphIterator it = vg->find(v, radius);
    dotfile << "digraph {\n";

    /*
    ref += "\tsubgraph cluster_0 {\n";
    ref += "\t\tlable=\"reference\";\n";
    */

    while(!it.done())
    {
      uint64_t node_id = (*it)->vertex_id();
			//uint64_t sample_id = (*it)->s_info(0).sample_id();

      //Todo: Get outgoing nodes
      for (const auto v : vg->topology.out_neighbors((*it)->vertex_id()))
      {
        const VariantGraphVertex out_node = vg->vertex_list.vertex(node_id);
        uint64_t out_node_id = out_node.vertex_id();
        dotfile << "\t" + std::to_string(node_id) + "->"
                << std::to_string(node_id) + ";\n";
      }

      ++it;
    }

    dotfile << "}";
    dotfile.close();
    return true;
  } // DotGraph

}

#endif
