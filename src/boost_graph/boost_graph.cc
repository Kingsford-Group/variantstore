#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_utility.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/config.hpp>
#include <iostream>
#include <stdlib.h>
#include <unordered_set>
#include <vector>
#include <openssl/rand.h>
#include <iterator>
#include <algorithm>
#include <unordered_map>

using namespace boost;


int main ( int argc, char *argv[] )
{

  /* copied from main.cc implementation */
  if (argc < 2) {
		fprintf(stderr, "Please specify the log of the number of slots in the CQF.\n");
		exit(1);
	}

	uint64_t qbits = atoi(argv[1]);
	uint64_t nslots = (1ULL << qbits);
	uint64_t nvals = 10*nslots/100;

	// Generate random values, stored in vals
	uint32_t *vals = (uint32_t*)malloc(nvals*sizeof(vals[0]));
	RAND_bytes((unsigned char *)vals, sizeof(*vals) * nvals);
	for (uint32_t i = 0; i < nvals; i++) {
		vals[i] = (1 * vals[i]) % UINT32_MAX;
	}

	srand(time(NULL)); //random seed
  /* copied */

  // Define graph type
  // Check Vertex and Edge Type ???
  typedef property<vertex_index_t, uint32_t> VertexProperty;
  typedef adjacency_list <setS, vecS, directedS, VertexProperty> graph_t;
  typedef graph_traits<graph_t>::vertex_descriptor Vertex;

  graph_t g;
  // Map from self-infined node to graph vertex index
  std::unordered_map<uint32_t, Vertex> id_vertex;
  std::unordered_map<Vertex, uint32_t> vertex_id;

  /* Randomly create adjacency list and store into BoostGraph */
  std::unordered_map<uint32_t, std::unordered_set<uint32_t>> adj_list;

	for (uint32_t i = 0; i < nvals; i++)
  {
		uint32_t source = vals[i]; // randomly select vertex id
		uint32_t nedges = rand() % 4;	//number of edges

    if (id_vertex.find(source) == id_vertex.end())
    {
      id_vertex[source] = add_vertex(source, g);
      vertex_id[id_vertex[source]] = source;
    }

    std::unordered_set<uint32_t> targets; // create targets for key

		for (uint32_t j = 0; j < nedges; j++)
    {
			uint32_t t = vals[rand() % nvals];

      if (t == source) {continue;}

      targets.insert(t);

      if (id_vertex.find(t) == id_vertex.end())
      {
        id_vertex[t] = add_vertex(t, g);
        vertex_id[id_vertex[t]] = t;
      }

      add_edge(id_vertex[source], id_vertex[t], g);
		}

    // Update targets for souce
    if (adj_list.find(source) == adj_list.end())
      adj_list[source] = targets;
    else
		  adj_list[source].insert(targets.begin(), targets.end());
	}

	std::cout << "Num vertices: " << num_vertices(g) << std::endl;
	std::cout << "Num edges: " << num_edges(g) << std::endl;


  // Correctness Check
  // If g contains adj_list
  auto print = [](const uint32_t& n) { std::cout << " " << n; };

  for (auto it = adj_list.begin(); it != adj_list.end(); ++it)
  {
    graph_t::adjacency_iterator vit, vend;
    std::unordered_set<uint32_t> out_v;
    for (boost::tie(vit, vend) = adjacent_vertices(id_vertex[it->first], g); vit != vend; vit++)
      out_v.insert(vertex_id[*vit]);

    if (it->second != out_v)
    {
      std::cout << "Correctness test failed!" << std::endl;
      std::cout << it->first << std::endl;
      std::cout << "Outcoming verteces: " << std::endl;
      std::for_each (out_v.begin(), out_v.end(), print);
      std::cout << "\n";
      std::cout << "Values in map: " << std::endl;
      std::for_each (it->second.begin(), it->second.end(), print);

			return EXIT_FAILURE;
    }
	}

  // If adj_list contains g
  adjacency_list<>::vertex_iterator vs, ve;

  for (boost::tie(vs, ve) = vertices(g); vs != ve; ++vs)
  {
    graph_t::adjacency_iterator vit, vend;
    std::unordered_set<uint32_t> out_v;

    for (boost::tie(vit, vend) = adjacent_vertices(*vs, g); vit != vend; vit++)
      out_v.insert(vertex_id[*vit]);

    if (adj_list[vertex_id[*vs]] != out_v)
    {
      std::cout << "Correctness test failed!" << std::endl;
      std::cout << vertex_id[*vs] << std::endl;
      std::cout << "Outcoming verteces: " << std::endl;
      std::for_each (out_v.begin(), out_v.end(), print);
      std::cout << "\n";
      std::cout << "Values in map: " << std::endl;
      std::for_each (adj_list[vertex_id[*vs]].begin(), adj_list[vertex_id[*vs]].end(), print);
			return EXIT_FAILURE;
    }
  }

  return EXIT_SUCCESS;
}
