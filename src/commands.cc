/*
 * ============================================================================
 *
 *       Filename:  construct.cc
 *
 *         Author:  Prashant Pandey (), ppandey2@cs.cmu.edu
 										Yinjie Gao, yinjieg@andrew.cmu.edu
 *   Organization:  Carnegie Mellon University
 *
 * ============================================================================
 */

#include <stdlib.h>
#include <tuple>
#include "spdlog/spdlog.h"

#include "query.h"
#include "index.h"
#include "variant_graph.h"
#include "progopts.h"


using namespace variantdb;

/*
 * ===  FUNCTION  =============================================================
 *         Name:  main
 *  Description:
 * ============================================================================
 */
	int
construct_main (ConstructOpts &opts)
{
	console->info("Creating variant graph");
	VariantGraph vg(opts.ref, opts.vcf, opts.prefix, READ_INDEX_ONLY);

	console->info("Graph stats:");
	console->info("Chromosome: {} #Vertices: {} #Edges: {} Seq length: {}",
								vg.get_chr(), vg.get_num_vertices() , vg.get_num_edges(),
								vg.get_seq_length());
	console->info("Serializing variant graph to disk");
	vg.serialize();

	console->info("Creating Index");
	Index idx(&vg);
	console->info("Serializing index to disk");
	idx.serialize(opts.prefix);

#if 0
	console->info("Loading variant graph");
	VariantGraph file_vg(opts.prefix, READ_COMPLETE_GRAPH);
	console->info("Graph stats:");
	console->info("Chromosome: {} #Vertices: {} #Edges: {} Seq length: {}",
								file_vg.get_chr(), file_vg.get_num_vertices() , file_vg.get_num_edges(),
								file_vg.get_seq_length());
#endif

	return EXIT_SUCCESS;
}				/* ----------  end of function main  ---------- */



std::vector<std::tuple<uint64_t, uint64_t>> read_regions (std::string region)
{
	std::vector<std::tuple<uint64_t, uint64_t>> regions;

	auto pos = region.find(',');
	while (true) {
    std::string token = region.substr(0, pos);
		auto pos2 = token.find(':');
    uint64_t beg = std::stoi(token.substr(0, pos2));
		uint64_t end = std::stoi(token.substr(pos2+1));
		regions.push_back (std::make_tuple(beg, end));

		if (pos == std::string::npos)
			break;

		region = region.substr(pos + 1);
		pos = region.find(',');
	}
	return regions;
}

	int
query_main ( QueryOpts& opts )
{
	console->info("Loading Index ...");
	Index idx(opts.prefix);
	console->info("Loading variant graph ...");
	VariantGraph vg(opts.prefix, READ_INDEX_ONLY);
	console->info("Graph stats:");
	console->info("Chromosome: {} #Vertices: {} #Edges: {} Seq length: {}",
								vg.get_chr(), vg.get_num_vertices() , vg.get_num_edges(),
								vg.get_seq_length());


	if (opts.type == 1) {
		console->info("Get variants in ref coordinate ...");
		std::vector<std::tuple<uint64_t, uint64_t>> regions = read_regions(opts.region);

		for (auto it = regions.begin(); it != regions.end(); it++) {
			console->info("Query region {}:{}", std::get<0>(*it), std::get<1>(*it));
			get_var_in_ref(&vg, &idx, std::get<0>(*it), std::get<1>(*it), true);
		}


	}

	// if (opts.type == 2) {
	// 	console->info("Get sample's variants in ref coordinate ...");
	// 	std::vector <Variant> vars = get_sample_var_in_ref(&vg, &idx, opts.begin, opts.end, opts.sample_name);
	// 	for (auto it = vars.begin(); it != vars.end(); it++) {
	// 		print_var(&(*it));
	// 	}
	// }
	//
	// if (opts.type == 3) {
	// 	console->info("Get the number of variants in sample coordinate ...");
	// 	std::vector <Variant> vars = get_sample_var_in_sample(&vg, &idx, opts.begin, opts.end, opts.sample_name);
	// 	// for (auto it = vars.begin(); it != vars.end(); it++) {
	// 	// 	print_var(&(*it));
	// 	// }
	// 	std::cout << vars.size();
	// }
	// if (opts.type == 4) {
	// 	console->info("Get sample's sequence in sample coordinate ...");
	// 	std::string s= query_sample_from_sample(&vg, &idx, opts.begin, opts.end, opts.sample_name);
	// 	std::cout << s;
	// }
	// if (opts.type == 5) {
	// 	console->info("Return closest mutation in ref coordinate ...");
	// 	Variant var;
	// 	closest_var(&vg, &idx, opts.begin, var);
	// 	print_var(&var);
	// }


	return EXIT_SUCCESS;
}				/* ----------  end of function main  ---------- */
