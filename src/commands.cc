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
#include "sys/time.h"

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
    uint64_t beg = 0;
		uint64_t end = 0;

		if (pos2 == std::string::npos) {
			beg = std::stoi(token);
		}
		else {
			end = std::stoi(token.substr(pos2+1));
			beg = std::stoi(token.substr(0, pos2));
		}

		regions.push_back (std::make_tuple(beg, end));

		if (pos == std::string::npos)
			break;

		region = region.substr(pos + 1);
		pos = region.find(',');
	}
	std::sort(regions.begin(), regions.end());
	return regions;
}

	int
query_main ( QueryOpts& opts )
{
	console->info("Loading Index ...");
	Index idx(opts.prefix);
	console->info("Loading variant graph ...");
	enum READ_TYPE mode;
	if (opts.mode == 0)
	{
		mode = READ_INDEX_ONLY;
		console->info("Read index only ..");
	}
	else
	{
		mode = READ_COMPLETE_GRAPH;
		console->info("Read complete graph ..");
	}


	VariantGraph vg(opts.prefix, mode);
	console->info("Graph stats:");
	console->info("Chromosome: {} #Vertices: {} #Edges: {} Seq length: {}",
								vg.get_chr(), vg.get_num_vertices() , vg.get_num_edges(),
								vg.get_seq_length());


	if (opts.type == 1) {
		console->info("1. Get sample's sequence in ref coordinate ...");
	}
	if (opts.type == 2) {
		console->info("2. Get sample's sequence in sample's coordinate ...");
	}
	if (opts.type == 3) {
		console->info("3. Return closest mutation in ref coordinate. ...");
	}
	if (opts.type == 4) {
		console->info("4. Get sample's variants in ref coordinate ...");
	}
	if (opts.type == 5) {
		console->info("5. Get sample's variants in sample coordinate ...");
	}
	if (opts.type == 6) {
		console->info("6. Get variants in ref coordinate. ...");
	}

	std::vector<std::tuple<uint64_t, uint64_t>> regions = read_regions(opts.region);

	uint32_t query_num=0;
	struct timeval start, end;
	struct timezone tzp;                                                    gettimeofday(&start, &tzp);


	for (auto it = regions.begin(); it != regions.end(); it++)
	{
		uint32_t max = 600;
		uint32_t id = rand() % max + 1;
		std::string sample_id = vg.get_sample_name(id);
		opts.sample_name = sample_id;

		console->debug("Query region {}:{}", std::get<0>(*it), std::get<1>(*it));
		if (opts.type == 1) {
			query_sample_from_ref(&vg, &idx, std::get<0>(*it), std::get<1>(*it), opts.sample_name, opts.verbose, opts.outfile);
		}
		if (opts.type == 2) {
			query_sample_from_sample(&vg, &idx, std::get<0>(*it), std::get<1>(*it), opts.sample_name, opts.verbose, opts.outfile);
		}
		if (opts.type == 3) {
			Variant var;
			closest_var(&vg, &idx, std::get<0>(*it), var, opts.verbose, opts.outfile);
		}
		if (opts.type == 4) {
			get_sample_var_in_sample(&vg, &idx, std::get<0>(*it), std::get<1>(*it), opts.sample_name, opts.verbose, opts.outfile);
		}
		if (opts.type == 5) {
			get_sample_var_in_ref(&vg, &idx, std::get<0>(*it), std::get<1>(*it), opts.sample_name, opts.verbose, opts.outfile);
		}
		if (opts.type == 6) {
			get_var_in_ref(&vg, &idx, std::get<0>(*it), std::get<1>(*it), opts.verbose, opts.outfile);
		}

		query_num += 1;
		if ((query_num == 10) || (query_num == 100) || (query_num == 1000))
		{
			gettimeofday(&end, &tzp);
			std::string dsc = "Query" ;
			dsc.append(std::to_string(query_num));
			dsc.append(": ");
			print_time_elapsed(dsc, &start, &end);
		}
	}
	gettimeofday(&end, &tzp);
	std::string dsc = "Query" ;
	dsc.append(std::to_string(query_num));
	dsc.append(": ");
	if (opts.type == 6)
		dsc.append("(query_var_in_ref) ");
	print_time_elapsed(dsc, &start, &end);

	return EXIT_SUCCESS;
}				/* ----------  end of function main  ---------- */
