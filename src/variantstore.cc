/*
 * ============================================================================
 *
 *       Filename:  variantstore.cc
 *
 *         Author:  Prashant Pandey (), ppandey2@cs.cmu.edu
 *   Organization:  Carnegie Mellon University
 *
 * ============================================================================
 */

#include <iostream>
#include <string>
#include <stdlib.h>

#include "spdlog/spdlog.h"

#include "progopts.h"
#include "clipp.h"

template <typename T>
void explore_options_verbose(T& res) {
  if(res.any_error()) { std::cerr << "error\n"; }

  //aggregated errors
  if(res.unmapped_args_count()) { std::cerr << "error unmapped args count\n"; /* ... */ }
  if(res.any_bad_repeat()) { std::cerr << "error bad repeat \n"; /* ... */ }
  if(res.any_blocked())    { std::cerr << "error blocked \n"; /* ... */ }
  if(res.any_conflict())   { std::cerr << "error conflict\n"; /* ... */ }

  for(const auto& m : res.missing()) {
    std::cerr << "missing " << m.param() << " after index " << m.after_index() << '\n';
  }

  //per-argument mapping
  for(const auto& m : res) {
    std::cerr << m.index() << ": " << m.arg() << " -> " << m.param();
    std::cerr << " repeat #" << m.repeat();
    if(m.blocked()) std::cerr << " blocked";
    if(m.conflict()) std::cerr << " conflict";
    std::cerr<< '\n';
  }
}

int construct_main(ConstructOpts& construct_opt);
int query_main(QueryOpts& query_opt);

std::shared_ptr<spdlog::logger> console;

/*
 * ===  FUNCTION  =============================================================
 *         Name:  main
 *  Description:
 * ============================================================================
 */
	int
main ( int argc, char *argv[] ) {
  using namespace clipp;
	enum class mode {construct, query, help};
  mode selected = mode::help;

	ConstructOpts construct_opt;
	QueryOpts query_opt;

	console = spdlog::default_logger();
#ifdef DEBUG_MODE
	console->set_level(spdlog::level::debug);
#endif

	auto construct_mode = (
									command("construct").set(selected, mode::construct),
									required("-r","--reference") & value("reference-file", construct_opt.ref) %
									"reference file",
									required("-v","--vcf") & value("vcf-file", construct_opt.vcf) %
									"variant call format file",
									required("-p","--output-prefix") & value(
																												 "output-prefix",
												construct_opt.prefix) %
									"output directory"
						 );

	auto ensure_sample_name = [](const QueryOpts& opts) -> bool {
		if (opts.type == 3 || opts.type == 7)
			return true;
		else if (opts.sample_name == "")
			return false;
		return false;
	};

	auto query_mode = (
									command("query").set(selected, mode::query),
									required("-p","--output-prefix") & value(
																												 "output-prefix",
												query_opt.prefix) %
									"output directory",
									required("-t","--type") & value("query-type", query_opt.type) %
									"Types of query. \n \
                  1. Return closest variant in reference coordinates. \n \
                  2. Get sample's sequence in reference coordinates. \n \
                  3. Get sample's sequence in sample coordinates. \n \
                  4. Get sample's variants in reference coordinates. \n \
                  5. Get sample's variants in sample coordinates. \n \
                  6. Get variants in reference coordinates. \n \
                  7. Return samples with a given mutation.",
                  required("-r","--region") & value("region", query_opt.region) %
                  "region in format <start>:<end>, regions separated by ','",
                  required("-m","--mode") & value("mode", query_opt.mode) %
                  "READ_INDEX_ONLY: 0, READ_COMPLETE_GRAPH:1",
                  option("-o","--output_file") & value("outfile", query_opt.outfile) %
                  "output_file",
									option("-s","--sample-name") & value("sample-name",
																											 query_opt.sample_name)
										%
									"sample name",
                  option("-a","--alt-seq") & value("alt-seq",
																											 query_opt.alt) %
									"alternative sequences, separated by ','",
                  option("-b","--ref-seq") & value("ref-seq",
																											 query_opt.ref) %
									"reference sequences, separated by ','",
                  option("-v","--verbose").set(query_opt.verbose, true) %
									"print vcf"
						 );

	auto cli = ((construct_mode | query_mode |
							 command("help").set(selected, mode::help)),
							option("-v", "--version").call([]{std::cout <<
																						 "version 0.1\n\n";}).doc("show version"));

  assert(construct_mode.flags_are_prefix_free());

  decltype(parse(argc, argv, cli)) res;
  try {
    res = parse(argc, argv, cli);
  } catch (std::exception& e) {
		std::cout << "\n\nParsing command line failed with exception: " <<
			e.what() << "\n";
    std::cout << "\n\n";
    std::cout << make_man_page(cli, "variantstore");
    return 1;
  }

  if(res) {
		switch(selected) {
			case mode::construct: construct_main(construct_opt); break;
			case mode::query:
				try {
					ensure_sample_name(query_opt);
				} catch (std::exception& e) {
					std::cout << "\n\nParsing command line failed with exception: " <<
						e.what() << "\n";
					std::cout << "\n\n";
					std::cout << make_man_page(cli, "variantstore");
					return 1;
				}// required("-b","--begin") & value("begin", query_opt.begin) %
									// "starting position",
									// required("-e","--end") & value("end", query_opt.end) %
									// "ending position",
				query_main(query_opt); break;
			case mode::help:  break;
		}
  } else {
    auto b = res.begin();
    auto e = res.end();
    if (std::distance(b,e) > 0) {
      if (b->arg() == "construct") {
        std::cout << make_man_page(construct_mode, "variantstore");
      } else if (b->arg() == "query") {
        std::cout << make_man_page(query_mode, "variantstore");
      }else {
        std::cout << "There is no command \"" << b->arg() << "\"\n";
        std::cout << usage_lines(cli, "variantstore") << '\n';
      }
    } else {
      std::cout << usage_lines(cli, "variantstore") << '\n';
    }
  }

	return EXIT_SUCCESS;
}				/* ----------  end of function main  ---------- */
