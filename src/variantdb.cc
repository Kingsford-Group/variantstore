/*
 * ============================================================================
 *
 *       Filename:  variantdb.cc
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

	auto query_mode = (
									command("query").set(selected, mode::query),
									required("-p","--output-prefix") & value(
																												 "output-prefix",
												construct_opt.prefix) %
									"output directory",
									required("-t","--type") & value("query-type", query_opt.type) %
									"types of query. \n \
									1. Get variants in ref coordinate. \n \
									2. Get the number of variants in sample coordinate. \n \
									3. Get sample's sequence in sample coordinate \n \
									4. Return closest mutation in ref coordinate.",
									required("-b","--begin") & value("begin", query_opt.begin) %
									"starting position",
									required("-e","--end") & value("end", query_opt.end) %
									"ending position",
									option("-s","--sample-name") & value("sample-name", query_opt.sample_name) %
									"sample name"
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
    std::cout << make_man_page(cli, "variantdb");
    return 1;
  }

  if(res) {
		switch(selected) {
			case mode::construct: construct_main(construct_opt); break;
			case mode::query: query_main(query_opt); break;
			case mode::help:  break;
		}
  } else {
    auto b = res.begin();
    auto e = res.end();
    if (std::distance(b,e) > 0) {
      if (b->arg() == "construct") {
        std::cout << make_man_page(construct_mode, "variantdb");
      } else if (b->arg() == "query") {
        std::cout << make_man_page(query_mode, "variantdb");
      }else {
        std::cout << "There is no command \"" << b->arg() << "\"\n";
        std::cout << usage_lines(cli, "variantdb") << '\n';
      }
    } else {
      std::cout << usage_lines(cli, "variantdb") << '\n';
    }
  }

	return EXIT_SUCCESS;
}				/* ----------  end of function main  ---------- */
