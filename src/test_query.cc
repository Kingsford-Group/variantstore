#include <iostream>
#include <fstream>
#include <vector>
#include <time.h>
#include <string>

#include "variant_graph.h"
#include "dot_graph.h"
#include "query.h"

using namespace variantstore;
std::shared_ptr<spdlog::logger> console;

void print_var (Variant *var) {
	PRINT("At pos: " << var->var_pos << ", ref: " << var->ref);
	for (auto a=var->alts.begin(); a!=var->alts.end(); a++) {
		std::cout << "alt: " << *a << ", samples: ";
		for (auto s=var->alt_sample_map[*a].begin();
							s != var-> alt_sample_map[*a].end(); s++) {
			std::cout << (*s) << ",";
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;
	return;
}

int main ( int argc, char *argv[] ) {
	if (argc < 2) {
		fprintf(stderr, "Please specify the reference fasta file and vcf file.\n");
		exit(1);
	}

	//Construct vg & index
	console = spdlog::default_logger();
	PRINT("Creating variant graph");
	std::string ref_file(argv[1]);
	PRINT("Creating variant graph");
	std::string vcf_file(argv[2]);
	PRINT("Creating variant graph");
	PRINT("Creating variant graph");
	VariantGraph vg(ref_file, vcf_file, "./ser", READ_COMPLETE_GRAPH);
	PRINT("Creating variant graph");

	PRINT("Creating Index");
	Index idx(&vg);
	createDotGraph(&vg, "./ser");
	/*
  PRINT("TEST get_prev_vertex_with_sample");
  int pos = 1;
  std::string sample_id;
  PRINT("Please specify position: ");
  std::cin >> pos;
  PRINT("Please specify sample id: ");
  std::cin.ignore();
  std::getline(std::cin, sample_id);
  while (pos != -1)
  {
    uint64_t ref_pos = 0;
    uint64_t sample_pos = 0;
    PRINT("Query...");
    get_prev_vertex_with_sample(&vg, &idx, pos, sample_id, ref_pos, sample_pos);
    PRINT("Please specify position: ");
    std::cin >> pos;
    PRINT("Please specify sample id: ");
    std::cin.ignore();
    std::getline(std::cin, sample_id);
  }

  PRINT("TEST query_sample_from_ref");
  int pos_x = 1;
  int pos_y = 1;

  while (pos_x != -1)
  {
    std::string sample_id;
    PRINT("Please specify position_x: ");
    std::cin >> pos_x;
    PRINT("Please specify position_y: ");
    std::cin >> pos_y;
    PRINT("Please specify sample id: ");
    std::cin.ignore();
    std::getline(std::cin, sample_id);
    PRINT("Query...");
    std::string seq = query_sample_from_ref(&vg, &idx, pos_x, pos_y, sample_id);
    PRINT("Sample " << sample_id << " has sequence " << seq
          << " in reference coordinate.");
  }

	PRINT("TEST query_sample_from_sample");
  int pos_x = 1;
  int pos_y = 1;

  while (pos_x != -1)
  {
    std::string sample_id;
    PRINT("Please specify position_x: ");
    std::cin >> pos_x;
    PRINT("Please specify position_y: ");
    std::cin >> pos_y;
    PRINT("Please specify sample id: ");
    std::cin.ignore();
    std::getline(std::cin, sample_id);
    PRINT("Query...");
    std::string seq = query_sample_from_sample(&vg, &idx, pos_x, pos_y, sample_id);
    PRINT("Sample " << sample_id << " has sequence " << seq
          << " in reference coordinate.");
  }


	// Correctness test for query_sample_from_sample
	PRINT("Correctness test for query_sample_from_sample");
	if (argc < 4) {
		fprintf(stderr, "Please specify the reference fasta file, vcf file, sample_id and ground truth.\n");
		exit(1);
	}

	std::string sample_id(argv[3]);
	std::string sample_seq_file(argv[4]);
	std::string chr;
	std::string sample_seq;
	read_fasta(sample_seq_file, chr, sample_seq);
	//PRINT("Sample sequence: " << sample_seq);

	for (uint64_t pos_x=1; pos_x <= sample_seq.length(); pos_x++) {
		for (uint64_t pos_y=pos_x+1; pos_y <= sample_seq.length(); pos_y++) {

			std::string query_seq = query_sample_from_sample( &vg, &idx, pos_x,
																												pos_y, sample_id);
			std::string true_seq = sample_seq.substr(pos_x-1, pos_y - pos_x);

			if (query_seq != true_seq) {
				ERROR("Queried position: [" << pos_x << ", " << pos_y << "]\n "
							<< "Queried seq: " << query_seq
							<< "True seq: " << true_seq);
				exit(1);
			}

			PRINT("Queried position: [" << pos_x << ", " << pos_y << "]\n ");
		}
	}

	PRINT("Correctness test successful!");


	// Efficiency test
	if (argc < 3) {
		fprintf(stderr, "Please specify the reference fasta file, vcf file, and sample_id.\n");
		exit(1);
	}

	std::string sample_id(argv[3]);
	PRINT("Efficiency test for query_sample_from_sample");
	auto start = std::chrono::high_resolution_clock::now();
	uint64_t num_queried = 0;

	for (uint64_t pos_x=1; pos_x <= vg.get_ref_length(); pos_x++) {
		for (uint64_t pos_y=pos_x+1; pos_y <= vg.get_ref_length(); pos_y++) {
			query_sample_from_sample( &vg, &idx, pos_x, pos_y, sample_id);
			num_queried++;

			if (num_queried % 500 == 0) {
				auto stop = std::chrono::high_resolution_clock::now();
				auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start);
				auto avg_time = duration.count() / 500;
				PRINT("The average time taken for query is: " << avg_time	<< "seconds");
			}
		}
	}

	auto stop = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start);
	auto avg_time = duration.count() / num_queried;
	PRINT("The average time taken for query is: "	<< avg_time << " seconds");
	return 0;


	// test next_variant_in_ref
	PRINT("TEST get_prev_vertex_with_sample");
  int pos = 1;
  std::string sample_id;
  PRINT("Please specify position: ");
  std::cin >> pos;
  PRINT("Please specify sample id: ");
	std::cin.ignore();
  std::getline(std::cin, sample_id);
  while (pos != -1)
  {
		Variant var;
    PRINT("Query...");
    if (!next_variant_in_ref(&vg, &idx, pos, var))
		{
			PRINT("Not found!");
		} else {
			print_var(&var);
		}

    PRINT("Please specify position: ");
    std::cin >> pos;
    PRINT("Please specify sample id: ");
		std::cin.ignore();
    std::getline(std::cin, sample_id);
  }


	// test get_sample_var_in_sample
	PRINT("TEST get_sample_var_in_sample");
  int pos_x = 1;
  int pos_y = 1;

  while (pos_x != -1)
  {
    std::string sample_id;
    PRINT("Please specify position_x: ");
    std::cin >> pos_x;
    PRINT("Please specify position_y: ");
    std::cin >> pos_y;
    PRINT("Please specify sample id: ");
    std::cin.ignore();
    std::getline(std::cin, sample_id);
    PRINT("Query...");
		std::vector <Variant> vars;
    vars = get_sample_var_in_sample(&vg, &idx, pos_x, pos_y, sample_id);
    for (auto it = vars.begin(); it != vars.end(); it++) {
			print_var(&(*it));
		}
  }

	// test next_variant_in_ref
	int pos = 1;
	while (pos != -1)
	{
		PRINT("TEST next_variant_in_ref");
		PRINT("Please specify position: ");
		std::cin >> pos;
		PRINT("Query...");
		Variant var;
		next_variant_in_ref (&vg, &idx, pos, var);
		print_var(&var);
	}


	pos = 1;
	while (pos != -1)
	{
		PRINT("TEST closest_var");
		PRINT("Please specify position: ");
		std::cin >> pos;
		PRINT("Query...");
		Variant var;
		closest_var (&vg, &idx, pos, var);
		print_var(&var);
	}

	// test get_sample_var_in_sample
	PRINT("TEST get_sample_var_in_ref");
	int pos_x = 1;
	int pos_y = 1;

	while (pos_x != -1)
	{
		std::string sample_id;
		PRINT("Please specify position_x: ");
		std::cin >> pos_x;
		PRINT("Please specify position_y: ");
		std::cin >> pos_y;
		PRINT("Please specify sample id: ");
		std::cin.ignore();
		std::getline(std::cin, sample_id);
		PRINT("Query...");
		std::vector <Variant> vars;
		vars = get_sample_var_in_ref(&vg, &idx, pos_x, pos_y, sample_id);
		for (auto it = vars.begin(); it != vars.end(); it++) {
			print_var(&(*it));
		}
	}
*/


	PRINT("TEST samples_has_var");
	int pos = 1;
	while (pos != -1)
	{
		std::string ref;
		std::string alt;
		PRINT("Please specify position: ");
		std::cin >> pos;
		PRINT("Please specify ref: ");
		std::cin.ignore();
		std::getline(std::cin, ref);
		PRINT("Please specify alt: ");
		//std::cin.ignore();
		std::getline(std::cin, alt);
		PRINT("Query...");
		std::vector <std::string> samples;
		samples = samples_has_var (&vg, &idx,pos, ref, alt);
		PRINT("Samples has var: ");
		for (auto it = samples.begin(); it != samples.end(); it++) {
			std::cout << (*it) << " ";
		}
		PRINT();
	}

}
