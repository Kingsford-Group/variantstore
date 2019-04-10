/*
 * ============================================================================
 *
 *         Author:  Prashant Pandey (), ppandey2@cs.cmu.edu
 *   Organization:  Stony Brook University
 *
 * ============================================================================
 */

#include <iostream>
#include <sstream>
#include <fstream>

#include <vector>
#include <set>
#include <unordered_set>
#include <bitset>
#include <cassert>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <stdint.h>
#include <inttypes.h>
#include <unistd.h>

#include "util.h"

namespace variantdb {
	float cal_time_elapsed(struct timeval* start, struct timeval* end)
	{
		struct timeval elapsed;
		if (start->tv_usec > end->tv_usec) {
			end->tv_usec += 1000000;
			end->tv_sec--;
		}
		elapsed.tv_usec = end->tv_usec - start->tv_usec;
		elapsed.tv_sec = end->tv_sec - start->tv_sec;
		return (elapsed.tv_sec * 1000000 + elapsed.tv_usec)/1000000.f;
	}
	
	void print_time_elapsed(std::string desc, struct timeval* start, struct
													timeval* end)
	{
		struct timeval elapsed;
		if (start->tv_usec > end->tv_usec) {
			end->tv_usec += 1000000;
			end->tv_sec--;
		}
		elapsed.tv_usec = end->tv_usec - start->tv_usec;
		elapsed.tv_sec = end->tv_sec - start->tv_sec;
		float time_elapsed = (elapsed.tv_sec * 1000000 + elapsed.tv_usec)/1000000.f;
		std::cout << desc << "Total Time Elapsed: " << std::to_string(time_elapsed) << 
			"seconds" << std::endl;
	}

	void read_fasta(const std::string fasta_file, std::string& chr, std::string&
									ref) {
		std::ifstream stream(fasta_file);
		if (!stream.good()) {
			std::cerr << "Failed to open input fasta file: " << fasta_file << '\n';
			abort();
		}
		bool found_ref = false;;
		std::string line;
		while(getline(stream, line)) {
			if (line.at(0) == '>') {
				if (found_ref) {
					std::cerr << "Found multiple references in the fasta file" << "\n";
					abort();
				}
				std::stringstream line_stream(line);
				getline(line_stream, chr, ' ');
				chr = chr.substr(1);
				found_ref = true;
			} else {
				ref.append(line);
			}
		}
	}
}
