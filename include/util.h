/*
 * ============================================================================
 *
 *         Author:  Prashant Pandey (), ppandey2@cs.cmu.edu
 *   Organization:  Carnegie Mellon University
 *
 * ============================================================================
 */

#ifndef _UTIL_H_
#define _UTIL_H_

#include <iostream>
#include <cstring>
#include <vector>
#include <cassert>
#include <fstream>
#include <unordered_map>

#ifdef DEBUG_MODE
#define PRINT_DEBUG 1
#else
#define PRINT_DEBUG 0
#endif

#define DEBUG(x) do { \
	if (PRINT_DEBUG) { std::cerr << x << std::endl; } \
} while (0)

#define ERROR(x) do { \
	{ std::cerr << x << std::endl; } \
} while (0)

#define PRINT(x) do { \
	{ std::cout << x << std::endl; } \
} while (0)

namespace variantdb {
	enum DNA_MAP {A, C, T, G, N};  // A=0, C=1, T=2, G=3, N=4
	char map_int(uint8_t base);
	uint8_t map_base(char base);
	float cal_time_elapsed(struct timeval* start, struct timeval* end);
	/* Print elapsed time using the start and end timeval */
	void print_time_elapsed(std::string desc, struct timeval* start, struct
													timeval* end);
	// Read the fasta file and sets the value of the chr and ref string.
	// TODO Assumes there's only one ref in the fasta file. May need to handle
	// the case with multiple references in one fasta file.
	void read_fasta(const std::string fasta_file, std::string& chr, std::string&
									ref);

}
#endif
