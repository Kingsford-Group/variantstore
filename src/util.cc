/*
 * ============================================================================
 *
 *         Author:  Prashant Pandey (), ppandey@cs.stonybrook.edu
 *   Organization:  Stony Brook University
 *
 * ============================================================================
 */

#include <iostream>
#include <algorithm>
#include <cstring>
#include <vector>
#include <set>
#include <unordered_set>
#include <bitset>
#include <cassert>

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <stdint.h>
#include <inttypes.h>
#include <unistd.h>

#include <stdexcept>
#include <signal.h>
#include <netinet/in.h>
#include <sys/errno.h>
#include <sys/socket.h>
#include <sys/types.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>

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

}
