/*
 * ============================================================================
 *
 *       Filename:  rank_select.h
 *
 *         Author:  Prashant Pandey (), ppandey2@cs.cmu.edu
 *   Organization:  Carnegie Mellon University
 *
 * ============================================================================
 */

#ifndef _RANK_SELECT_H_
#define _RANK_SELECT_H_

#include <sys/types.h>
#include <stdlib.h>
#include <stdint.h>


int popcnt(uint64_t val);

uint64_t bitselect(uint64_t val, int rank);

uint64_t word_rank(uint64_t val);

uint64_t word_select(uint64_t val, int rank);

#endif //_RANK_SELECT_H_
