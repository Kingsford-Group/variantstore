#!/bin/bash

for j in {1..1}
  		do
	      ./scripts/bm_vs_query_analysis.sh vs_v1/chr22 51304566 $1 0 HG01872 >> query_mem.txt;
done

