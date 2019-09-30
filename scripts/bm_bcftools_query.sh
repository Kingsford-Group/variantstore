#!/bin/bash
# Author: Yinjie Gao, yinjieg@andrew.cmu.edu
# THis is the query benchmark using bcftools
# equivalent to variantdb function -- get_var_in_ref

# vcf file path:
vcf="$1"
num="$2"
len="$4"
pct="$3"
query_len=$(($len/$pct))
echo "Query length: $query_len"
# Generate 10 random starting positions
pos=$(echo $(shuf -i 16050000-52420579 -n $num))
# samples=$(echo $(shuf -n 10 $sample_file))
# echo $samples16860

r_bcf=''

for i in $pos
do
  j=$(($i+$query_len))
  ./bcftools filter -r "22:${i}-${j}" $vcf > temp.txt
  r_bcf="${r_bcf}22:${i}-${j},"
done

r_bcf=${r_bcf::-1}
