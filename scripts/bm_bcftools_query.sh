#!/bin/sh
# Author: Yinjie Gao, yinjieg@andrew.cmu.edu
# THis is the query benchmark using bcftools
# equivalent to variantdb function -- get_sample_var_in_ref

# vcf file path:
vcf="$1"
sample_file="$2"
len="$3"
query_len=$(($len/10))

# Generate 10 random starting positions
pos=$(echo $(shuf -i 0-$len -n 10))
echo $pos

samples=$(echo $(shuf -n 10 $sample_file))
echo $samples

start=`date +%s`
for i in $pos
do
  j=$(($i+$query_len))

  for s in $samples
  do
    echo "$i, $j, $s"
    bcftools filter -r 22:$i-$j $vcf -s $s > temp.vcf
  done
done

end=`date +%s`
runtime=$(((end-start)/100))
echo "Average runtime: $runtime"
