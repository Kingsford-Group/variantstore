#!/bin/sh

# vcfpath:
vcf="ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"

# Generate 10 random positions
a=$(echo $(shuf -i 0-50000000 -n 10))
echo $a

start=`date +%s`
for i in $a
do
  j=$(($i+5000000))
  echo "$i, $j"
  bcftools filter -r 22:$i-$j $vcf > temp.vcf
done

end=`date +%s`
runtime=$(((end-start)/10))
echo "Average runtime: $runtime"
