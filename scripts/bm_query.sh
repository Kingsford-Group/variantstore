#!/bin/bash
num="$1"
pct="$2"
mode="0" #READ_INDEX_ONLY
# mode="1" #READ_COMPLETE_GRAPH
len="52420579"
file="../../ppandey/variantstore/vs_v4/chr22"
echo "Benchmark for $num queries of query length 1/$pct"
# cd /mnt/disk34/user/yinjieg/bcftools
# echo "[bcftools]"
# /usr/bin/time ../variantstore/scripts/bm_bcftools_query.sh ../../ppandey/vcfdata/ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz $num $pct $len

# cd /mnt/disk34/user/yinjieg/variantstore
echo "[variantstore]"
/usr/bin/time ./scripts/bm_vs_query.sh $file $num $pct $len 1 $mode
/usr/bin/time ./scripts/bm_vs_query.sh $file $num $pct $len 2 $mode
/usr/bin/time ./scripts/bm_vs_query.sh $file $num $pct $len 3 $mode
/usr/bin/time ./scripts/bm_vs_query.sh $file $num $pct $len 4 $mode
/usr/bin/time ./scripts/bm_vs_query.sh $file $num $pct $len 5 $mode
/usr/bin/time ./scripts/bm_vs_query.sh $file $num $pct $len 6 $mode
