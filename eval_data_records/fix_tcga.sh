i#!/bin/sh

dir=$1
echo $dir
for file in $(find /dbgap/ckingsf_A018251/ppandey2/tcga-new/tcga-data/$dir -name "*.vcf.gz"); do
        #echo $file;
        basename=$(basename "${file}")
        #echo $basename;
        bcftools annotate -x ^FORMAT/GT $file > /dbgap/ckingsf_A018251/ppandey2/tcga-new/tcga-chrom-bcf/$dir/$basename.vcf
 done;
