#!/bin/sh

dir=$1

#set -x
echo $dir
#Combine chrom 
	for i in {1..22}; do
		echo chr$i; 
		bcftools merge --force-sample -r chr$i -m all $(cat $dir"_first.lst") -o /dbgap/ckingsf_A018251/ppandey2/tcga-new/tcga-merged/$dir/chr$i"_first.vcf";
		bcftools merge --force-sample -r chr$i -m all $(cat $dir"_second.lst") -o /dbgap/ckingsf_A018251/ppandey2/tcga-new/tcga-merged/$dir/chr$i"_second.vcf";
	done;

#For chrom x and y
	echo chrX; 
	bcftools merge --force-sample -r chrX -m all $(cat $dir"_first.lst") -o /dbgap/ckingsf_A018251/ppandey2/tcga-new/tcga-merged/$dir/chrX_first.vcf;
	bcftools merge --force-sample -r chrX -m all $(cat $dir"_second.lst") -o /dbgap/ckingsf_A018251/ppandey2/tcga-new/tcga-merged/$dir/chrX_second.vcf;
	echo chrY; 
	bcftools merge --force-sample -r chrY -m all $(cat $dir"_first.lst") -o /dbgap/ckingsf_A018251/ppandey2/tcga-new/tcga-merged/$dir/chrY_first.vcf;
	bcftools merge --force-sample -r chrY -m all $(cat $dir"_second.lst") -o /dbgap/ckingsf_A018251/ppandey2/tcga-new/tcga-merged/$dir/chrY_second.vcf;

#Compress and index new combined chrom
	for file in $(find /dbgap/ckingsf_A018251/ppandey2/tcga-new/tcga-merged/$dir -name "*.vcf"); do ./compress_index_vcf.sh $file; done
