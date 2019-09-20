#!/bin/sh

dir=$1

#set -x
echo $dir
#Combine chrom
        for i in {1..22}; do
                echo chr$i;
                bcftools merge --force-sample -r chr$i -m all $(find /dbgap/ckingsf_A018251/ppandey2/tcga-new/tcga-chrom-bcf/$dir -name "*.vcf.gz") -o /dbgap/ckingsf_A018251/ppandey2/tcga-new/tcga-merged/$dir/chr$i.vcf;
        done;

#For chrom x and y
        echo chrX;
        bcftools merge --force-sample -r chrX -m all $(find /dbgap/ckingsf_A018251/ppandey2/tcga-new/tcga-chrom-bcf/$dir -name "*.vcf.gz") -o /dbgap/ckingsf_A018251/ppandey2/tcga-new/tcga-merged/$dir/chrX.vcf;
        echo chrY;
        bcftools merge --force-sample -r chrY -m all $(find /dbgap/ckingsf_A018251/ppandey2/tcga-new/tcga-chrom-bcf/$dir -name "*.vcf.gz") -o /dbgap/ckingsf_A018251/ppandey2/tcga-new/tcga-merged/$dir/chrY.vcf;

#Compress and index new combined chrom
        for file in $(find /dbgap/ckingsf_A018251/ppandey2/tcga-new/tcga-merged/$dir -name "*.vcf"); do ./compress_index_vcf.sh $file; done

