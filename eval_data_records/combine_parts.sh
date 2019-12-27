#!/bin/sh

dir=$1

for i in {2..22}; do
	echo chr$i;
	bcftools merge --force-sample -r chr$i -m all /dbgap/ckingsf_A018251/ppandey2/tcga-new/tcga-merged/$dir/chr$i"_first.vcf.gz"  /dbgap/ckingsf_A018251/ppandey2/tcga-new/tcga-merged/$dir/chr$i"_second.vcf.gz" -o /dbgap/ckingsf_A018251/ppandey2/tcga-new/tcga-merged/$dir/chr$i.vcf;
	./compress_index_vcf.sh /dbgap/ckingsf_A018251/ppandey2/tcga-new/tcga-merged/$dir/chr$i.vcf;
done;

echo chrX; 
bcftools merge --force-sample -r chrX -m all /dbgap/ckingsf_A018251/ppandey2/tcga-new/tcga-merged/$dir/chrX_first.vcf.gz  /dbgap/ckingsf_A018251/ppandey2/tcga-new/tcga-merged/$dir/chrX_second.vcf.gz -o /dbgap/ckingsf_A018251/ppandey2/tcga-new/tcga-merged/$dir/chrX.vcf;
./compress_index_vcf.sh /dbgap/ckingsf_A018251/ppandey2/tcga-new/tcga-merged/$dir/chrX.vcf;

echo chrY; 
bcftools merge --force-sample -r chrY -m all /dbgap/ckingsf_A018251/ppandey2/tcga-new/tcga-merged/$dir/chrY_first.vcf.gz  /dbgap/ckingsf_A018251/ppandey2/tcga-new/tcga-merged/$dir/chrY_second.vcf.gz -o /dbgap/ckingsf_A018251/ppandey2/tcga-new/tcga-merged/$dir/chrY.vcf;
./compress_index_vcf.sh /dbgap/ckingsf_A018251/ppandey2/tcga-new/tcga-merged/$dir/chrY.vcf;

#for file in $(find /dbgap/ckingsf_A018251/ppandey2/tcga-new/tcga-merged/$dir -name "*.vcf"); do ./compress_index_vcf.sh $file; done

