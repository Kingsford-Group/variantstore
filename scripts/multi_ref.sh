#!/bin/bash

#(seq 1 22; echo X; echo Y) | parallel -j 24 "minimap2 -ax asm5 -t 24 GRCh37/Homo_sapiens.GRCh37.dna.chromosome.{}.fa GRCh38/Homo_sapiens.GRCh38.dna.chromosome.{}.fa > alignment/align{}.sam"

#(seq 1 22; echo X; echo Y) | parallel -j 24 "samtools view -S -b alignment/align{}.sam > alignment/align{}.bam"

#(seq 1 22; echo X; echo Y) | parallel -j 24 "samtools sort -o alignment/align-sort{}.bam alignment/align{}.bam"

#(seq 1 22; echo X; echo Y) | parallel -j 24 "freebayes -f GRCh37/Homo_sapiens.GRCh37.dna.chromosome.{}.fa --ploidy 1 alignment/align-sort{}.bam > vcf/align{}.vcf"

# for file in $(find vcf/ -name "*.vcf"); do bgzip $file; done
# for file in $(find vcf/ -name "*.vcf.gz"); do tabix $file; done
# (seq 1 22; echo X; echo Y) | parallel -j 24 "tabix -f ALL.chr{}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"

# bcftools merge --force-samples -r X -m all ../ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chrX.phase3_shapeit2_mvncall_integrated_v1b.20130502.genotypes.vcf.gz vcf/alignX.vcf.gz > vcf/mergedX.vcf

# bcftools merge --force-samples -r Y -m all ../ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chrY.phase3_integrated_v2a.20130502.genotypes.vcf.gz vcf/alignY.vcf.gz > vcf/mergedY.vcf

#(seq 1 22;) | parallel -j 22 "bcftools merge --force-samples -r {} -m all ../ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr{}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz vcf/align{}.vcf.gz > vcf/merged{}.vcf"

for i in {1..22}
do
  bcftools merge --force-samples -r $i -m all ../ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr$i.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz vcf/align$i.vcf.gz > vcf/merged$i.vcf
  bgzip vcf/merged$i.vcf
done
