#! /usr/bin/env bash 

# Load conda environment with vcftools
source ~/mambaforge/etc/profile.d/conda.sh # Or path to where your conda is
conda activate genomics_env

# depth indv/locus
vcftools --vcf /scratch/phyletica/distichus/alignments/results/vcf/all-samples-merged-SNPs.recode.vcf --out filtering/raw_snps --depth
vcftools --vcf /scratch/phyletica/distichus/alignments/results/vcf/all-samples-merged-SNPs.recode.vcf --out filtering/raw_snps --site-mean-depth
vcftools --vcf /scratch/phyletica/distichus/alignments/results/vcf/all-samples-merged-SNPs.recode.vcf --out filtering/raw_snps --geno-depth

# missing data indv/locus
vcftools --vcf /scratch/phyletica/distichus/alignments/results/vcf/all-samples-merged-SNPs.recode.vcf --out filtering/raw_snps --missing-indv
vcftools --vcf /scratch/phyletica/distichus/alignments/results/vcf/all-samples-merged-SNPs.recode.vcf --out filtering/raw_snps --missing-site

# allele freq/indv freq buden
vcftools --vcf /scratch/phyletica/distichus/alignments/results/vcf/all-samples-merged-SNPs.recode.vcf --out filtering/raw_snps --indv-freq-burden
vcftools --vcf /scratch/phyletica/distichus/alignments/results/vcf/all-samples-merged-SNPs.recode.vcf --out filtering/raw_snps --freq2
vcftools --vcf /scratch/phyletica/distichus/alignments/results/vcf/all-samples-merged-SNPs.recode.vcf --out filtering/raw_snps --singletons
vcftools --vcf /scratch/phyletica/distichus/alignments/results/vcf/all-samples-merged-SNPs.recode.vcf --out filtering/raw_snps --012

# heterozygosity per individual
vcftools --vcf /scratch/phyletica/distichus/alignments/results/vcf/all-samples-merged-SNPs.recode.vcf --out filtering/raw_snps --het

# SNP call quality
vcftools --vcf /scratch/phyletica/distichus/alignments/results/vcf/all-samples-merged-SNPs.recode.vcf --out filtering/raw_snps --site-quality