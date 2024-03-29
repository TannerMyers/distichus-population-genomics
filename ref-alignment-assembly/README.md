
# ***Anolis distichus*** ddRADseq Project

## Reference genome alignment, variant calling, and variant filtering

This file outlines workflow of ddRADseq loci assembly via alignment to a reference genome. I demultplexed individual reads using Stacks' `process_radtags` function before mapping them to the *Anolis distichus* assembly with BWA MEM, then .sam files are converted to .bam and sorted with samtools, assessed for mapping quality with Qualimap, and, finally, SNPs are called using bcftools. Final dataset filtering is completed using vcftools.

**Setting up the environment**

I strongly recommend installing packages through a package manager (e.g., Anaconda, Miniconda, or Mamba). I used Mamba to install and manage necessary software, including:
    - BWA 0.7.17
    - samtools 1.15
    - bcftools 1.15.1
    - vcftools 0.1.16
    - PLINK 1.90

**Alignment**

First, I used the script `bwa_index.sh` to index the *distichus* genome (*AnoDis1.0.fasta.gz*) before aligning raw reads.

```
genome_fa=/home/tcm0036/distichus-ddRAD/genome/Anolis_distichus/AnoDis1.0.fasta.gz

bwa index $genome_fa
```

Next, I used the script `bwa-alignment-individual-jobs.sh` to align demultiplexed reads for each individual with `bwa mem`, convert .sam to .bam, sort, and index .bam files.

```
popmap=/home/tcm0036/distichus-ddRAD/info/distichus-popmap.tsv
# Isolate individual IDs in population map to use as looping variables
samples=`awk 'BEGIN {OFS = FS} {print $1}' $popmap`

for ID in $samples;
        do
        sbatch --job-name=$ID-alignment --partition=jro0014_amd --cpus-per-task=8 --mem=80G --time=4-24:00:00 --wrap="

        # Load modules 
        module load bwa
        module load samtools

        # Assign variables 

        # # Path to Anolis carolinensis genome fasta
        # genome=/home/tcm0036/distichus-ddRAD/genome/Anolis_carolinensis.AnoCar2.0.dna.toplevel.fa
        # # Path to including prefix of database generated by bwa index
        # bwa_db=/home/tcm0036/distichus-ddRAD/genome/bwa/anocar

        # Path to Anolis distichus genome fasta & BWA database
        genome=/home/tcm0036/distichus-ddRAD/genome/Anolis_distichus/AnoDis1.0.fasta.gz
        bwa_db=/home/tcm0036/distichus-ddRAD/genome/Anolis_distichus/AnoDis1.0.fasta.gz

        fq1=/home/tcm0036/distichus-ddRAD/samples/$ID.1.fq.gz
        fq2=/home/tcm0036/distichus-ddRAD/samples/$ID.2.fq.gz
        sam=/scratch/tcm0036/distichus-ddRAD/alignments/AnoDist/sam/$ID.aligned.sam
        bam=/scratch/tcm0036/distichus-ddRAD/alignments/AnoDist/bam/$ID.aligned.bam
        sorted_bam=/scratch/tcm0036/distichus-ddRAD/alignments/AnoDist/bam/$ID.aligned.sorted.bam

        bwa mem -t 8 \$bwa_db \$fq1 \$fq2 > \$sam
        samtools view -b \$sam > \$bam
        samtools sort -o \$sorted_bam \$bam
        samtools index \$sorted_bam

        echo $ID
        "
done
```

**Variant calling**

I used the script `bcftools-variant-calling.sh` to call variants.

```
genome=/home/tcm0036/distichus-ddRAD/genome/Anolis_distichus/AnoDis1.0.fasta
popmap=/home/tcm0036/distichus-ddRAD/info/popmap-bcftools.tsv
outputfile=/scratch/tcm0036/distichus-ddRAD/alignments/AnoDist/bcf/variants.bcf
tmp=mapped-files.tmp.list

ls /scratch/tcm0036/distichus-ddRAD/alignments/AnoDist/bam/*sorted.bam > $tmp 

bcftools mpileup -Ou --annotate AD,DP --fasta-ref $genome --bam-list $tmp | \
bcftools call --variants-only --threads $SLURM_CPUS_PER_TASK --multiallelic-caller --group-samples $popmap -Ob -o $outputfile
```

**VCF filtering**

I filtered variants based on:
    - Quality
    - Minimum depth of coverage
    - Maximum depth of coverage
    - Minor allele count
    - Missingness across loci (call rate)
    - Missingness within individuals

Initial Stats:
10,035,425 variants, non-filtered, across 203 specimens.

Before filtering for missing data:
529,597 variants were left across 203 specimens.

Filter 1: call rate > 50% and > 90% missing data
165,433 SNPs retained and specimens 4404 and 6081 were dropped. 201 specimens remaining after this step.

Filter 2: call rate > 60% and > 70% missing data
133,260 SNPs retained and specimens 4403 and 866 were dropped. 199 specimens remaining after this step.

Filter 3: call rate > 70% and > 50% missing data
100,687 SNPs retained and specimens 6072, 6388, and 6780 were dropped. 196 specimens remaining after this step..

Filter 4: call rate > 80% and > 30% missing data
70,854 SNPs retained and specimen 867 was dropped. 195 specimens remaining after this step.

Filter 5: call rate > 90% and > 25% missing data
37,822 SNPs retained and no specimens were dropped (195 specimens).

Filter 6: call rate > 95%
20,104 SNPs retained (195 specimens).

Filter 7: call rate > 100%
1,843 SNPs retained (195 specimens).

**Linkage Pruning**

I applied the same parameters for LD pruning on three datasets with varying levels of missingness: 80%, 90%, and 100%.

- window size of 50 kb
- step size of 100 bp
- r^2 of 0.1

**80%**

```
plink
 --allow-extra-chr
  --double-id
  --indep-pairwise 50 100 0.1
  --out distichus_geno80ind30_filtered_LD_50kb_100_0.1
  --set-missing-var-ids @:#
  --vcf ../../filtered/distichus_variants_renamed_minQ20minDP10maxDP291mac3geno80ind30_FIL-4.recode.vcf
```

Now, prune SNPs that don't meet thresholds:

```
plink
  --allow-extra-chr
  --double-id
  --extract distichus_geno80ind30_filtered_LD_50kb_100_0.1.prune.in
  --make-bed
  --out distichus_geno80_filtered_LD_50kb_100_0.1_pruned
  --recode vcf
  --set-missing-var-ids @:#
  --vcf ../../filtered/distichus_variants_renamed_minQ20minDP10maxDP291mac3geno80ind30_FIL-4.recode.vcf

```
This leaves 10160 SNPs

**90%**

```
plink
  --allow-extra-chr
  --double-id
  --indep-pairwise 50 100 0.1
  --out distichus_geno90ind25_filtered_LD_50kb_100_0.1
  --set-missing-var-ids @:#
  --vcf ../../filtered/distichus_variants_renamed_minQ20minDP10maxDP291mac3geno90ind25_FIL-4.recode.vcf
```

Now, prune SNPs that don't meet thresholds:

```
plink 
  --allow-extra-chr
  --double-id
  --extract distichus_geno90ind25_filtered_LD_50kb_100_0.1.prune.in
  --make-bed
  --out distichus_geno90_filtered_LD_50kb_100_0.1_pruned
  --recode vcf
  --set-missing-var-ids @:#
  --vcf ../../filtered/distichus_variants_renamed_minQ20minDP10maxDP291mac3geno90ind25_FIL-4.recode.vcf
```
This leaves 5886 SNPs

**100%**

```
plink
--allow-extra-chr
  --double-id
  --indep-pairwise 50 100 0.1
  --out distichus_geno100_filtered_LD_50kb_100_0.1
  --set-missing-var-ids @:#
  --vcf ../../filtered/distichus_variants_renamed_minQ20minDP10maxDP291mac3geno100_FIL-4.recode.vcf
```

Now, prune SNPs that don't meet thresholds:

```
plink 
  --allow-extra-chr
  --double-id
  --extract distichus_geno100_filtered_LD_50kb_100_0.1.prune.in
  --make-bed
  --out distichus_geno100_filtered_LD_50kb_100_0.1_pruned
  --recode vcf
  --set-missing-var-ids @:#
  --vcf 
../../filtered/distichus_variants_renamed_minQ20minDP10maxDP291mac3geno100_FIL-4.recode.vcf
```
This leaves 842 SNPs
