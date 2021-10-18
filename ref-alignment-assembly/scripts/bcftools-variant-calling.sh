#SBATCH --job-name variant-calling
#SBATCH --mail-type=ALL
#SBATCH --mail-user=tcm0036@auburn.edu
#SBATCH --time=30-00:00:00
#SBATCH --cpus-per-task 16 
#SBATCH --mem 60G
#SBATCH --partition jro0014_amd

genome=/home/tcm0036/distichus-ddRAD/genome/Anolis_carolinensis.AnoCar2.0.dna.toplevel.fa
popmap=/home/tcm0036/distichus-ddRAD/info/distichus-popmap.tsv
outputfile=/scratch/tcm0036/distichus-ddRAD/alignment/results/bcf/variants.bcf
tmp=mapped-files.tmp.list

ls /scratch/tcm0036/distichus-ddRAD/alignment/results/bam/*.bam > $tmp 

module load bcftools

bcftools mpileup -Ou -a AD --fasta-ref $genome --bam-list $tmp | \
bcftools call --variants-only --threads $SLURM_CPUS_PER_TASK --multiallelic-caller --group-samples $popmap -Ob -o $outputfile

rm $tmp
