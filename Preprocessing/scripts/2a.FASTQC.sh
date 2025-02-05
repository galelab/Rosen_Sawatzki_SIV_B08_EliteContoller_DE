#!/bin/sh
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 16
#SBATCH -t 5800-12
#SBATCH --partition HoldingPen
#SBATCH -e slurm-%j.err
#SBATCH -o slurm-%j.out
#SBATCH -J Watkins07_fastqc


PATH_FASTQ="../trimgalore_results/"
samples=("$PATH_FASTQ"*.fq.gz)

mkdir trimgalore_fastqc_results

for sample in  ${samples[*]}
do
	srun -c 16 /vol01/ngs_tools/FastQC/fastqc "$sample" --noextract -t 16 -o ./trimgalore_fastqc_results/
	wait
done
