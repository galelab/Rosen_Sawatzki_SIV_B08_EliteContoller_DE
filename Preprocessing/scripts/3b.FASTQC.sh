#!/bin/sh
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 16
#SBATCH -t 5800-12
#SBATCH --partition HoldingPen
#SBATCH -w austin
#SBATCH -e slurm-%j.err
#SBATCH -o slurm-%j.out
#SBATCH -J qc


PATH_FASTQ="./nohmrRNA_noglobin/"
samples=("$PATH_FASTQ"*.fastq.*.gz)

for sample in  ${samples[*]}
do
	srun -c 16 /vol01/ngs_tools/FastQC/fastqc "$sample" --noextract -t 16 -o ./nohmrRNA_noglobin/nohmrRNA_noglobin_fastqc_results/
	wait
done
