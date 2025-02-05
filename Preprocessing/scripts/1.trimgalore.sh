#!/bin/sh
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 24
#SBATCH -t 5800-12
#SBATCH --partition HoldingPen
#SBATCH -w mendel
#SBATCH -e slurm-%j.err
#SBATCH -o slurm-%j.out
#SBATCH -J Watkins07_trim

#AUTHOR: Leanne Whitmore ed KS
#Trim reads with trim galore

module use /share/tools/Modules/modulefiles-shared
module load TrimGalore/0.6.10
trim_galore --version

RAW_SOURCEDIR='../final-fastq-set/'
R1_files=("$RAW_SOURCEDIR"*R1*.fastq.gz)
R2_files=("$RAW_SOURCEDIR"*R2*.fastq.gz)
echo "Number samples = ${#R1_files[*]} R1 files"
echo "Number samples = ${#R2_files[*]} R2 files"

counter=0

while [ "$counter" -lt  ${#R1_files[*]} ]
do

     echo "Counter variable $counter"

     # Trim Galore
     ###PARAMETERS###
     #-q: Trims low quality ends from reads (defualt phred score is 20)
     #--phred33: instructs cutadapt to use ASCII+33 quality scores as Phred scores (this is the defualt)
     #--fastqc: run fastqc 
     #--output_dir: directory to put trimmed fastq files
     #no adapter sequence or option specified which means it runs autodetect (looks for Illumina universal, Nextera transposase or Illumina small RNA adapter sequences)
     if [ "$counter" -lt ${#R1_files[*]} ]
     then
        echo "Trim Galore running..."
        echo "${R1_files[$counter]}"
        echo "${R2_files[$counter]}"
        srun -c 4 trim_galore -q 20 --phred33 --gzip --cores 4 \
             --paired ${R1_files[$counter]} ${R2_files[$counter]} --output_dir trimgalore_results & 
     	counter=$((counter+1))
    fi
    if [ "$counter" -lt ${#R1_files[*]} ]
    then
        echo "Trim Galore running..."
        echo "${R1_files[$counter]}"
        echo "${R2_files[$counter]}"
	    srun -c 4 trim_galore -q 20 --phred33 --gzip --cores 4 \
            --paired ${R1_files[$counter]} ${R2_files[$counter]} --output_dir trimgalore_results &
        counter=$((counter+1))
    fi
    if [ "$counter" -lt ${#R1_files[*]} ]
    then
        echo "Trim Galore running..."
        echo "${R1_files[$counter]}"
        echo "${R2_files[$counter]}"
	    srun -c 4 trim_galore -q 20 --phred33 --gzip --cores 4 \
             --paired ${R1_files[$counter]} ${R2_files[$counter]} --output_dir trimgalore_results &
        counter=$((counter+1))
    fi
    if [ "$counter" -lt ${#R1_files[*]} ]
    then
        echo "Trim Galore running..."
        echo "${R1_files[$counter]}"
        echo "${R2_files[$counter]}"
	    srun -c 4 trim_galore -q 20 --phred33 --gzip --cores 4 \
            --paired ${R1_files[$counter]} ${R2_files[$counter]} --output_dir trimgalore_results &
        counter=$((counter+1))
    fi
    if [ "$counter" -lt ${#R1_files[*]} ]
    then
        echo "Trim Galore running..."
        echo "${R1_files[$counter]}"
        echo "${R2_files[$counter]}"
	    srun -c 4 trim_galore -q 20  --phred33 --gzip --cores 4 \
            --paired ${R1_files[$counter]} ${R2_files[$counter]} --output_dir trimgalore_results &
    	counter=$((counter+1))
    fi
    if [ "$counter" -lt ${#R1_files[*]} ]
    then
        echo "Trim Galore running..."
        echo "${R1_files[$counter]}"
        echo "${R2_files[$counter]}"
        srun -c 4 trim_galore -q 20 --phred33 --gzip --cores 4 \
            --paired ${R1_files[$counter]} ${R2_files[$counter]} --output_dir trimgalore_results &
    fi
    wait
    counter=$((counter+1))
done