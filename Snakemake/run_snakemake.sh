#!/bin/bash -l
#SBATCH -A b2011098
#SBATCH -p core
#SBATCH -n 16 
#SBATCH -t 80:00:00
#SBATCH -J Snakemake
#SBATCH --mail-type=All
#SBATCH --mail-user olofsson.anna12@gmail.com

module load bioinfo-tools
module load bowtie2
module load samtools
module load hmmer/3.1b1-gcc
module use /proj/b2013006/sw/modules
module load snakemake
module load htseq/0.6.1p1

snakemake -j 16
