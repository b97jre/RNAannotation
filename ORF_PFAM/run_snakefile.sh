#!/bin/bash -l
#SBATCH -A b2011098
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 2:00:00
#SBATCH -J Snakemake

module load bioinfo-tools
module load hmmer/3.1b1-gcc
module use /proj/b2013006/sw/modules
module load snakemake
module load htseq/0.6.1p1

snakemake -j 2 join_ORF_PFAM
snakemake --dag | dot -Tpdf > snake_graphic.pdf
snakemake --dag | dot -Tpng > snake_graphic.png
