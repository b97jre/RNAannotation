#!/bin/bash -l
#SBATCH -A b2014162
#SBATCH -p core
#SBATCH -n 8
#SBATCH -t 2:00:00
#SBATCH -J Snakemake
#SBATCH --mail-type=All
#SBATCH --mail-user johan.reimegard@scilifelab.se

module load bioinfo-tools
module load bowtie2
module load samtools
module use /proj/b2013006/sw/modules
module load snakemake
module load htseq/0.6.1p1 #special module made by Marcel in order to avoid "python conflict"

snakemake -j 8
snakemake --dag | dot -Tpdf > snake_graphic.pdf
snakemake --dag | dot -Tpng > snake_graphic.png
