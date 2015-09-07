Different code for how to run program
=====================================

Everythin here is specified how to run on uppmax. Lets start there and push for generality later.

ORF finder
==========

In: reference file (nt)

Out: reference file (aa)


Code



Mapping of reads
================

in: reads + reference file
out: counts per contig per reads file ()

code:

module load bioinfo-tools
module load bowtie2/2.2.3
module load samtools

bowtie2  --threads 16 -x /gulo/proj_nobackup/b2011098/private/deNovoAnnotation/references/bowtie2/human.sample123.trinity -1 /gulo/proj_nobackup/b2011098/private/deNovoAnnotation/reads/human/7_111116_AD0341ACXX_137_5_index5_1.fastq -2 /gu
lo/proj_nobackup/b2011098/private/deNovoAnnotation/reads/human/7_111116_AD0341ACXX_137_5_index5_2.fastq --un-conc /gulo/proj_nobackup/b2011098/private/deNovoAnnotation/bowtie2/human/7_111116_AD0341ACXX_137_5_index5_/7_111116_AD0341ACXX_13
7_5_index5__human.sample123.trinity.sam.noHit.fastq -S /gulo/proj_nobackup/b2011098/private/deNovoAnnotation/bowtie2/human/7_111116_AD0341ACXX_137_5_index5_/7_111116_AD0341ACXX_137_5_index5__human.sample123.trinity.sam


samtools view -bSh -o /gulo/proj_nobackup/b2011098/private/deNovoAnnotation/bowtie2/human/7_111116_AD0341ACXX_137_5_index5_/7_111116_AD0341ACXX_137_5_index5__human.sample123.trinity.bam /gulo/proj_nobackup/b2011098/private/deNovoAnnotation/bowtie2/human/7_111116_AD0341ACXX_137_5_index5_/7_111116_AD0341ACXX_137_5_index5__human.sample123.trinity.sam
samtools sort /gulo/proj_nobackup/b2011098/private/deNovoAnnotation/bowtie2/human/7_111116_AD0341ACXX_137_5_index5_/7_111116_AD0341ACXX_137_5_index5__human.sample123.trinity.bam /gulo/proj_nobackup/b2011098/private/deNovoAnnotation/bowtie2/human/7_111116_AD0341ACXX_137_5_index5_/7_111116_AD0341ACXX_137_5_index5__human.sample123.trinity.sorted 
samtools index /gulo/proj_nobackup/b2011098/private/deNovoAnnotation/bowtie2/human/7_111116_AD0341ACXX_137_5_index5_/7_111116_AD0341ACXX_137_5_index5__human.sample123.trinity.sorted.bam
samtools flagstat /gulo/proj_nobackup/b2011098/private/deNovoAnnotation/bowtie2/human/7_111116_AD0341ACXX_137_5_index5_/7_111116_AD0341ACXX_137_5_index5__human.sample123.trinity.sorted.bam >/gulo/proj_nobackup/b2011098/private/deNovoAnnot
ation/bowtie2/human/7_111116_AD0341ACXX_137_5_index5_/7_111116_AD0341ACXX_137_5_index5__human.sample123.trinity.sorted.bam.flagstat
samtools idxstats ....








