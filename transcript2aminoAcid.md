#Different code for how to run program
This is a github created for scripts that can help out in analysing de novo assembled transcriptomes.

For each script or file that is generated please write a nice markdown

https://help.github.com/articles/github-flavored-markdown/



##ORF finder

**In:**

test.fa (nt)

**Out:** 

test.info (table)
test.peptide.fa (aa) 
test.ORFs.fa  (cds)
test.523prime.fa (transcript 5 to 3 prime)

**Code:**

```
    java -jar /glob/johanr/bin/HTStools.jar -p sequencehandling orfs -i test.fa
```

##Mapping of reads

**in:** 

  reads + reference file

**out:**
  
  counts per contig per reads file ()

code:

```
  module load bioinfo-tools
  module load bowtie2/2.2.3
  module load samtools

  bowtie2  --threads 16 -x {REFfile} -1 {ReadsDir}/{sample}_1.fastq -2 {ReadsDir}/{sample}_2.fastq --un-conc {ResultDir}/{sample}.noHit.fastq -S {ResultDir}/{sample}.sam


  samtools view -bSh -o {ResultDir}/{sample}.bam {ResultDir}/{sample}.sam

  samtools sort {ResultDir}/{sample}.bam {ResultDir}/{sample}.sorted 

  samtools index {ResultDir}/{sample}.sorted.bam

  samtools flagstat {ResultDir}/{sample}.sorted.bam > {ResultDir}/{sample}.sorted.flagstat

  samtools idxstats {ResultDir}/{sample}.sorted.bam > {ResultDir}/{sample}.sorted.idxstats

```

##Rfam using Infernall

Information on how to run the latest version of Infernall is found here. 

http://infernal.janelia.org/

Especially the book chapter that goes through every step is found here.

http://selab.janelia.org/publications/Nawrocki13/Nawrocki13-preprint.pdf




##Pfam using HMMer

  **in:** 
  
  fastaFile, RFAM_database

  **out:**
  
  fastaFile.PfamAB.hmm.hmmer 


```

module load bioinfo-tools
module load hmmer/3.1b1-gcc

hmmscan --tblout test_0.PfamAB.hmm.hmmer /glob/johanr/references/Pfam/PfamAB.hmm test.fa
```




