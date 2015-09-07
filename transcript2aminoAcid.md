#Different code for how to run program



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

  **in:** 
  
  fastaFile, RFAM_database

  **out:**
  
  fastaFile.PfamAB.hmm.hmmer 


```

module load bioinfo-tools
module load hmmer/3.1b1-gcc

hmmscan --tblout test_0.PfamAB.hmm.hmmer /glob/johanr/references/Pfam/PfamAB.hmm test.fa
```

##Pfam using HMMer

  **in:** 
  
  fastaFile, RFAM_database

  **out:**






