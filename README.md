This is a github created for scripts that can help out in analysing de novo assembled transcriptomes.

For each script or file that is generated please write a nice markdown

https://help.github.com/articles/github-flavored-markdown/


This is the steps to run the RNA classifier.

The intent in the end is to run it on a set of RNA reads with meta information but right now we are trying to set it  up so that it takes a contigfile.

To get all the information to run the classifier we are using multiple other software. These should be installed before you can run the classifier. 

The program requires some steps. 
1. Create a table with the read files that you will use to identify the expression of the reads.This can be automatically done using a script. For more info see below. 
2. Run a Script that sets up the entirepath structure and makes sure that all the path names are correct.
3. Start the snakemake file that will run the pipeline from a set of reads and a file of contigs to get the data that you want and then categorize the contigs into coding, non-coding and other RNAs.




##Create read table
Will create a table with all the files that is below the **/path/to/reads/** that contains the suffix **fastq**. For paired end reads the flag **-suffix** has to be set and should contain the end of the files that seperates the files with the pairs of reads. 

**Code:**

```
    java -jar /glob/johanr/bin/SnakeMakeSetup.jar -p readTable -i /path/to/reads/ -suffix fastq -sep 1.fastq 2.fastq -o readsTable.tab.txt
```

**Out:**
*readsTable.tab.txt*

|SampleName|	forward|	reverse|
|----------|---------|---------|
|reads1	|/Users/johanreimegard/Vetenskap/Data/deNovoAnnotation/reads/reads1.1.fastq	|/Users/johanreimegard/Vetenskap/Data/deNovoAnnotation/reads/reads1.2.fastq|
|reads2	|/Users/johanreimegard/Vetenskap/Data/deNovoAnnotation/reads/reads2.1.fastq	|/Users/johanreimegard/Vetenskap/Data/deNovoAnnotation/reads/reads2.2.fastq|


##Create the snakemake configfile **pipeline.conf**

**in:**
ParameterFile.txt

Key|Value|Comment
---|-----|-------
#Reads info||
ReadsSuffix|.fastq|End of reads file in below reads directory 
readsTable|/path/to/reads.table.txt|Table with reads info|
#ReferenceInfo||
ReferenceFile|/path/to/refFile.fa|Reference file
ReferenceSuffix|.fa|End of reference file (most likely .fa or .fasta)
#Rfaminfo||
RfamLib|/path/to/Rfam.cm|RFAM library that will be used to find potential ncRNAs
cmsearchPath|/path/to/cmsearch|location of cmsearch program
#ProteinInfo||
ORFfinderPath|/path/to/HTStools.jar|Jar file for the ORF path finder
PfamLib|/path/to/Pfam-A.hmm"|PFAM library that will be used to find potential protein domains
hmmerPath|/path/to/hmmscan|Location of hmmer program 

**code:**

```
java -jar /glob/johanr/bin/SnakeMakeSetup.jar -p RNACLASSIFIER -wd /path/to/where/you/will/run/analysis -configFile /path/to/ParameterFile.txt
```


**out:**
1. **pipeline.conf**

2. Correct folder strutcture with soft links to files 
testrun
testrun/pipeline.conf
testrun/reads
testrun/reads/reads1.1.fastq
testrun/reads/reads1.2.fastq
testrun/reads/reads2.1.fastq
testrun/reads/reads2.2.fastq
testrun/reference
testrun/reference/fasta
testrun/reference/fasta/contigs.fa
testrun/reference/fasta/tmp
testrun/reference/fasta/tmp/contigs_0.fa
..


3. Contig file split up in smaller subset files of ~ 100 000 nt per file.

newly created folder structure that will contain the information that is needed 