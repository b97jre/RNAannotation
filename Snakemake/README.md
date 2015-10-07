

#Steps to run the RNA predictor pipeline

To get all the information to run the classifier we are using multiple other software. These should be installed before you can run the pipeline.
The program requires some steps decribed below. 

If you do not have the folder structure specified below, run the make_dir_structure.sh script.

Before proceeding you need to have the the snakemake file, configuration file and the input fasta
file in the main_dir directory.

Steps 4-7 are not needed if working into the UPPMAX cluster. In that case the needed modules are going to be loaded from the run_snakemake.sh script.

1. Specify input sequence fasta file and main directory details in the configuration file.
2. Specify the reads' directory and details in configuration file.
3. Download the Rfam.cm and Pfam-A.hmm files and specify their paths in the configuration file.
4. Download the infernal software and specify the cmsearch directory in the configuration file.
5. Make sure to have HTStools.jar needed for the splitter command and the ORF predictor.
6. Make sure to have hmmer/3.1b1-gcc software needed for the Pfam predictor.
7. Download/load bowtie2 and samtools needed to run the expression predictor.
8. Run the splitter command for the ncRNA predictor.
9. Start the snakemake that will run the pipeline (recommended to put it as a batch job).




##Splitter command in UPPMAX
Will split the input fasta file per every 100000 nucleotides. The -n flag along with a number will allow you to specify the number of nucleotides in each split. 


**Code:**

```
    java -jar /glob/johanr/bin/HTStools.jar -p sequenceHandling splitSize -i inputfile.fasta -suffix fasta
```

**Out:**
*tmp/inputfile_i.fasta, i=increment number of subset split*

## Folder Structure

  * main_dir
     * IDX
     * BowtieRef
     * intermediate_files
        * orf_pfam
        * expression

You can automatically create the above folder structure by running the make_dir_structure.sh script provided.

##Description of fields in pipeline.conf


Parameter|Value|Comment
---|-----|-------
#GENERAL INFO||
SUFFIX|.fasta or .fa|End of input sequence file in below main directory 
PREFIX||Name of input file without suffix|
MAIN_DIR|/path/to/main_directory/| The directory which snakemake will run
#EXPRESSION LEVEL||
READS_DIR|/path/to/reads.fastq/|The path to the reads directory
BowtieRef|/path/to/bowtieRef| The path to the folder which creates the bowtie2 index and also names the ref.files
reads|os.path.basename(f) for f in glob.glob('/path/to/reads/*.fastq') | The path to the reads and will only take the basename
READs| os.path.splitext(k)[0] for k in reads | same as above and will give you sample_1.fastq
samplenames|os.path.splitext(k)[0] for k in READSs| Takes the sample names without the prefix
SAMPLE|L[:-2] for L in samplenames| Removes the last 2 characters from the sample name to remove the read direction for pair end reads.
IDX_OUT|/path/to/idxstats/folder/ |Folder for all the idxstats files created
EXP_OUT|/path/to/intermediate/files/| Folder for the other expression intermediate outputs
#ORF + PFAM||
PFAM_A|/path/to/Pfam-A.hmm|PFAM library that will be used to find potential protein domains
ORF_PFAM_OUT|/path/to/result|Folder for all the outputs from PFAM and ORF
ORF_ORFs|/path/to/orf/PREFIX.ORFs.SUFFIX|Output from ORF 
ORF_PEP|/path/to/orf/PREFIX.peptide.SUFFIX|Output from ORF 
ORF_INF|/path/to/orf/PREFIX.info|Output from ORF
ORF_TEMP_TABLE|/path/to/intermediate/ORF_temp_table.txt|Temp orf table which will be updated according to PFAM out
TEMP1|/path/to/intermediate/temp1.txt|Temporary files for the ORF table update
TEMP2|/path/to/intermediate/temp2.txt|Temporary files for the ORF table update
ORF_TABLE|/path/to/intermediate/ORF_table.txt| Final ORF table
PFAM_RAW|/path/to/intermediate/PfamA.hmm.hmmer| Raw PFAM data (result)
PFAM_TABLE|/path/to/intermediate/PFAM_table.txt| Top hits for each contig from PFAM raw data
ORF_PF_JOIN|/path/to/intermediate/ORF_PFAM_join.txt| The joined table from ORF and PFAM result
#ncRNA||
CMSEARCH|/path/to/infernal|path to infernal cmsearch function
RFAM_CMS|/path/to/RFAM/cmfile| path to the RFAM.cm file
SAMPLES|os.path.splitext(f)[0] for f in glob.glob('tmp/*')| Python command to get all the filenames of the split files created by the splitter command in the tmp/ folder
NCRNA_TABLE|/path/to/intermediate/ncRNA_Results.txt| Path to the ncRNA result file
#GATHERING OUTPUTS||
ORF_PF_EL_JOIN|/path/to/ORF_PF_EL_join.txt| Path to the join between expression table and ORF+Pfam table
TRAIN_TABLE||Name of the final result table with the join results from all the different parts


