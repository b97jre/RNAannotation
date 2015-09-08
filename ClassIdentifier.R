# 
library(ggplot2)

setwd("/Users/johanreimegard/Vetenskap/Data/deNovoAnnotation")


Coding = read.table("out.psl", sep = "\t", header =  FALSE,stringsAsFactors=FALSE)
Coding$kind = "mRNA"

ncRNA = read.table("ncRNA.psl", sep = "\t", header =  FALSE,stringsAsFactors=FALSE)
ncRNA$kind = "ncRNA"

BlastOutput = rbind (Coding,ncRNA) 
head(BlastOutput)
tail(BlastOutput)
BlastOutputOrdered = BlastOutput[with(BlastOutput, order(V2, -V12)), ]
head(BlastOutputOrdered)
BlastOutputOrderedSingle = BlastOutputOrdered[ !duplicated(BlastOutputOrdered$V2), ]
dim(BlastOutputOrderedSingle)
head(BlastOutputOrderedSingle)

Contigs = read.table("Contigs.txt", sep = "\t", header =  FALSE,stringsAsFactors=FALSE)
Contigs$kind = "other"
Contigs = Contigs[with(Contigs, order(V1)), ]

Contigs$kind[Contigs$V1 %in% BlastOutputOrderedSingle$V2  ] = BlastOutputOrderedSingle$kind
head(Contigs, n = 100)

write.table(Contigs, sep = "\t", quote = FALSE, file = "contig.classes.txt", row.names = FALSE,col.names = FALSE)

ORFinfo =  read.table("Human_all_table.txt", sep = " ", header =  TRUE,stringsAsFactors=FALSE)
head(ORFinfo)
ORFinfo = ORFinfo[with(ORFinfo, order(Name)), ]

dim(ORFinfo)
dim(Contigs)
ORFinfo$Class = Contigs$kind
colnames(ORFinfo)

mRNA <- ggplot(ORFinfo, aes(x = pValue))
mRNA + geom_density(aes(fill=factor(Class)), size=0.5)+facet_grid(Class ~ .) 
ggsave("ORF_p_valueDistribution.pdf")

ncRNA <- ggplot(ORFinfo, aes(x = (Rfam_E.value)))
ncRNA + geom_density(aes(fill=factor(Class)), size=0.5)+facet_grid(Class ~ .) 
ggsave("RfamDistributoin.pdf")

IPS <- ggplot(ORFinfo, aes(x = (IPS_E.value)))
IPS + geom_density(aes(fill=factor(Class)), size=0.5)+facet_grid(Class ~ .) 
ggsave("InterProScanDistribution.pdf")

expressionLevel <- ggplot(ORFinfo, aes(x = log(SC_human_index1.idxstats,base = 10)))
expressionLevel + geom_density(aes(fill=factor(Class)), size=0.5)+facet_grid(Class ~ .) 
ggsave("ExpressionDistribution.pdf")


