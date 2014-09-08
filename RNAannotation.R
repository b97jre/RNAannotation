# TODO: Add comment
# 
# Author: johanreimegard
###############################################################################



# script to get the distribution size of ORFs in sequences

#install.packages("ggplot2") 
library("ggplot2")
#install.packages("reshape2") 
library("reshape2")

#install.packages("gridExtra")
library(gridExtra)
#install.packages("hexbin")
library(hexbin)

library(plyr)


cbPalette <- c(  "#E69F00", "#56B4E9", "#009E73", "#0072B2", "#F0E442", "#D55E00", "#CC79A7","#999999")

transcriptLength=200
dataDir=getwd()
lengthCutoff= 0.95
lengthCutoffFine = 0.5
breaksFine =100
breaks = 250
cutoff = 0.05

##datadir
dataDir = '/Users/johanreimegard/Vetenskap/Data/butterfly/deNovo'

## datafiles

GeneratedFileName = 'generated_0.info'
trueFileName = 'Butterfly.trinity.info'
positiveFileName = 'Homo_sapiens.GRCh37.74.cdna.all.info'

main <- function(dataDir=getwd()){
  
  EVDvalues <- getEVD(dataDir,GeneratedFileName1)
  
  info <- plotDeNovoInfo3(trueFileName,GeneratedFileName,positiveFileName,transcriptLength=100, cutoff=0.15 )
  write.classes(info, dataDir,trueFileName)
  
}

write.accepted <- function(sampleInfo ,outFileName , dataDir = getwd() , cutoff = 0.05 ,sample = "All", kind = "All"){
  
  if(sample != "All"){
    sampleInfo = sampleInfo[sampleInfo$Sample == sample, ]
  }
  
  if(kind != "All"){
    sampleInfo = sampleInfo[sampleInfo$kind == kind, ]
  }
  sampleInfo = sampleInfo[sampleInfo$pValue < cutoff, ]
  
  fileName = paste(outFileName, sample, kind, cutoff,"info", sep = ".")
  fileName = gsub(" ","_",fileName)
  
  write.table(sampleInfo,file= paste ( dataDir,fileName, sep= "/"), quote=FALSE,sep="\t",row.names=FALSE,col.names=TRUE)
  
  
  
  
}




plotDeNovoInfo2 <- function(trueFileName,GeneratedFileName,transcriptLength=200, dataDir=getwd(),cutoff=0.05,seq2=c(0,100,5000)){
  
  #plotHistograms(dataDir,GeneratedFileName1)
  plotHistogramsFromFile(dataDir,trueFileName)
  plotHistogramsFromFile(dataDir,GeneratedFileName)
  bp1 <- getDistribution(dataDir,GeneratedFileName)
  EVDvalues <- getEVD(dataDir,GeneratedFileName)
  #plotHistograms(dataDir,trueFileName,EVDvalues)
  
  names <- list()  	
  names[1] <- trueFileName
  names[2] <- GeneratedFileName
  deNovoInfo
  deNovoInfo <- calculatePValues(dataDir,trueFileName,EVDvalues,transcriptLength)
  GeneratedInfo <- calculatePValues(dataDir,GeneratedFileName,EVDvalues,transcriptLength)
  
  plotHistograms2(dataDir,list(deNovoInfo, GeneratedInfo),names)
  
  all <- splitTranscripts(dataDir,info,names,cutoff)
  distInfo <- plotAll(all,seq2=seq2)
  dev.off()
  
  ggplot(data=distInfo, aes(x=sample, y=count, fill=class)) +
    geom_bar(stat="identity")+
    theme(axis.text.x = element_text(colour = 'black', angle = 90))
  
  ggsave("ORF_classification2.pdf")  
  
  plotHistograms2(dataDir,info,names)
  return (all)
}


plotDeNovoInfo3 <- function(trueFileName,GeneratedFileName,positiveFileName,transcriptLength=200, dataDir=getwd(),lengthCutoff= 0.95,lengthCutoffFine = 0.5,breaksFine =100, breaks = 250 ,cutoff = 0.05){
  
  #plotHistograms(dataDir,GeneratedFileName1)
  EVDvalues <- getEVD(dataDir,GeneratedFileName)
  
  names <- list()  	
  names[1] <- trueFileName
  names[2] <- GeneratedFileName
  names[3] <- positiveFileName
  
  
  deNovoInfo <- calculatePValues(dataDir,trueFileName,EVDvalues,transcriptLength)
  deNovoInfo$Sample = "New Sample"
  
  GeneratedInfo <- calculatePValues(dataDir,GeneratedFileName,EVDvalues,transcriptLength)
  GeneratedInfo$Sample = "Negative Control"
  
  positiveInfo <- calculatePValues(dataDir,positiveFileName,EVDvalues,transcriptLength)
  positiveInfo$Sample = "Positive Controll"
  
  
  sampleInfo = rbind(deNovoInfo, GeneratedInfo,positiveInfo)
  sortedLength= sort(sampleInfo$SeqLength)
  maxLength = round_any(sortedLength[length(sortedLength)*lengthCutoff],1000,f=ceiling)
  maxLength2 = round_any(sortedLength[length(sortedLength)*lengthCutoffFine],1000,f=ceiling)
  sampleInfo$binFine = round_any(sampleInfo$SeqLength, breaksFine)
  sampleInfo$bin = round_any(sampleInfo$SeqLength, breaks)
  
  sampleInfo$decision = "accepted"
  sampleInfo$decision[which(sampleInfo$pValue >= cutoff)] = "rejected"
  sampleInfo$TUTRlength= sampleInfo$SeqLength-sampleInfo$stop
  sampleInfo$fraction = sampleInfo$ORFlength/sampleInfo$SeqLength
  sampleInfo$FUTRfraction = sampleInfo$start/sampleInfo$SeqLength
  sampleInfo$TUTRfraction = sampleInfo$TUTRlength/sampleInfo$SeqLength
  sampleInfo$GCfraction = sampleInfo$GC_count/sampleInfo$SeqLength
  sampleInfo$Sample = factor(sampleInfo$Sample, levels=c("Negative Control","New Sample","Positive Controll") ) 
  
  printOverview(sampleInfo,names)
  fileName = "test"
  
  samples = attributes(factor(sampleInfo$Sample))$levels
  
  for(i in 1:length(samples)){
    pdf(paste(paste(fileName, samples[i], sep = "."),"accepted.mRNAinfo.pdf", sep = "." ))
    plotSampleInfo(subset(sampleInfo,SeqLength < maxLength), samples[i],"accepted")
    dev.off()
    
    pdf(paste(paste(fileName, samples[i], sep = "."),"rejected.mRNAinfo.pdf", sep = "." ))
    plotSampleInfo(subset(sampleInfo,SeqLength < maxLength), samples[i],"rejected")
    dev.off()
    
    pdf(paste(paste(fileName, samples[i], sep = "."),"GC.mRNAinfo.pdf", sep = "." ))
    plotGCInfo(subset(sampleInfo,SeqLength < maxLength), samples[i])
    dev.off()
  }
  
  
  return (sampleInfo)
}






printGraphicalRepresentation <- function(){
  
  
  
  
  FULL = data.frame(fraction = c(0.05,0.25,0.75,0.95), kind = c("FULL","FULL","FULL","FULL"), size = c(1,2,1,1))
  ORF = data.frame(fraction = c(0.30,0.70), kind = c("ORF","ORF"), size = c(2,2))
  FUTR = data.frame(fraction = c(0.05,0.25,0.55), kind = c("5UTR","5UTR","5UTR"), size = c(1,2,1))
  TUTR = data.frame(fraction = c(0.45,0.75,0.95), kind = c("3UTR","3UTR","3UTR"), size = c(2,1,1))
  
  
  df = rbind (ORF,FULL,FUTR,TUTR)
  desired_order = sort(levels(as.factor(df$kind)))  
  df$kind  = factor(df$kind, levels=desired_order )  
  
  p <- ggplot(df, aes(x=fraction, y=kind, group=kind, color = kind))
  graphical_plot = p + geom_line(aes(size = size)) +
    expand_limits( x=c(-0.1,1.1))+ 
    scale_color_manual(values=cbPalette)+
    scale_size(guide = "none")+
    labs(title = "Graphical representation of different kind of RNAs", kind="Kind of RNA", x="fraction of a full length mRNA")+ 
    scale_x_continuous(breaks=c(0, 0.25, 0.5, 0.75, 1), 
                       labels=c("5'end", "AUG", 
                                "  ", "Stop", "3' end"))
  
  graphical_plot
  return(graphical_plot)
}


printOverview <- function(sampleInfo,names){
  
  
  
  
  graphicalPlot = printGraphicalRepresentation()
  
  m <- ggplot(sampleInfo, aes(x = SeqLength))
  SeqLenghPlot = m + stat_ecdf()+
    facet_grid( Sample ~ . )+ 
    coord_cartesian(xlim = c(0, 2000))+
    theme(text = element_text(size=8))+
    labs(title = "Cumulative distribution of RNA seq lesngth", y="Cumulative distribution", x="Sequence length")
  
  
  m <- ggplot(sampleInfo, aes(x = ORFlength))
  ORFLenghPlot = m + stat_ecdf()+
    facet_grid(Sample ~ . )+ coord_cartesian(xlim = c(0, 2000))+
    theme(text = element_text(size=8))+
    labs(title = "Cumulative distribution of ORF length", y="Cumulative distribution", x="ORF length")
  
  m <- ggplot(sampleInfo, aes(x = pValue))
  pValuePlot = m + geom_density(adjust=1/10)+
    scale_color_manual(values=cbPalette)+
    scale_size(guide = "none")+
    facet_grid(Sample ~ .  , scales = "free_y")+
    theme(text = element_text(size=8))+
    labs(title = "pValue distribution", y="Density", x="p-Value")
  
  
  
  levels(sampleInfo$kind)
  desired_order = sort(levels(as.factor(sampleInfo$kind)))  
  sampleInfo$kind  = factor(sampleInfo$kind, levels=desired_order )  
  
  countDist = qplot(kind, data = sampleInfo, geom="bar", fill = kind) +
    theme(text = element_text(size=8),legend.position = "left")+
    scale_fill_manual(values=cbPalette,breaks = sort(levels(sampleInfo$kind),decreasing = TRUE))+
    facet_grid(decision ~ Sample , scales = "free_y")+
    labs(title = "Distribution of RNAs", y="count")+theme(axis.text.x = element_text(colour = 'black', angle = 90))
  
  
  
  pdf("Overview.pdf")
  vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(8, 3)))
  print(graphicalPlot, vp = vplayout(1:2, 1:3))
  print(SeqLenghPlot, vp = vplayout(3:5, 1))
  print(ORFLenghPlot, vp = vplayout(3:5, 2))
  print(pValuePlot, vp = vplayout(3:5, 3))
  print(countDist, vp = vplayout(6:8, 1:3))
  
  dev.off()
  
  
}


printInfo <- function(fileName, sampleInfo,upperCutOff = 0.95, cutoff = 0.05,breaks = 250){  
  
  sortedLength= sort(sampleInfo$SeqLength)
  maxLength = round_any(sortedLength[length(sortedLength)*0.95],1000,f=ceiling)
  maxLength2 = round_any(sortedLength[length(sortedLength)*0.50],1000,f=ceiling)
  sampleInfo$roundedSeqLength = round_any(sampleInfoTrimmed$SeqLength, 100)
  sampleInfo$roundedSeqLength2 = round_any(sampleInfoTrimmed$SeqLength, breaks)
  
  #de novo info 
  
  
  
  printOverview(sampleInfo)
  
  
  
  m <- ggplot(subset(sampleInfoTrimmed,kind=="FULL"), aes(factor(roundedSeqLength2),fraction))
  futrPlot = m + geom_boxplot(aes(fill = factor(kind), color=factor(decision)))+
    theme(axis.text.x = element_text(colour = 'black', angle = 90))+
    facet_grid(Sample ~ kind )+
    scale_fill_manual(values=cbPalette[3])+
    scale_colour_manual(values=c("#000000","#999999"))+
    theme(text = element_text(size=8))+
    labs(title = "Boxplot of 5'UTR fraction of FULL RNAs", y="ORF length/Seq length", x="Seq length" )
  futrPlot
  
  
  
  
  
  m <- ggplot(subset(sampleInfo,kind != "ORF"), aes(factor(roundedSeqLength2),fraction))
  m + geom_boxplot(aes(fill = factor(decision)))+
    theme(axis.text.x = element_text(colour = 'black', angle = 90))+
    facet_grid(kind~Sample)+
    labs(title = "Boxplot of ORF fraction of full length", y="ORF length/Seq length", x="Seq length" )
  ggsave(paste(fileName,"ORF_fraction.pdf",sep="."))
  
  
  
  
  dist <- ggplot(sampleInfoTrimmed, aes(SeqLength, ..density.., colour = Sample, linetype = decision))
  dist + 
    geom_freqpoly(binwidth = 50)+
    facet_grid(kind~., scales="free_y",margins=TRUE)+
    labs(title = "Size distribution of sequences", x="Seq length")
  
  ggsave(paste(fileName,"SizeDistribution.pdf",sep="."))
  
  FUTRinfo <- sampleInfoTrimmed[sampleInfoTrimmed$kind =="FULL", ]
  FUTRinfo$kind <- as.factor(as.character(FUTRinfo$kind))
  m <- ggplot(FUTRinfo, aes(factor(roundedSeqLength2),FUTRfraction))
  m + geom_boxplot(aes(fill = factor(kind)))+
    theme(axis.text.x = element_text(colour = 'black', angle = 90))+
    facet_grid(decision~Sample)+
    labs(title = "Boxplot of 5'UTR fraction of full length RNAs", y="5'UTR length/Seq length", x="Seq length" )
  
  ggsave(paste(fileName,"FivePrimeFractionDistribution.pdf",sep="."))
  
  TUTRinfo <- sampleInfoTrimmed[sampleInfoTrimmed$kind =="FULL", ]
  TUTRinfo$kind <- as.factor(as.character(TUTRinfo$kind))
  m <- ggplot(TUTRinfo, aes(factor(roundedSeqLength2),TUTRfraction))
  m + geom_boxplot(aes(fill = factor(kind)))+
    theme(axis.text.x = element_text(colour = 'black', angle = 90))+
    facet_grid(decision~Sample)+
    labs(title = "Boxplot of 3'UTR fraction of full length RNAs", y="3'UTR length/Seq length", x="Seq length" )
  
  ggsave(paste(fileName,"ThreePrimeFractionDistribution.pdf",sep="."))
  
  
  
  dist <- ggplot(sampleInfoTrimmed, aes(pValue,..density.., colour=kind))
  dist +
    geom_freqpoly(binwidth = 0.01)+
    facet_grid(decision~Sample, scales="free_y")+
    labs(title = "pValue Distribution", x="p-value" )
  ggsave(paste(fileName,"PValueDistribution",sep="."))
  
  samples = attributes(factor(Allinfo$Sample))$levels
  
  for(i in 1:length(samples)){
    pdf(paste(paste(fileName, samples[i], sep = "."),"accepted.mRNAinfo.pdf", sep = "." ))
    plotSampleInfo(Allinfo, samples[i],"accepted")
    dev.off()
    
    pdf(paste(paste(fileName, samples[i], sep = "."),"rejected.mRNAinfo.pdf", sep = "." ))
    plotSampleInfo(Allinfo, samples[i],"rejected")
    dev.off()
    
    pdf(paste(paste(fileName, samples[i], sep = "."),"GC.mRNAinfo.pdf", sep = "." ))
    plotGCInfo(Allinfo, samples[i])
    dev.off()
  }
  dist <- ggplot(sampleInfoTrimmed, aes(pValue,..density..))
  
  gp2 = dist +
    geom_freqpoly(binwidth = 0.005)+
    facet_grid(Sample ~ ., scales="free_y")+
    labs(title = "pValue Distribution", x="p-value" )+theme(axis.text.x = element_text(colour = 'black', angle = 90))
  
  
  gp3 <- ggplot(sampleInfoTrimmed, aes(Sample,fill=kind)) + geom_bar()+facet_grid(decision ~ . )+theme(axis.text.x = element_text(colour = 'black', angle = 90))
  gp1 <- ggplot(sampleInfoTrimmed, aes(Sample,fill=kind)) + geom_bar()+facet_wrap(~ decision)
  
  grid.arrange(gp1,gp2,nrow=1)
  
  
}  



plotSampleInfo <- function(Allinfo, sampleName = "de Novo", verdict = "accepted"){
  
  
  ORFinfo <- ggplot(subset(subset(subset(Allinfo,kind != "ORF"), Sample==sampleName), decision == verdict), aes(factor(bin),fraction))
  ORFplot = ORFinfo + geom_boxplot(aes(fill = factor(kind)))+
    theme(axis.text.x = element_text(colour = 'black', angle = 90))+
    facet_grid(kind~.)+
    scale_fill_manual(values=cbPalette)+
    theme(text = element_text(size=8), legend.position="none")+
    labs(title = "Boxplot of ORF fraction of full length", y="ORF length/Seq length", x="Seq length" )
  
  m <- ggplot(subset(subset(subset(Allinfo, Sample==sampleName), decision == verdict),kind=="FULL"), aes(factor(bin),FUTRfraction))
  futrPlot = m + geom_boxplot(aes(fill = factor(kind)))+
    theme(axis.text.x = element_text(colour = 'black', angle = 90))+
    facet_grid(decision~.)+
    scale_fill_manual(values=cbPalette[3])+
    theme(text = element_text(size=8), legend.position="none")+
    labs(title = "Boxplot of 5'UTR fraction of FULL RNAs", y="5'UTR length/Seq length", x="Seq length" )
  
  m <- ggplot(subset(subset(subset(Allinfo, Sample==sampleName), decision == verdict),kind=="FULL"), aes(factor(bin),TUTRfraction))
  tutrPlot = m + geom_boxplot(aes(fill = kind))+
    theme(axis.text.x = element_text(colour = 'black', angle = 90))+
    facet_grid(decision~.)+
    scale_fill_manual(values=cbPalette[3])+   
    theme(text = element_text(size=8), legend.position="none")+
    labs(title = "Boxplot of 3'UTR fraction of FULL RNAs", y="3'UTR length/Seq length", x="Seq length" )
  
  ORFDistPlot =   qplot(ORFlength, data = subset(subset(Allinfo, Sample==sampleName), decision == verdict), geom = "freqpoly", binwidth = 50, colour = kind)+
    theme(text = element_text(size=8), legend.position="none")+
    scale_colour_manual(values=cbPalette)+
    coord_cartesian(xlim = c(0, maxLength))+
    labs(title = paste("Size distribution of ORFs in ", " RNAs", sep = verdict), y="count", x="Seq length" )
  
  sizeDistPlot =   qplot(SeqLength, data = subset(subset(Allinfo, Sample==sampleName), decision == verdict), geom = "freqpoly", binwidth = 50, colour = kind)+
    theme(text = element_text(size=8), legend.position="none")+
    scale_colour_manual(values=cbPalette)+
    coord_cartesian(xlim = c(0, maxLength))+
    labs(title =  paste("Size distribution of ", " RNAs", sep = verdict), y="count", x="Seq length" )
  
  
  countDist = qplot(kind, data = subset(subset(Allinfo, Sample==sampleName), decision == verdict), geom="bar", fill = kind) +
    theme(text = element_text(size=8),legend.position = "left")+
    scale_fill_manual(values=cbPalette)+
    labs(title = paste("Distribution of ", " RNAs", sep = verdict), y="count")+theme(axis.text.x = element_text(colour = 'black', angle = 90))
  
  
  vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(15, 6)))
  print(grid.text (paste("Sample:", sampleName,sep = " ") , vp = vplayout(1,2)))
  print(grid.text (paste("Filter:", verdict,sep = " ") , vp = vplayout(1, 5)))
  print(sizeDistPlot, vp = vplayout(2:4, 1:4))
  print(ORFDistPlot, vp = vplayout(5:7, 1:4))
  print(countDist, vp = vplayout(2:7, 5:6))
  print(ORFplot, vp = vplayout(8:15, 1:3))
  print(tutrPlot, vp = vplayout(8:11, 4:6))
  print(futrPlot, vp = vplayout(12:15, 4:6))
  
  
  
}

plotGCInfo <- function(Allinfo, sampleName = "de Novo"){
  
  m <- ggplot(subset(Allinfo, Sample==sampleName), aes(factor(bin),GCfraction))
  GCplot = m + geom_boxplot(aes(color = factor(decision)))+
    theme(axis.text.x = element_text(colour = 'black', angle = 90))+
    #    facet_grid(kind~decision.)+
    labs(title = "Boxplot of GC percentage distribtuion", y="% GC", x="Seq length" )
  
  d <- ggplot(subset(Allinfo, Sample==sampleName), aes(GCfraction,SeqLength))
  
  GC_Length_Plot = d + stat_binhex()+facet_grid(decision~kind)
  
  
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(7, 1)))
  print(grid.text (paste("Sample:", sampleName,sep = " ") , vp = vplayout(1,1)))
  print(GCplot, vp = vplayout(2:3, 1))
  print(GC_Length_Plot, vp = vplayout(4:7, 1))
  
}


vplayout <- function(x, y) viewport(layout.pos.row=x, layout.pos.col=y)



calculateDeNovoInfo <- function(trueFileName,GeneratedFileName,positiveFileName,transcriptLength=200, dataDir=getwd()){
  
  
  #plotHistograms(dataDir,GeneratedFileName1)
  plotHistogramsFromFile(dataDir,trueFileName)
  plotHistogramsFromFile(dataDir,GeneratedFileName)
  bp1 <- getDistribution(dataDir,GeneratedFileName)
  EVDvalues <- getEVD(dataDir,GeneratedFileName)
  #plotHistograms(dataDir,trueFileName,EVDvalues)
  
  
  names <- list()		
  names[1] <- trueFileName
  names[2] <- GeneratedFileName
  names[3] <- positiveFileName
  
  info <- list()
  info[1] <- list(calculatePValues(dataDir,trueFileName,EVDvalues,transcriptLength))
  info[2] <- list(calculatePValues(dataDir,GeneratedFileName,EVDvalues,transcriptLength))
  info[3] <- list(calculatePValues(dataDir,positiveFileName,EVDvalues,transcriptLength))
  
  return (info)
}

plotErrorDistribution <- function(stats,fileName="generated",seq2=c(0,100,4500)){
  matrix <- do.call(rbind,stats)
  matrix.sd <- apply(matrix,2,sd)
  matrix.mean <- apply(matrix,2,mean)
  classes2 <-c(seq(seq2[1],seq2[3],by=seq2[2]),4600)
  df <- data.frame(matrix.mean,matrix.sd,classes2)
  
  pdf(paste(fileName,'ORFlength.barplot.pdf',sep="."))
  ggplot(df, aes(classes2,matrix.mean*classes2))+geom_point()+geom_errorbar(aes(ymin=((matrix.mean-matrix.sd)*classes2),ymax=((matrix.mean+matrix.sd)*classes2)))+ylab("95 percentile mean")+xlab("length of transcript")+theme(title = expression("95 percent quantile length distriubtion of ORFs in generated sequences (n=5)"))
  dev.off()
  
  
  pdf(paste(fileName,'ORFlength.barplot.relative.pdf',sep="."))
  ggplot(df, aes(classes2,matrix.mean))+geom_point()+geom_errorbar(aes(ymin=((matrix.mean-matrix.sd)),ymax=((matrix.mean+matrix.sd))))+ylab("95 percentile mean")+xlab("length of transcript")+theme(title = expression("95 percent quantile length distriubtion of ORFs in generated sequences (n=5)"))
  dev.off()
}


splitTranscripts <- function(dataDir,sampleInfo,cutoff=0.05){
  
  ORFsrejected <- sampleInfo[(sampleInfo$kind == 'ORF'  & sampleInfo$pValue >= cutoff ),]
  ORFsaccepted <- sampleInfo[(sampleInfo$kind == 'ORF'  & sampleInfo$pValue<cutoff ),]
  FUTRsgood    <- sampleInfo[(sampleInfo$kind == '5UTR' & sampleInfo$pValue<cutoff ),]
  FUTRsbad <- sampleInfo[(sampleInfo$kind == '5UTR' &sampleInfo$pValue>=cutoff ),]
  TUTRsgood <- sampleInfo[(sampleInfo$kind == '3UTR' &sampleInfo$pValue<cutoff ),]
  TUTRsbad <- sampleInfo[(sampleInfo$kind == '3UTR'  &sampleInfo$pValue>=cutoff),]
  accepted<- sampleInfo[(sampleInfo$pValue<cutoff &sampleInfo$kind == 'FULL' ),]
  rejected<- sampleInfo[(sampleInfo$pValue>=cutoff &sampleInfo$kind =='FULL' ),]
  
  all <- list()
  all$accepted_ORFs<-ORFsaccepted
  all$accepted_5UTRs <-FUTRsgood
  all$accepted_3UTRs <-TUTRsgood
  all$accepted_mRNAs <-accepted
  all$rejected_ORFs<-ORFsrejected
  all$rejected_5UTRs <-FUTRsbad
  all$rejected_3UTRs <-TUTRsbad
  all$rejected_mRNAs <-rejected
  
  return(all)
  
}


splitDeNovoTranscripts <- function(dataDir,deNovoFrame,cutoff=0.05){
  
  all <- list()
  all$accepted_All<-deNovoFrame[deNovoFrame$pValue < cutoff,]
  all$accepted_ORFs <-  deNovoFrame[(deNovoFrame[7] == 'ORF'  & deNovoFrame$pValue<cutoff  ),]
  all$accepted_5UTRs <- deNovoFrame[(deNovoFrame[7] == '5UTR' & deNovoFrame$pValue < cutoff ),]
  all$accepted_3UTRs <- deNovoFrame[(deNovoFrame[7] == '3UTR' & deNovoFrame$pValue < cutoff ),]
  all$accepted_mRNAs <- deNovoFrame[(deNovoFrame[7] == 'FULL' & deNovoFrame$pValue < cutoff ),]
  all$rejected_All<-deNovoFrame[deNovoFrame$pValue >= cutoff,]
  all$rejected_ORFs <-  deNovoFrame[(deNovoFrame[7] == 'ORF'  & deNovoFrame$pValue >= cutoff ),]
  all$rejected_5UTRs <- deNovoFrame[(deNovoFrame[7] == '5UTR' & deNovoFrame$pValue >= cutoff ),]
  all$rejected_3UTRs <- deNovoFrame[(deNovoFrame[7] == '3UTR' & deNovoFrame$pValue >= cutoff),]
  all$rejected_mRNAs <- deNovoFrame[(deNovoFrame[7] =='FULL'  & deNovoFrame$pValue >= cutoff ),]
  
  return(all)
  
}






plotHistograms2 <- function(dataDir,dist,names){
  pdf(paste(dataDir,'ORF.pValueDistribution.EVD.pdf',sep="/"))
  par(mfrow=c(2,3))
  for(i in seq(1,length(dist), by=1)){
    hist(dist[[i]]$pValue, freq=FALSE, breaks=40,main=names[i])
  }
  dev.off()
}

plotDeNovoHistograms <- function(dataDir,dist,DeNovoDist){
  pdf(paste(dataDir,'ORF.pValueDistribution.EVD.pdf',sep="/"))
  par(mfrow=c(3,5))
  # plotAll
  hist(dist[[1]]$pValue, freq=FALSE, breaks=40,main = "Negative" )
  hist(dist[[1]]$pValue, freq=FALSE, breaks=40,main=  "de Novo"  )
  hist(dist[[1]]$pValue, freq=FALSE, breaks=40,main=  "Positive" )
  plot.new()
  plot.new()
  
  names[i] <- colnames(dist)
  for(i in seq(1,length(dist), by=1)){
    hist(dist[[i]]$pValue, freq=FALSE, breaks=40,main=names[i])
  }
  dev.off()
}

loadInfoFile <- function(fileName,dir=getwd()){
  data.file=file.path(dir,fileName)
  
  transcripts <- read.table(data.file,sep="\t",header = TRUE)
  transcripts <- transcripts[with(transcripts, order(SeqLength,ORFlength)),]
  
  return (transcripts)
}


loadBlastFile <- function(fileName,dir=getwd()){
  data.file=file.path(dir,fileName)
  
  blastTab <- read.table(data.file,sep="\t",header=TRUE)
  
  return (blastTab)
}


plotSizeDistributionFromFile <- function(dataDir,fileName,seq2=c(0,100,4500)){
  data.file=file.path(dataDir,fileName)
  
  transcripts <- read.table(data.file,sep="\t",header = TRUE)
  transcripts <- transcripts[with(transcripts, order(SeqLength,ORFlength)),]
  
  j = 1;
  distributions = list()
  classes = list()
  length = list()
  
  for(i in seq(seq2[1],seq2[3],by=seq2[2])){
    test.ordered.group <- transcripts[transcripts[[4]] < i+100 & transcripts[[4]]>i,]
    
    distributions[j] = list(test.ordered.group[[3]]/test.ordered.group[[4]])
    length[j] <- paste(i,(i+100),sep="-")
    classes[j] = list(table(test.ordered.group[[7]]))
    j=j+1
  }
  test.ordered.group <- transcripts[transcripts[[4]] > (seq2[3]+seq2[2]),]
  distributions[j] = list(test.ordered.group[[3]]/test.ordered.group[[4]])
  length[j] <- paste((seq2[3]+seq2[2]),max(transcripts[[4]]),sep="-")
  classes[j] = list(table(test.ordered.group[[7]]))
  j=j+1
  
  matrix <- do.call(rbind,classes)
  total = rowSums(matrix)-4
  
  pdf(paste(fileName,'sizeDistribution.barplot.pdf',sep="."))
  par(mfrow=c(2,1))
  
  barplot(t(matrix), names=length)
  barplot(total, names=length)
  
  dev.off()
}


getSizeHistogramFromFile <-function(fileName,dataDir=getwd(),seq2=c(0,300,10000),length=300,rowstart){
  
  data.file=file.path(dataDir,fileName)
  
  transcripts <- read.table(data.file,sep="\t",header = TRUE)
  transcriptsLengths <- transcripts[transcripts[,4]>length,4]
  histBreaks = c(seq(seq2[1],seq2[3],by=seq2[2]),max(c(transcriptsLengths,seq2[3]+seq2[2])))
  histInfo <- hist(transcriptsLengths,histBreaks)
  
  df = data.frame(row.names=histInfo$mids)
  df$mids = histInfo$mids
  df$counts <-histInfo$counts
  df$density <-histInfo$density
  df$name = fileName
  row.names(df)<-seq(rowstart+1,rowstart+length(df$name))	
  return (df)
}

getSizeHistogramsFromFile <-function(fileNames,dataDir=getwd(),seq2=c(0,300,10000),length=300){
  
  
  histInfo <- getSizeHistogramFromFile(fileNames[1],dataDir,seq2,length,0)
  for(i in 2:length(fileNames)){
    temp <- getSizeHistogramFromFile(fileNames[i],dataDir,seq2,length,length(histInfo$name))
    histInfo <- rbind(histInfo,temp)
  }
  return(histInfo)
  
}



plotHistogramsFromFile <- function(dataDir,fileName,seq2=c(0,100,4500)){
  data.file=file.path(dataDir,fileName)
  
  transcripts <- read.table(data.file,sep="\t",header = TRUE)
  transcripts <- transcripts[with(transcripts, order(SeqLength,ORFlength)),]
  
  j = 1;
  distributions = list()
  classes = list()
  length = list()
  
  for(i in seq(seq2[1],seq2[3],by=seq2[2])){
    test.ordered.group <- transcripts[transcripts[[4]] < i+100 & transcripts[[4]]>i,]
    
    distributions[j] = list(test.ordered.group[[3]]/test.ordered.group[[4]])
    length[j] <- paste(i,(i+100),sep="-")
    classes[j] = list(table(test.ordered.group[[7]]))
    j=j+1
  }
  test.ordered.group <- transcripts[transcripts[[4]] > (seq2[3]+seq2[2]),]
  distributions[j] = list(test.ordered.group[[3]]/test.ordered.group[[4]])
  length[j] <- paste((seq2[3]+seq2[2]),max(transcripts[[4]]),sep="-")
  classes[j] = list(table(test.ordered.group[[7]]))
  j=j+1
  
  
  matrix <- do.call(rbind,classes)
  pdf(paste(fileName,'distribution.boxplot.pdf',sep="."))
  par(mfrow=c(2,1))
  
  boxplot(distributions, names=length)
  bp <-barplot(t(matrix), names=length)
  gcol = gray.colors(4)
  
  legend("topright",legend=colnames(matrix), fill = gcol)
  dev.off()
}

plotAll <- function (all,seq2=c(0,100,3000)){
  plotRNAfeatures(all,seq2)
  plotClassCount(all,seq2)
}

plotRNAfeatures<- function(all,seq2=c(0,100,3000)){
  classes<-names(all)
  samples<-names(all[[1]])
  
  for(i in seq(1,length(samples))){
    DistributionInfos <- list()
    for(j in seq(1,length(classes))){
      fileName <-paste(samples[i],classes[j],sep="_")
      fileName <-paste(fileName,"distribution.pdf",sep=".")
      DistributionInfos[classes[j]] <- list(plotHistograms(fileName,seq2,all[[j]][[i]]))
    }
    fileName <-paste(fileName,"summary.distribution.pdf",sep=".")
    plotSampleDistributionSummary(DistributionInfos)
    
  }
}


plotClassCount <- function(fileName,DistributionInfos){
  classes <- names(DistributionInfos)
  bin <- DistributionInfos[[classes[1]]]$ORF
  ORFmean <- DistributionInfos[[classes[1]]]$ORF$stats[3,]
  TUTRmean <- DistributionInfos[[classes[1]]]$TUTR$stats[3,]
  FUTRmean <- DistributionInfos[[classes[1]]]$FUTR$stats[3,]
  count <- DistributionInfos[[classes[1]]]$FUTR$n
  classNames <- rep(classes[1],length(DistributionInfos[[classes[1]]]$FUTR$n))
  
  xlabels  = DistributionInfos[[classes[1]]]$ORF$names
  
  for(j in seq(2,length(classes))){
    ORFmean <- cbind(ORFmean,DistributionInfos[[classes[j]]]$ORF$stats[3,])
    TUTRmean <- cbind(TUTRmean,DistributionInfos[[classes[j]]]$TUTR$stats[3,])
    FUTRmean <- cbind(FUTRmean,DistributionInfos[[classes[j]]]$FUTR$stats[3,])
    count <- cbind(count,DistributionInfos[[classes[j]]]$FUTR$n)
    classNames <- cbind(classNames,rep(classes[j],length(DistributionInfos[[classes[j]]]$FUTR$n)))
    
  }
  
  DI <- cbind(as.numeric(ORFmean), as.numeric(TUTRmean),as.numeric(FUTRmean),as.numeric(count),as.factor(classNames))
  colnames(DI) <- c("ORF","3UTR", "5UTR","count","classNames")
  
  
  
  
}


plotClassCount <- function(all,seq2=c(0,100,3000)){
  classes<-names(all)
  samples<-names(all[[1]])
  
  distribution <-matrix(0,nrow=length(samples),ncol=length(classes))
  for(i in seq(1,length(samples))){
    for(j in seq(1,length(classes))){
      distribution[i,j]=length(all[[j]][[i]][[1]])
    }
  }
  relativeDistribution = distribution/rowSums(distribution)
  
  dist2 <-matrix(0,nrow=length(samples)*length(classes),ncol=2)
  count2 <-matrix(0,nrow=length(samples)*length(classes),ncol=1)
  for(i in seq(1,length(samples))){
    for(j in seq(1,length(classes))){
      dist2[(i-1)*length(classes)+j,1]=samples[i]
      dist2[(i-1)*length(classes)+j,2]=classes[j]
      count2[(i-1)*length(classes)+j,1]=relativeDistribution[i,j]
    }
  }
  distribution2 <- as.data.frame(dist2)
  count3 <- as.data.frame(count2)
  dist3 <- cbind(count3,distribution2)
  colnames(dist3)<-c("count","sample","class")
  dist3 <- dist3[with(dist3, order(sample,class)),]
  
  ggplot(data=dist3, aes(x=sample, y=count, fill=class)) +
    geom_bar(stat="identity")+
    theme(axis.text.x = element_text(colour = 'black', angle = 90))
  ggsave("ORF_classification.pdf")
  
  return(dist3)
  
  
  
}


plotHistograms <- function(fileName,seq2=c(0,100,3000),transcripts){
  j = 1;
  distributions = list()
  TUTR=list()
  FUTR=list()
  classes = list()
  length = list()
  
  for(i in seq(seq2[1],seq2[3],by=seq2[2])){
    test.ordered.group <- transcripts[transcripts[[4]] < i+seq2[2] & transcripts[[4]]>i,]
    
    distributions[j] = list(test.ordered.group[[3]]/test.ordered.group[[4]])
    TUTR[j] = list((test.ordered.group[[4]]-test.ordered.group[[6]]))
    FUTR[j] = list(test.ordered.group[[5]])
    
    length[j] <- paste(i,(i+seq2[2]),sep="-")
    classes[j] = list(table(test.ordered.group[[7]]))
    j=j+1
  }
  test.ordered.group <- transcripts[transcripts[[4]] > (seq2[3]+seq2[2]),]
  #distributions[j] = list(test.ordered.group[[3]]/test.ordered.group[[4]])
  #TUTR[j] = list((test.ordered.group[[4]]-test.ordered.group[[6]])/test.ordered.group[[4]])
  #FUTR[j] = list(test.ordered.group[[5]]/test.ordered.group[[4]])
  
  #	length[j] <- paste((seq2[3]+seq2[2]),max(transcripts[[4]]),sep="-")
  classes[j] = list(table(test.ordered.group[[7]]))
  j=j+1
  
  matrix <- do.call(rbind,classes)
  pdf(paste(fileName,'distribution.boxplot.pdf',sep="."))
  par(mfrow=c(4,1))
  
  ORF <-boxplot(distributions, names=length, main="Boxplot of the fraction of ORF length compared to the transcript length dependent on transtript lengths",xlab="Transcript length",ylab="Fraction")
  
  FUTR <-boxplot(FUTR, names=length,ylim=c(0,2000), main="Boxplot of average 5' UTR lengths dependent on transcript length",xlab="Transcript length",ylab="5' UTR length")
  TUTR <- boxplot(TUTR, names=length,ylim=c(0,2000), main="Boxplot of average 3' UTR lengths dependent on transcript length",xlab="Transcript length",ylab="3' UTR length")
  LengthDistribution <- hist(transcripts[[4]],breaks=c(seq(seq2[1],seq2[3],by=seq2[2]),max(transcripts[[4]])),xlim=c(0,(seq2[3]-(seq2[2]*3))),main="Histogram transcript length",freq=TRUE,xlab="Transcript length")
  dev.off()
  DistributionInfo <- list()
  DistributionInfo$ORF <- ORF
  DistributionInfo$FUTR <- FUTR
  DistributionInfo$TUTR <- TUTR
  DistributionInfo$LengthDistribution <- LengthDistribution
  
  return(DistributionInfo)
  
  
  
  
  
}



calculatePValues <- function(dataDir,fileName,EVD,TranscriptLength=400){
  data.file=file.path(dataDir,fileName)
  
  transcripts <- read.table(data.file,sep="\t",header = TRUE)
  transcripts <- transcripts[with(transcripts, order(SeqLength,ORFlength)),]
  transcripts <- transcripts[transcripts[[4]] >TranscriptLength,]
  j = 1;
  #	transcripts$distributions <- apply(transcripts[ ,3:4],1,function(row) row[[1]]/row[[2]])
  distributions <- transcripts[[3]]/transcripts[[4]]
  transcripts$distributions <- distributions
  slope <- EVD[3]+transcripts[[4]]*EVD[4]
  intercept <- EVD[1]+transcripts[[4]]*EVD[2]
  Theta <- -1/slope
  Xi <- Theta*intercept
  pValue <-1 - exp(-exp(-((distributions-Xi)/Theta)))
  transcripts$pValue <- pValue
  return(transcripts)
}




splitORFs <- function(boxplot,dataDir,fileName){
  data.file=file.path(dataDir,fileName)
  
  transcripts <- read.table(data.file,sep="\t",header = TRUE)
  transcripts <- transcripts[with(transcripts, order(SeqLength,ORFlength)),]
  transcripts$GoB =transcripts[[7]]
  
  
  pointer = 1
  cutoff = 400
  j = 1
  
  
  DF$newcol <- apply(DF,1,function(row) mean(vec[ row[1] : row[2] ] ))
  if(transcripts[[3]][j]/transcripts[[4]][j]>boxplot$stats[5,pointer]){
    transcripts$GoB[[j]] =1
  }else{
    transcripts$GoB[[j]] =0
  }
  print(j)
  j = j+1
}
categoizeORFS <- function(row, boxplot){
  newValue <- c(0)
  
  if(row[3]/row[4]>boxplot$stats[5,row[4]/100]){
    
    
    
  }
  
  
  pointer = pointer+1
  cutoff = cutoff + 100
  
  while(cutoff < 8501 & j < length(transcripts[[1]])){
    while(transcripts[[4]][j] <cutoff & j < length(transcripts[[1]])){
      if(transcripts[[3]][j]/transcripts[[4]][j]>boxplot$stats[5,pointer]){
        transcripts$GoB[j] =1
      }else{
        transcripts$GoB[j] =0
      }
      j = j+1
    }
    pointer = pointer+1
    cutoff = cutoff + 100
  }
  while( j < length(transcripts[[1]])){
    if(transcripts[[3]][j]/transcripts[[4]][j]>boxplot$stats[5,pointer]){
      transcripts$GoB[j] =1
    }else{
      transcripts$GoB[j] =0
    }
    j = j+1
  }
  
  transcripts
  
}

getDistribution<- function(dataDir,fileName,seq2=c(0,100,4500)){
  data.file=file.path(dataDir,fileName)
  
  transcripts <- read.table(data.file,sep="\t",header = TRUE)
  transcripts <- transcripts[with(transcripts, order(SeqLength,ORFlength)),]
  
  j = 1;
  distributions = list()
  lengths2 = list()
  classes = list()
  length = list()
  info = list()
  pdf(paste(fileName,"ExtremeValueDistribution.pdf",sep="."))
  par(mfrow=c(2,2))
  
  PLOT = TRUE
  for(i in seq(seq2[1],seq2[3],by=seq2[2])){
    test.ordered.group <- transcripts[transcripts[[4]] < i+25 & transcripts[[4]]>i-25,]
    distributions[j] = list(sort(test.ordered.group[[3]]/test.ordered.group[[4]]))
    l1 <- length(distributions[[j]])
    if(l1 >500){
      l2 <- c(seq(1/l1,1-1/l1,by=1/l1),1-1/l1)
      dLOne <-distributions[[j]]<1
      d1 <- distributions[[j]]
      d2<- d1[dLOne]
      ll = log(-log(l2))
      ll2 <- ll[1:length(d2)]
      linModel2 <- lm(ll2~sort(d2))
      if(PLOT){
        plot(sort(d2),ll2)
        PLOT=FALSE
      }
      summ <- summary(linModel2)
      info[j] <- list(c(linModel2[[1]][[1]],linModel2[[1]][[2]],summ[[9]]))
    }
    else{
      info[j] <-list(c(0,0,0))
    }
    length[j] <- i
    classes[j] = list(table(test.ordered.group[[7]]))
    j=j+1
  }
  
  
  
  infoMatrix <- matrix(unlist(info),3,length(unlist(info))/3)
  keep <- infoMatrix[1,]>0
  intercept <- infoMatrix[1,][keep]
  slope <- infoMatrix[2,][keep]
  Theta <- -1/slope
  Xi <- Theta*intercept
  
  Rsquare <- infoMatrix[3,][keep]
  TranscriptLength <- unlist(length)[keep]
  lmIntercept <- lm(intercept~TranscriptLength)
  lmSlope <- lm(slope~TranscriptLength)
  
  
  
  plot(TranscriptLength,intercept)
  abline(lmIntercept[[1]][[1]], lmIntercept[[1]][[2]])
  
  plot(TranscriptLength,slope)
  abline(lmSlope[[1]][[1]], lmSlope[[1]][[2]])
  plot(TranscriptLength,Rsquare)
  dev.off()
  
  test.ordered.group <- transcripts[transcripts[[4]]>(seq2[3]+seq2[2]),]
  distributions[j] = list(test.ordered.group[[3]]/test.ordered.group[[4]])
  length[j] <- paste((seq2[3]+seq2[2]),max(transcripts[[4]]),sep="-")
  classes[j] = list(table(test.ordered.group[[7]]))
  j=j+1
  
  boxplotInfo <- boxplot(distributions, names=length,plot=F)
  
  return(boxplotInfo)
  
  
}

getEVD <- function(dataDir,fileName,seq2=c(0,100,4500)){
  data.file=file.path(dataDir,fileName)
  
  transcripts <- read.table(data.file,sep="\t", header=TRUE)
  transcripts <- transcripts[with(transcripts, order(SeqLength,ORFlength)),]
  
  j = 1;
  distributions = list()
  lengths2 = list()
  classes = list()
  length = list()
  info = list()
  pdf(paste(fileName,"ExtremeValueDistribution.pdf",sep="."))
  par(mfrow=c(3,2))
  
  PLOT = TRUE
  for(i in seq(seq2[1],seq2[3],by=50)){
    test.ordered.group <- transcripts[transcripts[[4]] < i+25 & transcripts[[4]]>i-25,]
    distributions[j] = list(sort(test.ordered.group[[3]]/test.ordered.group[[4]]))
    l1 <- length(distributions[[j]])
    if(l1 >500){
      l2 <- c(seq(1/l1,1-1/l1,by=1/l1),1-1/l1)
      dLOne <-distributions[[j]]<1
      d1 <- distributions[[j]]
      d2<- d1[dLOne]
      l3<-l2[1:length(d2)]
      ll = log(-log(l2))
      ll2 <- ll[1:length(d2)]
      linModel2 <- lm(ll2~sort(d2))
      if(i==500){
        plot(density(distributions[[j]]))
        plot(sort(d2),l3)
        plot(sort(d2),ll2)
        PLOT=FALSE
      }
      summ <- summary(linModel2)
      info[j] <- list(c(linModel2[[1]][[1]],linModel2[[1]][[2]],summ[[9]]))
    }
    else{
      info[j] <-list(c(0,0,0))
    }
    length[j] <- i
    classes[j] = list(table(test.ordered.group[[7]]))
    j=j+1
  }
  
  
  
  infoMatrix <- matrix(unlist(info),3,length(unlist(info))/3)
  keep <- infoMatrix[1,]>0
  intercept <- infoMatrix[1,][keep]
  slope <- infoMatrix[2,][keep]
  Theta <- -1/slope
  Xi <- Theta*intercept
  
  Rsquare <- infoMatrix[3,][keep]
  TranscriptLength <- unlist(length)[keep]
  lmIntercept <- lm(intercept~TranscriptLength)
  lmSlope <- lm(slope~TranscriptLength)
  
  plot(TranscriptLength,Rsquare)
  
  
  plot(TranscriptLength,intercept)
  abline(lmIntercept[[1]][[1]], lmIntercept[[1]][[2]])
  
  plot(TranscriptLength,slope)
  abline(lmSlope[[1]][[1]], lmSlope[[1]][[2]])
  dev.off()
  
  test.ordered.group <- transcripts[transcripts[[4]]>(seq2[3]+seq2[2]),]
  distributions[j] = list(test.ordered.group[[3]]/test.ordered.group[[4]])
  length[j] <- paste((seq2[3]+seq2[2]),max(transcripts[[4]]),sep="-")
  classes[j] = list(table(test.ordered.group[[7]]))
  j=j+1
  info <- c(II=lmIntercept[[1]][[1]], IS=lmIntercept[[1]][[2]],SI=lmSlope[[1]][[1]],SS=lmSlope[[1]][[2]] )
  
  
}




plotSensitivity <- function(all,blastInfo,pointer){
  #t2 <- read.table(file='arabidopsis.ATH_cDNA.blast.bestCoverage',header=FALSE)
  #t1 <- read.table(file='arabidopsis.oases.blast.bestCoverage',header=FALSE)
  blastInfo$Info = "unsorted"
  for(i in seq(1:length(all))){
    subset <- all[[i]][[pointer]][[1]]	
    blastInfo$Info[is.element(as.character(blastInfo[[1]]),sub(">","",subset))]=names(all)[i]
  }
}





addPfamInfo <- function(fileName,dir=getwd(),infoDF){
  
  data.file=file.path(dir,fileName)
  pfamInfo <- read.table(data.file, header=FALSE)
  pfamInfo2 <- pfamInfo[with(pfamInfo,order(V13)),]
  pfamInfoBest <- pfamInfo2[!duplicated(pfamInfo2[[1]]),]
  pfamInfoBest[[1]] = gsub(pfamInfoBest[[1]],pattern="_forward",replacement="")
  pfamInfoBest[[1]] = gsub(pfamInfoBest[[1]],pattern="_reverse",replacement="")
  infoDF$pfam = 1
  infoDF$pfam[match(pfamInfoBest[[1]],infoDF$V1)]=pfamInfoBest[[13]]
  
  
}



plotClassification <-function(infoDF,fileName = "classification.pdf"){
  pdf(fileName)
  
  
  dev.off()
}


addClassInfo <- function(infoDF,cutoff=-3){
  infoDF$classification="NaN"
  infoDF$classification[(infoDF[7] == 0 )]="ORF"
  infoDF$classification[(infoDF[7] == 1 ) & infoDF$class_Score<cutoff]="accepted 5UTR"
  infoDF$classification[(infoDF[7] == 2 ) & infoDF$class_Score<cutoff]="accepted 3UTR"
  infoDF$classification[(infoDF[7] == 3 ) & infoDF$class_Score<cutoff]="accepted mRNA"
  infoDF$classification[(infoDF[7] == 1 ) & infoDF$class_Score>=cutoff]="rejected 5UTR"
  infoDF$classification[(infoDF[7] == 2 ) & infoDF$class_Score>=cutoff]="rejected 3UTR"
  infoDF$classification[(infoDF[7] == 3 ) & infoDF$class_Score>=cutoff]="rejected mRNA"
  
  infoDF$classification2="NaN"
  infoDF$classification2[(infoDF[7] == 0 )]="ORF"
  infoDF$classification2[(infoDF[7] == 1 ) & infoDF$class_Score<cutoff]="accepted"
  infoDF$classification2[(infoDF[7] == 2 ) & infoDF$class_Score<cutoff]="accepted"
  infoDF$classification2[(infoDF[7] == 3 ) & infoDF$class_Score<cutoff]="accepted"
  infoDF$classification2[(infoDF[7] == 1 ) & infoDF$class_Score>=cutoff]="rejected"
  infoDF$classification2[(infoDF[7] == 2 ) & infoDF$class_Score>=cutoff]="rejected"
  infoDF$classification2[(infoDF[7] == 3 ) & infoDF$class_Score>=cutoff]="rejected"
  infoDF$classification2 <- as.factor(infoDF$classification2)
  infoDF$classification <- as.factor(infoDF$classification)
  
  return (infoDF)
}



addClassInfopValue <- function(infoDF,cutoff=0.05){
  infoDF$classification="NaN"
  infoDF$classification[(infoDF[7] == 0 )]="ORF"
  infoDF$classification[(infoDF[7] == 1 ) & infoDF$pValue<cutoff]="accepted 5UTR"
  infoDF$classification[(infoDF[7] == 2 ) & infoDF$pValue<cutoff]="accepted 3UTR"
  infoDF$classification[(infoDF[7] == 3 ) & infoDF$pValue<cutoff]="accepted mRNA"
  infoDF$classification[(infoDF[7] == 1 ) & infoDF$pValue>=cutoff]="rejected 5UTR"
  infoDF$classification[(infoDF[7] == 2 ) & infoDF$pValue>=cutoff]="rejected 3UTR"
  infoDF$classification[(infoDF[7] == 3 ) & infoDF$pValue>=cutoff]="rejected mRNA"
  
  infoDF$classification2="NaN"
  infoDF$classification2[(infoDF[7] == 0 )]="ORF"
  infoDF$classification2[(infoDF[7] == 1 ) & infoDF$pValue<cutoff]="accepted"
  infoDF$classification2[(infoDF[7] == 2 ) & infoDF$pValue<cutoff]="accepted"
  infoDF$classification2[(infoDF[7] == 3 ) & infoDF$pValue<cutoff]="accepted"
  infoDF$classification2[(infoDF[7] == 1 ) & infoDF$pValue>=cutoff]="rejected"
  infoDF$classification2[(infoDF[7] == 2 ) & infoDF$pValue>=cutoff]="rejected"
  infoDF$classification2[(infoDF[7] == 3 ) & infoDF$pValue>=cutoff]="rejected"
  infoDF$classification2 <- as.factor(infoDF$classification2)
  infoDF$classification <- as.factor(infoDF$classification)
  
  return (infoDF)
}



addClassInfopfam <- function(infoDF,cutoff=0.00001){
  infoDF$classification="NaN"
  infoDF$classification[(infoDF[7] == 0 )]="ORF"
  infoDF$classification[(infoDF[7] == 1 ) & infoDF$pfam<cutoff]="accepted 5UTR"
  infoDF$classification[(infoDF[7] == 2 ) & infoDF$pfam<cutoff]="accepted 3UTR"
  infoDF$classification[(infoDF[7] == 3 ) & infoDF$pfam<cutoff]="accepted mRNA"
  infoDF$classification[(infoDF[7] == 1 ) & infoDF$pfam>=cutoff]="rejected 5UTR"
  infoDF$classification[(infoDF[7] == 2 ) & infoDF$pfam>=cutoff]="rejected 3UTR"
  infoDF$classification[(infoDF[7] == 3 ) & infoDF$pfam>=cutoff]="rejected mRNA"
  
  infoDF$classification2="NaN"
  infoDF$classification2[(infoDF[7] == 0 )]="ORF"
  infoDF$classification2[(infoDF[7] == 1 ) & infoDF$pfam<cutoff]="accepted"
  infoDF$classification2[(infoDF[7] == 2 ) & infoDF$pfam<cutoff]="accepted"
  infoDF$classification2[(infoDF[7] == 3 ) & infoDF$pfam<cutoff]="accepted"
  infoDF$classification2[(infoDF[7] == 1 ) & infoDF$pfam>=cutoff]="rejected"
  infoDF$classification2[(infoDF[7] == 2 ) & infoDF$pfam>=cutoff]="rejected"
  infoDF$classification2[(infoDF[7] == 3 ) & infoDF$pfam>=cutoff]="rejected"
  infoDF$classification2 <- as.factor(infoDF$classification2)
  infoDF$classification <- as.factor(infoDF$classification)
  
  return (infoDF)
}


addBlastInfo <- function(fileName,dir=getwd(),infoDF){
  data.file=file.path(dir,fileName)
  blastInfo <- read.table(data.file, header=TRUE)
  infoDF$blast=blastInfo[[3]][match(infoDF$V1,blastInfo[[1]])]
  infoDF$blast=blastInfo[[3]][match(infoDF$V1,blastInfo[[1]])]
  infoDF$blast=blastInfo[[3]][match(infoDF$V1,blastInfo[[1]])]
  infoDF$blast=blastInfo[[3]][match(infoDF$V1,blastInfo[[1]])]
  infoDF$blast=blastInfo[[3]][match(infoDF$V1,blastInfo[[1]])]
  
}

addBlastInfo <- function(fileName,dir=getwd(),infoDF){
  data.file=file.path(dir,fileName)
  blastInfo <- read.table(data.file, header=FALSE)
  infoDF$blast = 0
  infoDF$blast[match(blastInfo[[1]],infoDF$V1)]=pfamInfoBest[[3]]
  return (infoDF)
}

printDistributions <- function(infoDF,fileName,dir=getwd()){
  fileName1 <- paste(fileName,"pValue.pdf",sep="_")
  pdf(fileName1)
  infoDF = addClassInfopValue(infoDF)
  infoDF_noORF<- infoDF[infoDF$classification!="ORF" ,]
  ggplot(infoDF_noORF,aes(blast))+geom_bar() +facet_wrap(~ classification)
  dev.off()
  
  fileName1 <- paste(fileName,"pfam.pdf",sep="_")
  pdf(fileName1)
  infoDF = addClassInfopfam(infoDF)
  infoDF_noORF<- infoDF[infoDF$classification!="ORF" ,]
  ggplot(infoDF_noORF,aes(blast))+geom_bar() +facet_wrap(~ classification)
  dev.off()
  
  fileName1 <- paste(fileName,"both.pdf",sep="_")
  pdf(fileName1)
  infoDF = addClassInfo(infoDF)
  infoDF_noORF<- infoDF[infoDF$classification!="ORF" ,]
  ggplot(infoDF_noORF,aes(blast))+geom_bar() +facet_wrap(~ classification)
  dev.off()
  
}







displayDistribution <-function(infoQuery,infoHit,blastInfo,Name){
  forwardORF <- matrix(ncol=100,nrow=nlevels(infoQuery$classification),data=0)
  reverseORF <- matrix(ncol=100,nrow=nlevels(infoQuery$classification),data=0)
  forward5UTR <- matrix(ncol=100,nrow=nlevels(infoQuery$classification),data=0)
  reverse5UTR <- matrix(ncol=100,nrow=nlevels(infoQuery$classification),data=0)
  forward3UTR <- matrix(ncol=100,nrow=nlevels(infoQuery$classification),data=0)
  reverse3UTR <- matrix(ncol=100,nrow=nlevels(infoQuery$classification),data=0)
  classes <-  levels(infoQuery$classification)
  
  
  blastInfo <- blastInfo[blastInfo$BitScore>0,]
  
  for(i in seq(1,nlevels(infoQuery$classification))){
    infoQuerytemp <- infoQuery[infoQuery['classification']==classes[i],]
    blastInfoTemp <- blastInfo[match(infoQuerytemp[[1]],blastInfo$Query),]
    infoHitTemp<-infoHit[match(blastInfoTemp$Hit,infoHit[[1]]),]
    stepNr = 10;
    printNr = 10;
    for(j in seq(1,length(blastInfoTemp[[1]]))){
      if(!is.na(blastInfoTemp$queryStart[[j]])){
        Qstart = blastInfoTemp$queryStart[[j]]
        Qstop = blastInfoTemp$queryStop[[j]]
        Hstart = blastInfoTemp$hitStart[[j]]
        Hstop = blastInfoTemp$hitStop[[j]]
        QI =infoQuerytemp[infoQuerytemp$V1==blastInfoTemp$Query[[j]],]
        HI =infoHitTemp[match(blastInfoTemp$Hit[[j]],infoHitTemp$V1),]
        if(!is.na(HI$V2)){
          #print(blastInfoTemp[j,])
          #print(HI)
          #print(QI)
          if(j > printNr){	
            print(j) 
            printNr =printNr+stepNr
          }
          
          Adir=1
          if(Hstart>Hstop){ #if reverse alignment change start and stop location
            Adir=-1
            temp = Hstart
            Hstart =Hstop
            Hstop = temp
          }
          if(HI$V2==-1){ #ORF other direction 
            temp=HI$V4-Hstop
            Hstop = HI$V4-Hstart
            Hstart = temp
          }
          
          # find the start
          if(Hstart< HI$V5){
            FUTRstart = round(Hstart*100/HI$V5)
            ORFstart=0
            TUTRstart=0							
          }else if(Hstart<HI$V6){
            FUTRstart = 100
            ORFstart =  round((Hstart-HI$V5)*100/(HI$V6-HI$V5))
            TUTRstart=0							
          }else{
            FUTRstart = 100
            ORFstart =  100
            TUTRstart=round((Hstart-HI$V6)*100/(HI$V4-HI$V6))							
            
          }
          # find the stop
          if(Hstop< HI$V5){
            FUTRstop = round(Hstop*100/HI$V5)
            ORFstop=0
            TUTRstop=0							
          }else if(Hstop<HI$V6){
            FUTRstop = 100
            ORFstop =  round((Hstop-HI$V5)*100/(HI$V6-HI$V5))
            TUTRstop=0							
          }else{ 
            FUTRstop = 100
            ORFstop =  100
            if(HI$V4!=HI$V6){
              TUTRstop=round((Hstop-HI$V6)*100/(HI$V4-HI$V6))
            }else{
              TUTRstop=0
            }
            
          }
          if((HI$V2*QI$V2*Adir ==1)){
            if(FUTRstart!=FUTRstop){
              forward5UTR[i,FUTRstart:FUTRstop]=forward5UTR[i,FUTRstart:FUTRstop]+1
            }
            if(ORFstart!=ORFstop){
              forwardORF[i,ORFstart:ORFstop]=forwardORF[i,ORFstart:ORFstop]+1
            }
            if(TUTRstart!=TUTRstop){
              forward3UTR[i,TUTRstart:TUTRstop]=forward3UTR[i,TUTRstart:TUTRstop]+1
            }
          }else{
            if(FUTRstart!=FUTRstop){
              reverse5UTR[i,FUTRstart:FUTRstop]=reverse5UTR[i,FUTRstart:FUTRstop]+1
            }
            if(ORFstart!=ORFstop){
              reverseORF[i,ORFstart:ORFstop]=reverseORF[i,ORFstart:ORFstop]+1
            }
            if(TUTRstart!=TUTRstop){
              reverse3UTR[i,TUTRstart:TUTRstop]=reverse3UTR[i,TUTRstart:TUTRstop]+1
            }
          }
        }
      }
    }
  }
  
  
  for(i in seq(1,nlevels(infoQuery$classification))){
    numberOfGenes[i] = length(infoQuery[infoQuery['classification']==classes[i],][[1]])
  }
  
  
  forwardORF2 <- forwardORF/ numberOfGenes
  reverseORF2 <- reverseORF/ numberOfGenes
  forward5UTR2 <- forward5UTR/ numberOfGenes
  reverse5UTR2 <- reverse5UTR/ numberOfGenes
  forward3UTR2 <- forward3UTR/ numberOfGenes
  reverse3UTR2<-reverse3UTR/ numberOfGenes
  
  forward <- cbind(forward5UTR2,forwardORF2)
  forward <- cbind(forward,forward3UTR2)
  
  reverse <- cbind(reverse5UTR2,reverseORF2)
  reverse <- cbind(reverse,reverse3UTR2)
  
  pdf("Arabidopsis_Coverage")
  
  matplot(x,t(forward[1:4,]), type = "l",lty=1,ylim=c(0,1))
  matlines(x,t(reverse[1:4,]), type = "l",lty=2)
  legendInfoAccepted = paste(classes,"correct direction",sep=" ")[1:4]
  legendInfoAccepted = cbind(legendInfoAccepted,paste(classes,"wrong direction",sep=" ")[1:4])
  legend("topleft",legendInfoAccepted,lty=c(1,1,1,1,2,2,2,2),col=seq(1,4),cex=0.7)
  
  matplot(x,t(forward[5:7,]), type = "l",lty=1,ylim=c(0,1))
  matlines(x,t(reverse[5:7,]), type = "l",lty=2)
  legendInfoRejected = paste(classes,"correct direction",sep=" ")[5:7]
  legendInfoRejected = cbind(legendInfoRejected,paste(classes,"wrong direction",sep=" ")[5:7])
  legend("topleft",legendInfoRejected,lty=c(1,1,1,2,2,2),col=seq(1,3),cex=0.7)	
  
  dev.off()
  
  distributions = list()	
  distributions$forward5UTR <- forward5UTR
  distributions$forward3UTR <- forward3UTR
  distributions$forwardORF <- forward3UTR
  
  
  
}

