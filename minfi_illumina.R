#!/usr/bin/env Rscript

source("http://bioconductor.org/biocLite.R")
biocLite("minfi")
biocLite("IlluminaHumanMethylation450kanno.ilmn12.hg19")
biocLite("IlluminaHumanMethylation450kmanifest")
biocLite("IlluminaHumanMethylationEPICmanifest")
biocLite("IlluminaHumanMethylationEPICanno.ilm10b2.hg19")
biocLite("limma")
library("minfi")
library("IlluminaHumanMethylationEPICmanifest")
library("IlluminaHumanMethylationEPICanno.ilm10b2.hg19")
library("IlluminaHumanMethylation450kmanifest")
library("IlluminaHumanMethylation450kanno.ilmn12.hg19")
library("Biobase")
library("RColorBrewer")
library("limma")


###Reading input IDAT SampleSheet###
baseDir <- getwd()
outDir <- getwd()
targets <-  read.metharray.sheet(baseDir)

###Read IDAT intensity data as RGChannelSet (Red/Green Channel set)
RGSet <- read.metharray.exp(targets=targets)
phenoData <- pData(RGSet)
manifest <- getManifest(RGSet)

#####QC and plots####
qcReportpdf = qcReport(RGSet,sampNames = targets$Sample_Name,sampGroups=targets$Sample_Group,pdf="qcReport.pdf")
qcReportpdf

Mset.illumina <- preprocessIllumina(RGSet, bg.correct=TRUE, normalize="controls") ##After Norm
qc <- getQC(Mset.illumina)
png(file=paste(outDir,"/QC.png", sep=""),width=2048,height=2048,pointsize=50)
plotQC(qc)
dev.off()

###Removing poor quality samples###
detP <- detectionP(RGSet)
colnames(detP) <- RGSet$Sample_Group
# examine mean detection p-values across all samples to identify any failed samples
pal <- brewer.pal(8,"Dark2")
par(mfrow=c(1,2))
png(file=paste(outDir,"/pvalue_sample_filter.png", sep=""),width=2048,height=2048,pointsize=50)
barplot(colMeans(detP), col=pal[factor(targets$Sample_Group)], las=2,
        cex.names=0.5, cex.axis=0.5, names.arg=targets$Sample_Name, ylab="Mean Detection p-values")
abline(h=0.05,col="red")
legend("topleft", legend=levels(factor(targets$Sample_Group)), fill=pal,
       bg="white")
dev.off()

# remove poor quality samples
keep <- colMeans(detP) < 0.05
RGSet <- RGSet[,keep]
RGSet
# remove poor quality samples from targets data
targets <- targets[keep,]
detP <- detP[,keep]
dim(detP)
mSetSq <- preprocessQuantile(RGSet)

# ensure probes are in the same order in the mSetSq and detP objects
detP <- detP[match(featureNames(mSetSq),rownames(detP)),]
# remove any probes that have failed in one or more samples
keep <- rowSums(detP < 0.01) == ncol(mSetSq)
table(keep)
gset.funnorm <- mSetSq[keep,]
dim(gset.funnorm)

###Normalisation###
gset.funnorm <- addSnpInfo(gset.funnorm) ##add the genomic ranges info to gset
gset.funnorm <- dropLociWithSnps(gset.funnorm,snps=c("SBE", "CpG"), maf=0) ##drop the loci which has snps

###Get annotation###
annot = getAnnotation(gset.funnorm)

###Remove sex probes###
sex_probes = annot$Name[annot$chr %in% c("chrX", "chrY")]
gset.funnorm = gset.funnorm[!(rownames(gset.funnorm) %in% sex_probes),]

###Get beta values to obtain the data matrix###
gset.funnorm.beta <- getBeta(gset.funnorm)
colnames(gset.funnorm.beta) <- gset.funnorm$Sample_Group
write.csv(gset.funnorm.beta,file="beta.csv")

#MDSplots
pal <- brewer.pal(8,"Dark2")
png(file=paste(outDir,"/MDSplot.png", sep=""),width=2048,height=2048,pointsize=40)
plotMDS(gset.funnorm.beta, top=1000, gene.selection="common",col=pal[factor(targets$Sample_Group)],pch=16)
legend("topright", legend=levels(factor(targets$Sample_Group)), text.col=pal,bg="white", cex=0.9)
dev.off()
write.csv(gset.funnorm.beta,file="beta.csv") ##write to a csv the beta value datamatrix

###Differential Methylation###
condition <- pData(gset.funnorm)$Sample_Group #provide the sample groups for diff.methylation (if all groups)

dmp <- dmpFinder(gset.funnorm.beta, pheno=condition, type="categorical") #run dmp on the data
dmp <- cbind(dmp, ID=rownames(dmp))
dmp_annot_combined <- cbind(annot[row.names(dmp),],dmp) #combine the dmp pvalue data with the annotation.
write.csv(dmp_annot_combined,file="differential_methylation_withAnnotation.csv")

###Plot the heatmap###
library(gplots)
cell_colors = colorRampPalette( c("#010F57", "#010F57", "#FAFAFA", "#B21212", "#B21212") )(300)
f <- factor(targets$Sample_Group)
png(file=paste(outDir,"/Heatmap(Top 100).png", sep=""),width=2050,height=2048,pointsize=50)
heatmap.2(gset.funnorm.beta[row.names(dmp[1:100,]),],trace = 'none',key.title="Methylation",labCol = targets$Sample_Name, cexCol = 0.5, scale = 'row',col = cell_colors,cexRow = 0.52,key.xlab = "BetaValue",main = "Heatmap (Top 100)",ColSideColors = pal[factor(targets$Sample_Group)])
legend("topright", legend=levels(f), col=pal[factor(levels(f))], pch=15,cex = 0.45)
dev.off()


###Plot Boxplot and Violin plots of global methylation###
library(reshape2)
library(ggplot2)
df <- melt(gset.funnorm.beta,value.name = "Methylation")
colnames(df) <- c("Probes","Samples","Methylation")
png(file=paste(outDir,"/GlobalMethylation_Boxplot.png", sep=""),width=800,height=800,pointsize=50)
ggplot(aes(y=Methylation,x=Samples,fill=Samples),data = df)+geom_boxplot()
dev.off()

png(file=paste(outDir,"/GlobalMethylation_Violinplot.png", sep=""),width=800,height=800,pointsize=50)
ggplot(aes(y=Methylation,x=Samples,fill=Samples),data = df)+geom_violin()
dev.off()