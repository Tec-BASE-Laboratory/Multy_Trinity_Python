## Tec BASE Differential Expression Analysis
### Step 1: Load package libraries required for data import and analysis--------------------------------------
library(DESeq2)
library(tximport)
library(ggplot2)
#-------------------------------------------------------------------------------------------------------------
### Step 2: Import Salmon derived quantification files and sample information-----------------------------------
#Set the dir where we will be working, ir is important you have the files you want to analyse in this directory 
dir <- "D:/Emili/Desktop/Prueba_Deseq/tximp/"

#you may also need to set your working directory on the 'Session' tab
#Remember that in the 'samples.txt'file you must list all the files names exactly as is in the working directory
#You also have to provide a condition describing the treatment on your sample
samples <-read.table(file.path(dir,"samples.txt"), header=TRUE) #txt with file names as ref

#Gather the files given on you txt file, if the names differ this will carry on more errors later on
files <- file.path(dir, samples$Sample)
names(files) <- paste(samples$Sample)

#---------Step 3: create experimental design objects for DEA--------------------------------------------------------
#txi will allow us to import the files in a way DESEQ can work, we use the function lengthScaledTPM to normalize the abundance of the transcripts given its lenght
txi <- tximport(files, type="salmon", txOut= TRUE, countsFromAbundance = "lengthScaledTPM",  dropInfReps= TRUE)
#we run the DESeqDataSetFromTximport preliminary to the real DEA
dds <- DESeqDataSetFromTximport(txi, colData= samples, design = ~Condition)

#------------------Step 4: We run properly the analysis----------------------------------------------

#We run the main Differential Expression Analysis
ddsDE <- DESeq(dds) #Main differential express analysis
#Extracts results table with log2 fold changes, p values and adjusted p values.
res <- results(ddsDE) 
#Adjusting the results with an specific P value
res05 <- results(ddsDE, alpha=0.05)
#Log fold change shrinkage for visualization and ranking
resLFC <- lfcShrink(ddsDE, coef="Condition_FSN2_vs_Control", type="apeglm")
#We can order our results table by the smallest p value
resOrdered <- res[order(res$pvalue),]

#---------------------Standar Plots and Visualizations--------------------------------------------

#Dispersion plot for diagnostic
plotDispEsts(ddsDE)
#shows the log2 fold changes attributable to a given variable over the mean of normalized counts for all the samples in the DESeqDataSet
plotMA(ddsDE, ylim= c(-5,5))
#plot with the shrunken value log2 fold changes, which remove the noise associated with 
#log2 fold changes from low count genes without requiring arbitrary filtering thresholds.
plotMA(resLFC, ylim=c(-2,2))
#plot the counts of an specific gene in this case the one with the lowest p value
plotCounts(ddsDE, gene=which.min(res$padj), intgroup="Condition")
#the one with the max p value available
plotCounts(ddsDE, gene=which.max(res$padj), intgroup="Condition") 
#once annotated we can perform this analysis for our genes of interest
#plots the standard deviation of the transformed data, across samples, against the mean, using the shifted logarithm transformation
ntd <- normTransform(ddsDE)
library("vsn")
meanSdPlot(assay(ntd))

#-------------------------Heat map of the data (preliminary)---------------------------------------

#standard heat map obtain with the data given to the program
library("pheatmap")
select <- order(rowMeans(counts(ddsDE,normalized=TRUE)),
                decreasing=TRUE)[1:20]
df <- as.data.frame(colData(ddsDE)[,c("Condition","Sample")])
pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)


#standard heat map obtain with the data given to the program
library("pheatmap")
select <- order(rowMeans(counts(ddsDE,normalized=TRUE)),
                decreasing=TRUE)[1:20]
df <- as.data.frame(colData(ddsDE)[,c("Sample","Clusters")])
pheatmap(assay(ntd)[select,], cluster_rows=TRUE, show_rownames=FALSE,
         cluster_cols=TRUE, annotation_col=df)


#----------------------------------Volcano plot-------------------------------------

#plotting the heat map first with the -log10 of the pvalue and then with a theme
p <- ggplot(data=DeseqWithNames, aes(x=log2FoldChange, y=-log10(pvalue))) + geom_point()
p <- ggplot(data=DeseqWithNames, aes(x=log2FoldChange, y=-log10(pvalue))) + geom_point() + theme_minimal()

# Since we want to indicate the ones that are diff express we create a new column
DeseqWithNames$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
DeseqWithNames$diffexpressed[DeseqWithNames$log2FoldChange > 0.6 & DeseqWithNames$pvalue < 0.05] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
DeseqWithNames$diffexpressed[DeseqWithNames$log2FoldChange < -0.6 & DeseqWithNames$pvalue < 0.05] <- "DOWN"

# Re-plot but this time color the points with "diffexpressed"
p <- ggplot(data=DeseqWithNames, aes(x=log2FoldChange, y=-log10(pvalue), col=diffexpressed)) + geom_point() + theme_minimal()

# Add lines to the plot
p2 <- p + geom_vline(xintercept=c(-0.6, 0.6), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")


# 1. by default, it is assigned to the categories in an alphabetical order):
p3 <- p2 + scale_color_manual(values=c("blue", "black", "red"))

# 2. to automate a bit: ceate a named vector: the values are the colors to be used, the names are the categories they will be assigned to:
mycolors <- c("blue", "red", "black")
names(mycolors) <- c("DOWN", "UP", "NO")
p3 <- p2 + scale_colour_manual(values = mycolors)

#We create a new column containing the -log values of the pvalues
DeseqWithNames$maj12 <- -log(DeseqWithNames$pvalue)

# This column will contain labels that will be used on the graph
DeseqWithNames$delabel <- "NO"
DeseqWithNames$delabel[DeseqWithNames$maj12 > 12] <- "HIT"

#We concatenate the names here if the genes is diff expressed
DeseqWithNames$delabel[DeseqWithNames$delabel == "HIT"] <- DeseqWithNames$conc.x[DeseqWithNames$delabel == "HIT"]


#Ploting without the labels
ggplot(data=DeseqWithNames, aes(x=log2FoldChange, y=-log10(pvalue), col=diffexpressed)) + 
  geom_point() + 
  theme_minimal() 

#Ploting with the labels
ggplot(data=DeseqWithNames, aes(x=log2FoldChange, y=-log10(pvalue), col=diffexpressed, label=delabel)) + 
  geom_point() + 
  theme_minimal() +
  geom_text()


