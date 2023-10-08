# if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
# BiocManager::install("DESeq2")
# BiocManager::install("EnhancedVolcano")
# BiocManager::install("pheatmap")
# BiocManager::install("RColorBrewer")

# Load requier libraries
library("DESeq2")
library("EnhancedVolcano")
library("pheatmap")
library("RColorBrewer")

#Set working directory
## Change this PATH !!!!!!!
path="/your_path/htseq"
setwd(path)

#Create variable with file names with count 
files<-c("SRR2936836.htseq-counts.tsv","SRR2936837.htseq-counts.tsv","SRR2936838.htseq-counts.tsv","SRR2936839.htseq-counts.tsv",
"SRR2936840.htseq-counts.tsv","SRR2936841.htseq-counts.tsv","SRR2936842.htseq-counts.tsv","SRR2936843.htseq-counts.tsv",
"SRR2936844.htseq-counts.tsv","SRR2936845.htseq-counts.tsv","SRR2936846.htseq-counts.tsv","SRR2936847.htseq-counts.tsv",
"SRR2936848.htseq-counts.tsv","SRR2936849.htseq-counts.tsv","SRR2936850.htseq-counts.tsv","SRR2936851.htseq-counts.tsv",
"SRR2936852.htseq-counts.tsv","SRR2936853.htseq-counts.tsv","SRR2936854.htseq-counts.tsv","SRR2936855.htseq-counts.tsv",
"SRR2936856.htseq-counts.tsv","SRR2936857.htseq-counts.tsv","SRR2936858.htseq-counts.tsv","SRR2936859.htseq-counts.tsv",
"SRR2936860.htseq-counts.tsv","SRR2936861.htseq-counts.tsv","SRR2936862.htseq-counts.tsv","SRR2936863.htseq-counts.tsv",
"SRR2936864.htseq-counts.tsv","SRR2936865.htseq-counts.tsv")

#Define study group
group=c(rep("ctrl", each=16), rep("knoc", each=14))
table <- data.frame(sampleName = files, fileName = files, condition = group)
table

#Create Deseq2 data set
ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = table, directory = path, design= ~condition)

#Remove low counts
keep <- rowSums(counts(ddsHTSeq)) >= 5
ddsHTSeq <- ddsHTSeq[keep,]

#Run DESeq2
dds <- DESeq(ddsHTSeq)
res<-results( dds, contrast=c("condition","knoc","ctrl"))
res<-res[order(res$padj),]

#Print top results
head(res)
de<-as.data.frame(res)

#Export results to CSV
write.csv(de,file='de.csv')
write.csv(counts(dds, normalized=T),file='normalized_counts.csv')
write.csv(counts(dds, normalized=F),file='raw_counts.csv')

#Plot MA
plotMA(dds,ylim=c(-10,10),main='DE pAdjValue < 0.05')
dev.copy(png,'DE MAplot.png')
dev.off()

# Draw sample to sample heatmaps
vsd <- vst(dds, blind=FALSE)
rld <- rlog(dds, blind=FALSE)
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$condition, vsd$type)
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
df <- as.data.frame(colData(dds)[,c("condition")])
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         annotate=df)
dev.copy(png,'sample-to-sample.png')
dev.off()

#Plot PCA
plotPCA(vsd, intgroup=c("condition"))
dev.copy(png,'pca.png')
dev.off()

# Draw top genes heatmap
topVarGenes <- head( order( rowVars( assay(rld) ), decreasing=TRUE ), 15 )
pheatmap( assay(rld)[ topVarGenes, ], scale="row", 
          trace="none", dendrogram="column", 
          col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255),
          annotate_col=df )

dev.copy(png,'15-gene-distance.png')
dev.off()

count_de<-dim(de[de[5]<0.05,])[1]
topVarGenes <- head( order( rowVars( assay(rld) ), decreasing=TRUE ),count_de )
pheatmap(assay(rld)[ topVarGenes, ], scale="row", 
         trace="none", dendrogram="column", 
         col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255),
         annotate_col=df )

dev.copy(png,'all-gene-distance.png')
dev.off()

# Draw volcano plot
EnhancedVolcano(res,
                lab = rownames(res),
                x = "log2FoldChange",
                y = "pvalue")

dev.copy(png,'volcano_plot.png')
dev.off()