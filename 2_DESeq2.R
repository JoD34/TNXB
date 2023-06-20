# Import libraries ----
library(DESeq2)
library(EnhancedVolcano)
library(flow)
library(ggplot2)
library(ggrepel)
library(plotly)
library(pheatmap)
library(RColorBrewer)
library(ReportingTools)
library(scales)
library(stringr)
library(tidyverse)
library(vsn)

# functions ----
getNames <- function(list){
  name <- names(list) %>% 
    str_extract("LF_WT|1_2") %>% 
    strsplit("_") %>% 
    unlist() %>% 
    sapply(function(elem){
      elem <- ifelse(grepl("LF|WT", elem),
        ifelse(grepl("LF",elem),"Laminar Flow", "Wild Type"), 
        paste0("Sample ", elem))
    })
}

# Pdf data collection
pdf("graphs_tnxb.pdf", width = 11, height = 8.5, onefile = T)

# Make the design table ----
design.table <- data.frame(
  condition = str_extract(col.names, "LF|WT"),
  sample = str_extract(col.names, "(1|2)$")
) %>% 
  mutate(
    condition = factor(condition, levels=c("WT", "LF")),
    sample = factor(sample)) %>% 
  'rownames<-'(col.names)

# DESeq2 data generation ----
dds <- DESeqDataSetFromTximport(txi.kallisto, 
                                colData = design.table, 
                                design = ~ condition + sample) %>% 
  DESeq()

# Implement contrast between conditions + correct for NAs ----
graph_list <- list(
  con_LF_WT = results(dds,
                      contrast = c("condition", "LF", "WT")) %>% 
    na.omit(),
  sam_1_2 = results(dds,
                    contrast = c("sample", "1", "2")) %>%
    na.omit())

# Organize data for functional analysis ----
contrast_list <- list()
GSEA_data_list <- list()
GO_list <- list()

for(i in seq_along(graph_list)){
  
  # Extract list name
  listName <- names(graph_list[i]) %>% 
    str_extract("LF_WT|1_2") %>% 
    ifelse(. == "1_2", paste0("_", .), .)
  subList <- graph_list[[i]]
  
  # Filtering for p-adj < 0.01 & log2(FC) > 1.5
  filtered <- subList[subList$padj < 0.01 & 
                        abs(subList$log2FoldChange) > log2(1.5), ]
  # Ordered by p-adj
  ordered <- filtered[order(filtered$padj), ]
  
  # Get Go data
  GO_list <- append(GO_list, list(ordered))
  names(GO_list)[i] <- listName
  
  # Get data for GSEA analysis
  GSEA_data_list <- append(GSEA_data_list, 
                           list(ordered["log2FoldChange"]))
  names(GSEA_data_list)[i] <- listName
  
  # Seperate the upregulated from the downregulated
  decisionSet <- ordered$log2FoldChange > 0 
  myList <- list(
    pos = ordered[decisionSet > 0, ],
    neg = ordered[!decisionSet > 0, ]
  )

  # Get data for GO analysis
  contrast_list <- append(contrast_list, list(myList))
  names(contrast_list)[i] <- listName
}

# Volcano plots ----
for(i in seq_along(graph_list)){
  
  names <- getNames(graph_list[i])
  value <- graph_list[[i]]
  myLabels <- value[order(log(value$padj)), ][1:20, ]
  
  myPlot <- ggplot(as_tibble(value)) +
    aes(x = log2FoldChange, 
        y = -log10(padj)) +
    geom_point(aes(color = log2FoldChange > 0),
               show.legend = F, 
               shape = 1, 
               size = 2) +
    geom_vline(xintercept = 0,
               linetype = "dashed",
               alpha = 0.5) +
    geom_text_repel(data = as_tibble(myLabels),
                    label = rownames(myLabels),
                    size = 3,
                    box.padding = 0.5,
                    point.padding = 0.2) +
    labs(title = paste0("Volcano plot: ", names[[2]], 
                        " in relation to ", names[[1]]),
         x = expression(log[2]("Fold Change")),
         y = expression(-log[10]("adjusted pValue")),
         caption = paste("Produced on", Sys.time())) +
    theme_bw()
  
  y <- diff(layer_scales(myPlot)$y$range$range)/2
  
  myPlot <- myPlot +
    annotate("text",
             label = "Downregulated", 
             alpha = 0.5, 
             angle = 90, 
             x = -0.5, 
             y = y,
             size = 3) +
    annotate("text",
             label = "Upregulated", 
             alpha = 0.5, 
             angle = -90, 
             x = 0.5, 
             y = y,
             size = 3)
    print(myPlot)
}
# Apply log fold change shrinkage -> lower noises from log2(FC)----
ma_data <- list(
con_LF_WT.noNa.filt = lfcShrink(dds,
                                coef=c("condition_LF_vs_WT"), 
                                type = "apeglm") %>% 
  na.omit(),
sam_1_2.noNA.filt = lfcShrink(dds,
                                 coef=c("sample_2_vs_1"), 
                                 type = "apeglm") %>% 
  na.omit()
)
for (i in seq_along(ma_data)){
  names <- getNames(ma_data[i])
  value <- ma_data[[i]]
  
  myMA <- plotMA(value, colSig= "royalblue4", alpha = 0.01, 
         main=paste0("MAplot: ", names[[2]], " in relation to ", names[[1]]),
         ylab="log (FC)") 
  print(myMA)
}

# Effects of transformations on the variance ----

# vst -> variance stabilizing transformation
vsd <- vst(dds, blind=F)
meanSD <- meanSdPlot(assay(vsd))
print(meanSD)

# PCA analysis
pcaData <- plotPCA(vsd, 
                   intgroup=c("condition", "sample"),
                   returnData=T,
                   ntop=50000)
percentVar <- round(100 * attr(pcaData, "percentVar"))
namesPCA <- attr(pcaData,"names")

# PCA plot
pca <- ggplot(pcaData) +
  aes(PC1, PC2, color=condition, shape=sample) +
  geom_point(size=3) +
  labs(title = paste("Principal Component Analysis (PCA) - TNF\u03B1 data"),
       x = paste0(namesPCA[1], ": ", percentVar[1], " %"),
       y = paste0(namesPCA[2], ": ", percentVar[2], " %"),
       caption = paste("Produced on", Sys.time())) +
  coord_fixed() +
  theme_bw()

print(pca)


# Heatmap of the count matrix
select <- order(rowMeans(counts(dds,normalized=TRUE)),decreasing=TRUE)[1:200]
df <- as.data.frame(colData(dds))

heat <- pheatmap(assay(vsd)[select,], cluster_rows=T, show_rownames=FALSE,
         cluster_cols=F, annotation_col=df)
print(heat)

# Get sample-to-sample distance
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)

# Generate heatmap of the distance matrix 
rownames(sampleDistMatrix) <- paste(vsd$condition, vsd$sample, sep=":")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette(rev(brewer.pal(name="BuPu", n=9)))(255)
dist_mat <- pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors,
         main="Sample-to-sample Distance visualization")

print(dist_mat)

# Print datas in PDF
grDevices::dev.off()

