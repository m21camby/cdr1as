.libPaths(c("/data/rajewsky/home/skim/R/usr_lib_Seurat/", "/data/rajewsky/home/skim/R/usr_lib_v4/"))
library(limma)
library(DESeq2)
library(ggplot2)
library(dplyr)
library(ggrepel)

DGEm <- readRDS(file = "/data/rajewsky/projects/cdr1as_ko_snRNA/codes_github/cdr1as/DGEm/WTN_DGEm.rda")
Counts <- t(as.matrix(DGEm))

batch <- c("B", "M", "M","N", "B", "M", "M", "N", "M", "M", "M", "M", "M", "M", "M", "M")
mys <- colnames(DGEm)
batch2 = factor(c("Rep1", "Rep2", "Rep3","Rep4",
                  "Rep1", "Rep2", "Rep3","Rep4", 
                  "Rep5", "Rep6", "Rep7","Rep8",
                  "Rep5", "Rep6", "Rep7","Rep8"))

cleanNames <- function(x){ 
  return(substr(x,1,nchar(x)-2))
}

isKO <- function(x){
  return(substr(x,1,2))
}

isOE <- function(x){
  if (nchar(x)>2){
    return("OE")
  }
  else {
    return("NE") # normal expression
  }
}

getLineID <- function(x){
  return(substr(x,nchar(x),nchar(x)))
}

mys <- rownames(Counts)
batch2 <- sapply(mys, getLineID)
myg <- sapply(mys, cleanNames)
myg1 <- sapply(mys, isKO)
myg2 <- sapply(myg, isOE)

myp <- data.frame(id=mys, genotype=myg1, cond=myg2, batch=batch2, batchBuch=batch, group=myg)
mmCorrected <- model.matrix(~batchBuch+group, myp)

dds <- DESeqDataSetFromMatrix(t(Counts), colData = myp, mmCorrected)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

###################
# PCA plot
###################

pca_plot <- function(pc = p1, group = group, x_pc = "32", y_pc ="19"){
  ggplot(pc ,aes(x=PC1,y=PC2, color=group)) +  geom_point(size = 5, aes(shape = group)) + 
    xlab(paste0("PC1: ", x_pc, "% variance")) + 
    ylab(paste0("PC2: ", y_pc, "% variance")) +
    #  geom_text(aes(label=name),hjust=.5, vjust=2, size = 3) +
    theme(axis.text = element_text(face = "bold", size = 15, color = "black"),
          axis.title = element_text(face ="bold", size = 15, color = "black"),
          axis.line = element_line(colour="black"),
          axis.ticks = element_line(),
          panel.grid.major = element_line(color="white"),
          panel.background = element_rect(fill = "white", colour = "white"),
          panel.grid.minor = element_blank(),
          legend.text = element_text(size = 15),
          legend.title = element_blank(),
          legend.key = element_rect(colour = "transparent", fill = "white")) +
    scale_color_manual(values = c("#FF0066", "#880E4F", "#0099FF", "#0D47A1")) + 
    scale_shape_manual(values=c(19, 18, 17, 15))
  
}

# original PCA
vsd <- varianceStabilizingTransformation(dds)
vsd2 <- vsd
plotPCA(vsd, intgroup = "group")
p1 <- plotPCA(vsd, intgroup = "group", returnData = TRUE)

p1$group <- factor(p1$group, levels = c("WT", "WTm7o", "KO", "KOm7o"))
g1 <- pca_plot(pc = p1, group = group, x_pc = "51", y_pc ="21")
g1

ggsave(filename = "/data/rajewsky/projects/cdr1as_ko_snRNA/codes_github/cdr1as/Figures/WTN_PCA_orig.pdf",
       plot = g1,
       scale = 1, width = 6, height = 5, units = "in", device = cairo_pdf,
       dpi = 300)

# ------------------------- #
# site correction
# ------------------------- #
vsd <- vst(t(Counts[, keep]), blind = F) 

CountsNew <- vsd
mm2 <- model.matrix(~0+genotype+genotype:cond, myp)  
r <- removeBatchEffect(CountsNew, batch=batch, design=mm2)
assay(vsd2) <- r
plotPCA(vsd2, intgroup = "group")
p1 <- plotPCA(vsd2, intgroup = "group", returnData = TRUE)
p1$group <- factor(p1$group, levels = c("WT", "WTm7o", "KO", "KOm7o"))
g1 <- pca_plot(pc = p1, group = group, x_pc = "33", y_pc ="18")
g1
ggsave(filename = "/data/rajewsky/projects/cdr1as_ko_snRNA/codes_github/cdr1as/Figures/WTN_PCA_site_corr.pdf",
       plot = g1,
       scale = 1, width = 6, height = 5, units = "in", device = cairo_pdf,
       dpi = 300)

# ------------------------- #
# nested correction
# ------------------------- #
vsd <- vst(t(Counts[, keep]), blind = F) 
CountsNew <- vsd
vsd <- varianceStabilizingTransformation(dds)

r <- removeBatchEffect(CountsNew, batch=batch2, design=mm2)
vsd2 <- vsd
assay(vsd2) <- r
plotPCA(vsd2, intgroup = "group")

p1 <- plotPCA(vsd2, intgroup = "group", returnData = TRUE)
p1$group <- factor(p1$group, levels = c("WT", "WTm7o", "KO", "KOm7o"))
g1 <- pca_plot(pc = p1, group = group, x_pc = "51", y_pc ="20")
g1
ggsave(filename = "/data/rajewsky/projects/cdr1as_ko_snRNA/codes_github/cdr1as/Figures/WTN_PCA_nested.pdf",
       plot = g1,
       scale = 1, width = 6, height = 5, units = "in", device = cairo_pdf,
       dpi = 300)

# ------------------------- #
# batch & nested correction
# ------------------------- #
vsd <- vst(t(Counts[, keep]), blind = F) 
CountsNew <- vsd
vsd <- varianceStabilizingTransformation(dds3)
r <- removeBatchEffect(CountsNew, batch=batch2, batch2=batch, design=mm2)
vsd2 <- vsd
assay(vsd2) <- r
plotPCA(vsd2, intgroup = "group")
p1 <- plotPCA(vsd2, intgroup = "group", returnData = TRUE)
p1$group <- factor(p1$group, levels = c("WT", "WTm7o", "KO", "KOm7o"))
g1 <- pca_plot(pc = p1, group = group, x_pc = "42", y_pc ="22")
g1
ggsave(filename = "/data/rajewsky/projects/cdr1as_ko_snRNA/codes_github/cdr1as/Figures/WTN_PCA_batch_corr_nested.pdf",
       plot = g1,
       scale = 1, width = 6, height = 5, units = "in", device = cairo_pdf,
       dpi = 300)

# save batch & nested corrected matrix for WGCNA and RF
r.matrix <- r %>% as.matrix()
r.df <- r %>% as.data.frame()
saveRDS(r.df, file = "/data/rajewsky/projects/cdr1as_ko_snRNA/codes_github/cdr1as/DGEm/WTN_batch_corr_nested_DGEm.rda")


# ------------------------- #
# batch & nested correction DESeq2
# ------------------------- #
mm <- model.matrix(~0 + batchBuch + genotype + genotype:batch + genotype:cond, myp)
wc <- which(colnames(mm)=="genotypeWT:batch4" | colnames(mm)=="genotypeWT:batch3")
mmCorrected <- mm[,-wc]
dds <- DESeqDataSetFromMatrix(DGEm, colData = myp, mmCorrected)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds <- DESeq(dds)
resultsNames(dds)

DGEm_normalized <- counts(dds, normalized=T) %>% as.data.frame
DGEm_normalized_WT <- DGEm_normalized[,c(1:4)] %>% rowMeans

# WTm7oe vs WT
res2 <- results(dds, list( c("genotypeWT.condOE","genotypeWT") ))
res2.df <- as.data.frame(res2)
res2.df$gene <- rownames(res2.df)
res2.df$WT <- DGEm_normalized_WT
res2.df$gene[res2.df$gene == "1700020I14Rik"] <- "Cyrano"
res2.df$gene[res2.df$gene == "C230004F18Rik"] <- "Cdr1os"
saveRDS(res2.df, file = "/data/rajewsky/projects/cdr1as_ko_snRNA/codes_github/cdr1as/Results/DESeq2_WTN_batch_nested_corrected_res2.rda")

# KOm7oe vs KO
res3 <- results(dds, name=c("genotypeKO.condOE"))
res3.df <- as.data.frame(res3)
res3.df$gene <- rownames(res3.df)
res3.df$WT <- DGEm_normalized_WT
res3.df$gene[res3.df$gene == "1700020I14Rik"] <- "Cyrano"
res3.df$gene[res3.df$gene == "C230004F18Rik"] <- "Cdr1os"
saveRDS(res3.df, file = "/data/rajewsky/projects/cdr1as_ko_snRNA/codes_github/cdr1as/Results/DESeq2_WTN_batch_nested_corrected_res3.rda")

# ------------------------- #
# nested correction DESeq2
# ------------------------- #
myp <- data.frame(id=mys, genotype=factor(myg1, levels=c("WT", "KO")), cond=myg2, batch=batch2, group=myg)
mm <- model.matrix(~genotype+genotype:batch+genotype:cond, myp)

dds <- DESeqDataSetFromMatrix(DGEm, colData = myp, mm)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds <- DESeq(dds)
resultsNames(dds)

# WTm7oe vs WT
res2 <- results(dds, name = c("genotypeWT.condOE"))
res2.df <- as.data.frame(res2)
res2.df$gene <- rownames(res2.df)
res2.df$WT <- DGEm_normalized_WT
res2.df$gene[res2.df$gene == "1700020I14Rik"] <- "Cyrano"
res2.df$gene[res2.df$gene == "C230004F18Rik"] <- "Cdr1os"
saveRDS(res2.df, file = "/data/rajewsky/projects/cdr1as_ko_snRNA/codes_github/cdr1as/Results/DESeq2_WTN_nested_corrected_res2.rda")

# KOm7oe vs KO
res3 <- results(dds, name=c("genotypeKO.condOE"))
res3.df <- as.data.frame(res3)
res3.df$gene <- rownames(res3.df)
res3.df$WT <- DGEm_normalized_WT
res3.df$gene[res3.df$gene == "1700020I14Rik"] <- "Cyrano"
res3.df$gene[res3.df$gene == "C230004F18Rik"] <- "Cdr1os"
saveRDS(res3.df, file = "/data/rajewsky/projects/cdr1as_ko_snRNA/codes_github/cdr1as/Results/DESeq2_WTN_nested_corrected_res3.rda")


# ------------------------- #
# batch corrected DESeq2
# ------------------------- #

dds <- DESeqDataSetFromMatrix(DGEm, colData = myp, design = ~ 0 + group + batchBuch)
keep <- rowSums(counts(dds)) >= 10
table(keep)
dds <- dds[keep,]
dds$group <- relevel(dds$group, ref = "WT")
dds <- DESeq(dds)
resultsNames(dds)

# DGEm_normalized <- counts(dds, normalized=T) %>% as.data.frame
# DGEm_normalized_WT <- DGEm_normalized[,c(1:3)] %>% rowMeans

# KO vs WT
res1 <- results(dds, contrast=c("group","KO", "WT"), test="Wald")
res1.df <- as.data.frame(res1)
res1.df$gene <- rownames(res1.df)
res1.df$WT <- DGEm_normalized_WT
res1.df$gene[res1.df$gene == "1700020I14Rik"] <- "Cyrano"
res1.df$gene[res1.df$gene == "C230004F18Rik"] <- "Cdr1os"
saveRDS(res1.df, file = "/data/rajewsky/projects/cdr1as_ko_snRNA/codes_github/cdr1as/Results/DESeq2_WTN_batch_corrected_res1.rda")

# WTm7oe vs WT
res2 <- results(dds, contrast=c("group","WTm7o", "WT"), test="Wald")
res2.df <- as.data.frame(res2)
res2.df$gene <- rownames(res2.df)
res2.df$WT <- DGEm_normalized_WT
res2.df$gene[res2.df$gene == "1700020I14Rik"] <- "Cyrano"
res2.df$gene[res2.df$gene == "C230004F18Rik"] <- "Cdr1os"
saveRDS(res2.df, file = "/data/rajewsky/projects/cdr1as_ko_snRNA/codes_github/cdr1as/Results/DESeq2_WTN_batch_corrected_res2.rda")

# KOm7oe vs KO
res3 <- results(dds, contrast=c("group","KOm7o", "KO"), test="Wald")
res3.df <- as.data.frame(res3)
res3.df$gene <- rownames(res3.df)
res3.df$WT <- DGEm_normalized_WT
res3.df$gene[res3.df$gene == "1700020I14Rik"] <- "Cyrano"
res3.df$gene[res3.df$gene == "C230004F18Rik"] <- "Cdr1os"
saveRDS(res3.df, file = "/data/rajewsky/projects/cdr1as_ko_snRNA/codes_github/cdr1as/Results/DESeq2_WTN_batch_corrected_res3.rda")


# KOm7oe vs WT
res4 <- results(dds, contrast=c("group","KOm7o", "WT"), test="Wald")
res4.df <- as.data.frame(res4)
res4.df$gene <- rownames(res4.df)
res4.df$WT <- DGEm_normalized_WT
res4.df$gene[res4.df$gene == "1700020I14Rik"] <- "Cyrano"
res4.df$gene[res4.df$gene == "C230004F18Rik"] <- "Cdr1os"
saveRDS(res4.df, file = "/data/rajewsky/projects/cdr1as_ko_snRNA/codes_github/cdr1as/Results/DESeq2_WTN_batch_corrected_res4.rda")

target_MA_Plot <- function(DF, labels){
  ggplot(DF, aes(x = WT, y = log2FoldChange)) + 
    geom_point(data = DF[which(DF$padj > 0.05), ], aes(x = WT, y = log2FoldChange), alpha = 0.3, color = "darkgrey", shape = 1, na.rm = TRUE) + 
    geom_point(data = DF[which(DF$padj < 0.05 & DF$log2FoldChange > 0.5), ], aes(x = WT, y = log2FoldChange), color = "darkred", shape = 1, size = 2.5) + 
    geom_point(data = DF[which(DF$padj < 0.05 & DF$log2FoldChange < -0.5), ], aes(x = WT, y = log2FoldChange), color = "navy", shape = 1, size = 2.5) + 
    theme(axis.title = element_text(face = "bold", size = 12), 
          axis.text = element_text(face = "bold", size = 12, color = "black"), 
          panel.background = element_blank(), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          axis.line = element_line(colour = "black"), 
          legend.position = c(0.9, 0.8), 
          legend.title = element_blank(),
          legend.key = element_blank()) + 
    xlab("Mean of Normalized Counts (WT)") + 
    ylab(paste0(labels)) +
    scale_x_continuous(trans='log10') +
    geom_hline(yintercept = 0, size = 0.3) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "darkorange", size = 0.2) + 
    geom_hline(yintercept = -1, linetype = "dashed", color = "darkorange", size = 0.2) 
}

target_genes <- c("Cyrano", "Cdr1os")


g1 <- target_MA_Plot(res1.df, labels = "log2FC (WTN vs FullKO)")  +
  geom_point(data = res1.df[which(res1.df$gene %in% target_genes & res1.df$padj < 0.05), ], aes(x = WT, y = log2FoldChange), color = "darkorange", size = 3) + 
  geom_point(data = res1.df[which(res1.df$gene %in% target_genes & res1.df$padj > 0.05), ], aes(x = WT, y = log2FoldChange), color = "darkgreen", size = 3) + 
  coord_cartesian(ylim = c(-5,5)) + 
  geom_text_repel(data = subset(res1.df, gene %in% c(target_genes[1])),aes(label = gene),
                  size = 4, box.padding = unit(0.35, "lines"), point.padding = unit(0.3, "lines"),
                  nudge_y = 0.2) + 
  geom_text_repel(data = subset(res1.df, gene %in% c(target_genes[2])),aes(label = gene),
                  size = 4, box.padding = unit(0.35, "lines"), point.padding = unit(0.3, "lines"),
                  nudge_y =  0.4)
g1

ggsave(filename = "/data/rajewsky/projects/cdr1as_ko_snRNA/codes_github/cdr1as/Figures/DESeq2_WTN_res1.pdf",
       plot = g1,
       scale = 1, width = 5, height = 4, units = "in", device = cairo_pdf,
       dpi = 300)

g2 <- target_MA_Plot(res2.df, labels = "log2FC (WTN vs WTN-miR-7oe)")  +
  geom_point(data = res2.df[which(res2.df$gene %in% target_genes & res2.df$padj < 0.05), ], aes(x = WT, y = log2FoldChange), color = "darkorange", size = 3) + 
  geom_point(data = res2.df[which(res2.df$gene %in% target_genes & res2.df$padj > 0.05), ], aes(x = WT, y = log2FoldChange), color = "darkgreen", size = 3) + 
  coord_cartesian(ylim = c(-5,5)) + 
  geom_text_repel(data = subset(res2.df, gene %in% c(target_genes[1])),aes(label = gene),
                  size = 4, box.padding = unit(0.35, "lines"), point.padding = unit(0.3, "lines"),
                  nudge_y = - 0.4) + 
  geom_text_repel(data = subset(res2.df, gene %in% c(target_genes[2])),aes(label = gene),
                  size = 4, box.padding = unit(0.35, "lines"), point.padding = unit(0.3, "lines"),
                  nudge_y =  0.6)
g2

ggsave(filename = "/data/rajewsky/projects/cdr1as_ko_snRNA/codes_github/cdr1as/Figures/DESeq2_WTN_res2.pdf",
       plot = g2,
       scale = 1, width = 5, height = 4, units = "in", device = cairo_pdf,
       dpi = 300)

g3 <- target_MA_Plot(res3.df, labels = "log2FC (FullKO vs FullKO-miR-7oe)")  +
  geom_point(data = res3.df[which(res3.df$gene %in% target_genes & res3.df$padj < 0.05), ], aes(x = WT, y = log2FoldChange), color = "darkorange", size = 3) + 
  geom_point(data = res3.df[which(res3.df$gene %in% target_genes & res3.df$padj > 0.05), ], aes(x = WT, y = log2FoldChange), color = "darkgreen", size = 3) + 
  coord_cartesian(ylim = c(-5,5)) + 
  geom_text_repel(data = subset(res3.df, gene %in% c(target_genes[1])),aes(label = gene),
                  size = 4, box.padding = unit(0.35, "lines"), point.padding = unit(0.3, "lines"),
                  nudge_y = - 0.4) + 
  geom_text_repel(data = subset(res3.df, gene %in% c(target_genes[2])),aes(label = gene),
                  size = 4, box.padding = unit(0.35, "lines"), point.padding = unit(0.3, "lines"),
                  nudge_y =  0.4)
g3

ggsave(filename = "/data/rajewsky/projects/cdr1as_ko_snRNA/codes_github/cdr1as/Figures/DESeq2_WTN_res3.pdf",
       plot = g3,
       scale = 1, width = 5, height = 4, units = "in", device = cairo_pdf,
       dpi = 300)

g4 <- target_MA_Plot(res4.df, labels = "log2FC (WTN vs FullKO-miR-7oe)")  +
  geom_point(data = res4.df[which(res4.df$gene %in% target_genes & res4.df$padj < 0.05), ], aes(x = WT, y = log2FoldChange), color = "darkorange", size = 3) + 
  geom_point(data = res4.df[which(res4.df$gene %in% target_genes & res4.df$padj > 0.05), ], aes(x = WT, y = log2FoldChange), color = "darkgreen", size = 3) + 
  coord_cartesian(ylim = c(-5,5)) + 
  geom_text_repel(data = subset(res4.df, gene %in% c(target_genes[1])),aes(label = gene),
                  size = 4, box.padding = unit(0.35, "lines"), point.padding = unit(0.3, "lines"),
                  nudge_x = 0.2, nudge_y = -0.5) + 
  geom_text_repel(data = subset(res4.df, gene %in% c(target_genes[2])),aes(label = gene),
                  size = 4, box.padding = unit(0.35, "lines"), point.padding = unit(0.3, "lines"),
                  nudge_x = 0.2,nudge_y =  0.5)
g4

ggsave(filename = "/data/rajewsky/projects/cdr1as_ko_snRNA/codes_github/cdr1as/Figures/DESeq2_WTN_res4.pdf",
       plot = g4,
       scale = 1, width = 5, height = 4, units = "in", device = cairo_pdf,
       dpi = 300)


# ------------------------- #
# nested only MA plots
# ------------------------- #

res2.df <- readRDS(file = "/data/rajewsky/projects/cdr1as_ko_snRNA/codes_github/cdr1as/Results/DESeq2_WTN_nested_corrected_res2.rda")

g2 <- target_MA_Plot(res2.df, labels = "log2FC (WTN vs WTN-miR-7oe)")  +
  geom_point(data = res2.df[which(res2.df$gene %in% target_genes & res2.df$padj < 0.05), ], aes(x = WT, y = log2FoldChange), color = "darkorange", size = 3) + 
  geom_point(data = res2.df[which(res2.df$gene %in% target_genes & res2.df$padj > 0.05), ], aes(x = WT, y = log2FoldChange), color = "darkgreen", size = 3) + 
  coord_cartesian(ylim = c(-5,5)) + 
  geom_text_repel(data = subset(res2.df, gene %in% c(target_genes[1])),aes(label = gene),
                  size = 4, box.padding = unit(0.35, "lines"), point.padding = unit(0.3, "lines"),
                  nudge_y = - 0.4) + 
  geom_text_repel(data = subset(res2.df, gene %in% c(target_genes[2])),aes(label = gene),
                  size = 4, box.padding = unit(0.35, "lines"), point.padding = unit(0.3, "lines"),
                  nudge_y =  - 0.6)
g2

ggsave(filename = "/data/rajewsky/projects/cdr1as_ko_snRNA/codes_github/cdr1as/Figures/DESeq2_WTN_res2_nested.pdf",
       plot = g2,
       scale = 1, width = 5, height = 4, units = "in", device = cairo_pdf,
       dpi = 300)

res3.df <- readRDS(file = "/data/rajewsky/projects/cdr1as_ko_snRNA/codes_github/cdr1as/Results/DESeq2_WTN_nested_corrected_res3.rda")

g3 <- target_MA_Plot(res3.df, labels = "log2FC (FullKO vs FullKO-miR-7oe)")  +
  geom_point(data = res3.df[which(res3.df$gene %in% target_genes & res3.df$padj < 0.05), ], aes(x = WT, y = log2FoldChange), color = "darkorange", size = 3) + 
  geom_point(data = res3.df[which(res3.df$gene %in% target_genes & res3.df$padj > 0.05), ], aes(x = WT, y = log2FoldChange), color = "darkgreen", size = 3) + 
  coord_cartesian(ylim = c(-5,5)) + 
  geom_text_repel(data = subset(res3.df, gene %in% c(target_genes[1])),aes(label = gene),
                  size = 4, box.padding = unit(0.35, "lines"), point.padding = unit(0.3, "lines"),
                  nudge_y = - 0.4) + 
  geom_text_repel(data = subset(res3.df, gene %in% c(target_genes[2])),aes(label = gene),
                  size = 4, box.padding = unit(0.35, "lines"), point.padding = unit(0.3, "lines"),
                  nudge_y =  0.4)
g3

ggsave(filename = "/data/rajewsky/projects/cdr1as_ko_snRNA/codes_github/cdr1as/Figures/DESeq2_WTN_res3_nested.pdf",
       plot = g3,
       scale = 1, width = 5, height = 4, units = "in", device = cairo_pdf,
       dpi = 300)
