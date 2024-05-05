
# This file is for analyzing 5'SSKO data to validate FullKO
.libPaths(c("/data/rajewsky/home/skim/R/usr_lib_Seurat/", "/data/rajewsky/home/skim/R/usr_lib_v4/"))
#.libPaths()
unloadNamespace("mgcv")
unloadNamespace("Matrix")
library(ggplot2)
library(plotly)
library(gridExtra)
library(ggrepel)
library(dplyr)
library(cowplot)
library(DESeq2)
library(limma)

SS_DGEm <- readRDS(file = "/data/rajewsky/projects/cdr1as_ko_snRNA/codes_github/cdr1as/DGEm/WTJ_5SSKO_DGEm.rda")




group <- c("WTJ", "WTJ", "WTJ", 
           "WTJm7oe", "WTJm7oe", "WTJm7oe", 
           "SSKO","SSKO","SSKO",
           "SSKOm7oe","SSKOm7oe","SSKOm7oe")

genotype <- c("WTJ", "WTJ", "WTJ", 
              "WTJ", "WTJ", "WTJ", 
              "SSKO","SSKO","SSKO",
              "SSKO","SSKO","SSKO")


cond <- c("0","0", "0",
          "oe","oe","oe",
          "0","0","0",
          "oe","oe","oe")

individual.nested <- sprintf("%d", c(1,2,3,
                                     1,2,3,
                                     1,2,3,
                                     1,2,3))


sample_info <- data.frame(exp = factor(group, levels=c("WTJ", "SSKO", "WTJm7oe", "SSKOm7oe")), 
                          genotype= factor(genotype, levels=c("WTJ", "SSKO")), 
                          cond= factor(cond, levels=c("0", "oe")), 
                          individual=factor(individual.nested, levels=c("1", "2", "3")))

m1 <- model.matrix(~0+ genotype + genotype:individual + genotype:cond, sample_info)
v <- DESeqDataSetFromMatrix(SS_DGEm, sample_info, m1)
keep <- rowSums(counts(v)) >= 10
v <- v[keep,]

# run DESeq2
ds <- DESeq(v)
resultsNames(ds)

# calculated WT counts
DGEm_normalized <- counts(ds, normalized=T) %>% as.data.frame
DGEm_normalized_WT <- DGEm_normalized[,c(1:3)] %>% rowMeans
saveRDS(DGEm_normalized, file = "/data/rajewsky/projects/cdr1as_ko_snRNA/codes_github/cdr1as/DGEm/WTJ_nested_normalized_DGEm.rda")



# ---------------------------------- #
# PCA
# ---------------------------------- #

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
          legend.title = element_text(size = 15),
          legend.key = element_rect(colour = "transparent", fill = "white")) +
    scale_color_manual(values = c("#FFCC00", "#FF6600", "#99FF00", "#006600")) + 
    scale_shape_manual(values=c(19, 18, 17, 15))
  
}

Counts <- t(as.matrix(SS_DGEm))
vsd <- varianceStabilizingTransformation(ds)
vsd2 <- vsd

vsd$group <- group
p1 <- plotPCA(vsd, intgroup = "group", returnData = TRUE)
p1$group <- factor(p1$group, levels = c("WTJ", "WTJm7oe", "SSKO", "SSKOm7oe"))


g1 <- pca_plot(pc = p1, group = group, x_pc = "82", y_pc ="8")
g1
ggsave(filename = "/data/rajewsky/projects/cdr1as_ko_snRNA/codes_github/cdr1as/Figures/WTJ_PCA_orig.pdf",
       plot = g1,
       scale = 1, width = 6, height = 5, units = "in", device = cairo_pdf,
       dpi = 300)

# nested corrected PCA

vsd <- vst(t(Counts[, keep]), blind = F) 
CountsNew <- vsd
vsd <- varianceStabilizingTransformation(ds)
r <- removeBatchEffect(CountsNew, batch=individual.nested, design=m1)
vsd2 <- vsd
assay(vsd2) <- r
# save batch & nested corrected matrix for WGCNA and RF
r.matrix <- r %>% as.matrix()
r.df <- r %>% as.data.frame()
saveRDS(r.df, file = "/data/rajewsky/projects/cdr1as_ko_snRNA/codes_github/cdr1as/DGEm/WTJ_nested_DGEm.rda")



vsd2$group <- group
p2 <- plotPCA(vsd2, intgroup = "group", returnData = TRUE)
p2$group <- factor(p2$group, levels = c("WTJ", "WTJm7oe", "SSKO", "SSKOm7oe"))

g2 <- pca_plot(pc = p2, group = group, x_pc = "78", y_pc ="12")
g2
ggsave(filename = "/data/rajewsky/projects/cdr1as_ko_snRNA/codes_github/cdr1as/Figures/WTJ_PCA_nested.pdf",
       plot = g2,
       scale = 1, width = 6, height = 5, units = "in", device = cairo_pdf,
       dpi = 300)

# ---------------------------------- #
# DESeq2 results
# ---------------------------------- #
# WTm7oe vs WT

result2 <- results(ds, name="genotypeWTJ.condoe")
result2.df <- as.data.frame(result2)
result2.df$gene <- rownames(result2.df)
result2.df$WT <- DGEm_normalized_WT
df2 <- mcols(result2) %>% as.data.frame()
result2.df$gene[result2.df$gene == "1700020I14Rik"] <- "Cyrano"
result2.df$gene[result2.df$gene == "C230004F18Rik"] <- "Cdr1os"
saveRDS(result2.df, file = "/data/rajewsky/projects/cdr1as_ko_snRNA/codes_github/cdr1as/Results/DESeq2_WTJ_5SS_nested_corrected_res2.rda")
result2.df <- readRDS(file = "/data/rajewsky/projects/cdr1as_ko_snRNA/codes_github/cdr1as/Results/DESeq2_WTJ_5SS_nested_corrected_res2.rda")
write.csv(result2.df, file = "/data/rajewsky/projects/cdr1as_ko_snRNA/codes_github/cdr1as/Results/DESeq2_WTJ_5SS_nested_corrected_res2.csv")






# KOm7oe vs KO
result3 <- results(ds, name="genotypeSSKO.condoe")
result3.df <- as.data.frame(result3)
result3.df$gene <- rownames(result3.df)
result3.df$WT <- DGEm_normalized_WT
result3.df$gene[result3.df$gene == "1700020I14Rik"] <- "Cyrano"
result3.df$gene[result3.df$gene == "C230004F18Rik"] <- "Cdr1os"
saveRDS(result3.df, file = "/data/rajewsky/projects/cdr1as_ko_snRNA/codes_github/cdr1as/Results/DESeq2_WTJ_5SS_nested_corrected_res3.rda")
result3.df <- readRDS(file = "/data/rajewsky/projects/cdr1as_ko_snRNA/codes_github/cdr1as/Results/DESeq2_WTJ_5SS_nested_corrected_res3.rda")
write.csv(result3.df, file = "/data/rajewsky/projects/cdr1as_ko_snRNA/codes_github/cdr1as/Results/DESeq2_WTJ_5SS_nested_corrected_res3.csv")




# ---------------------------------- #
# DESeq2 results w/o nested corrected
# ---------------------------------- #

m2 <- model.matrix(~ exp, sample_info)
v <- DESeqDataSetFromMatrix(SS_DGEm, sample_info, m2)
keep <- rowSums(counts(v)) >= 10
v <- v[keep,]

# run DESeq2
ds <- DESeq(v)
resultsNames(ds)


# KO vs WT
result1 <- results(ds, name= "expSSKO")
result1.df <- as.data.frame(result1)
result1.df$gene <- rownames(result1.df)
result1.df$WT <- DGEm_normalized_WT
result1.df$gene[result1.df$gene == "1700020I14Rik"] <- "Cyrano"
result1.df$gene[result1.df$gene == "C230004F18Rik"] <- "Cdr1os"
# check reference level
df <- mcols(result1) %>% as.data.frame()
saveRDS(result1.df, file = "/data/rajewsky/projects/cdr1as_ko_snRNA/codes_github/cdr1as/Results/DESeq2_WTJ_5SS_res1.rda")
result1.df <- readRDS(file = "/data/rajewsky/projects/cdr1as_ko_snRNA/codes_github/cdr1as/Results/DESeq2_WTJ_5SS_res1.rda")
write.csv(result1.df, file = "/data/rajewsky/projects/cdr1as_ko_snRNA/codes_github/cdr1as/Results/DESeq2_WTJ_5SS_res1.csv")

# WT vs KOm7oe
result4 <- results(ds, name = "expSSKOm7oe")
result4.df <- as.data.frame(result4)
result4.df$gene <- rownames(result4.df)
result4.df$WT <- DGEm_normalized_WT
result4.df$gene[result4.df$gene == "1700020I14Rik"] <- "Cyrano"
result4.df$gene[result4.df$gene == "C230004F18Rik"] <- "Cdr1os"
saveRDS(result4.df, file = "/data/rajewsky/projects/cdr1as_ko_snRNA/codes_github/cdr1as/Results/DESeq2_WTJ_5SS_res4.rda")
result4.df <- readRDS(file = "/data/rajewsky/projects/cdr1as_ko_snRNA/codes_github/cdr1as/Results/DESeq2_WTJ_5SS_res4.rda")
write.csv(result4.df, file = "/data/rajewsky/projects/cdr1as_ko_snRNA/codes_github/cdr1as/Results/DESeq2_WTJ_5SS_res4.csv")

# ---------------------------------- #
# MA plots
# ---------------------------------- #


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


g1 <- target_MA_Plot(result1.df, labels = "log2FC (WTJ vs 5'SSKO)")  +
  geom_point(data = result1.df[which(result1.df$gene %in% target_genes & result1.df$padj < 0.05), ], aes(x = WT, y = log2FoldChange), color = "darkorange", size = 3) +
  geom_point(data = result1.df[which(result1.df$gene %in% target_genes & result1.df$padj > 0.05), ], aes(x = WT, y = log2FoldChange), color = "darkgreen", size = 3) +
  coord_cartesian(ylim = c(-5,5)) +
  geom_text_repel(data = subset(result1.df, gene %in% c(target_genes[1])),aes(label = gene),
                  size = 4, box.padding = unit(0.35, "lines"), point.padding = unit(0.3, "lines"),
                  nudge_x = 1) +
  geom_text_repel(data = subset(result1.df, gene %in% c(target_genes[2])),aes(label = gene),
                  size = 4, box.padding = unit(0.35, "lines"), point.padding = unit(0.3, "lines"),
                  nudge_x =  -1)
g1

ggsave(filename = "/data/rajewsky/projects/cdr1as_ko_snRNA/codes_github/cdr1as/Figures/DESeq2_WTJ_res1.pdf",
       plot = g1,
       scale = 1, width = 5, height = 4, units = "in", device = cairo_pdf,
       dpi = 300)

g2 <- target_MA_Plot(result2.df, labels = "log2FC (WTJ vs WTJ-miR-7oe)")  +
  geom_point(data = result2.df[which(result2.df$gene %in% target_genes & result2.df$padj < 0.05), ], aes(x = WT, y = log2FoldChange), color = "darkorange", size = 3) +
  geom_point(data = result2.df[which(result2.df$gene %in% target_genes & result2.df$padj > 0.05), ], aes(x = WT, y = log2FoldChange), color = "darkgreen", size = 3) +
  coord_cartesian(ylim = c(-5,5)) +
  geom_text_repel(data = subset(result2.df, gene %in% c(target_genes[1])),aes(label = gene),
                  size = 4, box.padding = unit(0.35, "lines"), point.padding = unit(0.3, "lines"),
                  nudge_y = - 0.2) +
  geom_text_repel(data = subset(result2.df, gene %in% c(target_genes[2])),aes(label = gene),
                  size = 4, box.padding = unit(0.35, "lines"), point.padding = unit(0.3, "lines"),
                  nudge_y =  0.2)
g2

ggsave(filename = "/data/rajewsky/projects/cdr1as_ko_snRNA/codes_github/cdr1as/Figures/DESeq2_WTJ_res2.pdf",
       plot = g2,
       scale = 1, width = 5, height = 4, units = "in", device = cairo_pdf,
       dpi = 300)

g3 <- target_MA_Plot(result3.df, labels = "log2FC (5'SSKO vs 5'SSKO-miR-7oe)")  +
  geom_point(data = result3.df[which(result3.df$gene %in% target_genes & result3.df$padj < 0.05), ], aes(x = WT, y = log2FoldChange), color = "darkorange", size = 3) +
  geom_point(data = result3.df[which(result3.df$gene %in% target_genes & result3.df$padj > 0.05), ], aes(x = WT, y = log2FoldChange), color = "darkgreen", size = 3) +
  coord_cartesian(ylim = c(-5,5)) +
  geom_text_repel(data = subset(result3.df, gene %in% c(target_genes[1])),aes(label = gene),
                  size = 4, box.padding = unit(0.35, "lines"), point.padding = unit(0.3, "lines"),
                  nudge_x = 1) +
  geom_text_repel(data = subset(result3.df, gene %in% c(target_genes[2])),aes(label = gene),
                  size = 4, box.padding = unit(0.35, "lines"), point.padding = unit(0.3, "lines"),
                  nudge_y =  1.2)
g3

ggsave(filename = "/data/rajewsky/projects/cdr1as_ko_snRNA/codes_github/cdr1as/Figures/DESeq2_WTJ_res3.pdf",
       plot = g3,
       scale = 1, width = 5, height = 4, units = "in", device = cairo_pdf,
       dpi = 300)

g4 <- target_MA_Plot(result4.df, labels = "log2FC (WTJ vs 5'SSKO-miR-7oe)")  +
  geom_point(data = result4.df[which(result4.df$gene %in% target_genes & result4.df$padj < 0.05), ], aes(x = WT, y = log2FoldChange), color = "darkorange", size = 3) +
  geom_point(data = result4.df[which(result4.df$gene %in% target_genes & result4.df$padj > 0.05), ], aes(x = WT, y = log2FoldChange), color = "darkgreen", size = 3) +
  coord_cartesian(ylim = c(-5,5)) +
  geom_text_repel(data = subset(result4.df, gene %in% c(target_genes[1])),aes(label = gene),
                  size = 4, box.padding = unit(0.35, "lines"), point.padding = unit(0.3, "lines"),
                  nudge_y = 0.1) +
  geom_text_repel(data = subset(result4.df, gene %in% c(target_genes[2])),aes(label = gene),
                  size = 4, box.padding = unit(0.35, "lines"), point.padding = unit(0.3, "lines"),
                  nudge_y =  0.5, nudge_x = -1)
g4

ggsave(filename = "/data/rajewsky/projects/cdr1as_ko_snRNA/codes_github/cdr1as/Figures/DESeq2_WTJ_res4.pdf",
       plot = g4,
       scale = 1, width = 5, height = 4, units = "in", device = cairo_pdf,
       dpi = 300)
