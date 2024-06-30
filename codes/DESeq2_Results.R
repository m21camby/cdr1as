.libPaths(c("/data/rajewsky/home/skim/R/usr_lib_Seurat/", "/data/rajewsky/home/skim/R/usr_lib_v4/"))

library(DESeq2)
library(ggplot2)
library(dplyr)
library(Cairo)
library(ggplot2)
library(genefilter)
library(paletteer)
library(ggrepel)

res1.df <- readRDS(file = "/data/rajewsky/home/skim/Neuronal_Activity_Cledi/201911_m7_overexpression/publication_codes_figures/DE_Analysis_DESeq2_res1_KO_WT.rda")
res1.df$gene[res1.df$gene == "1700020I14Rik"] <- "Cyrano"
res1.df$gene[res1.df$gene == "C230004F18Rik"] <- "Cdr1os"

result2.df <- readRDS(file = "/data/rajewsky/home/skim/Neuronal_Activity_Cledi/201911_m7_overexpression/publication_codes_figures/DE_Analysis_DESeq2_nested_result2.rda")
result2.df$gene[result2.df$gene == "1700020I14Rik"] <- "Cyrano"
result2.df$gene[result2.df$gene == "C230004F18Rik"] <- "Cdr1os"

result3.df <- readRDS(file = "/data/rajewsky/home/skim/Neuronal_Activity_Cledi/201911_m7_overexpression/publication_codes_figures/DE_Analysis_DESeq2_nested_result3.rda")
result3.df$gene[result3.df$gene == "1700020I14Rik"] <- "Cyrano"
result3.df$gene[result3.df$gene == "C230004F18Rik"] <- "Cdr1os"

res4.df <- readRDS(file = "/data/rajewsky/home/skim/Neuronal_Activity_Cledi/201911_m7_overexpression/publication_codes_figures/DE_Analysis_DESeq2_res4_KOm7ove_WT.rda")
res4.df$gene[res4.df$gene == "1700020I14Rik"] <- "Cyrano"
res4.df$gene[res4.df$gene == "C230004F18Rik"] <- "Cdr1os"

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

# KO vs WT
g1 <- target_MA_Plot(res1.df, labels = "log2FC (WT vs KO)")  +
  geom_point(data = res1.df[which(res1.df$gene %in% target_genes & res1.df$padj < 0.05), ], aes(x = WT, y = log2FoldChange), color = "darkorange", size = 3) + 
  geom_point(data = res1.df[which(res1.df$gene %in% target_genes & res1.df$padj > 0.05), ], aes(x = WT, y = log2FoldChange), color = "darkgreen", size = 3) + 
  coord_cartesian(ylim = c(-5,5)) + 
  geom_text_repel(data = subset(res1.df, gene %in% c(target_genes[1])),aes(label = gene),
                  size = 4, box.padding = unit(0.35, "lines"), point.padding = unit(0.3, "lines"),
                  nudge_y = 0.2) + 
  geom_text_repel(data = subset(res1.df, gene %in% c(target_genes[2])),aes(label = gene),
                  size = 4, box.padding = unit(0.35, "lines"), point.padding = unit(0.3, "lines"),
                  nudge_y =  0.4) + 
  geom_text_repel(data = subset(res1.df, gene %in% c(target_genes[3])),aes(label = gene),
                  size = 4, box.padding = unit(0.35, "lines"), point.padding = unit(0.3, "lines"),
                  nudge_y =  0.4)
g1
ggsave(filename = "/data/rajewsky/projects/cdr1as_ko_snRNA/codes/DESeq2_Results/DESeq2_Results_MA_KO_WT.pdf",
       plot = g1,
       scale = 1, width = 5, height = 4, units = "in", device = cairo_pdf,
       dpi = 300)

# WT miR-7 vs WT
g1 <- target_MA_Plot(result2.df, labels = "log2FC (WT vs WT-miR-7-oe)")  +
  geom_point(data = result2.df[which(result2.df$gene %in% target_genes & result2.df$padj < 0.05), ], aes(x = WT, y = log2FoldChange), color = "darkorange", size = 3) + 
  geom_point(data = result2.df[which(result2.df$gene %in% target_genes & result2.df$padj > 0.05), ], aes(x = WT, y = log2FoldChange), color = "darkgreen", size = 3) + 
  coord_cartesian(ylim = c(-5,5)) + 
  geom_text_repel(data = subset(result2.df, gene %in% c(target_genes[1])),aes(label = gene),
                  size = 4, box.padding = unit(0.35, "lines"), point.padding = unit(0.3, "lines"),
                  nudge_x = 0.2, nudge_y = -0.7) + 
  geom_text_repel(data = subset(result2.df, gene %in% c(target_genes[2])),aes(label = gene),
                  size = 4, box.padding = unit(0.35, "lines"), point.padding = unit(0.3, "lines"),
                  nudge_x = 0.1, nudge_y =  0.8) + 
  geom_text_repel(data = subset(result2.df, gene %in% c(target_genes[3])),aes(label = gene),
                  size = 4, box.padding = unit(0.35, "lines"), point.padding = unit(0.3, "lines"),
                  nudge_y =  -0.4)
g1
ggsave(filename = "/data/rajewsky/projects/cdr1as_ko_snRNA/codes/DESeq2_Results/DESeq2_Results_MA_WTm7oe_WT.pdf",
       plot = g1,
       scale = 1, width = 5, height = 4, units = "in", device = cairo_pdf,
       dpi = 300)


# KO miR-7 vs KO
g1 <- target_MA_Plot(result3.df, labels = "log2FC (KO vs KO-miR-7-oe)")  +
  geom_point(data = result3.df[which(result3.df$gene %in% target_genes & result3.df$padj < 0.05), ], aes(x = WT, y = log2FoldChange), color = "darkorange", size = 3) + 
  geom_point(data = result3.df[which(result3.df$gene %in% target_genes & result3.df$padj > 0.05), ], aes(x = WT, y = log2FoldChange), color = "darkgreen", size = 3) + 
  coord_cartesian(ylim = c(-5,5)) + 
  geom_text_repel(data = subset(result3.df, gene %in% c(target_genes[1])),aes(label = gene),
                  size = 4, box.padding = unit(0.35, "lines"), point.padding = unit(0.3, "lines"),
                  nudge_x = 0, nudge_y = -0.7) + 
  geom_text_repel(data = subset(result3.df, gene %in% c(target_genes[2])),aes(label = gene),
                  size = 4, box.padding = unit(0.35, "lines"), point.padding = unit(0.3, "lines"),
                  nudge_x = 0, nudge_y =  0.8) + 
  geom_text_repel(data = subset(result3.df, gene %in% c(target_genes[3])),aes(label = gene),
                  size = 4, box.padding = unit(0.35, "lines"), point.padding = unit(0.3, "lines"),
                  nudge_y =  -0.4)
g1
ggsave(filename = "/data/rajewsky/projects/cdr1as_ko_snRNA/codes/DESeq2_Results/DESeq2_Results_MA_KOm7oe_KO.pdf",
       plot = g1,
       scale = 1, width = 5, height = 4, units = "in", device = cairo_pdf,
       dpi = 300)

# KO miR-7 vs WT
g1 <- target_MA_Plot(res4.df, labels = "log2FC (WT vs KO-miR-7-oe)")  +
  geom_point(data = res4.df[which(res4.df$gene %in% target_genes & res4.df$padj < 0.05), ], aes(x = WT, y = log2FoldChange), color = "darkorange", size = 3) + 
  geom_point(data = res4.df[which(res4.df$gene %in% target_genes & res4.df$padj > 0.05), ], aes(x = WT, y = log2FoldChange), color = "darkgreen", size = 3) + 
  coord_cartesian(ylim = c(-5,5)) + 
  geom_text_repel(data = subset(res4.df, gene %in% c(target_genes[1])),aes(label = gene),
                  size = 4, box.padding = unit(0.35, "lines"), point.padding = unit(0.3, "lines"),
                  nudge_x = 0.2, nudge_y = -0.7) + 
  geom_text_repel(data = subset(res4.df, gene %in% c(target_genes[2])),aes(label = gene),
                  size = 4, box.padding = unit(0.35, "lines"), point.padding = unit(0.3, "lines"),
                  nudge_x = 0.2, nudge_y =  0.8) + 
  geom_text_repel(data = subset(res4.df, gene %in% c(target_genes[3])),aes(label = gene),
                  size = 4, box.padding = unit(0.35, "lines"), point.padding = unit(0.3, "lines"),
                  nudge_y =  -0.4)
g1
ggsave(filename = "/data/rajewsky/projects/cdr1as_ko_snRNA/codes/DESeq2_Results/DESeq2_Results_MA_KOm7oe_WT.pdf",
       plot = g1,
       scale = 1, width = 5, height = 4, units = "in", device = cairo_pdf,
       dpi = 300)

