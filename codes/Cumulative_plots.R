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

res1.df <- readRDS(file = "/data/rajewsky/home/skim/Neuronal_Activity_Cledi/201911_m7_overexpression/publication_codes_figures/DE_Analysis_DESeq2_res1_KO_WT.rda")

res4.df <- readRDS(file = "/data/rajewsky/home/skim/Neuronal_Activity_Cledi/201911_m7_overexpression/publication_codes_figures/DE_Analysis_DESeq2_res4_KOm7ove_WT.rda")

result2.df <- readRDS(file = "/data/rajewsky/home/skim/Neuronal_Activity_Cledi/201911_m7_overexpression/publication_codes_figures/DE_Analysis_DESeq2_nested_result2.rda")

result3.df <- readRDS(file = "/data/rajewsky/home/skim/Neuronal_Activity_Cledi/201911_m7_overexpression/publication_codes_figures/DE_Analysis_DESeq2_nested_result3.rda")

# normal
# up-genes
set1 <- res1.df[which(res1.df$log2FoldChange > 0.5 &  res1.df$padj < 0.05), ]
set2 <- result2.df[which(result2.df$log2FoldChange > 0.5 &  result2.df$padj < 0.05), ]
set3 <- result3.df[which(result3.df$log2FoldChange > 0.5 &  result3.df$padj < 0.05), ]
set4 <- res4.df[which(res4.df$log2FoldChange > 0.5 &  res4.df$padj < 0.05), ] 


# down-genes
set5 <- res1.df[which(res1.df$log2FoldChange < -0.5 &  res1.df$padj < 0.05), ]
set6 <- result2.df[which(result2.df$log2FoldChange < -0.5 &  result2.df$padj < 0.05), ]
set7 <- result3.df[which(result3.df$log2FoldChange < -0.5 &  result3.df$padj < 0.05), ]
set8 <- res4.df[which(res4.df$log2FoldChange < -0.5 &  res4.df$padj < 0.05), ] 

genes <- c(set1$gene, set2$gene, set3$gene, set4$gene,
           set5$gene, set6$gene, set7$gene, set8$gene)

# miR-7 target genes 
mir7target.df <- read.table("/data/rajewsky/home/skim/Neuronal_Activity_Cledi/201911_m7_overexpression/R_Scripts/external_files/TargetScan7.2__miR-7-5p.predicted_targets_lists_without_header.txt")

# miR-7 targe genes from mirdb
mir7mirdb.df <- read.csv("/data/rajewsky/home/skim/Neuronal_Activity_Cledi/201911_m7_overexpression/R_Scripts/20220423_mirdb_mir_7_extract.csv", header = TRUE)

# miR-122 target genes 
mir122target.df <- read.table("/data/rajewsky/home/skim/Neuronal_Activity_Cledi/201911_m7_overexpression/R_Scripts/external_files/TargetScan7.2__miR-122-5p.predicted_targets_without_header.txt")

# miR-122 targe genes from mirdb
mir122mirdb.df <- read.csv("/data/rajewsky/home/skim/Neuronal_Activity_Cledi/201911_m7_overexpression/R_Scripts/20220423_mirdb_mir_122_extract.csv", header = TRUE)

# miR-7 target genes (in both DB)
mir7_genes <- intersect(mir7target.df$V1, mir7mirdb.df$X4)

# miR_122 target genes (in both DB)
mir122_genes <- intersect(mir122target.df$V1, mir122mirdb.df$X4)

module_df2 <- read.csv(file = "/data/rajewsky/home/skim/Neuronal_Activity_Cledi/201911_m7_overexpression/publication_codes_figures/WGCNA_Analysis/WGCNA_Analysis_batch_corrected_1350_module2.csv")

# colnames(mir7target.df) <- "gene_id"
# mir7target_module <- left_join(mir7target.df, module_df2, by = "gene_id")

###### function for plot
cumulative_plot <- function(DF, target_labels, labels){
  ggplot(DF, aes(x = log2FoldChange, color = is.target)) +
    stat_ecdf() +
    theme_classic() +
    ylab("Cumulative fraction") +
    xlim(-1, 1) +
    scale_color_manual(values = c("black", "red"),
                       breaks = c(F, T),
                       labels = c(paste0("nontargets (", sum(!DF$is.target), ")"),
                                  paste0(target_labels, sum(DF$is.target), ")")),
                       name   = "") +
    theme(axis.line.y=element_line(size=.5),
          axis.line.x=element_line(size=.5),
          axis.text.x=element_text(size=18,  family="Helvetica", color="black"),
          axis.text.y=element_text(size=18,  family="Helvetica", color="black"),
          axis.title.x=element_text(size=20, family="Helvetica", color="black"),
          axis.title.y=element_text(size=20, family="Helvetica", color="black"),
          legend.text = element_text(size=12, family="Helvetica", color="black"),
          legend.position = c(0.01, 1),
          legend.justification = c(0, 1),
          aspect.ratio=1,
          plot.title=element_text(size=15, face="bold", family="Helvetica", hjust=0.5)) +
    annotate("text", x = 0.6, y = 0.5, label = labels, colour = "#003366", face = "bold", size = 5, family="Helvetica")
}

#################
# WT vs KO 
#################
# miR-7
res1.df$is.target <- rownames(res1.df) %in% mir7_genes
wilcox.test(log2FoldChange ~ is.target, data = res1.df, alt = "two.sided")

c1 <- cumulative_plot(res1.df, "miR-7-5p targets (", "p-value = 0.06342") + ggtitle("WT vs KO")

ggsave(filename = "/data/rajewsky/projects/cdr1as_ko_snRNA/codes/Cumulative_plots/Cumulative_plots_miR_7_KO_WT.pdf",
       plot = c1,
       scale = 1, width = 7, height = 5.5, units = "in", device = cairo_pdf,
       dpi = 300)

# miR-122
res1.df$is.target <- rownames(res1.df) %in% mir122_genes
wilcox.test(log2FoldChange ~ is.target, data = res1.df, alt = "two.sided")

c1 <- cumulative_plot(res1.df, "miR-122-5p targets (", "p-value = 0.8849") + ggtitle("WT vs KO")

ggsave(filename = "/data/rajewsky/projects/cdr1as_ko_snRNA/codes/Cumulative_plots/Cumulative_plots_miR_122_KO_WT.pdf",
       plot = c1,
       scale = 1, width = 7, height = 5.5, units = "in", device = cairo_pdf,
       dpi = 300)

# miR-7 in blue module
res1.df$is.target <- rownames(res1.df) %in% mir7_genes[mir7_genes %in% module_df2[module_df2$color_check %in% "blue",]$gene_id]
wilcox.test(log2FoldChange ~ is.target, data = res1.df, alt = "two.sided")
c1 <- cumulative_plot(res1.df, "miR-7-5p (blue) targets (", "p-value = 3.783e-05") + ggtitle("WT vs KO")

ggsave(filename = "/data/rajewsky/projects/cdr1as_ko_snRNA/codes/Cumulative_plots/Cumulative_plots_miR_7_blue_KO_WT.pdf",
       plot = c1,
       scale = 1, width = 7, height = 5.5, units = "in", device = cairo_pdf,
       dpi = 300)

# miR-7 in non-blue module
res1.df$is.target <- rownames(res1.df) %in% mir7_genes[mir7_genes %in% module_df2[!module_df2$color_check %in% "blue",]$gene_id]
wilcox.test(log2FoldChange ~ is.target, data = res1.df, alt = "two.sided")
c1 <- cumulative_plot(res1.df, "miR-7-5p (non-blue) targets (", "p-value = 0.04749") + ggtitle("WT vs KO")

ggsave(filename = "/data/rajewsky/projects/cdr1as_ko_snRNA/codes/Cumulative_plots/Cumulative_plots_miR_7_non_blue_KO_WT.pdf",
       plot = c1,
       scale = 1, width = 7, height = 5.5, units = "in", device = cairo_pdf,
       dpi = 300)


#################
# WTm7oe vs WT 
#################
# miR-7
result2.df$is.target <- rownames(result2.df) %in% mir7_genes
wilcox.test(log2FoldChange ~ is.target, data = result2.df, alt = "two.sided")

c1 <- cumulative_plot(result2.df, "miR-7-5p targets (", "p-value < 2.2e-16") + ggtitle("WT vs WTm7oe")

ggsave(filename = "/data/rajewsky/projects/cdr1as_ko_snRNA/codes/Cumulative_plots/Cumulative_plots_miR_7_WTm7oe_WT.pdf",
       plot = c1,
       scale = 1, width = 7, height = 5.5, units = "in", device = cairo_pdf,
       dpi = 300)

# miR-122
result2.df$is.target <- rownames(result2.df) %in% mir122_genes
wilcox.test(log2FoldChange ~ is.target, data = result2.df, alt = "two.sided")

c1 <- cumulative_plot(result2.df, "miR-122-5p targets (", "p-value = 0.3514") + ggtitle("WT vs WTm7oe")

ggsave(filename = "/data/rajewsky/projects/cdr1as_ko_snRNA/codes/Cumulative_plots/Cumulative_plots_miR_122_WTm7oe_WT.pdf",
       plot = c1,
       scale = 1, width = 7, height = 5.5, units = "in", device = cairo_pdf,
       dpi = 300)

# miR-7 in blue module
result2.df$is.target <- rownames(result2.df) %in% mir7_genes[mir7_genes %in% module_df2[module_df2$color_check %in% "blue",]$gene_id]
wilcox.test(log2FoldChange ~ is.target, data = result2.df, alt = "two.sided")

c1 <- cumulative_plot(result2.df, "miR-7-5p (blue) targets (", "p-value < 2.2e-16") + ggtitle("WT vs WTm7oe")

ggsave(filename = "/data/rajewsky/projects/cdr1as_ko_snRNA/codes/Cumulative_plots/Cumulative_plots_miR_7_blue_WTm7oe_WT.pdf",
       plot = c1,
       scale = 1, width = 7, height = 5.5, units = "in", device = cairo_pdf,
       dpi = 300)

# miR-7 in non-blue module
result2.df$is.target <- rownames(result2.df) %in% mir7_genes[mir7_genes %in% module_df2[!module_df2$color_check %in% "blue",]$gene_id]
wilcox.test(log2FoldChange ~ is.target, data = result2.df, alt = "two.sided")

c1 <- cumulative_plot(result2.df, "miR-7-5p (non-blue) targets (", "p-value = 0.0001022") + ggtitle("WT vs WTm7oe")

ggsave(filename = "/data/rajewsky/projects/cdr1as_ko_snRNA/codes/Cumulative_plots/Cumulative_plots_miR_7_non_blue_WTm7oe_WT.pdf",
       plot = c1,
       scale = 1, width = 7, height = 5.5, units = "in", device = cairo_pdf,
       dpi = 300)

#################
# KOm7oe vs KO 
#################
# miR-7
result3.df$is.target <- rownames(result3.df) %in% mir7_genes
wilcox.test(log2FoldChange ~ is.target, data = result3.df, alt = "two.sided")

c1 <- cumulative_plot(result3.df, "miR-7-5p targets (", "p-value < 2.2e-16") + ggtitle("KO vs KOm7oe")

ggsave(filename = "/data/rajewsky/projects/cdr1as_ko_snRNA/codes/Cumulative_plots/Cumulative_plots_miR_7_KOm7oe_KO.pdf",
       plot = c1,
       scale = 1, width = 7, height = 5.5, units = "in", device = cairo_pdf,
       dpi = 300)

# miR-122
result3.df$is.target <- rownames(result3.df) %in% mir122_genes
wilcox.test(log2FoldChange ~ is.target, data = result3.df, alt = "two.sided")

c1 <- cumulative_plot(result3.df, "miR-122-5p targets (", "p-value = 0.6354") + ggtitle("KO vs KOm7oe")

ggsave(filename = "/data/rajewsky/projects/cdr1as_ko_snRNA/codes/Cumulative_plots/Cumulative_plots_miR_122_KOm7oe_KO.pdf",
       plot = c1,
       scale = 1, width = 7, height = 5.5, units = "in", device = cairo_pdf,
       dpi = 300)

# miR-7 in blue module
result3.df$is.target <- rownames(result3.df) %in% mir7_genes[mir7_genes %in% module_df2[module_df2$color_check %in% "blue",]$gene_id]
wilcox.test(log2FoldChange ~ is.target, data = result3.df, alt = "two.sided")

c1 <- cumulative_plot(result3.df, "miR-7-5p (blue) targets (", "p-value < 2.2e-16") + ggtitle("KO vs KOm7oe")

ggsave(filename = "/data/rajewsky/projects/cdr1as_ko_snRNA/codes/Cumulative_plots/Cumulative_plots_miR_7_blue_KOm7oe_KO.pdf",
       plot = c1,
       scale = 1, width = 7, height = 5.5, units = "in", device = cairo_pdf,
       dpi = 300)

# miR-7 in non-blue module
result3.df$is.target <- rownames(result3.df) %in% mir7_genes[mir7_genes %in% module_df2[!module_df2$color_check %in% "blue",]$gene_id]
wilcox.test(log2FoldChange ~ is.target, data = result3.df, alt = "two.sided")

c1 <- cumulative_plot(result3.df, "miR-7-5p (non-blue) targets (", "p-value = 1.41e-08") + ggtitle("KO vs KOm7oe")

ggsave(filename = "/data/rajewsky/projects/cdr1as_ko_snRNA/codes/Cumulative_plots/Cumulative_plots_miR_7_non_blue_KOm7oe_KO.pdf",
       plot = c1,
       scale = 1, width = 7, height = 5.5, units = "in", device = cairo_pdf,
       dpi = 300)

#################
# KOm7oe vs WT 
#################
# miR-7
res4.df$is.target <- rownames(res4.df) %in% mir7_genes
wilcox.test(log2FoldChange ~ is.target, data = res4.df, alt = "two.sided")

c1 <- cumulative_plot(res4.df, "miR-7-5p targets (", "p-value < 2.2e-16") + ggtitle("WT vs KOm7oe")

ggsave(filename = "/data/rajewsky/projects/cdr1as_ko_snRNA/codes/Cumulative_plots/Cumulative_plots_miR_7_KOm7oe_WT.pdf",
       plot = c1,
       scale = 1, width = 7, height = 5.5, units = "in", device = cairo_pdf,
       dpi = 300)

# miR-122
res4.df$is.target <- rownames(res4.df) %in% mir122_genes
wilcox.test(log2FoldChange ~ is.target, data = res4.df, alt = "two.sided")

c1 <- cumulative_plot(res4.df, "miR-122-5p targets (", "p-value = 0.3168") + ggtitle("WT vs KOm7oe")

ggsave(filename = "/data/rajewsky/projects/cdr1as_ko_snRNA/codes/Cumulative_plots/Cumulative_plots_miR_122_KOm7oe_WT.pdf",
       plot = c1,
       scale = 1, width = 7, height = 5.5, units = "in", device = cairo_pdf,
       dpi = 300)

# miR-7 in blue module
res4.df$is.target <- rownames(res4.df) %in% mir7_genes[mir7_genes %in% module_df2[module_df2$color_check %in% "blue",]$gene_id]
wilcox.test(log2FoldChange ~ is.target, data = res4.df, alt = "two.sided")

c1 <- cumulative_plot(res4.df, "miR-7-5p (blue) targets (", "p-value < 2.2e-16") + ggtitle("WT vs KOm7oe")

ggsave(filename = "/data/rajewsky/projects/cdr1as_ko_snRNA/codes/Cumulative_plots/Cumulative_plots_miR_7_blue_KOm7oe_WT.pdf",
       plot = c1,
       scale = 1, width = 7, height = 5.5, units = "in", device = cairo_pdf,
       dpi = 300)

# miR-7 in non-blue module
res4.df$is.target <- rownames(res4.df) %in% mir7_genes[mir7_genes %in% module_df2[!module_df2$color_check %in% "blue",]$gene_id]
wilcox.test(log2FoldChange ~ is.target, data = res4.df, alt = "two.sided")

c1 <- cumulative_plot(res4.df, "miR-7-5p (non-blue) targets (", "p-value = 8.789e-09") + ggtitle("WT vs KOm7oe")

ggsave(filename = "/data/rajewsky/projects/cdr1as_ko_snRNA/codes/Cumulative_plots/Cumulative_plots_miR_7_non_blue_KOm7oe_WT.pdf",
       plot = c1,
       scale = 1, width = 7, height = 5.5, units = "in", device = cairo_pdf,
       dpi = 300)
