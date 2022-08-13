
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


# ------------------------- #
# data loading
# ------------------------- #

result1.df <- readRDS(file = "/data/rajewsky/projects/cdr1as_ko_snRNA/codes_github/cdr1as/Results/DESeq2_WTJ_5SS_res1.rda")
result2.df <- readRDS(file = "/data/rajewsky/projects/cdr1as_ko_snRNA/codes_github/cdr1as/Results/DESeq2_WTJ_5SS_nested_corrected_res2.rda")
result3.df <- readRDS(file = "/data/rajewsky/projects/cdr1as_ko_snRNA/codes_github/cdr1as/Results/DESeq2_WTJ_5SS_nested_corrected_res3.rda")
result4.df <- readRDS(file = "/data/rajewsky/projects/cdr1as_ko_snRNA/codes_github/cdr1as/Results/DESeq2_WTJ_5SS_res4.rda")

# miR-7 target genes 
mir7target.df <- read.table("/data/rajewsky/projects/cdr1as_ko_snRNA/codes_github/cdr1as/ext_files/TargetScan7.2__miR-7-5p.predicted_targets_lists_without_header.txt")
# miR-7 targe genes from mirdb
mir7mirdb.df <- read.csv("/data/rajewsky/projects/cdr1as_ko_snRNA/codes_github/cdr1as/ext_files//20220423_mirdb_mir_7_extract.csv", header = TRUE)

# miR-7 target genes 
mir122target.df <- read.table("/data/rajewsky/projects/cdr1as_ko_snRNA/codes_github/cdr1as/ext_files/TargetScan7.2__miR-122-5p.predicted_targets_without_header.txt")
# miR-7 targe genes from mirdb
mir122mirdb.df <- read.csv("/data/rajewsky/projects/cdr1as_ko_snRNA/codes_github/cdr1as/ext_files/20220423_mirdb_mir_122_extract.csv", header = TRUE)

# miR-7 target genes (in both DB)
mir7_genes <- intersect(mir7target.df$V1, mir7mirdb.df$X4)
mir122_genes <- intersect(mir122target.df$V1, mir122mirdb.df$X4)


###### function for plot
cumulative_plot <- function(DF, target_labels, labels){
  ggplot(DF, aes(x = log2FoldChange, color = is.target)) +
    stat_ecdf() +
    theme_classic() +
    ylab("Cumulative fraction") +
    xlim(-2, 2) +
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
    annotate("text", x = 1, y = 0.5, label = labels, colour = "#003366", face = "bold", size = 5, family="Helvetica")
}

#################
# WT vs KO miR-7
#################
result1.df$is.target <- rownames(result1.df) %in% mir7_genes
wilcox.test(log2FoldChange ~ is.target, data = result1.df, alt = "two.sided")

c1 <- cumulative_plot(result1.df, "miR-7-5p targets (", "p-value = 1.795e-10") + ggtitle("WTJ vs 5'SSKO")

ggsave(filename = "/data/rajewsky/projects/cdr1as_ko_snRNA/codes_github/cdr1as/Figures/cumulative_plots_miR_7_5SSKO_WTJ.pdf",
       plot = c1,
       scale = 1, width = 7, height = 5.5, units = "in", device = cairo_pdf,
       dpi = 300)

#################
# WT vs KO miR-122
#################
result1.df$is.target <- rownames(result1.df) %in% mir122_genes
wilcox.test(log2FoldChange ~ is.target, data = result1.df, alt = "two.sided")

c1 <- cumulative_plot(result1.df, "miR-122-5p targets (", "p-value = 0.03493") + ggtitle("WTJ vs 5'SSKO")

ggsave(filename = "/data/rajewsky/projects/cdr1as_ko_snRNA/codes_github/cdr1as/Figures/cumulative_plots_miR_122_5SSKO_WTJ.pdf",
       plot = c1,
       scale = 1, width = 7, height = 5.5, units = "in", device = cairo_pdf,
       dpi = 300)

#################
# WT vs WTm7oe miR-7
#################
result2.df$is.target <- rownames(result2.df) %in% mir7_genes
wilcox.test(log2FoldChange ~ is.target, data = result2.df, alt = "two.sided")

c1 <- cumulative_plot(result2.df, "miR-7-5p targets (", "p-value < 2.2e-16") + ggtitle("WTJ vs WTJm7oe")

ggsave(filename = "/data/rajewsky/projects/cdr1as_ko_snRNA/codes_github/cdr1as/Figures/cumulative_plots_miR_7_WTJm7oe_WTJ.pdf",
       plot = c1,
       scale = 1, width = 7, height = 5.5, units = "in", device = cairo_pdf,
       dpi = 300)

#################
# WT vs WTm7oe miR-122
#################
result2.df$is.target <- rownames(result2.df) %in% mir122_genes
wilcox.test(log2FoldChange ~ is.target, data = result2.df, alt = "two.sided")

c1 <- cumulative_plot(result2.df, "miR-122-5p targets (", "p-value = 0.2113") + ggtitle("WTJ vs WTJm7oe")

ggsave(filename = "/data/rajewsky/projects/cdr1as_ko_snRNA/codes_github/cdr1as/Figures/cumulative_plots_miR_122_WTJm7oe_WTJ.pdf",
       plot = c1,
       scale = 1, width = 7, height = 5.5, units = "in", device = cairo_pdf,
       dpi = 300)

#################
# KO vs KOm7oe miR-7
#################
result3.df$is.target <- rownames(result3.df) %in% mir7_genes
wilcox.test(log2FoldChange ~ is.target, data = result3.df, alt = "two.sided")

c1 <- cumulative_plot(result3.df, "miR-7-5p targets (", "p-value = 0.02037") + ggtitle("5'SSKO vs 5'SSKOm7oe")

ggsave(filename = "/data/rajewsky/projects/cdr1as_ko_snRNA/codes_github/cdr1as/Figures/cumulative_plots_miR_7_5SSKOm7oe_5SSKO.pdf",
       plot = c1,
       scale = 1, width = 7, height = 5.5, units = "in", device = cairo_pdf,
       dpi = 300)

#################
# KO vs KOm7oe miR-122
#################
result3.df$is.target <- rownames(result3.df) %in% mir122_genes
wilcox.test(log2FoldChange ~ is.target, data = result3.df, alt = "two.sided")

c1 <- cumulative_plot(result3.df, "miR-122-5p targets (", "p-value = 0.1409") + ggtitle("5'SSKO vs 5'SSKOm7oe")

ggsave(filename = "/data/rajewsky/projects/cdr1as_ko_snRNA/codes_github/cdr1as/Figures/cumulative_plots_miR_122_5SSKOm7oe_5SSKO.pdf",
       plot = c1,
       scale = 1, width = 7, height = 5.5, units = "in", device = cairo_pdf,
       dpi = 300)

#################
# WT vs KOm7oe miR-7
#################
result4.df$is.target <- rownames(result4.df) %in% mir7_genes
wilcox.test(log2FoldChange ~ is.target, data = result4.df, alt = "two.sided")

c1 <- cumulative_plot(result4.df, "miR-7-5p targets (", "p-value  = 2.144e-08") + ggtitle("WTJ vs 5'SSKOm7oe")

ggsave(filename = "/data/rajewsky/projects/cdr1as_ko_snRNA/codes_github/cdr1as/Figures/cumulative_plots_miR_7_5SSKOm7oe_WTJ.pdf",
       plot = c1,
       scale = 1, width = 7, height = 5.5, units = "in", device = cairo_pdf,
       dpi = 300)

#################
# WT vs KOm7oe miR-122
#################
result4.df$is.target <- rownames(result4.df) %in% mir122_genes
wilcox.test(log2FoldChange ~ is.target, data = result4.df, alt = "two.sided")

c1 <- cumulative_plot(result4.df, "miR-122-5p targets (", "p-value = 0.3379") + ggtitle("WTJ vs 5'SSKOm7oe")

ggsave(filename = "/data/rajewsky/projects/cdr1as_ko_snRNA/codes_github/cdr1as/Figures/cumulative_plots_miR_122_5SSKOm7oe_WTJ.pdf",
       plot = c1,
       scale = 1, width = 7, height = 5.5, units = "in", device = cairo_pdf,
       dpi = 300)
