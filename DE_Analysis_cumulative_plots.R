
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
# batch corrected data
res1.df <- readRDS(file = "/data/rajewsky/projects/cdr1as_ko_snRNA/codes_github/cdr1as/Results/DESeq2_WTN_batch_corrected_res1.rda")
res2.df <- readRDS(file = "/data/rajewsky/projects/cdr1as_ko_snRNA/codes_github/cdr1as/Results/DESeq2_WTN_batch_corrected_res2.rda")
res3.df <- readRDS(file = "/data/rajewsky/projects/cdr1as_ko_snRNA/codes_github/cdr1as/Results/DESeq2_WTN_batch_corrected_res3.rda")
res4.df <- readRDS(file = "/data/rajewsky/projects/cdr1as_ko_snRNA/codes_github/cdr1as/Results/DESeq2_WTN_batch_corrected_res4.rda")

# nested & batch corrected data
nested_res2.df <- readRDS(file = "/data/rajewsky/projects/cdr1as_ko_snRNA/codes_github/cdr1as/Results/DESeq2_WTN_batch_nested_corrected_res2.rda")
nested_res3.df <- readRDS(file = "/data/rajewsky/projects/cdr1as_ko_snRNA/codes_github/cdr1as/Results/DESeq2_WTN_batch_nested_corrected_res3.rda")



# miR-7 target genes 
mir7target.df <- read.table("/data/rajewsky/projects/cdr1as_ko_snRNA/codes_github/cdr1as/ext_files/TargetScan7.2__miR-7-5p.predicted_targets_lists_without_header.txt")
# miR-7 targe genes from mirdb
mir7mirdb.df <- read.csv("/data/rajewsky/projects/cdr1as_ko_snRNA/codes_github/cdr1as/ext_files//20220423_mirdb_mir_7_extract.csv", header = TRUE)

# miR-7 target genes 
mir122target.df <- read.table("/data/rajewsky/projects/cdr1as_ko_snRNA/codes_github/cdr1as/ext_files/TargetScan7.2__miR-122-5p.predicted_targets_without_header.txt")
# miR-7 targe genes from mirdb
mir122mirdb.df <- read.csv("/data/rajewsky/projects/cdr1as_ko_snRNA/codes_github/cdr1as/ext_files/20220423_mirdb_mir_122_extract.csv", header = TRUE)

# miR-7 target genes (in both DB)
mir7_genes <- intersect(mir7target.df$gene_id, mir7mirdb.df$X4)
mir122_genes <- intersect(mir122target.df$V1, mir122mirdb.df$X4)


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
# WT vs KO miR-7
#################
res1.df$is.target <- rownames(res1.df) %in% mir7_genes
wilcox.test(log2FoldChange ~ is.target, data = res1.df, alt = "two.sided")

c1 <- cumulative_plot(res1.df, "miR-7-5p targets (", "p-value = 0.7805") + ggtitle("WTN vs FullKO")

ggsave(filename = "/data/rajewsky/projects/cdr1as_ko_snRNA/codes_github/cdr1as/Figures/cumulative_plots_miR_7_KO_WTN.pdf",
       plot = c1,
       scale = 1, width = 7, height = 5.5, units = "in", device = cairo_pdf,
       dpi = 300)


#################
# WT vs KO miR-122
#################
res1.df$is.target <- rownames(res1.df) %in% mir122_genes
wilcox.test(log2FoldChange ~ is.target, data = res1.df, alt = "two.sided")

c1 <- cumulative_plot(res1.df, "miR-122-5p targets (", "p-value = 0.9384") + ggtitle("WTN vs FullKO")

ggsave(filename = "/data/rajewsky/projects/cdr1as_ko_snRNA/codes_github/cdr1as/Figures/cumulative_plots_miR_122_KO_WTN.pdf",
       plot = c1,
       scale = 1, width = 7, height = 5.5, units = "in", device = cairo_pdf,
       dpi = 300)

#################
# WT vs WTm7oe miR-7
#################
res2.df$is.target <- rownames(res2.df) %in% mir7_genes
wilcox.test(log2FoldChange ~ is.target, data = res2.df, alt = "two.sided")

c1 <- cumulative_plot(res2.df, "miR-7-5p targets (", "p-value < 2.2e-16") + ggtitle("WTN vs WTNm7oe")

ggsave(filename = "/data/rajewsky/projects/cdr1as_ko_snRNA/codes_github/cdr1as/Figures/cumulative_plots_miR_7_WTNm7oe_WTN.pdf",
       plot = c1,
       scale = 1, width = 7, height = 5.5, units = "in", device = cairo_pdf,
       dpi = 300)

#################
# WT vs WTm7oe miR-122
#################
res2.df$is.target <- rownames(res2.df) %in% mir122_genes
wilcox.test(log2FoldChange ~ is.target, data = res2.df, alt = "two.sided")

c1 <- cumulative_plot(res2.df, "miR-122-5p targets (", "p-value = 0.7937") + ggtitle("WTN vs WTNm7oe")

ggsave(filename = "/data/rajewsky/projects/cdr1as_ko_snRNA/codes_github/cdr1as/Figures/cumulative_plots_miR_122_WTNm7oe_WTN.pdf",
       plot = c1,
       scale = 1, width = 7, height = 5.5, units = "in", device = cairo_pdf,
       dpi = 300)


#################
# KO vs KOm7oe miR-7
#################
res3.df$is.target <- rownames(res3.df) %in% mir7_genes
wilcox.test(log2FoldChange ~ is.target, data = res3.df, alt = "two.sided")

c1 <- cumulative_plot(res3.df, "miR-7-5p targets (", "p-value < 2.2e-16") + ggtitle("FullKO vs KOm7oe")

ggsave(filename = "/data/rajewsky/projects/cdr1as_ko_snRNA/codes_github/cdr1as/Figures/cumulative_plots_miR_7_KOm7oe_KO.pdf",
       plot = c1,
       scale = 1, width = 7, height = 5.5, units = "in", device = cairo_pdf,
       dpi = 300)

#################
# KO vs KOm7oe miR-122
#################
res3.df$is.target <- rownames(res3.df) %in% mir122_genes
wilcox.test(log2FoldChange ~ is.target, data = res3.df, alt = "two.sided")

c1 <- cumulative_plot(res3.df, "miR-122-5p targets (", "p-value = 0.6405") + ggtitle("FullKO vs KOm7oe")

ggsave(filename = "/data/rajewsky/projects/cdr1as_ko_snRNA/codes_github/cdr1as/Figures/cumulative_plots_miR_122_KOm7oe_KO.pdf",
       plot = c1,
       scale = 1, width = 7, height = 5.5, units = "in", device = cairo_pdf,
       dpi = 300)


#################
# WT vs KOm7oe miR-7
#################
res4.df$is.target <- rownames(res4.df) %in% mir7_genes
wilcox.test(log2FoldChange ~ is.target, data = res4.df, alt = "two.sided")

c1 <- cumulative_plot(res4.df, "miR-7-5p targets (", "p-value < 2.2e-16") + ggtitle("WTN vs KOm7oe")

ggsave(filename = "/data/rajewsky/projects/cdr1as_ko_snRNA/codes_github/cdr1as/Figures/cumulative_plots_miR_7_KOm7oe_WTN.pdf",
       plot = c1,
       scale = 1, width = 7, height = 5.5, units = "in", device = cairo_pdf,
       dpi = 300)

#################
# WT vs KOm7oe miR-122
#################
res4.df$is.target <- rownames(res4.df) %in% mir122_genes
wilcox.test(log2FoldChange ~ is.target, data = res4.df, alt = "two.sided")

c1 <- cumulative_plot(res4.df, "miR-122-5p targets (", "p-value = 0.4317") + ggtitle("WTN vs KOm7oe")

ggsave(filename = "/data/rajewsky/projects/cdr1as_ko_snRNA/codes_github/cdr1as/Figures/cumulative_plots_miR_122_KOm7oe_WTN.pdf",
       plot = c1,
       scale = 1, width = 7, height = 5.5, units = "in", device = cairo_pdf,
       dpi = 300)


#----------------------------#
# batch and nested corrected
#----------------------------#

#################
# WT vs WTm7oe miR-7
#################
nested_res2.df$is.target <- rownames(nested_res2.df) %in% mir7_genes
wilcox.test(log2FoldChange ~ is.target, data = nested_res2.df, alt = "two.sided")

c1 <- cumulative_plot(nested_res2.df, "miR-7-5p targets (", "p-value < 2.2e-16") + ggtitle("WTN vs WTNm7oe")

ggsave(filename = "/data/rajewsky/projects/cdr1as_ko_snRNA/codes_github/cdr1as/Figures/cumulative_plots_miR_7_WTNm7oe_WTN_nested.pdf",
       plot = c1,
       scale = 1, width = 7, height = 5.5, units = "in", device = cairo_pdf,
       dpi = 300)

#################
# WT vs WTm7oe miR-122
#################
nested_res2.df$is.target <- rownames(nested_res2.df) %in% mir122_genes
wilcox.test(log2FoldChange ~ is.target, data = nested_res2.df, alt = "two.sided")

c1 <- cumulative_plot(nested_res2.df, "miR-122-5p targets (", "p-value = 0.484") + ggtitle("WTN vs WTNm7oe")

ggsave(filename = "/data/rajewsky/projects/cdr1as_ko_snRNA/codes_github/cdr1as/Figures/cumulative_plots_miR_122_WTNm7oe_WTN_nested.pdf",
       plot = c1,
       scale = 1, width = 7, height = 5.5, units = "in", device = cairo_pdf,
       dpi = 300)


#################
# KO vs KOm7oe miR-7
#################
nested_res3.df$is.target <- rownames(nested_res3.df) %in% mir7_genes
wilcox.test(log2FoldChange ~ is.target, data = nested_res3.df, alt = "two.sided")

c1 <- cumulative_plot(nested_res3.df, "miR-7-5p targets (", "p-value < 2.2e-16") + ggtitle("FullKO vs KOm7oe")

ggsave(filename = "/data/rajewsky/projects/cdr1as_ko_snRNA/codes_github/cdr1as/Figures/cumulative_plots_miR_7_KOm7oe_KO_nested.pdf",
       plot = c1,
       scale = 1, width = 7, height = 5.5, units = "in", device = cairo_pdf,
       dpi = 300)

#################
# KO vs KOm7oe miR-122
#################
nested_res3.df$is.target <- rownames(nested_res3.df) %in% mir122_genes
wilcox.test(log2FoldChange ~ is.target, data = nested_res3.df, alt = "two.sided")

c1 <- cumulative_plot(nested_res3.df, "miR-122-5p targets (", "p-value = 0.6541") + ggtitle("FullKO vs KOm7oe")

ggsave(filename = "/data/rajewsky/projects/cdr1as_ko_snRNA/codes_github/cdr1as/Figures/cumulative_plots_miR_122_KOm7oe_KO_nested.pdf",
       plot = c1,
       scale = 1, width = 7, height = 5.5, units = "in", device = cairo_pdf,
       dpi = 300)












