
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
# batch OR nested corrected data
res1.df <- readRDS(file = "/data/rajewsky/projects/cdr1as_ko_snRNA/codes_github/cdr1as/Results/DESeq2_WTN_batch_corrected_res1.rda")
res2.df <- readRDS(file = "/data/rajewsky/projects/cdr1as_ko_snRNA/codes_github/cdr1as/Results/DESeq2_WTN_nested_corrected_res2.rda")
res3.df <- readRDS(file = "/data/rajewsky/projects/cdr1as_ko_snRNA/codes_github/cdr1as/Results/DESeq2_WTN_nested_corrected_res3.rda")
res4.df <- readRDS(file = "/data/rajewsky/projects/cdr1as_ko_snRNA/codes_github/cdr1as/Results/DESeq2_WTN_batch_corrected_res4.rda")



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

cumulative_combine_plot <- function(DF, target_labels1, target_labels2, labels){
  ggplot(DF, aes(x = log2FoldChange, color = is.target2)) +
    stat_ecdf() +
    theme_classic() + 
    ylab("Cumulative fraction") +
    xlim(-1, 1) +
    scale_color_manual(values = c("#FF9900", "#66CCCC", "darkgray"),
                       breaks = c(F, T, "NA"),
                       labels = c(paste0(target_labels1, " (","269", ")"),
                                  paste0(target_labels2, " (","269", ")"),
                                  paste0("nontargets (", sum(!DF$is.target), ")"))) + 
    theme(axis.line.y=element_line(size=.5),
          axis.line.x=element_line(size=.5),
          axis.text.x=element_text(size=18,  family="Helvetica", color="black"),
          axis.text.y=element_text(size=18,  family="Helvetica", color="black"),
          axis.title.x=element_text(size=20, family="Helvetica", color="black"),
          axis.title.y=element_text(size=20, family="Helvetica", color="black"),
          legend.title = element_blank(),
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

c1 <- cumulative_plot(res2.df, "miR-122-5p targets (", "p-value = 0.6777") + ggtitle("WTN vs WTNm7oe")

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

c1 <- cumulative_plot(res3.df, "miR-122-5p targets (", "p-value = 0.6541") + ggtitle("FullKO vs KOm7oe")

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



#################
# All combined
#################
res1.df$is.target <- rownames(res1.df) %in% mir7_genes
res2.df$is.target <- rownames(res2.df) %in% mir7_genes
res3.df$is.target <- rownames(res3.df) %in% mir7_genes
res4.df$is.target <- rownames(res4.df) %in% mir7_genes

res1.df$exp <- "KO_vs_WT"
res2.df$exp <- "WTm7oe_vs_WT"
res3.df$exp <- "KOm7oe_vs_KO"
res4.df$exp <- "KOm7oe_vs_WT"

# KO_vs_WT & WTm7oe_vs_WT
res12 <- rbind(res1.df, res2.df[res2.df$is.target %in% TRUE, ])
res12$is.target2 <- ifelse(res12$exp %in% "WTm7oe_vs_WT", TRUE,
                           ifelse((res12$exp %in% "KO_vs_WT" & res12$is.target %in% TRUE), FALSE, NA))

wilcox.test(log2FoldChange ~ is.target2, data = res12, alt = "two.sided")

res23$is.target2 <- ifelse(res12$exp %in% "WTm7oe_vs_WT", TRUE,
                           ifelse((res12$exp %in% "KO_vs_WT" & res12$is.target %in% TRUE), FALSE, "NA"))

c1 <- cumulative_combine_plot(res12, "targets (KO vs WT)", "targets (WTm7oe vs WT)", "p-value < 2.2e-16")

ggsave(filename = "/data/rajewsky/projects/cdr1as_ko_snRNA/codes_github/cdr1as/Figures/cumulative_plots_miR_7_WTm7oe_WT_KO_WT.pdf",
       plot = c1,
       scale = 1, width = 7, height = 5.5, units = "in", device = cairo_pdf,
       dpi = 300)



# KOm7oe_vs_KO & WTm7oe_vs_WT
res23 <- rbind(res2.df, res3.df[res3.df$is.target %in% TRUE, ])
res23$is.target2 <- ifelse(res23$exp %in% "KOm7oe_vs_KO", TRUE,
                           ifelse((res23$exp %in% "WTm7oe_vs_WT" & res23$is.target %in% TRUE), FALSE, NA))

wilcox.test(log2FoldChange ~ is.target2, data = res23, alt = "two.sided")

res23$is.target2 <- ifelse(res23$exp %in% "KOm7oe_vs_KO", TRUE,
                           ifelse((res23$exp %in% "WTm7oe_vs_WT" & res23$is.target %in% TRUE), FALSE, "NA"))

c1 <- cumulative_combine_plot(res23, "targets (WTm7oe vs WT)", "targets (KOm7oe vs KO)", "p-value = 3.459e-13")

ggsave(filename = "/data/rajewsky/projects/cdr1as_ko_snRNA/codes_github/cdr1as/Figures/cumulative_plots_miR_7_KOm7oe_KO_WTm7oe_WT.pdf",
       plot = c1,
       scale = 1, width = 7, height = 5.5, units = "in", device = cairo_pdf,
       dpi = 300)

# KOm7oe_vs_WT & WTm7oe_vs_WT
res24 <- rbind(res2.df, res4.df[res4.df$is.target %in% TRUE, ])
res24$is.target2 <- ifelse(res24$exp %in% "KOm7oe_vs_WT", TRUE,
                           ifelse((res24$exp %in% "WTm7oe_vs_WT" & res24$is.target %in% TRUE), FALSE, NA))

wilcox.test(log2FoldChange ~ is.target2, data = res24, alt = "two.sided")

res24$is.target2 <- ifelse(res24$exp %in% "KOm7oe_vs_WT", TRUE,
                           ifelse((res24$exp %in% "WTm7oe_vs_WT" & res24$is.target %in% TRUE), FALSE, "NA"))

c1 <- cumulative_combine_plot(res24, "targets (WTm7oe vs WT)", "targets (KOm7oe vs WT)", "p-value = 1.112e-06")

ggsave(filename = "/data/rajewsky/projects/cdr1as_ko_snRNA/codes_github/cdr1as/Figures/cumulative_plots_miR_7_KOm7oe_WT_WTm7oe_WT.pdf",
       plot = c1,
       scale = 1, width = 7, height = 5.5, units = "in", device = cairo_pdf,
       dpi = 300)


# KOm7oe_vs_WT & KOm7oe_vs_KO
res34 <- rbind(res3.df, res4.df[res4.df$is.target %in% TRUE, ])
res34$is.target2 <- ifelse(res34$exp %in% "KOm7oe_vs_WT", TRUE,
                           ifelse((res34$exp %in% "KOm7oe_vs_KO" & res34$is.target %in% TRUE), FALSE, NA))

wilcox.test(log2FoldChange ~ is.target2, data = res34, alt = "two.sided")

res34$is.target2 <- ifelse(res34$exp %in% "KOm7oe_vs_WT", TRUE,
                           ifelse((res34$exp %in% "KOm7oe_vs_KO" & res34$is.target %in% TRUE), FALSE, "NA"))

c1 <- cumulative_combine_plot(res34, "targets (KOm7oe vs KO)", "targets (KOm7oe vs WT)", "p-value = 0.1281")

ggsave(filename = "/data/rajewsky/projects/cdr1as_ko_snRNA/codes_github/cdr1as/Figures/cumulative_plots_miR_7_KOm7oe_WT_KOm7oe_KO.pdf",
       plot = c1,
       scale = 1, width = 7, height = 5.5, units = "in", device = cairo_pdf,
       dpi = 300)


res1234 <- rbind(res1.df, res2.df[res2.df$is.target %in% TRUE, ], res3.df[res3.df$is.target %in% TRUE, ], res4.df[res4.df$is.target %in% TRUE, ])

res1234$is.target2 <- ifelse((res1234$exp %in% "KO_vs_WT" & res1234$is.target %in% FALSE), "miR-7_nontargets",
                           ifelse((res1234$exp %in% "KO_vs_WT" & res1234$is.target %in% TRUE), "miR-7_KO_vs_WT", 
                                  ifelse((res1234$exp %in% "WTm7oe_vs_WT" & res1234$is.target %in% TRUE), "miR-7_WTm7oe_vs_WT",
                                         ifelse((res1234$exp %in% "KOm7oe_vs_KO" & res1234$is.target %in% TRUE), "miR-7_KOm7oe_vs_KO",
                                                ifelse((res1234$exp %in% "KOm7oe_vs_WT" & res1234$is.target %in% TRUE), "miR-7_KOm7oe_vs_WT", "NA")))))


g1 <- ggplot(res1234, aes(x = log2FoldChange, color = is.target2)) +
  stat_ecdf() +
  theme_classic() + 
  ylab("Cumulative fraction") +
  xlim(-1, 1) + 
  scale_color_manual(values = c("darkgray", "#FF9900", "purple",  "darkred", "#66CCCC"),
                                breaks = c("miR-7_nontargets", "miR-7_KO_vs_WT", "miR-7_WTm7oe_vs_WT","miR-7_KOm7oe_vs_KO","miR-7_KOm7oe_vs_WT")) + 
  theme(axis.line.y=element_line(size=.5),
        axis.line.x=element_line(size=.5),
        axis.text.x=element_text(size=18,  family="Helvetica", color="black"),
        axis.text.y=element_text(size=18,  family="Helvetica", color="black"),
        axis.title.x=element_text(size=20, family="Helvetica", color="black"),
        axis.title.y=element_text(size=20, family="Helvetica", color="black"),
        legend.title = element_blank(),
        legend.text = element_text(size=12, family="Helvetica", color="black"),
        legend.position = c(0.55, 0.5),
        legend.justification = c(0, 1),
        aspect.ratio=1,
        plot.title=element_text(size=15, face="bold", family="Helvetica", hjust=0.5)) 

ggsave(filename = "/data/rajewsky/projects/cdr1as_ko_snRNA/codes_github/cdr1as/Figures/cumulative_plots_all.pdf",
       plot = g1,
       scale = 1, width = 7, height = 5.5, units = "in", device = cairo_pdf,
       dpi = 300)

g1 <- ggplot(res1234, aes(x = log2FoldChange, color = is.target2)) +
  stat_ecdf() +
  theme_classic() + 
  ylab("Cumulative fraction") +
  xlim(-1, 1) + 
  scale_color_manual(values = c("darkgray", "#FF9900", "purple",  "darkred", "#66CCCC"),
                     breaks = c("miR-7_nontargets", "miR-7_KO_vs_WT", "miR-7_WTm7oe_vs_WT","miR-7_KOm7oe_vs_KO","miR-7_KOm7oe_vs_WT")) + 
  theme(axis.line.y=element_line(size=.5),
        axis.line.x=element_line(size=.5),
        axis.text.x=element_text(size=18,  family="Helvetica", color="black"),
        axis.text.y=element_text(size=18,  family="Helvetica", color="black"),
        axis.title.x=element_text(size=20, family="Helvetica", color="black"),
        axis.title.y=element_text(size=20, family="Helvetica", color="black"),
        legend.title = element_blank(),
        legend.text = element_text(size=12, family="Helvetica", color="black"),
        legend.position = c(0.55, 0.5),
        legend.justification = c(0, 1),
        aspect.ratio=1,
        plot.title=element_text(size=15, face="bold", family="Helvetica", hjust=0.5)) +
        annotate("text", x = -0.64, y = 0.95, label = "1 vs 2: p-value < 2.2e-16", colour = "#003366", face = "bold", size = 4, family="Helvetica") + 
        annotate("text", x = -0.6, y = 0.9, label = "2 vs 3: p-value = 3.459e-13", colour = "#003366", face = "bold", size = 4, family="Helvetica") +
        annotate("text", x = -0.6, y = 0.85, label = "2 vs 4: p-value = 1.112e-06", colour = "#003366", face = "bold", size = 4, family="Helvetica") +
        annotate("text", x = -0.65, y = 0.8, label = "3 vs 4: p-value = 0.1281", colour = "#003366", face = "bold", size = 4, family="Helvetica") 


ggsave(filename = "/data/rajewsky/projects/cdr1as_ko_snRNA/codes_github/cdr1as/Figures/cumulative_plots_all_p_value.pdf",
       plot = g1,
       scale = 1, width = 7, height = 5.5, units = "in", device = cairo_pdf,
       dpi = 300)


