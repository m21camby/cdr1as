
library(VennDiagram)
library(ggvenn)


#########################
# DE data load
#########################
res1.df <- readRDS(file = "/data/rajewsky/projects/cdr1as_ko_snRNA/codes_github/cdr1as/Results/DESeq2_WTN_batch_corrected_res1.rda")
res2.df <- readRDS(file = "/data/rajewsky/projects/cdr1as_ko_snRNA/codes_github/cdr1as/Results/DESeq2_WTN_nested_corrected_res2.rda")
res3.df <- readRDS(file = "/data/rajewsky/projects/cdr1as_ko_snRNA/codes_github/cdr1as/Results/DESeq2_WTN_nested_corrected_res3.rda")
res4.df <- readRDS(file = "/data/rajewsky/projects/cdr1as_ko_snRNA/codes_github/cdr1as/Results/DESeq2_WTN_batch_corrected_res4.rda")

# up-genes
set1 <- res1.df[which(res1.df$log2FoldChange > 0.5 &  res1.df$padj < 0.05), ]$gene
set2 <- res2.df[which(res2.df$log2FoldChange > 0.5 &  res2.df$padj < 0.05), ]$gene
set3 <- res3.df[which(res3.df$log2FoldChange > 0.5 &  res3.df$padj < 0.05), ]$gene
set4 <- res4.df[which(res4.df$log2FoldChange > 0.5 &  res4.df$padj < 0.05), ]$gene 

genes1 <- c(set1, set2, set3, set4) %>% unique()

# down-genes
set5 <- res1.df[which(res1.df$log2FoldChange < -0.5 &  res1.df$padj < 0.05), ]$gene
set6 <- res2.df[which(res2.df$log2FoldChange < -0.5 &  res2.df$padj < 0.05), ]$gene
set7 <- res3.df[which(res3.df$log2FoldChange < -0.5 &  res3.df$padj < 0.05), ]$gene
set8 <- res4.df[which(res4.df$log2FoldChange < -0.5 &  res4.df$padj < 0.05), ]$gene 

genes2 <- c(set5, set6, set7, set8) %>% unique()

all_genes <- c(genes1, genes2) %>% unique()
# genes exist in both up- and down-
intersect(genes1, genes2)
# "Tgm2"

# ---------------------- #
# up-regulated
# ---------------------- #

x <- list(
  KO_vs_WT = res1.df[which(res1.df$log2FoldChange > 0.5 &  res1.df$padj < 0.05), ]$gene,
  WTm7oe_vs_WT = res2.df[which(res2.df$log2FoldChange > 0.5 &  res2.df$padj < 0.05), ]$gene,
  KOm7oe_vs_KO = res3.df[which(res3.df$log2FoldChange > 0.5 &  res3.df$padj < 0.05), ]$gene,
  KOm7oe_vs_WT = res4.df[which(res4.df$log2FoldChange > 0.5 &  res4.df$padj < 0.05), ]$gene
)

g1 <- ggvenn(
  x, 
  stroke_linetype = 1, 
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
  stroke_size = 1, set_name_size = 9, text_size = 15
) + 
  theme(plot.margin=unit(c(0,0,0,0),"cm"))

g1


ggsave(filename = "/data/rajewsky/projects/cdr1as_ko_snRNA/codes_github/cdr1as/Figures/venn_diagram_up_all.pdf",
       plot = g1,
       scale = 1, width = 18, height = 12, units = "in", device = cairo_pdf,
       dpi = 300)


g1 <- ggvenn(
  x, 
  stroke_linetype = 1, 
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
  stroke_size = 1, set_name_size = 9, text_size = 15, show_percentage = FALSE
) + 
  theme(plot.margin=unit(c(0,0,0,0),"cm"))

g1

ggsave(filename = "/data/rajewsky/projects/cdr1as_ko_snRNA/codes_github/cdr1as/Figures/venn_diagram_up_all_no_per.pdf",
       plot = g1,
       scale = 1, width = 9, height = 6, units = "in", device = cairo_pdf,
       dpi = 300)

# ---------------------- #
# down-regulated
# ---------------------- #

x <- list(
  KO_vs_WT = res1.df[which(res1.df$log2FoldChange < -0.5 &  res1.df$padj < 0.05), ]$gene,
  WTm7oe_vs_WT = res2.df[which(res2.df$log2FoldChange < -0.5 &  res2.df$padj < 0.05), ]$gene,
  KOm7oe_vs_KO = res3.df[which(res3.df$log2FoldChange < -0.5 &  res3.df$padj < 0.05), ]$gene,
  KOm7oe_vs_WT = res4.df[which(res4.df$log2FoldChange < -0.5 &  res4.df$padj < 0.05), ]$gene
)


g1 <- ggvenn(
  x, 
  stroke_linetype = 1, 
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
  stroke_size = 1, set_name_size = 9, text_size = 15
) + 
  theme(plot.margin=unit(c(0,0,0,0),"cm"))
g1

ggsave(filename = "/data/rajewsky/projects/cdr1as_ko_snRNA/codes_github/cdr1as/Figures/venn_diagram_down_all.pdf",
       plot = g1,
       scale = 1, width = 18, height = 12, units = "in", device = cairo_pdf,
       dpi = 300)

g1 <- ggvenn(
  x, 
  stroke_linetype = 1, 
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
  stroke_size = 1, set_name_size = 9, text_size = 15, show_percentage = FALSE
) + 
  theme(plot.margin=unit(c(0,0,0,0),"cm"))
g1

ggsave(filename = "/data/rajewsky/projects/cdr1as_ko_snRNA/codes_github/cdr1as/Figures/venn_diagram_down_all_no_per.pdf",
       plot = g1,
       scale = 1, width = 9, height = 6, units = "in", device = cairo_pdf,
       dpi = 300)



# co-up specific KOm7oe vs KO & WTm7oe vs WT 
list1 <- Reduce(intersect, list(set3, set2))
#list2 <- Reduce(intersect, list(set4, set1))
list <- list1[!list1 %in% c(set4, set1)]
list
# "Susd2"

# co-down WTm7oe vs WT & KOm7oe vs WT
list1 <- Reduce(intersect, list(set8, set6, set5))
list <- list1[!list1 %in% c(set7)]
list

# co-down specific KOm7oe vs KO & WTm7oe vs WT
list1 <- Reduce(intersect, list(set7, set6))
#list2 <- Reduce(intersect, list(set4, set1))
list <- list1[!list1 %in% c(set8, set5)]
list

