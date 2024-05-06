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

SS_DGEm <- readRDS(file = "/data/rajewsky/home/skim/Neuronal_Activity_Cledi/202207_m7_oex_5ssKO/publication_code/SSKO_Matrix.rda")



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

sample_info <- data.frame(exp = group, genotype=genotype, cond=cond, individual=individual.nested)


m1 <- model.matrix(~0+ genotype + genotype:individual + genotype:cond, sample_info)
v <- DESeqDataSetFromMatrix(SS_DGEm, sample_info, m1)
keep <- rowSums(counts(v)) >= 10
v <- v[keep,]

# run DESeq2
ds <- DESeq(v)
#resultsNames(ds)

# calculated WT counts
DGEm_normalized <- counts(ds, normalized=T) %>% as.data.frame
DGEm_normalized_WT <- DGEm_normalized[,c(1:3)] %>% rowMeans

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


# ggsave(filename = "/data/rajewsky/projects/cdr1as_ko_snRNA/codes_github/cdr1as/Figures/WTJ_PCA_orig.pdf",
#        plot = g1,
#        scale = 1, width = 6, height = 5, units = "in", device = cairo_pdf,
#        dpi = 300)



# nested corrected PCA

vsd <- vst(t(Counts[, keep]), blind = F) 
CountsNew <- vsd
vsd <- varianceStabilizingTransformation(ds)
r <- removeBatchEffect(CountsNew, batch=individual.nested, design=m1)
vsd2 <- vsd
assay(vsd2) <- r

vsd2$group <- group
p2 <- plotPCA(vsd2, intgroup = "group", returnData = TRUE)
p2$group <- factor(p2$group, levels = c("WTJ", "WTJm7oe", "SSKO", "SSKOm7oe"))

g2 <- pca_plot(pc = p2, group = group, x_pc = "82", y_pc ="8")
g2

# ggsave(filename = "/data/rajewsky/projects/cdr1as_ko_snRNA/codes_github/cdr1as/Figures/WTJ_PCA_nested.pdf",
#        plot = g2,
#        scale = 1, width = 6, height = 5, units = "in", device = cairo_pdf,
#        dpi = 300)


# ---------------------------------- #
# DESeq2 results
# ---------------------------------- #
# KO vs WT
result1 <- results(ds, contrast=list("genotypeSSKO","genotypeWTJ"))
result1.df <- as.data.frame(result1)
result1.df$gene <- rownames(result1.df)
result1.df$WT <- DGEm_normalized_WT
# check reference level
df <- mcols(result1) %>% as.data.frame()
saveRDS(result1.df, file = "/data/rajewsky/home/skim/Neuronal_Activity_Cledi/202207_m7_oex_5ssKO/publication_code/SSKO_DESeq_nested_analysis/SSKO_DESeq_res1_SSKO_WT.rda")

# WTm7oe vs WT
result2 <- results(ds, name="genotypeWTJ.condoe")
result2.df <- as.data.frame(result2)
result2.df$gene <- rownames(result2.df)
result2.df$WT <- DGEm_normalized_WT
df2 <- mcols(result2) %>% as.data.frame()
saveRDS(result2.df, file = "/data/rajewsky/home/skim/Neuronal_Activity_Cledi/202207_m7_oex_5ssKO/publication_code/SSKO_DESeq_nested_analysis/SSKO_DESeq_res2_WTJm7oe_WTJ.rda")


# KOm7oe vs KO
result3 <- results(ds, name="genotypeSSKO.condoe")
result3.df <- as.data.frame(result3)
result3.df$gene <- rownames(result3.df)
result3.df$WT <- DGEm_normalized_WT
saveRDS(result3.df, file = "/data/rajewsky/home/skim/Neuronal_Activity_Cledi/202207_m7_oex_5ssKO/publication_code/SSKO_DESeq_nested_analysis/SSKO_DESeq_res3_SSKOm7oe_SSKO.rda")

# KOm7oe vs WT
result4 <- results(ds, contrast=list(c("genotypeSSKO","genotypeSSKO.condoe"), c("genotypeWTJ")))
result4.df <- as.data.frame(result4)
result4.df$gene <- rownames(result4.df)
result4.df$WT <- DGEm_normalized_WT
saveRDS(result4.df, file = "/data/rajewsky/home/skim/Neuronal_Activity_Cledi/202207_m7_oex_5ssKO/publication_code/SSKO_DESeq_nested_analysis/SSKO_DESeq_res4_SSKOm7oe_WTJ.rda")

