# This file is for DESeq2 of 4 experiments (WT, WT with miR-7 oe, Cdr1as KO Cdr1as KO with miR-7 oe)
# In here, I tested 3 different methods of DESeq2. 1) Batch correction, 2) nested design, 3) both batch and nested design.


library(limma)
library(DESeq2)

DGEm_normalized <- readRDS(file = "/data/rajewsky/home/skim/Neuronal_Activity_Cledi/201911_m7_overexpression/publication_codes_figures/DGE_matrix_normalized_nested.rda")

DGEm <- readRDS(file = "/data/rajewsky/home/skim/Neuronal_Activity_Cledi/201911_m7_overexpression/publication_codes_figures/DGE_matrix.rda")
Counts <- t(as.matrix(DGEm))


 

batch <- c("B", "M", "M", "B", "M", "M", "M", "M", "M", "M", "M", "M", "M", "M")
mys <- colnames(DGEm)
batch2 = factor(c("Rep1", "Rep2", "Rep3",
                 "Rep1", "Rep2", "Rep3",
                 "Rep4", "Rep5", "Rep6", "Rep7",
                 "Rep4", "Rep5", "Rep6", "Rep7"))



# ------------------- #
# functions
# ------------------- #

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




dds3 <- DESeqDataSetFromMatrix(t(Counts), colData = myp, mmCorrected)
keep <- rowSums(counts(dds3)) >= 10
dds3 <- dds3[keep,]

###################
# PCA plot
###################

vsd <- vst(t(Counts[, keep]), blind = F) 
CountsNew <- vsd
#vsd <- vst(dds3, blind = F) 


# original PCA
vsd <- varianceStabilizingTransformation(dds3)
vsd2 <- vsd
mm2 <- model.matrix(~0+genotype+genotype:cond, myp)  

plotPCA(vsd, intgroup = "group")
p1 <- plotPCA(vsd, intgroup = "group", returnData = TRUE)

# ----------------- #
# PCA by exp
# ----------------- #

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
    scale_color_manual(values = c("#FF0066", "#880E4F", "#0099FF", "#0D47A1")) + 
    scale_shape_manual(values=c(19, 18, 17, 15))
  
}

p1$group <- factor(p1$group, levels = c("WT", "WTm7o", "KO", "KOm7o"))
g1 <- pca_plot(pc = p1, group = group, x_pc = "60", y_pc ="11")
g1
ggsave(filename = "/data/rajewsky/projects/cdr1as_ko_snRNA/codes/DE_Analysis_DESeq2_nested_corrected_data/PCA_orig.pdf",
       plot = g1,
       scale = 1, width = 6, height = 5, units = "in", device = cairo_pdf,
       dpi = 300)



# ------------------------- #
# site correction
# ------------------------- #

r <- removeBatchEffect(CountsNew, batch=batch, design=mm2)
assay(vsd2) <- r
plotPCA(vsd2, intgroup = "group")
p1 <- plotPCA(vsd2, intgroup = "group", returnData = TRUE)
p1$group <- factor(p1$group, levels = c("WT", "WTm7o", "KO", "KOm7o"))
g1 <- pca_plot(pc = p1, group = group, x_pc = "34", y_pc ="19")
g1
ggsave(filename = "/data/rajewsky/projects/cdr1as_ko_snRNA/codes/DE_Analysis_DESeq2_nested_corrected_data/PCA_batch_corr.pdf",
       plot = g1,
       scale = 1, width = 6, height = 5, units = "in", device = cairo_pdf,
       dpi = 300)

# ------------------------- #
# nested correction
# ------------------------- #
vsd <- vst(t(Counts[, keep]), blind = F) 
CountsNew <- vsd
vsd <- varianceStabilizingTransformation(dds3)

r <- removeBatchEffect(CountsNew, batch=batch2, design=mm2)
vsd2 <- vsd
assay(vsd2) <- r
plotPCA(vsd2, intgroup = "group")

p1 <- plotPCA(vsd2, intgroup = "group", returnData = TRUE)
p1$group <- factor(p1$group, levels = c("WT", "WTm7o", "KO", "KOm7o"))
g1 <- pca_plot(pc = p1, group = group, x_pc = "55", y_pc ="16")
g1
ggsave(filename = "/data/rajewsky/projects/cdr1as_ko_snRNA/codes/DE_Analysis_DESeq2_nested_corrected_data/PCA_nested_corr.pdf",
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
g1 <- pca_plot(pc = p1, group = group, x_pc = "45", y_pc ="21")
g1
ggsave(filename = "/data/rajewsky/projects/cdr1as_ko_snRNA/codes/DE_Analysis_DESeq2_nested_corrected_data/PCA_batch_nested_corr.pdf",
       plot = g1,
       scale = 1, width = 6, height = 5, units = "in", device = cairo_pdf,
       dpi = 300)


# save batch & nested corrected matrix for WGCNA and RF
r.matrix <- r %>% as.matrix()
r.df <- r %>% as.data.frame()
saveRDS(r.df, file = "/data/rajewsky/projects/cdr1as_ko_snRNA/codes/DE_Analysis_DESeq2_nested_corrected_data/DE_Analysis_DESeq2_nested_corrected_data_DGEm.rda")


# ------------------------- #
# batch & nested correction DESeq2
# ------------------------- #

mm <- model.matrix(~batchBuch+genotype+genotype:batch+genotype:cond, myp)
wc <- which(colnames(mm)=="genotypeWT:batch4" | colnames(mm)=="genotypeWT:batch3")
mmCorrected <- mm[,-wc]

decisions <- matrix(0, ncol(Counts), length(coef.of.interest))
colnames(decisions) <- coef.of.interest
rownames(decisions) <- colnames(Counts)

dds3 <- DESeqDataSetFromMatrix(t(Counts), colData = myp, mmCorrected)
keep <- rowSums(counts(dds3)) >= 10
dds3 <- dds3[keep,]
dds3 <- DESeq(dds3)
resultsNames(dds3)

DGEm_normalized <- counts(dds3, normalized=T) %>% as.data.frame
DGEm_normalized_WT <- DGEm_normalized[,c(1:3)] %>% rowMeans


# KO vs WT
res <- results(dds3, name=c("genotypeWT"), test="Wald")
res.df <- as.data.frame(res)
res.df$gene <- rownames(res.df)
res.df$WT <- DGEm_normalized_WT
res.df$log2FoldChange <- res.df$log2FoldChange * -1
saveRDS(res.df, file = "/data/rajewsky/cdr1as_ko_snRNA/codes/DE_Analysis_DESeq2_nested_corrected_data/DE_Analysis_DESeq2_nested_corrected_res1.rda")


# WTm7oe vs WT
res2 <- results(dds3, list( c("genotypeWT.condOE","genotypeWT") ))
res2.df <- as.data.frame(res2)
res2.df$gene <- rownames(res2.df)
res2.df$WT <- DGEm_normalized_WT
saveRDS(res2.df, file = "/data/rajewsky/cdr1as_ko_snRNA/codes/DE_Analysis_DESeq2_nested_corrected_data/DE_Analysis_DESeq2_nested_corrected_res2.rda")

# KOm7oe vs KO
res3 <- results(dds3, name=c("genotypeKO.condOE"))
res3.df <- as.data.frame(res3)
res3.df$gene <- rownames(res3.df)
res3.df$WT <- DGEm_normalized_WT
saveRDS(res3.df, file = "/data/rajewsky/cdr1as_ko_snRNA/codes/DE_Analysis_DESeq2_nested_corrected_data/DE_Analysis_DESeq2_nested_corrected_res3.rda")

# KOm7oe vs WT
res4 <- results(dds3, list( c("genotypeKO.condOE","genotypeWT") ))
res4.df <- as.data.frame(res4)
res4.df$gene <- rownames(res4.df)
res4.df$WT <- DGEm_normalized_WT
saveRDS(res4.df, file = "/data/rajewsky/cdr1as_ko_snRNA/codes/DE_Analysis_DESeq2_nested_corrected_data/DE_Analysis_DESeq2_nested_corrected_res4.rda")


# ------------------------- #
# batch correction DESeq2
# ------------------------- #

sample_info <- data.frame(condition = c("WT","WT","WT",
                                        "WTm7ove","WTm7ove","WTm7ove",
                                        "KO","KO","KO","KO",
                                        "KOm7ove","KOm7ove","KOm7ove","KOm7ove"), 
                          batch = c("Buch", "Mitte", "Mitte",
                                    "Buch", "Mitte","Mitte",
                                    "Mitte", "Mitte", "Mitte", "Mitte",
                                    "Mitte", "Mitte", "Mitte", "Mitte"), row.names = names(DGEm))

dds <- DESeqDataSetFromMatrix(DGEm, colData = sample_info, design = ~ condition + batch)

keep <- rowSums(counts(dds)) >= 10
table(keep)
dds <- dds[keep,]
dds$condition <- relevel(dds$condition, ref = "WT")
dds <- DESeq(dds)
resultsNames(dds)

DGEm_normalized <- counts(dds, normalized=T) %>% as.data.frame
DGEm_normalized_WT <- DGEm_normalized[,c(1:3)] %>% rowMeans

# KO vs WT
res <- results(dds, name=c("condition_KO_vs_WT"), test="Wald")
res.df <- as.data.frame(res)
res.df$gene <- rownames(res.df)
res.df$WT <- DGEm_normalized_WT
saveRDS(res.df, file = "/data/rajewsky/home/skim/Neuronal_Activity_Cledi/201911_m7_overexpression/publication_codes_figures/DE_Analysis_DESeq2_res1_KO_WT.rda")

# WTm7oe vs WT
res <- results(dds, name=c("condition_WTm7ove_vs_WT"), test="Wald")
res.df <- as.data.frame(res)
res.df$gene <- rownames(res.df)
res.df$WT <- DGEm_normalized_WT
saveRDS(res.df, file = "/data/rajewsky/home/skim/Neuronal_Activity_Cledi/201911_m7_overexpression/publication_codes_figures/DE_Analysis_DESeq2_res2_WTm7ove_WT.rda")

# KOm7oe vs KO
res <- results(dds, contrast=c("condition","KOm7ove", "KO"), test="Wald")
res.df <- as.data.frame(res)
res.df$gene <- rownames(res.df)
res.df$WT <- DGEm_normalized_WT
saveRDS(res.df, file = "/data/rajewsky/home/skim/Neuronal_Activity_Cledi/201911_m7_overexpression/publication_codes_figures/DE_Analysis_DESeq2_res3_KOm7ove_KO.rda")

# KOm7oe vs WT
res <- results(dds, contrast=c("condition","WTm7ove", "KO"), test="Wald")
res.df <- as.data.frame(res)
res.df$gene <- rownames(res.df)
res.df$WT <- DGEm_normalized_WT
saveRDS(res.df, file = "/data/rajewsky/home/skim/Neuronal_Activity_Cledi/201911_m7_overexpression/publication_codes_figures/DE_Analysis_DESeq2_res5_WTm7ove_KO.rda")



# ------------------------- #
# nested correction DESeq2
# ------------------------- #
sample_info <- DataFrame(id = names(DGEm),
                         genotype = factor(c("WT","WT","WT",
                                             "WT","WT","WT",
                                             "KO","KO","KO","KO",
                                             "KO","KO","KO","KO")), 
                         batch = factor(c("Rep1", "Rep2", "Rep3",
                                          "Rep1", "Rep2", "Rep3",
                                          "Rep4", "Rep5", "Rep6", "Rep7",
                                          "Rep4", "Rep5", "Rep6", "Rep7")), 
                         condition = factor(c("Ctrl","Ctrl","Ctrl",
                                              "m7ove","m7ove","m7ove",
                                              "Ctrl","Ctrl","Ctrl","Ctrl",
                                              "m7ove","m7ove","m7ove","m7ove")), row.names = names(DGEm))



sample_info$batch.n <- c(rep(2:4), rep(2:4), rep(1:4), rep(1:4)) %>% as.character()
sample_info$batch <- NULL
sample_info.df <- sample_info %>% as.data.frame()
sample_info$batch.n <- sample_info$batch.n %>% as.character()

m1 <- model.matrix(~0+ genotype + genotype:batch.n + genotype:condition, sample_info.df)

wc <- which(colnames(m1)=="genotypeWT:batch.n3") # WT has one sample less
m2 <- m1[,-wc]

v <- DESeqDataSetFromMatrix(DGEm, sample_info.df, m2)
keep <- rowSums(counts(v)) >= 10
v <- v[keep,]

# run DESeq2
ds <- DESeq(v)

#design(ds)
resultsNames(ds)

# calculated WT counts
DGEm_normalized <- counts(ds, normalized=T) %>% as.data.frame
DGEm_normalized_WT <- DGEm_normalized[,c(1:3)] %>% rowMeans

# KO vs WT
result1 <- results(ds, contrast=list("genotypeKO","genotypeWT"))
result1.df <- as.data.frame(result1)
result1.df$gene <- rownames(result1.df)
result1.df$WT <- DGEm_normalized_WT
saveRDS(result1.df, file = "/data/rajewsky/home/skim/Neuronal_Activity_Cledi/201911_m7_overexpression/publication_codes_figures/DE_Analysis_DESeq2_nested_result1.rda")

# WTm7oe vs WT
result2 <- results(ds, name="genotypeWT.conditionm7ove")
result2.df <- as.data.frame(result2)
result2.df$gene <- rownames(result2.df)
result2.df$WT <- DGEm_normalized_WT
saveRDS(result2.df, file = "/data/rajewsky/home/skim/Neuronal_Activity_Cledi/201911_m7_overexpression/publication_codes_figures/DE_Analysis_DESeq2_nested_result2.rda")

# KOm7oe vs KO
result3 <- results(ds, name="genotypeKO.conditionm7ove")
result3.df <- as.data.frame(result3)
result3.df$gene <- rownames(result3.df)
result3.df$WT <- DGEm_normalized_WT
saveRDS(result3.df, file = "/data/rajewsky/home/skim/Neuronal_Activity_Cledi/201911_m7_overexpression/publication_codes_figures/DE_Analysis_DESeq2_nested_result3.rda")

# KOm7oe vs WT
result4 <- results(ds, contrast=list(c("genotypeKO","genotypeKO.conditionm7ove"), c("genotypeWT")))
result4.df <- as.data.frame(result4)
result4.df$gene <- rownames(result4.df)
result4.df$WT <- DGEm_normalized_WT
saveRDS(result4.df, file = "/data/rajewsky/home/skim/Neuronal_Activity_Cledi/201911_m7_overexpression/publication_codes_figures/DE_Analysis_DESeq2_nested_result4.rda")

