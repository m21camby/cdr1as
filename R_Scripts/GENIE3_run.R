

# Two runs of GENIE3

#----------#
# 1st run
#----------#

# "targets" are the phenotype-related genes
# "regulators" are all other down- and upregulated genes by miR-7 overexpression

#----------#
# 2nd run
#----------#

# "targets" are potential regulators from run A that were indeed selected as true regulators and are not known miR-7 targets
# regulators  are known miR-7 targets that are downregulated by miR-7 overexpression.

#----------#
# up-date from Elisabeth
#----------#
# GENIE3 can assign targets & regulons
# Snca is miR-7 target and remove from phenotype
# move node for hierarchy
# keep 3 ~ 4 genes per node
# different color by DE test
# global run

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
library(pcaExplorer)
library(topGO)
library(paletteer)
library(scales)
library(stringr)
library(GENIE3)
library(readxl)

res1.df <- readRDS(file = "/data/rajewsky/home/skim/Neuronal_Activity_Cledi/201911_m7_overexpression/publication_codes_figures/DE_Analysis_DESeq2_res1_KO_WT.rda")

res2.df <- readRDS(file = "/data/rajewsky/home/skim/Neuronal_Activity_Cledi/201911_m7_overexpression/publication_codes_figures/DE_Analysis_DESeq2_res2_WTm7ove_WT.rda")

res3.df <- readRDS(file = "/data/rajewsky/home/skim/Neuronal_Activity_Cledi/201911_m7_overexpression/publication_codes_figures/DE_Analysis_DESeq2_res3_KOm7ove_KO.rda")

res4.df <- readRDS(file = "/data/rajewsky/home/skim/Neuronal_Activity_Cledi/201911_m7_overexpression/publication_codes_figures/DE_Analysis_DESeq2_res4_KOm7ove_WT.rda")

mir7target.df <- read.table("/data/rajewsky/home/skim/Neuronal_Activity_Cledi/201911_m7_overexpression/R_Scripts/external_files/TargetScan7.2__miR-7-5p.predicted_targets_lists_without_header.txt")

DGEm_normalized <- readRDS(file = "/data/rajewsky/home/skim/Neuronal_Activity_Cledi/201911_m7_overexpression/publication_codes_figures/DGE_matrix_normalized.rda")

mat.df <- readRDS(file = "/data/rajewsky/home/skim/Neuronal_Activity_Cledi/201911_m7_overexpression/publication_codes_figures/DE_Analysis_DESeq2_batch_corrected_matrix.rda")

cledi_pheno_genes <- read_excel("/data/rajewsky/home/skim/Neuronal_Activity_Cledi/201911_m7_overexpression/R_Scripts/external_files_2nd/biased_list_of_phenotype_related_genes.xlsx", col_names = FALSE)

##################################
# 1st: only miR-7 DE genes
##################################

# up-genes
set1 <- res1.df[which(res1.df$log2FoldChange > 0.75 &  res1.df$padj < 0.05), ]$gene
set2 <- res2.df[which(res2.df$log2FoldChange > 0.75 &  res2.df$padj < 0.05), ]$gene
set3 <- res3.df[which(res3.df$log2FoldChange > 0.75 &  res3.df$padj < 0.05), ]$gene
set4 <- res4.df[which(res4.df$log2FoldChange > 0.75 &  res4.df$padj < 0.05), ]$gene 

genes1 <- c(set2, set3, set4) %>% unique()

# down-genes
set1 <- res1.df[which(res1.df$log2FoldChange < -0.75 &  res1.df$padj < 0.05), ]$gene
set2 <- res2.df[which(res2.df$log2FoldChange < -0.75 &  res2.df$padj < 0.05), ]$gene
set3 <- res3.df[which(res3.df$log2FoldChange < -0.75 &  res3.df$padj < 0.05), ]$gene
set4 <- res4.df[which(res4.df$log2FoldChange < -0.75 &  res4.df$padj < 0.05), ]$gene 

genes2 <- c(set2, set3, set4) %>% unique()

# all DE genes
genes <- c(genes1, genes2)

# phenotypic genes
pheno_genes <- c("Zdhhc9", "Parp1", "Myrip", "Wipf2", "Phactr1", "Basp1", "Pfn2", "Snca", "Cplx1")

##################
# 1st Run GENIE3
##################

pheno_genes <- pheno_genes[!pheno_genes %in% mir7target.df$V1] 
genes1 <- c(genes, pheno_genes) 
#genes1 <- genes1[c(!genes1 %in% mir7target.df$V1)]
genes2 <- genes[c(!genes %in% pheno_genes & !genes %in% mir7target.df$V1)]


matrix_test <- DGEm_normalized[rownames(DGEm_normalized) %in% genes1, ] %>% as.matrix


weightMat_test <- GENIE3(matrix_test, regulators = genes2, targets = pheno_genes)
linkList_test <- getLinkList(weightMat_test)

saveRDS(linkList_test, file = "/data/rajewsky/home/skim/Neuronal_Activity_Cledi/201911_m7_overexpression/publication_codes_figures/GENIE3_run_1st_run_Elisabeth.rda")

# # subset by DE genes
# matrix_test <- DGEm_normalized[rownames(DGEm_normalized) %in% c(genes, pheno_genes), ] %>% as.matrix
# weightMat_test <- GENIE3(matrix_test)
# linkList_test <- getLinkList(weightMat_test)
# 
# # Filter by weight
# linkList_test2 <- linkList_test %>% subset(weight > 0.01 & targetGene %in% pheno_genes)
# saveRDS(linkList_test2, file = "/data/rajewsky/home/skim/Neuronal_Activity_Cledi/201911_m7_overexpression/publication_codes_figures/GENIE3_run_1st_run.rda")

##################
# 2nd Run GENIE3
##################

genes2 <- c(set2, set3, set4) %>% unique()
regul_genes <- genes2[genes2 %in% mir7target.df$V1]

target_genes <- linkList_test$regulatoryGene[!linkList_test$regulatoryGene %in% mir7target.df$V1] %>% as.character() %>% unique()

# subset by DE genes
matrix_test <- DGEm_normalized[rownames(DGEm_normalized) %in% c(regul_genes, target_genes), ] %>% as.matrix

weightMat_test <- GENIE3(matrix_test, regulators = regul_genes, targets = target_genes)
linkList_test3 <- getLinkList(weightMat_test)

saveRDS(linkList_test3, file = "/data/rajewsky/home/skim/Neuronal_Activity_Cledi/201911_m7_overexpression/publication_codes_figures/GENIE3_run_2nd_run_Elisabeth.rda")

# weightMat_test <- GENIE3(matrix_test)
# linkList_test3 <- getLinkList(weightMat_test)
# 
# linkList_test4 <- linkList_test3 %>% 
#   dplyr::filter(weight > 0.01) %>% 
#   dplyr::filter(targetGene %in% target_genes & regulatoryGene %in% regul_genes) %>% 
#   dplyr::filter(!targetGene %in% pheno_genes & !regulatoryGene %in% pheno_genes) %>%
#   dplyr::filter(targetGene %in% linkList_test2$regulatoryGene) %>% 
#   dplyr::filter(!regulatoryGene %in% linkList_test2$regulatoryGene)
# 
# saveRDS(linkList_test4, file = "/data/rajewsky/home/skim/Neuronal_Activity_Cledi/201911_m7_overexpression/publication_codes_figures/GENIE3_run_2nd_run.rda")



##################################
# 2nd: Including Cdr1as DE genes
##################################

# up-genes
set1 <- res1.df[which(res1.df$log2FoldChange > 0.75 &  res1.df$padj < 0.05), ]$gene
set2 <- res2.df[which(res2.df$log2FoldChange > 0.75 &  res2.df$padj < 0.05), ]$gene
set3 <- res3.df[which(res3.df$log2FoldChange > 0.75 &  res3.df$padj < 0.05), ]$gene
set4 <- res4.df[which(res4.df$log2FoldChange > 0.75 &  res4.df$padj < 0.05), ]$gene 

genes1 <- c(set1, set2, set3, set4) %>% unique()

# down-genes
set1 <- res1.df[which(res1.df$log2FoldChange < -0.75 &  res1.df$padj < 0.05), ]$gene
set2 <- res2.df[which(res2.df$log2FoldChange < -0.75 &  res2.df$padj < 0.05), ]$gene
set3 <- res3.df[which(res3.df$log2FoldChange < -0.75 &  res3.df$padj < 0.05), ]$gene
set4 <- res4.df[which(res4.df$log2FoldChange < -0.75 &  res4.df$padj < 0.05), ]$gene 

genes2 <- c(set1, set2, set3, set4) %>% unique()

# all DE genes
genes <- c(genes1, genes2) %>% unique()

# phenotypic genes
pheno_genes <- c("Zdhhc9", "Parp1", "Myrip", "Wipf2", "Phactr1", "Basp1", "Pfn2", "Snca", "Cplx1")


##################
# 1st Run GENIE3 for all
##################

pheno_genes <- pheno_genes[!pheno_genes %in% mir7target.df$V1] 
genes1 <- c(genes, pheno_genes) 
#genes1 <- genes1[c(!genes1 %in% mir7target.df$V1)]
genes2 <- genes[c(!genes %in% pheno_genes & !genes %in% mir7target.df$V1)]

matrix_test <- DGEm_normalized[rownames(DGEm_normalized) %in% genes1, ] %>% as.matrix

weightMat_test <- GENIE3(matrix_test, regulators = genes2, targets = pheno_genes)
linkList_test <- getLinkList(weightMat_test)

saveRDS(linkList_test, file = "/data/rajewsky/home/skim/Neuronal_Activity_Cledi/201911_m7_overexpression/publication_codes_figures/GENIE3_run_all_1st_run_Elisabeth.rda")

# matrix_test <- DGEm_normalized[rownames(DGEm_normalized) %in% c(genes, pheno_genes), ] %>% as.matrix
# 
# weightMat_test <- GENIE3(matrix_test)
# linkList_test <- getLinkList(weightMat_test)
# 
# # Filter by weight
# linkList_test2 <- linkList_test %>% subset(weight > 0.01 & targetGene %in% pheno_genes)
# 
# saveRDS(linkList_test2, file = "/data/rajewsky/home/skim/Neuronal_Activity_Cledi/201911_m7_overexpression/publication_codes_figures/GENIE3_run_all_1st_run.rda")

##################
# 2nd Run GENIE3 for all
##################

genes2 <- c(set2, set3, set4) %>% unique()
regul_genes <- genes2[genes2 %in% mir7target.df$V1]

target_genes <- linkList_test$regulatoryGene[!linkList_test$regulatoryGene %in% mir7target.df$V1] %>% as.character() %>% unique()

# subset by DE genes
matrix_test <- DGEm_normalized[rownames(DGEm_normalized) %in% c(regul_genes, target_genes), ] %>% as.matrix

weightMat_test <- GENIE3(matrix_test, regulators = regul_genes, targets = target_genes)
linkList_test3 <- getLinkList(weightMat_test)

saveRDS(linkList_test3, file = "/data/rajewsky/home/skim/Neuronal_Activity_Cledi/201911_m7_overexpression/publication_codes_figures/GENIE3_run_all_2nd_run_Elisabeth.rda")

# weightMat_test <- GENIE3(matrix_test)
# linkList_test3 <- getLinkList(weightMat_test)
# 
# linkList_test4 <- linkList_test3 %>% 
#   dplyr::filter(weight > 0.01) %>% 
#   dplyr::filter(targetGene %in% target_genes & regulatoryGene %in% regul_genes) %>% 
#   dplyr::filter(!targetGene %in% pheno_genes & !regulatoryGene %in% pheno_genes) %>%
#   dplyr::filter(targetGene %in% linkList_test2$regulatoryGene) %>% 
#   dplyr::filter(!regulatoryGene %in% linkList_test2$regulatoryGene)
# 
# saveRDS(linkList_test4, file = "/data/rajewsky/home/skim/Neuronal_Activity_Cledi/201911_m7_overexpression/publication_codes_figures/GENIE3_run_all_2nd_run.rda")



################################################
# 3rd: Cledi phenotype genes and miR-7 DE genes
################################################

# up-genes
set1 <- res1.df[which(res1.df$log2FoldChange > 0.75 &  res1.df$padj < 0.05), ]$gene
set2 <- res2.df[which(res2.df$log2FoldChange > 0.75 &  res2.df$padj < 0.05), ]$gene
set3 <- res3.df[which(res3.df$log2FoldChange > 0.75 &  res3.df$padj < 0.05), ]$gene
set4 <- res4.df[which(res4.df$log2FoldChange > 0.75 &  res4.df$padj < 0.05), ]$gene 

genes1 <- c(set2, set3, set4) %>% unique()

# down-genes
set1 <- res1.df[which(res1.df$log2FoldChange < -0.75 &  res1.df$padj < 0.05), ]$gene
set2 <- res2.df[which(res2.df$log2FoldChange < -0.75 &  res2.df$padj < 0.05), ]$gene
set3 <- res3.df[which(res3.df$log2FoldChange < -0.75 &  res3.df$padj < 0.05), ]$gene
set4 <- res4.df[which(res4.df$log2FoldChange < -0.75 &  res4.df$padj < 0.05), ]$gene 

genes2 <- c(set2, set3, set4) %>% unique()

# all DE genes
genes <- c(genes1, genes2)

# phenotypic genes

pheno_genes <- cledi_pheno_genes$...1

##################
# 1st Run GENIE3
##################

pheno_genes <- pheno_genes[!pheno_genes %in% mir7target.df$V1] 
genes1 <- c(genes, pheno_genes) 
#genes1 <- genes1[c(!genes1 %in% mir7target.df$V1)]
genes2 <- genes[c(!genes %in% pheno_genes & !genes %in% mir7target.df$V1)]

matrix_test <- DGEm_normalized[rownames(DGEm_normalized) %in% genes1, ] %>% as.matrix

weightMat_test <- GENIE3(matrix_test, regulators = genes2, targets = pheno_genes)
linkList_test <- getLinkList(weightMat_test)

saveRDS(linkList_test, file = "/data/rajewsky/home/skim/Neuronal_Activity_Cledi/201911_m7_overexpression/publication_codes_figures/GENIE3_run_cledi_pheno_1st_run_Elisabeth.rda")

# # subset by DE genes
# matrix_test <- DGEm_normalized[rownames(DGEm_normalized) %in% c(genes, pheno_genes), ] %>% as.matrix
# 
# weightMat_test <- GENIE3(matrix_test)
# linkList_test <- getLinkList(weightMat_test)
# 
# # Filter by weight
# linkList_test2 <- linkList_test %>% subset(weight > 0.01 & targetGene %in% pheno_genes)
# 
# saveRDS(linkList_test2, file = "/data/rajewsky/home/skim/Neuronal_Activity_Cledi/201911_m7_overexpression/publication_codes_figures/GENIE3_run_cledi_pheno_1st_run.rda")


##################
# 2nd Run GENIE3
##################
genes2 <- c(set2, set3, set4) %>% unique()
regul_genes <- genes2[genes2 %in% mir7target.df$V1]

target_genes <- linkList_test$regulatoryGene[!linkList_test$regulatoryGene %in% mir7target.df$V1] %>% as.character() %>% unique()

# subset by DE genes
matrix_test <- DGEm_normalized[rownames(DGEm_normalized) %in% c(regul_genes, target_genes), ] %>% as.matrix

weightMat_test <- GENIE3(matrix_test, regulators = regul_genes, targets = target_genes)
linkList_test3 <- getLinkList(weightMat_test)

saveRDS(linkList_test3, file = "/data/rajewsky/home/skim/Neuronal_Activity_Cledi/201911_m7_overexpression/publication_codes_figures/GENIE3_run_cledi_pheno_2nd_run_Elisabeth.rda")



# weightMat_test <- GENIE3(matrix_test)
# linkList_test3 <- getLinkList(weightMat_test)
# 
# linkList_test4 <- linkList_test3 %>% 
#   dplyr::filter(weight > 0.01) %>% 
#   dplyr::filter(targetGene %in% target_genes & regulatoryGene %in% regul_genes) %>% 
#   dplyr::filter(!targetGene %in% pheno_genes & !regulatoryGene %in% pheno_genes) %>%
#   dplyr::filter(targetGene %in% linkList_test2$regulatoryGene) %>% 
#   dplyr::filter(!regulatoryGene %in% linkList_test2$regulatoryGene)
# 
# saveRDS(linkList_test4, file = "/data/rajewsky/home/skim/Neuronal_Activity_Cledi/201911_m7_overexpression/publication_codes_figures/GENIE3_run_cledi_pheno_2nd_run.rda")

##################################
# 4th: Cledi phenotype genes and Including Cdr1as DE genes
##################################

# up-genes
set1 <- res1.df[which(res1.df$log2FoldChange > 0.75 &  res1.df$padj < 0.05), ]$gene
set2 <- res2.df[which(res2.df$log2FoldChange > 0.75 &  res2.df$padj < 0.05), ]$gene
set3 <- res3.df[which(res3.df$log2FoldChange > 0.75 &  res3.df$padj < 0.05), ]$gene
set4 <- res4.df[which(res4.df$log2FoldChange > 0.75 &  res4.df$padj < 0.05), ]$gene 

genes1 <- c(set1, set2, set3, set4) %>% unique()

# down-genes
set1 <- res1.df[which(res1.df$log2FoldChange < -0.75 &  res1.df$padj < 0.05), ]$gene
set2 <- res2.df[which(res2.df$log2FoldChange < -0.75 &  res2.df$padj < 0.05), ]$gene
set3 <- res3.df[which(res3.df$log2FoldChange < -0.75 &  res3.df$padj < 0.05), ]$gene
set4 <- res4.df[which(res4.df$log2FoldChange < -0.75 &  res4.df$padj < 0.05), ]$gene 

genes2 <- c(set1, set2, set3, set4) %>% unique()

# all DE genes
genes <- c(genes1, genes2) %>% unique()

# phenotypic genes

pheno_genes <- cledi_pheno_genes$...1

##################
# 1st Run GENIE3
##################

pheno_genes <- pheno_genes[!pheno_genes %in% mir7target.df$V1] 
genes1 <- c(genes, pheno_genes) 
#genes1 <- genes1[c(!genes1 %in% mir7target.df$V1)]
genes2 <- genes[c(!genes %in% pheno_genes & !genes %in% mir7target.df$V1)]

matrix_test <- DGEm_normalized[rownames(DGEm_normalized) %in% genes1, ] %>% as.matrix

weightMat_test <- GENIE3(matrix_test, regulators = genes2, targets = pheno_genes)
linkList_test <- getLinkList(weightMat_test)

saveRDS(linkList_test, file = "/data/rajewsky/home/skim/Neuronal_Activity_Cledi/201911_m7_overexpression/publication_codes_figures/GENIE3_run_cledi_pheno_all_1st_run_Elisabeth.rda")

# # subset by DE genes
# matrix_test <- DGEm_normalized[rownames(DGEm_normalized) %in% c(genes, pheno_genes), ] %>% as.matrix
# 
# weightMat_test <- GENIE3(matrix_test)
# linkList_test <- getLinkList(weightMat_test)
# 
# # Filter by weight
# linkList_test2 <- linkList_test %>% subset(weight > 0.01 & targetGene %in% pheno_genes)
# 
# saveRDS(linkList_test2, file = "/data/rajewsky/home/skim/Neuronal_Activity_Cledi/201911_m7_overexpression/publication_codes_figures/GENIE3_run_cledi_pheno_all_1st_run.rda")


##################
# 2nd Run GENIE3
##################

genes2 <- c(set2, set3, set4) %>% unique()

regul_genes <- genes2[genes2 %in% mir7target.df$V1]

target_genes <- linkList_test$regulatoryGene[!linkList_test$regulatoryGene %in% mir7target.df$V1] %>% as.character() %>% unique()

# subset by DE genes
matrix_test <- DGEm_normalized[rownames(DGEm_normalized) %in% c(regul_genes, target_genes), ] %>% as.matrix

weightMat_test <- GENIE3(matrix_test, regulators = regul_genes, targets = target_genes)
linkList_test3 <- getLinkList(weightMat_test)

saveRDS(linkList_test3, file = "/data/rajewsky/home/skim/Neuronal_Activity_Cledi/201911_m7_overexpression/publication_codes_figures/GENIE3_run_cledi_pheno_all_2nd_run_Elisabeth.rda")



# weightMat_test <- GENIE3(matrix_test)
# linkList_test3 <- getLinkList(weightMat_test)
# 
# linkList_test4 <- linkList_test3 %>% 
#   dplyr::filter(weight > 0.01) %>% 
#   dplyr::filter(targetGene %in% target_genes & regulatoryGene %in% regul_genes) %>% 
#   dplyr::filter(!targetGene %in% pheno_genes & !regulatoryGene %in% pheno_genes) %>%
#   dplyr::filter(targetGene %in% linkList_test2$regulatoryGene) %>% 
#   dplyr::filter(!regulatoryGene %in% linkList_test2$regulatoryGene)
# 
# saveRDS(linkList_test4, file = "/data/rajewsky/home/skim/Neuronal_Activity_Cledi/201911_m7_overexpression/publication_codes_figures/GENIE3_run_cledi_pheno_all_2nd_run.rda")


################################################
# 5th: Cledi phenotype genes only statistical significant and miR-7 DE genes
################################################

# up-genes
set1 <- res1.df[which(res1.df$log2FoldChange > 0.75 &  res1.df$padj < 0.05), ]$gene
set2 <- res2.df[which(res2.df$log2FoldChange > 0.75 &  res2.df$padj < 0.05), ]$gene
set3 <- res3.df[which(res3.df$log2FoldChange > 0.75 &  res3.df$padj < 0.05), ]$gene
set4 <- res4.df[which(res4.df$log2FoldChange > 0.75 &  res4.df$padj < 0.05), ]$gene 

genes1 <- c(set2, set3, set4) %>% unique()

# down-genes
set1 <- res1.df[which(res1.df$log2FoldChange < -0.75 &  res1.df$padj < 0.05), ]$gene
set2 <- res2.df[which(res2.df$log2FoldChange < -0.75 &  res2.df$padj < 0.05), ]$gene
set3 <- res3.df[which(res3.df$log2FoldChange < -0.75 &  res3.df$padj < 0.05), ]$gene
set4 <- res4.df[which(res4.df$log2FoldChange < -0.75 &  res4.df$padj < 0.05), ]$gene 

genes2 <- c(set2, set3, set4) %>% unique()

# all DE genes
genes <- c(genes1, genes2)

# phenotypic genes

pheno_genes <- cledi_pheno_genes$...1

# subset only statistical significant
pheno_genes <- pheno_genes[pheno_genes %in% genes]

##################
# 1st Run GENIE3
##################


pheno_genes <- pheno_genes[!pheno_genes %in% mir7target.df$V1] 
genes1 <- c(genes, pheno_genes) 
#genes1 <- genes1[c(!genes1 %in% mir7target.df$V1)]
genes2 <- genes[c(!genes %in% pheno_genes & !genes %in% mir7target.df$V1)]

matrix_test <- DGEm_normalized[rownames(DGEm_normalized) %in% genes1, ] %>% as.matrix

weightMat_test <- GENIE3(matrix_test, regulators = genes2, targets = pheno_genes)
linkList_test <- getLinkList(weightMat_test)

saveRDS(linkList_test, file = "/data/rajewsky/home/skim/Neuronal_Activity_Cledi/201911_m7_overexpression/publication_codes_figures/GENIE3_run_cledi_pheno_stat_1st_run_Elisabeth.rda")


# # subset by DE genes
# matrix_test <- DGEm_normalized[rownames(DGEm_normalized) %in% c(genes, pheno_genes), ] %>% as.matrix
# 
# weightMat_test <- GENIE3(matrix_test)
# linkList_test <- getLinkList(weightMat_test)
# 
# # Filter by weight
# linkList_test2 <- linkList_test %>% subset(weight > 0.01 & targetGene %in% pheno_genes)
# saveRDS(linkList_test2, file = "/data/rajewsky/home/skim/Neuronal_Activity_Cledi/201911_m7_overexpression/publication_codes_figures/GENIE3_run_cledi_pheno_stat_1st_run.rda")


##################
# 2nd Run GENIE3
##################
genes2 <- c(set2, set3, set4) %>% unique()

regul_genes <- genes2[genes2 %in% mir7target.df$V1]

target_genes <- linkList_test$regulatoryGene[!linkList_test$regulatoryGene %in% mir7target.df$V1] %>% as.character() %>% unique()

# subset by DE genes
matrix_test <- DGEm_normalized[rownames(DGEm_normalized) %in% c(regul_genes, target_genes), ] %>% as.matrix

weightMat_test <- GENIE3(matrix_test, regulators = regul_genes, targets = target_genes)
linkList_test3 <- getLinkList(weightMat_test)

saveRDS(linkList_test3, file = "/data/rajewsky/home/skim/Neuronal_Activity_Cledi/201911_m7_overexpression/publication_codes_figures/GENIE3_run_cledi_pheno_stat_2nd_run_Elisabeth.rda")

# weightMat_test <- GENIE3(matrix_test)
# linkList_test3 <- getLinkList(weightMat_test)
# 
# linkList_test4 <- linkList_test3 %>% 
#   dplyr::filter(weight > 0.01) %>% 
#   dplyr::filter(targetGene %in% target_genes & regulatoryGene %in% regul_genes) %>% 
#   dplyr::filter(!targetGene %in% pheno_genes & !regulatoryGene %in% pheno_genes) %>%
#   dplyr::filter(targetGene %in% linkList_test2$regulatoryGene) %>% 
#   dplyr::filter(!regulatoryGene %in% linkList_test2$regulatoryGene)
# 
# saveRDS(linkList_test4, file = "/data/rajewsky/home/skim/Neuronal_Activity_Cledi/201911_m7_overexpression/publication_codes_figures/GENIE3_run_cledi_pheno_stat_2nd_run.rda")

##################################
# 6th: Cledi phenotype genes only statistical significant and Including Cdr1as DE genes
##################################

# up-genes
set1 <- res1.df[which(res1.df$log2FoldChange > 0.75 &  res1.df$padj < 0.05), ]$gene
set2 <- res2.df[which(res2.df$log2FoldChange > 0.75 &  res2.df$padj < 0.05), ]$gene
set3 <- res3.df[which(res3.df$log2FoldChange > 0.75 &  res3.df$padj < 0.05), ]$gene
set4 <- res4.df[which(res4.df$log2FoldChange > 0.75 &  res4.df$padj < 0.05), ]$gene 

genes1 <- c(set1, set2, set3, set4) %>% unique()

# down-genes
set1 <- res1.df[which(res1.df$log2FoldChange < -0.75 &  res1.df$padj < 0.05), ]$gene
set2 <- res2.df[which(res2.df$log2FoldChange < -0.75 &  res2.df$padj < 0.05), ]$gene
set3 <- res3.df[which(res3.df$log2FoldChange < -0.75 &  res3.df$padj < 0.05), ]$gene
set4 <- res4.df[which(res4.df$log2FoldChange < -0.75 &  res4.df$padj < 0.05), ]$gene 

genes2 <- c(set1, set2, set3, set4) %>% unique()

# all DE genes
genes <- c(genes1, genes2) %>% unique()

# phenotypic genes

pheno_genes <- cledi_pheno_genes$...1

# subset only statistical significant
pheno_genes <- pheno_genes[pheno_genes %in% genes]

##################
# 1st Run GENIE3
##################

pheno_genes <- pheno_genes[!pheno_genes %in% mir7target.df$V1] 
genes1 <- c(genes, pheno_genes) 
#genes1 <- genes1[c(!genes1 %in% mir7target.df$V1)]
genes2 <- genes[c(!genes %in% pheno_genes & !genes %in% mir7target.df$V1)]

matrix_test <- DGEm_normalized[rownames(DGEm_normalized) %in% genes1, ] %>% as.matrix

weightMat_test <- GENIE3(matrix_test, regulators = genes2, targets = pheno_genes)
linkList_test <- getLinkList(weightMat_test)

saveRDS(linkList_test, file = "/data/rajewsky/home/skim/Neuronal_Activity_Cledi/201911_m7_overexpression/publication_codes_figures/GENIE3_run_cledi_pheno_stat_all_1st_run_Elisabeth.rda")


# # subset by DE genes
# matrix_test <- DGEm_normalized[rownames(DGEm_normalized) %in% c(genes, pheno_genes), ] %>% as.matrix
# 
# weightMat_test <- GENIE3(matrix_test)
# linkList_test <- getLinkList(weightMat_test)
# 
# # Filter by weight
# linkList_test2 <- linkList_test %>% subset(weight > 0.01 & targetGene %in% pheno_genes)
# 
# saveRDS(linkList_test2, file = "/data/rajewsky/home/skim/Neuronal_Activity_Cledi/201911_m7_overexpression/publication_codes_figures/GENIE3_run_cledi_pheno_stat_all_1st_run.rda")


##################
# 2nd Run GENIE3
##################

genes2 <- c(set2, set3, set4) %>% unique()

regul_genes <- genes2[genes2 %in% mir7target.df$V1]

target_genes <- linkList_test$regulatoryGene[!linkList_test$regulatoryGene %in% mir7target.df$V1] %>% as.character() %>% unique()

# subset by DE genes
matrix_test <- DGEm_normalized[rownames(DGEm_normalized) %in% c(regul_genes, target_genes), ] %>% as.matrix

weightMat_test <- GENIE3(matrix_test, regulators = regul_genes, targets = target_genes)
linkList_test3 <- getLinkList(weightMat_test)

saveRDS(linkList_test3, file = "/data/rajewsky/home/skim/Neuronal_Activity_Cledi/201911_m7_overexpression/publication_codes_figures/GENIE3_run_cledi_pheno_stat_all_2nd_run_Elisabeth.rda")


# weightMat_test <- GENIE3(matrix_test)
# linkList_test3 <- getLinkList(weightMat_test)
# 
# linkList_test4 <- linkList_test3 %>% 
#   dplyr::filter(weight > 0.01) %>% 
#   dplyr::filter(targetGene %in% target_genes & regulatoryGene %in% regul_genes) %>% 
#   dplyr::filter(!targetGene %in% pheno_genes & !regulatoryGene %in% pheno_genes) %>%
#   dplyr::filter(targetGene %in% linkList_test2$regulatoryGene) %>% 
#   dplyr::filter(!regulatoryGene %in% linkList_test2$regulatoryGene)
# 
# saveRDS(linkList_test4, file = "/data/rajewsky/home/skim/Neuronal_Activity_Cledi/201911_m7_overexpression/publication_codes_figures/GENIE3_run_cledi_pheno_stat_all_2nd_run.rda")
# 



##################################
# 7th: Cledi phenotype genes only statistical significant and Including Cdr1as DE genes and 3rd DE group (KO miR7 DE)
##################################

# phenotype genes from Cdr1as KO vs WT DE is empty so skip
# regul genes is Aldoc only in miR7 vs WT so error so skip
# Error message:
# Error in .checkArguments(exprMatrix = exprMatrix, regulators = regulators,  : Provide at least 2 potential regulators.


# up-genes
set1 <- res1.df[which(res1.df$log2FoldChange > 0.75 &  res1.df$padj < 0.05), ]$gene
set2 <- res2.df[which(res2.df$log2FoldChange > 0.75 &  res2.df$padj < 0.05), ]$gene
set3 <- res3.df[which(res3.df$log2FoldChange > 0.75 &  res3.df$padj < 0.05), ]$gene
set4 <- res4.df[which(res4.df$log2FoldChange > 0.75 &  res4.df$padj < 0.05), ]$gene 

genes1 <- c(set1, set2, set3, set4) %>% unique()

# down-genes
set5 <- res1.df[which(res1.df$log2FoldChange < -0.75 &  res1.df$padj < 0.05), ]$gene
set6 <- res2.df[which(res2.df$log2FoldChange < -0.75 &  res2.df$padj < 0.05), ]$gene
set7 <- res3.df[which(res3.df$log2FoldChange < -0.75 &  res3.df$padj < 0.05), ]$gene
set8 <- res4.df[which(res4.df$log2FoldChange < -0.75 &  res4.df$padj < 0.05), ]$gene 

genes2 <- c(set5, set6, set7, set8) %>% unique()

# all DE genes
genes <- c(genes1, genes2) %>% unique()

# phenotypic genes

pheno_genes <- cledi_pheno_genes$...1

# subset only statistical significant
pheno_genes <- pheno_genes[pheno_genes %in% genes]

##################
# 1st Run GENIE3
##################

pheno_genes <- pheno_genes[pheno_genes %in% c(set3, set7)]

pheno_genes <- pheno_genes[!pheno_genes %in% mir7target.df$V1] 

genes <- c(set3, set7) %>% unique()

genes1 <- c(genes, pheno_genes) 
#genes1 <- genes1[c(!genes1 %in% mir7target.df$V1)]
genes2 <- genes[c(!genes %in% pheno_genes & !genes %in% mir7target.df$V1)]

matrix_test <- DGEm_normalized[rownames(DGEm_normalized) %in% genes1, ] %>% as.matrix

weightMat_test <- GENIE3(matrix_test, regulators = genes2, targets = pheno_genes)
linkList_test <- getLinkList(weightMat_test)

saveRDS(linkList_test, file = "/data/rajewsky/home/skim/Neuronal_Activity_Cledi/201911_m7_overexpression/publication_codes_figures/GENIE3_run_cledi_pheno_stat_KO_miR7_all_1st_run_Elisabeth.rda")

##################
# 2nd Run GENIE3
##################

genes2 <- c(set7) %>% unique()

regul_genes <- genes2[genes2 %in% mir7target.df$V1]

target_genes <- linkList_test$regulatoryGene[!linkList_test$regulatoryGene %in% mir7target.df$V1] %>% as.character() %>% unique()

# subset by DE genes
matrix_test <- DGEm_normalized[rownames(DGEm_normalized) %in% c(regul_genes, target_genes), ] %>% as.matrix

weightMat_test <- GENIE3(matrix_test, regulators = regul_genes, targets = target_genes)
linkList_test3 <- getLinkList(weightMat_test)

saveRDS(linkList_test3, file = "/data/rajewsky/home/skim/Neuronal_Activity_Cledi/201911_m7_overexpression/publication_codes_figures/GENIE3_run_cledi_pheno_stat_KO_miR7_all_2nd_run_Elisabeth.rda")

##################################
# 8th: Cledi phenotype genes only statistical significant and Including Cdr1as DE genes and 4th DE group (KO miR7 vs WT DE)
##################################

# up-genes
set1 <- res1.df[which(res1.df$log2FoldChange > 0.75 &  res1.df$padj < 0.05), ]$gene
set2 <- res2.df[which(res2.df$log2FoldChange > 0.75 &  res2.df$padj < 0.05), ]$gene
set3 <- res3.df[which(res3.df$log2FoldChange > 0.75 &  res3.df$padj < 0.05), ]$gene
set4 <- res4.df[which(res4.df$log2FoldChange > 0.75 &  res4.df$padj < 0.05), ]$gene 

genes1 <- c(set1, set2, set3, set4) %>% unique()

# down-genes
set5 <- res1.df[which(res1.df$log2FoldChange < -0.75 &  res1.df$padj < 0.05), ]$gene
set6 <- res2.df[which(res2.df$log2FoldChange < -0.75 &  res2.df$padj < 0.05), ]$gene
set7 <- res3.df[which(res3.df$log2FoldChange < -0.75 &  res3.df$padj < 0.05), ]$gene
set8 <- res4.df[which(res4.df$log2FoldChange < -0.75 &  res4.df$padj < 0.05), ]$gene 

genes2 <- c(set5, set6, set7, set8) %>% unique()

# all DE genes
genes <- c(genes1, genes2) %>% unique()

# phenotypic genes

pheno_genes <- cledi_pheno_genes$...1

# subset only statistical significant
pheno_genes <- pheno_genes[pheno_genes %in% genes]

##################
# 1st Run GENIE3
##################

pheno_genes <- pheno_genes[pheno_genes %in% c(set4, set8)]

pheno_genes <- pheno_genes[!pheno_genes %in% mir7target.df$V1] 

genes <- c(set4, set8) %>% unique()

genes1 <- c(genes, pheno_genes) 
#genes1 <- genes1[c(!genes1 %in% mir7target.df$V1)]
genes2 <- genes[c(!genes %in% pheno_genes & !genes %in% mir7target.df$V1)]

matrix_test <- DGEm_normalized[rownames(DGEm_normalized) %in% genes1, ] %>% as.matrix

weightMat_test <- GENIE3(matrix_test, regulators = genes2, targets = pheno_genes)
linkList_test <- getLinkList(weightMat_test)

saveRDS(linkList_test, file = "/data/rajewsky/home/skim/Neuronal_Activity_Cledi/201911_m7_overexpression/publication_codes_figures/GENIE3_run_cledi_pheno_stat_KO_miR7_WT_all_1st_run_Elisabeth.rda")

##################
# 2nd Run GENIE3
##################

genes2 <- c(set8) %>% unique()

regul_genes <- genes2[genes2 %in% mir7target.df$V1]

target_genes <- linkList_test$regulatoryGene[!linkList_test$regulatoryGene %in% mir7target.df$V1] %>% as.character() %>% unique()

# subset by DE genes
matrix_test <- DGEm_normalized[rownames(DGEm_normalized) %in% c(regul_genes, target_genes), ] %>% as.matrix

weightMat_test <- GENIE3(matrix_test, regulators = regul_genes, targets = target_genes)
linkList_test3 <- getLinkList(weightMat_test)

saveRDS(linkList_test3, file = "/data/rajewsky/home/skim/Neuronal_Activity_Cledi/201911_m7_overexpression/publication_codes_figures/GENIE3_run_cledi_pheno_stat_KO_miR7_WT_all_2nd_run_Elisabeth.rda")

