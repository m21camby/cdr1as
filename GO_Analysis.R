
.libPaths(c("/data/rajewsky/home/skim/R/usr_lib_Seurat/", "/data/rajewsky/home/skim/R/usr_lib_v4/"))

library(DESeq2)
library(dplyr)
library(pcaExplorer)
library(topGO)
library(ggplot2)

res1.df <- readRDS(file = "/data/rajewsky/projects/cdr1as_ko_snRNA/codes_github/cdr1as/Results/DESeq2_WTN_batch_corrected_res1.rda")
res2.df <- readRDS(file = "/data/rajewsky/projects/cdr1as_ko_snRNA/codes_github/cdr1as/Results/DESeq2_WTN_nested_corrected_res2.rda")
res3.df <- readRDS(file = "/data/rajewsky/projects/cdr1as_ko_snRNA/codes_github/cdr1as/Results/DESeq2_WTN_nested_corrected_res3.rda")
res4.df <- readRDS(file = "/data/rajewsky/projects/cdr1as_ko_snRNA/codes_github/cdr1as/Results/DESeq2_WTN_batch_corrected_res4.rda")


# background genes
bg_ids <- rownames(res1.df)

#####################
# TopGO function
#####################
run_GO_all <- function(GO_targets = GO_targets1, name = "KO_WT_UP"){
  topgo_BP <- topGOtable(GO_targets, bg_ids,
                         ontology = "BP",
                         mapping = "org.Mm.eg.db",
                         geneID = "symbol")
  
  saveRDS(topgo_BP, file = paste0("/data/rajewsky/projects/cdr1as_ko_snRNA/codes_github/cdr1as/Results/GO_Analysis_", name,"_BP.rda"))
  
  ########
  # Up MF
  ########
  topgo_MF <- topGOtable(GO_targets, bg_ids,
                         ontology = "MF",
                         mapping = "org.Mm.eg.db",
                         geneID = "symbol")
  
  saveRDS(topgo_MF, file = paste0("/data/rajewsky/projects/cdr1as_ko_snRNA/codes_github/cdr1as/Results/GO_Analysis_", name,"_MF.rda"))
  
  ########
  # Up CC
  ########
  topgo_CC <- topGOtable(GO_targets, bg_ids,
                         ontology = "CC",
                         mapping = "org.Mm.eg.db",
                         geneID = "symbol")
  
  saveRDS(topgo_CC, file = paste0("/data/rajewsky/projects/cdr1as_ko_snRNA/codes_github/cdr1as/Results/GO_Analysis_", name,"_CC.rda"))
  
  
}




####################
# KO vs WT
####################

########
# GO Up 
########
GO_targets1 <- rownames(res1.df[which(res1.df$padj < .05 & res1.df$log2FoldChange > 0.75), ])

run_GO_all(GO_targets = GO_targets1, name = "KO_WT_UP")

###########
# GO Down 
###########
GO_targets1 <- rownames(res1.df[which(res1.df$padj < .05 & res1.df$log2FoldChange < -0.75), ])

run_GO_all(GO_targets = GO_targets1, name = "KO_WT_DOWN")



##########################
# WTmiR7ove vs WT
##########################
########
# GO Up 
########
GO_targets1 <- rownames(res2.df[which(res2.df$padj < .05 & res2.df$log2FoldChange > 0.75), ])

run_GO_all(GO_targets = GO_targets1, name = "WTm7oe_WT_UP")

###########
# GO Down 
###########
GO_targets1 <- rownames(res2.df[which(res2.df$padj < .05 & res2.df$log2FoldChange < -0.75), ])

run_GO_all(GO_targets = GO_targets1, name = "WTm7oe_WT_DOWN")


##########################
# KOmiR7ove vs KO
##########################

########
# GO Up 
########
GO_targets1 <- rownames(res3.df[which(res3.df$padj < .05 & res3.df$log2FoldChange > 0.75), ])

run_GO_all(GO_targets = GO_targets1, name = "KOm7oe_KO_UP")

###########
# GO Down 
###########
GO_targets1 <- rownames(res3.df[which(res3.df$padj < .05 & res3.df$log2FoldChange < -0.75), ])

run_GO_all(GO_targets = GO_targets1, name = "KOm7oe_KO_DOWN")


##########################
# KOmiR7ove vs WT
##########################

########
# GO Up 
########
GO_targets1 <- rownames(res4.df[which(res4.df$padj < .05 & res4.df$log2FoldChange > 0.75), ])

run_GO_all(GO_targets = GO_targets1, name = "KOm7oe_WT_UP")

###########
# GO Down 
###########
GO_targets1 <- rownames(res4.df[which(res4.df$padj < .05 & res4.df$log2FoldChange < -0.75), ])

run_GO_all(GO_targets = GO_targets1, name = "KOm7oe_WT_DOWN")





###########################
# KOm7oe vs KO up specific
###########################

set1 <- res1.df[which(res1.df$log2FoldChange > 0.75 &  res1.df$padj < 0.05), ]$gene
set2 <- res2.df[which(res2.df$log2FoldChange > 0.75 &  res2.df$padj < 0.05), ]$gene
set3 <- res3.df[which(res3.df$log2FoldChange > 0.75 &  res3.df$padj < 0.05), ]$gene
set4 <- res4.df[which(res4.df$log2FoldChange > 0.75 &  res4.df$padj < 0.05), ]$gene

genes <- Reduce(setdiff, list(set3,set1,set2, set4))

run_GO_all(GO_targets = genes, name = "KOm7ove_KO_UP_SPECIFIC")

###########################
# KOm7oe vs WT up specific
###########################

genes <- Reduce(setdiff, list(set4,set1,set2, set3))

run_GO_all(GO_targets = genes, name = "KOm7ove_WT_UP_SPECIFIC")

###########################
# KOm7oe vs KO down specific
###########################

set1 <- res1.df[which(res1.df$log2FoldChange < -0.75 &  res1.df$padj < 0.05), ]$gene
set2 <- res2.df[which(res2.df$log2FoldChange < -0.75 &  res2.df$padj < 0.05), ]$gene
set3 <- res3.df[which(res3.df$log2FoldChange < -0.75 &  res3.df$padj < 0.05), ]$gene
set4 <- res4.df[which(res4.df$log2FoldChange < -0.75 &  res4.df$padj < 0.05), ]$gene

genes <- Reduce(setdiff, list(set3,set1,set2, set4))

run_GO_all(GO_targets = genes, name = "KOm7ove_KO_DOWN_SPECIFIC")

###########################
# KOm7oe vs WT down specific
###########################

genes <- Reduce(setdiff, list(set4,set1,set2, set3))

run_GO_all(GO_targets = genes, name = "KOm7ove_WT_DOWN_SPECIFIC")


