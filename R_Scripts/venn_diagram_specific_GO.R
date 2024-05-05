
.libPaths(c("/data/rajewsky/home/skim/R/usr_lib_Seurat/", "/data/rajewsky/home/skim/R/usr_lib_v4/"))

library(DESeq2)
library(dplyr)
library(pcaExplorer)
library(topGO)
library(ggplot2)
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

UP_KO_KOm7oe_specific <- Reduce(setdiff, list(set3,set1,set2, set4)) 

UP_WT_KOm7oe_specific <- Reduce(setdiff, list(set4,set1,set2, set3)) 

DOWN_KO_KOm7oe_specific <- Reduce(setdiff, list(set7,set5,set6, set8)) 

DOWN_WT_KOm7oe_specific <- Reduce(setdiff, list(set8,set5,set6, set7)) 

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

########################
# GO Up KO vs KOm7oe 
########################
run_GO_all(GO_targets = UP_KO_KOm7oe_specific, name = "KOm7oe_KO_SPECIFIC_UP")

########################
# GO Down KO vs KOm7oe 
########################
run_GO_all(GO_targets = DOWN_KO_KOm7oe_specific, name = "KOm7oe_KO_SPECIFIC_DOWN")

########################
# GO Up WT vs KOm7oe 
########################
run_GO_all(GO_targets = UP_WT_KOm7oe_specific, name = "KOm7oe_WT_SPECIFIC_UP")

########################
# GO Down WT vs KOm7oe 
########################
run_GO_all(GO_targets = DOWN_WT_KOm7oe_specific, name = "KOm7oe_WT_SPECIFIC_DOWN")



