library(clusterProfiler)
library(topGO)
library(pcaExplorer)

res1.df <- readRDS(file = "/data/rajewsky/home/skim/Neuronal_Activity_Cledi/201911_m7_overexpression/publication_codes_figures/DE_Analysis_DESeq2_res1_KO_WT.rda")

res4.df <- readRDS(file = "/data/rajewsky/home/skim/Neuronal_Activity_Cledi/201911_m7_overexpression/publication_codes_figures/DE_Analysis_DESeq2_res4_KOm7ove_WT.rda")

result2.df <- readRDS(file = "/data/rajewsky/home/skim/Neuronal_Activity_Cledi/201911_m7_overexpression/publication_codes_figures/DE_Analysis_DESeq2_nested_result2.rda")

result3.df <- readRDS(file = "/data/rajewsky/home/skim/Neuronal_Activity_Cledi/201911_m7_overexpression/publication_codes_figures/DE_Analysis_DESeq2_nested_result3.rda")

set1 <- res1.df[which(res1.df$log2FoldChange > 0.5 &  res1.df$padj < 0.05), ]
set2 <- result2.df[which(result2.df$log2FoldChange > 0.5 &  result2.df$padj < 0.05), ]
set3 <- result3.df[which(result3.df$log2FoldChange > 0.5 &  result3.df$padj < 0.05), ]
set4 <- res4.df[which(res4.df$log2FoldChange > 0.5 &  res4.df$padj < 0.05), ] 


# down-genes
set5 <- res1.df[which(res1.df$log2FoldChange < -0.5 &  res1.df$padj < 0.05), ]
set6 <- result2.df[which(result2.df$log2FoldChange < -0.5 &  result2.df$padj < 0.05), ]
set7 <- result3.df[which(result3.df$log2FoldChange < -0.5 &  result3.df$padj < 0.05), ]
set8 <- res4.df[which(res4.df$log2FoldChange < -0.5 &  res4.df$padj < 0.05), ] 

bg_ids <- res4.df$gene  


Run_GO <- function(GO_targets1 = GO_targets2, exp = "KO_WT_Up"){
  
  topgo_BP <- topGOtable(GO_targets1, bg_ids,
                             ontology = "BP",
                             mapping = "org.Mm.eg.db",
                             geneID = "symbol")

  saveRDS(topgo_BP, file = paste0("/data/rajewsky/projects/cdr1as_ko_snRNA/codes/GO_analysis/", exp, "_GO_BP.rda"))

  topgo_MF <- topGOtable(GO_targets1, bg_ids,
                             ontology = "MF",
                             mapping = "org.Mm.eg.db",
                             geneID = "symbol")

  saveRDS(topgo_MF, file = paste0("/data/rajewsky/projects/cdr1as_ko_snRNA/codes/GO_analysis/", exp, "_GO_MF.rda"))

  topgo_CC <- topGOtable(GO_targets1, bg_ids,
                             ontology = "CC",
                             mapping = "org.Mm.eg.db",
                             geneID = "symbol")

  saveRDS(topgo_CC, file = paste0("/data/rajewsky/projects/cdr1as_ko_snRNA/codes/GO_analysis/", exp, "_GO_CC.rda"))

  
}

GO_targets1 <- c(set1$gene) %>% unique()
Run_GO(GO_targets1 = GO_targets1, exp = "KO_WT_Up")
  
GO_targets2 <- c(set5$gene) %>% unique()
Run_GO(GO_targets1 = GO_targets2, exp = "KO_WT_Down")

GO_targets1 <- c(set2$gene) %>% unique()
Run_GO(GO_targets1 = GO_targets1, exp = "WTm7oe_WT_Up")

GO_targets1 <- c(set6$gene) %>% unique()
Run_GO(GO_targets1 = GO_targets1, exp = "WTm7oe_WT_Down")

GO_targets1 <- c(set3$gene) %>% unique()
Run_GO(GO_targets1 = GO_targets1, exp = "KOm7oe_KO_Up")

GO_targets1 <- c(set7$gene) %>% unique()
Run_GO(GO_targets1 = GO_targets1, exp = "KOm7oe_KO_Down")

GO_targets1 <- c(set4$gene) %>% unique()
Run_GO(GO_targets1 = GO_targets1, exp = "KOm7oe_WT_Up")

GO_targets1 <- c(set8$gene) %>% unique()
Run_GO(GO_targets1 = GO_targets1, exp = "KOm7oe_WT_Down")


