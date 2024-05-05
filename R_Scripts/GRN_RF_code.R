.libPaths(c("/data/rajewsky/home/skim/R/usr_lib_Seurat/", "/data/rajewsky/home/skim/R/usr_lib_v4/"))
#.libPaths()
unloadNamespace("mgcv")
unloadNamespace("Matrix")


library(ggplot2)
library(gridExtra)
library(ggrepel)
library(dplyr)
library(cowplot)
library(paletteer)
library(stringr)
library(scales)
library(readxl)
library(randomForest)
library(pcaExplorer)
library(topGO)

cor_mat <- readRDS("/data/rajewsky/projects/cdr1as_ko_snRNA/codes_github/cdr1as/DGEm/WTN_batch_corr_nested_DGEm.rda")

# miR-7 target genes 
mir7target.df <- read.table("/data/rajewsky/projects/cdr1as_ko_snRNA/codes_github/cdr1as/ext_files/TargetScan7.2__miR-7-5p.predicted_targets_lists_without_header.txt")
# miR-7 targe genes from mirdb
mir7mirdb.df <- read.csv("/data/rajewsky/projects/cdr1as_ko_snRNA/codes_github/cdr1as/ext_files/20220423_mirdb_mir_7_extract.csv", header = TRUE)

# load pheno genes
cledi_pheno_genes <- read_excel("/data/rajewsky/projects/cdr1as_ko_snRNA/codes_github/cdr1as/ext_files/biased_list_of_phenotype_related_genes.xlsx", col_names = FALSE)

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

# list of DE genes by experiment
Cdr1as_DE_genes <- c(set1, set5)
WT_miR7_DE_genes <- c(set2, set6)
KO_miR7_DE_genes <- c(set3, set7)
KO_WT_miR7_DE_genes <- c(set4, set8)


computeRMSE <- function(pred, y){
  r <- sqrt(mean((pred-y)^2))
  return(r)
}

computePairCorrelation <- function(pred, y, ytrain){
  ytm <- matrix(rep(ytrain, length(y)), length(y), length(ytrain), byrow=T) # rep not necessary, but make intention clearer
  gt <- c(matrix(rep(y, length(ytrain)), length(y), length(ytrain), byrow=F)-ytm) # ground truth # c() makes column-wise concatenation 
  pr <- c(matrix(rep(pred, length(ytrain)), length(pred), length(ytrain), byrow=F)-ytm)
  r <- (pr %*% gt)/(norm(pr, type="2")*norm(gt, type="2")) # cosine similarity # −1 meaning exactly opposite, to 1 meaning exactly the same, V and aV are maximally similar for any constant a
  return(r)
}


# DE genes
DE_genes <- c(genes1, genes2)

# miR-7 target genes (in both DB)
mir7_genes <- intersect(mir7target.df$V1, mir7mirdb.df$X4)

# pheno type (17 genes)
# genes in DE
pheno_genes <- cledi_pheno_genes$...1[cledi_pheno_genes$...1 %in% DE_genes]
pheno_genes <- pheno_genes[!pheno_genes %in% mir7_genes]


cor_mat <- cor_mat[rownames(cor_mat) %in% DE_genes, ]

r <- t(as.matrix(cor_mat))

# save for https://lce.biohpc.swmed.edu/ glasso test
# write.csv(r, "/data/rajewsky/home/skim/Neuronal_Activity_Cledi/201911_m7_overexpression/publication_codes_figures/DE_Analysis_DESeq2_batch_corrected_matrix_1350_genes.csv", row.names = FALSE) 

ugroup <- c(1:16)
group <- c(1:16)

Run_RF <- function(interesting_gene = pheno_genes[1], outdir = "/data/rajewsky/projects/cdr1as_ko_snRNA/codes_github/cdr1as/Results/GRN_RF/", layer = "low", seed = 1){
  
  
  dir.create(paste0(outdir, "seed_", seed,"/", interesting_gene), recursive = TRUE)
  
  evaluationAll <- matrix(NA, length(ugroup), 2)
  colnames(evaluationAll) <- c("RMSE", "Cosine similarity")
  rownames(evaluationAll) <- c(1:16)
  
  xRaw <- r[, !colnames(r) %in% interesting_gene ]
  if(layer == "low"){
    print("low")
    xRaw <- xRaw[, !colnames(xRaw) %in% c(mir7_genes, pheno_genes) ]
  }
  if(layer == "high"){
    print("high")
    xRaw <- xRaw[, !colnames(xRaw) %in% c(middle_genes, pheno_genes) ]
  }
  
  yRaw <- r[, interesting_gene]
  
  
  for (i in 1:length(ugroup)){
    #print(i)
    x <- xRaw
    y <- yRaw
    
    wtest <- which(group %in% c(ugroup[i]))
    
    
    print("start training")
    set.seed(seed)
    model <- randomForest(x[-wtest, , drop=F], y[-wtest], ntree=1000, nodesize=2, importance=T, keep.forest=T)
    
    saveRDS(model, paste0(outdir,"seed_", seed,"/",interesting_gene,"/", interesting_gene ,"_" , i,".rda"))
    pred <- predict(model, x[wtest, , drop=F])
    evaluation <- c(computeRMSE(pred, y[wtest]), computePairCorrelation(pred, y[wtest], y[-wtest]))
    evaluationAll[i,] <- evaluation
    
  }
  saveRDS(evaluationAll, paste0(outdir,"seed_", seed,"/",interesting_gene,"/", interesting_gene ,"_" , "evaluation",".rda"))
  
}

for(i in c(1:10)){

Run_RF(interesting_gene = pheno_genes[1], outdir = "/data/rajewsky/projects/cdr1as_ko_snRNA/codes_github/cdr1as/Results/GRN_RF/low/", seed = i, layer = "low")

for(j in c(2:15)){
  print(j)
  Run_RF(interesting_gene = pheno_genes[j], outdir = "/data/rajewsky/projects/cdr1as_ko_snRNA/codes_github/cdr1as/Results/GRN_RF/low/",seed = i, layer = "low")
  
}
}





extract_top_interaction <- function(gene = "Cplx1", path = "/data/rajewsky/projects/cdr1as_ko_snRNA/codes_github/cdr1as/Results/GRN_RF/low/seed_1/", k = 10){
  
  model <- readRDS(paste0(path,gene, "/", gene, "_1.rda"))
  
  impt <- model$importance[,2, drop=F]
  
  for(i in c(2:16)){
    model <- readRDS(paste0(path,gene, "/", gene, "_", i,".rda"))
    
    impt2 <- model$importance[,2, drop=F]
    
    impt <- cbind(impt, impt2)
    
  }
  
  
  
  a <- apply(impt, 1, median)
  
  s <- sort(a, decreasing=T, na.last=T, index.return=T) 
  
  a2 <- a[s$ix[1:k]] %>% as.data.frame
  
  colnames(a2) <- "IncNodePurity"
  a2$target <- gene
  a2$regulator <- rownames(a2)
  
  return(a2)
}



for(i in c(1:10)){

tmp <- extract_top_interaction(gene = pheno_genes[1], path = paste0("/data/rajewsky/projects/cdr1as_ko_snRNA/codes_github/cdr1as/Results/GRN_RF/low/seed_", i, "/"), k = 20)
for(j in c(2:15)){
  tmp2 <- extract_top_interaction(gene = pheno_genes[j], path = paste0("/data/rajewsky/projects/cdr1as_ko_snRNA/codes_github/cdr1as/Results/GRN_RF/low/seed_",i, "/"), k = 20)
  tmp <- rbind(tmp, tmp2)

}



tmp3 <- table(tmp$regulator) %>% as.data.frame()
tmp4 <- tmp3[tmp3$Freq > 1, ]
final_tmp <- tmp %>% filter(regulator %in% tmp4$Var1)

saveRDS(final_tmp, paste0("/data/rajewsky/projects/cdr1as_ko_snRNA/codes_github/cdr1as/Results/GRN_RF/low/seed_",i,"/low_layer.rda"))

middle_genes <- tmp4$Var1 %>% as.character()



# High layer
Run_RF(interesting_gene = middle_genes[1], outdir = "/data/rajewsky/projects/cdr1as_ko_snRNA/codes_github/cdr1as/Results/GRN_RF/high/", seed = i, layer = "high")

for(j in c(2:length(middle_genes))){
  #print(i)
  Run_RF(interesting_gene = middle_genes[j], outdir = "/data/rajewsky/projects/cdr1as_ko_snRNA/codes_github/cdr1as/Results/GRN_RF/high/", seed = i, layer = "high")

}

high_tmp <- extract_top_interaction(gene = middle_genes[1], path = paste0("/data/rajewsky/projects/cdr1as_ko_snRNA/codes_github/cdr1as/Results/GRN_RF/high/seed_",i,"/"), k = 20)

for(j in c(2:length(middle_genes))){
  high_tmp2 <- extract_top_interaction(gene = middle_genes[j],path = paste0("/data/rajewsky/projects/cdr1as_ko_snRNA/codes_github/cdr1as/Results/GRN_RF/high/seed_", i,"/"), k = 20)
  high_tmp <- rbind(high_tmp, high_tmp2)

}


high_tmp3 <- table(high_tmp$regulator) %>% as.data.frame()
high_tmp4 <- high_tmp3[high_tmp3$Var1 %in% mir7_genes,  ]
high_tmp4 <- high_tmp4[high_tmp4$Freq > 1, ]

final_high_tmp <- high_tmp %>% filter(regulator %in% high_tmp4$Var1)
saveRDS(final_high_tmp, paste0("/data/rajewsky/projects/cdr1as_ko_snRNA/codes_github/cdr1as/Results/GRN_RF/high/seed_",i,"/high_layer.rda"))

}


# ---------------------- #
# for graph and GO
# ---------------------- #

# low layer 
low_count <- data.frame()

for(i in c(1:10)){
  final_low <- readRDS(paste0("/data/rajewsky/projects/cdr1as_ko_snRNA/codes_github/cdr1as/Results/GRN_RF/low/seed_", i,"/low_layer.rda"))
  final_low$seed <- i
  low_count <- rbind(low_count, final_low)
  
  
}

low_count2 <- low_count %>% dplyr::group_by(target, regulator)%>%
  dplyr::summarise(count = n())

# high layer
high_count <- data.frame()

for(i in c(1:10)){
  final_high <- readRDS(paste0("/data/rajewsky/projects/cdr1as_ko_snRNA/codes_github/cdr1as/Results/GRN_RF/high/seed_", i,"/high_layer.rda"))
  final_high$seed <- i
  high_count <- rbind(high_count, final_high)
  
  
}

high_count2 <- high_count %>% dplyr::group_by(target, regulator)%>%
  dplyr::summarise(count = n())

disconnetcted_genes <- setdiff(low_count2$regulator, high_count2$target)
low_count2 <- low_count2 %>% filter(!regulator %in% disconnetcted_genes)
input <- rbind(low_count2, high_count2)
input <- input[, c("regulator", "target" , "count")]
colnames(input) <- c("regulatoryGene", "targetGene", "weight")

input <- input[input$weight > 1, ]

input$weight <- input$weight / 10

saveRDS(input, "/data/rajewsky/projects/cdr1as_ko_snRNA/codes_github/cdr1as/Results/GRN_RF/input.rda")


# ---------------------- #
# GO for RF results
# ---------------------- #

res1.df <- readRDS(file = "/data/rajewsky/projects/cdr1as_ko_snRNA/codes_github/cdr1as/Results/DESeq2_WTN_batch_corrected_res1.rda")

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
  
  saveRDS(topgo_BP, file = paste0("/data/rajewsky/projects/cdr1as_ko_snRNA/codes_github/cdr1as/Results/GRN_RF/", name,"_BP.rda"))
  
  ########
  # Up MF
  ########
  topgo_MF <- topGOtable(GO_targets, bg_ids,
                         ontology = "MF",
                         mapping = "org.Mm.eg.db",
                         geneID = "symbol")
  
  saveRDS(topgo_MF, file = paste0("/data/rajewsky/projects/cdr1as_ko_snRNA/codes_github/cdr1as/Results/GRN_RF/", name,"_MF.rda"))
  
  ########
  # Up CC
  ########
  topgo_CC <- topGOtable(GO_targets, bg_ids,
                         ontology = "CC",
                         mapping = "org.Mm.eg.db",
                         geneID = "symbol")
  
  saveRDS(topgo_CC, file = paste0("/data/rajewsky/projects/cdr1as_ko_snRNA/codes_github/cdr1as/Results/GRN_RF/", name,"_CC.rda"))
  
  
}

GO_targets1 <- c(input$regulatoryGene, input$targetGene)

run_GO_all(GO_targets = GO_targets1, name = "GRN_RF")

