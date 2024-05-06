.libPaths(c("/data/rajewsky/home/skim/R/usr_lib_Seurat/", "/data/rajewsky/home/skim/R/usr_lib_v4/"))

library(DESeq2)
library(dplyr)
library(pcaExplorer)
library(topGO)
library(WGCNA)
library(ggplot2)
library(genefilter)
library(gplots)
library(igraph)
library(penalized)
library(visNetwork)
library(flashClust)

DGEm <- readRDS(file = "/data/rajewsky/home/skim/Neuronal_Activity_Cledi/201911_m7_overexpression/publication_codes_figures/DGE_matrix.rda")

sample_info <- data.frame(condition = c("WT","WT","WT",
                                        "WTm7ove","WTm7ove","WTm7ove",
                                        "KO","KO","KO","KO",
                                        "KOm7ove","KOm7ove","KOm7ove","KOm7ove"), 
                          batch = c("Buch", "Mitte", "Mitte",
                                    "Buch", "Mitte","Mitte",
                                    "Mitte", "Mitte", "Mitte", "Mitte",
                                    "Mitte", "Mitte", "Mitte", "Mitte"), row.names = names(DGEm))

dds <- DESeqDataSetFromMatrix(DGEm, colData = sample_info, design = ~ condition + batch)
#dds <- DESeqDataSetFromMatrix(DGEm, colData = sample_info, design = ~ 0 + condition + overexpression + batch)

keep <- rowSums(counts(dds)) >= 10
table(keep)
dds <- dds[keep,]
dds$condition <- relevel(dds$condition, ref = "WT")
dds <- DESeq(dds)

##################################
# acquire batch corrected matrix
##################################
vsd <- varianceStabilizingTransformation(dds)
mat <- assay(vsd)
mat2 <- limma::removeBatchEffect(mat, vsd$batch)
assay(vsd) <- mat2
counts_batch_corrected <- assay(vsd)


##################
# HVG top 40%
##################

genes.df <- counts_batch_corrected %>% as.matrix %>% varFilter(var.func=IQR, var.cutoff=0.6, filterByQuantile=TRUE) %>% as.data.frame

gene.names <- rownames(genes.df)
genes_t.df <- t(genes.df)

# write.csv(genes_t.df, file = "/data/rajewsky/home/skim/Neuronal_Activity_Cledi/201911_m7_overexpression/R_Scripts/20211109_WGCNA_batch_corrected_VST_final/20211109_WGCNA_batch_corrected_VST_final_batch_genes_t_df.csv")

powers <- c(c(1:30))
sft  <- pickSoftThreshold(genes_t.df, 
                          dataIsExpr = TRUE, 
                          powerVector = powers, 
                          corFnc = bicor, 
                          verbose = 5,
                          corOptions = list(use ='p'), 
                          networkType = "signed")

par(mfrow = c(1,2));
cex1 = 0.9;

plot(sft$fitIndices[, 1],
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)",
     ylab = "Scale Free Topology Model Fit, signed R^2",
     main = paste("Scale independence")
)
text(sft$fitIndices[, 1],
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers, cex = cex1, col = "red"
)
abline(h = 0.90, col = "red")
plot(sft$fitIndices[, 1],
     sft$fitIndices[, 5],
     xlab = "Soft Threshold (power)",
     ylab = "Mean Connectivity",
     type = "n",
     main = paste("Mean connectivity")
)
text(sft$fitIndices[, 1],
     sft$fitIndices[, 5],
     labels = powers,
     cex = cex1, col = "red")



softPower = 30
adjacency = adjacency(genes_t.df, power = softPower, type = "signed")
# dim(adjacency)
# adjacency[c(1:10),c(1:10)]

cor <- WGCNA::cor
netwk <- blockwiseModules(genes_t.df,

                          corType="bicor",

                          # == Adjacency Function ==
                          power = softPower,
                          networkType = "signed",

                          # == Tree and Block Options ==
                          deepSplit = 2,
                          pamRespectsDendro = F,
                          # detectCutHeight = 0.75,
                          minModuleSize = 30,
                          maxBlockSize = 10000,

                          # == Module Adjustments ==
                          reassignThreshold = 0,
                          mergeCutHeight = 0.25,
                          minKMEtoStay = 0,

                          # == TOM == Archive the run results in TOM file (saves time)
                          TOMType="signed",
                          saveTOMs = T,
                          saveTOMFileBase = "/data/rajewsky/home/skim/Neuronal_Activity_Cledi/201911_m7_overexpression/publication_codes_figures/WGCNA_Analysis_batch_corrected_VST_final_TOM",

                          # == Output Options
                          numericLabels = F,
                          verbose = 3)


save(netwk,file= "/data/rajewsky/home/skim/Neuronal_Activity_Cledi/201911_m7_overexpression/publication_codes_figures/WGCNA_Analysis_batch_corrected_VST_signed_blockmodule.RData")

color = netwk[[1]]

plotDendroAndColors(
  netwk$dendrograms[[1]],
  #mergedColors,
  #unmergedColors,
  color,
  "Module colors",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05 )


load(file= "/data/rajewsky/home/skim/Neuronal_Activity_Cledi/201911_m7_overexpression/publication_codes_figures/WGCNA_Analysis_batch_corrected_VST_final_TOM-block.1.RData")

dissTOM = 1-TOM %>% as.matrix
diag(dissTOM) <- NA
myheatcol = colorpanel(250,'red',"orange",'lemonchiffon')

TOMplot(dissTOM, netwk$dendrograms[[1]], color, col=myheatcol)

module_df2 <- data.frame(
  gene_id = names(netwk$colors),
  colors = labels2colors(netwk$colors),
  colors2 = netwk$unmergedColors,
  color_check = netwk[[1]]
)

# netwk[[1]] is unmergedColors
# labels2colors = mergedColors

#table(module_df2$colors)
#table(module_df2$colors2)
table(module_df2$color_check)

colorh = netwk[[1]]
tophub <- chooseTopHubInEachModule(
  genes_t.df, 
  colorh, 
  omitColors = "grey", 
  power = 12, 
  type = "signed")

tophub.df <- tophub %>% as.data.frame


#################################
# 1350 genes by batch corrected 
#################################

cor_mat <- readRDS("/data/rajewsky/projects/cdr1as_ko_snRNA/codes/DE_Analysis_DESeq2_nested_corrected_data/DE_Analysis_DESeq2_nested_corrected_data_DGEm.rda")

res1.df <- readRDS(file = "/data/rajewsky/home/skim/Neuronal_Activity_Cledi/201911_m7_overexpression/publication_codes_figures/DE_Analysis_DESeq2_res1_KO_WT.rda")

res4.df <- readRDS(file = "/data/rajewsky/home/skim/Neuronal_Activity_Cledi/201911_m7_overexpression/publication_codes_figures/DE_Analysis_DESeq2_res4_KOm7ove_WT.rda")

result2.df <- readRDS(file = "/data/rajewsky/home/skim/Neuronal_Activity_Cledi/201911_m7_overexpression/publication_codes_figures/DE_Analysis_DESeq2_nested_result2.rda")

result3.df <- readRDS(file = "/data/rajewsky/home/skim/Neuronal_Activity_Cledi/201911_m7_overexpression/publication_codes_figures/DE_Analysis_DESeq2_nested_result3.rda")

# normal
# up-genes
set1 <- res1.df[which(res1.df$log2FoldChange > 0.5 &  res1.df$padj < 0.05), ]
set2 <- result2.df[which(result2.df$log2FoldChange > 0.5 &  result2.df$padj < 0.05), ]
set3 <- result3.df[which(result3.df$log2FoldChange > 0.5 &  result3.df$padj < 0.05), ]
set4 <- res4.df[which(res4.df$log2FoldChange > 0.5 &  res4.df$padj < 0.05), ] 


# down-genes
set5 <- res1.df[which(res1.df$log2FoldChange < -0.5 &  res1.df$padj < 0.05), ]
set6 <- result2.df[which(result2.df$log2FoldChange < -0.5 &  result2.df$padj < 0.05), ]
set7 <- result3.df[which(result3.df$log2FoldChange < -0.5 &  result3.df$padj < 0.05), ]
set8 <- res4.df[which(res4.df$log2FoldChange < -0.5 &  res4.df$padj < 0.05), ] 

genes <- c(set1$gene, set2$gene, set3$gene, set4$gene,
           set5$gene, set6$gene, set7$gene, set8$gene)

cor_mat_sub <- cor_mat[rownames(cor_mat) %in% genes, ]


# WGCNA

genes_t.df <- t(cor_mat_sub)


powers <- c(c(1:30))
sft  <- pickSoftThreshold(genes_t.df, 
                          dataIsExpr = TRUE, 
                          powerVector = powers, 
                          corFnc = bicor, 
                          verbose = 5,
                          corOptions = list(use ='p'), 
                          networkType = "signed")

par(mfrow = c(1,2));
cex1 = 0.9;

plot(sft$fitIndices[, 1],
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)",
     ylab = "Scale Free Topology Model Fit, signed R^2",
     main = paste("Scale independence")
)

text(sft$fitIndices[, 1],
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers, cex = cex1, col = "red"
)
abline(h = 0.90, col = "red")
plot(sft$fitIndices[, 1],
     sft$fitIndices[, 5],
     xlab = "Soft Threshold (power)",
     ylab = "Mean Connectivity",
     type = "n",
     main = paste("Mean connectivity")
)
text(sft$fitIndices[, 1],
     sft$fitIndices[, 5],
     labels = powers,
     cex = cex1, col = "red")



softPower = 30
adjacency = adjacency(genes_t.df, power = softPower, type = "signed")
# dim(adjacency)
# adjacency[c(1:10),c(1:10)]

cor <- WGCNA::cor
netwk <- blockwiseModules(genes_t.df,
                          
                          corType="bicor",
                          
                          # == Adjacency Function ==
                          power = softPower,
                          networkType = "signed",
                          
                          # == Tree and Block Options ==
                          deepSplit = 2,
                          pamRespectsDendro = F,
                          # detectCutHeight = 0.75,
                          minModuleSize = 30,
                          maxBlockSize = 10000,
                          
                          # == Module Adjustments ==
                          reassignThreshold = 0,
                          mergeCutHeight = 0.2,
                          minKMEtoStay = 0,
                          
                          # == TOM == Archive the run results in TOM file (saves time)
                          TOMType="signed",
                          saveTOMs = T,
                          saveTOMFileBase = "/data/rajewsky/home/skim/Neuronal_Activity_Cledi/201911_m7_overexpression/publication_codes_figures/WGCNA_Analysis/WGCNA_Analysis_batch_corrected_VST_final_TOM_1350",
                          
                          # == Output Options
                          numericLabels = F,
                          verbose = 3)


save(netwk,file= "/data/rajewsky/home/skim/Neuronal_Activity_Cledi/201911_m7_overexpression/publication_codes_figures/WGCNA_Analysis/WGCNA_Analysis_batch_corrected_VST_signed_blockmodule_1350.RData")

load(file= "/data/rajewsky/home/skim/Neuronal_Activity_Cledi/201911_m7_overexpression/publication_codes_figures/WGCNA_Analysis/WGCNA_Analysis_batch_corrected_VST_signed_blockmodule_1350.RData")


color = netwk[[1]]

plotDendroAndColors(
  netwk$dendrograms[[1]],
  #mergedColors,
  #unmergedColors,
  color,
  "Module colors",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05 )


load(file= "/data/rajewsky/home/skim/Neuronal_Activity_Cledi/201911_m7_overexpression/publication_codes_figures/WGCNA_Analysis/WGCNA_Analysis_batch_corrected_VST_final_TOM_1350-block.1.RData")

dissTOM = 1-TOM %>% as.matrix
diag(dissTOM) <- NA
myheatcol = colorpanel(250,'red',"orange",'lemonchiffon')

TOMplot(dissTOM, netwk$dendrograms[[1]], color, col=myheatcol)

module_df2 <- data.frame(
  gene_id = names(netwk$colors),
  colors = labels2colors(netwk$colors),
  colors2 = netwk$unmergedColors,
  color_check = netwk[[1]],
  num = c(1:1350)
)

write.csv(module_df2, file = "/data/rajewsky/home/skim/Neuronal_Activity_Cledi/201911_m7_overexpression/publication_codes_figures/WGCNA_Analysis/WGCNA_Analysis_batch_corrected_1350_module2.csv")

# netwk[[1]] is unmergedColors
# labels2colors = mergedColors

#table(module_df2$colors)
#table(module_df2$colors2)
table(module_df2$color_check)

colorh = netwk[[1]]
tophub <- chooseTopHubInEachModule(
  genes_t.df, 
  colorh, 
  omitColors = "grey", 
  power = 12, 
  type = "signed")

tophub.df <- tophub %>% as.data.frame

saveRDS(tophub.df, file = "/data/rajewsky/home/skim/Neuronal_Activity_Cledi/201911_m7_overexpression/publication_codes_figures/WGCNA_Analysis/WGCNA_Analysis_tophub_df_1350.rda")

interesting_colors <- module_df2$color_check %>% unique


# extract number of edges

num_edge_extract <- function(color = "blue", edge_thre = 0.2){
  
  # extract length of genes in color
  color_num <- length(module_df2[module_df2$color_check %in% color,]$num)
  
  # extract TOM by color
  TOM2 <- subsetTOM(
    genes_t.df, 
    subset = module_df2[module_df2$color_check %in% color,]$num,
    corFnc = "bicor", 
    weights = NULL,  
    networkType = "signed", 
    power = 30, 
    verbose = 1, indent = 0)
  
  # TOM to adjacency matrix
  adj <- TOM2
  adj[adj > edge_thre] = 1
  adj[adj != 1] = 0
  
  # adjacency matrix to network
  network <- graph.adjacency(adj, mode = "directed")
  
  # extract edges per gene
  edge_matrix <- network[1] %>% as.data.frame()
  
  for(i in c(2:color_num)){
    edge_matrix2 <- network[i] %>% as.data.frame()
    edge_matrix <- cbind(edge_matrix, edge_matrix2)
  }
  
  colnames(edge_matrix) <- rownames(edge_matrix)
  
  final_output <- colSums(edge_matrix) %>% as.data.frame()
  
  colnames(final_output) <- "edges"
  
  final_output$gene <- rownames(final_output)
  
  return(final_output)
}


for(i in interesting_colors){
  print(i)
  # run
  final <- num_edge_extract(color = i, edge_thre = 0.1)
  saveRDS(final, file = paste0("/data/rajewsky/home/skim/Neuronal_Activity_Cledi/201911_m7_overexpression/publication_codes_figures/WGCNA_Analysis/WGCNA_Analysis_edges_",i ,"_1350.rda"))
}


blue <- readRDS(file = paste0("/data/rajewsky/home/skim/Neuronal_Activity_Cledi/201911_m7_overexpression/publication_codes_figures/WGCNA_Analysis/WGCNA_Analysis_edges_","blue" ,"_1350.rda"))

#################
# adjency matrix
#################

dissTOM <- 1-TOM %>% as.matrix
dissTOM.df <- dissTOM %>% as.data.frame()
colnames(dissTOM.df) <- names(netwk$colors)
rownames(dissTOM.df) <- names(netwk$colors)


softPower = 30
adjacency = adjacency(genes_t.df, power = softPower, type = "signed")

adjacency.df <- adjacency %>% as.data.frame()

saveRDS(adjacency.df, file = "/data/rajewsky/home/skim/Neuronal_Activity_Cledi/201911_m7_overexpression/publication_codes_figures/WGCNA_Analysis/WGCNA_Analysis_1350_adjacency_matrix.rda")


####################
# correlation plot
####################
# miR-7 target genes 
mir7target.df <- read.table("/data/rajewsky/home/skim/Neuronal_Activity_Cledi/201911_m7_overexpression/R_Scripts/external_files/TargetScan7.2__miR-7-5p.predicted_targets_lists_without_header.txt")

# miR-7 targe genes from mirdb
mir7mirdb.df <- read.csv("/data/rajewsky/home/skim/Neuronal_Activity_Cledi/201911_m7_overexpression/R_Scripts/20220423_mirdb_mir_7_extract.csv", header = TRUE)

# miR-122 target genes 
mir122target.df <- read.table("/data/rajewsky/home/skim/Neuronal_Activity_Cledi/201911_m7_overexpression/R_Scripts/external_files/TargetScan7.2__miR-122-5p.predicted_targets_without_header.txt")

# miR-122 targe genes from mirdb
mir122mirdb.df <- read.csv("/data/rajewsky/home/skim/Neuronal_Activity_Cledi/201911_m7_overexpression/R_Scripts/20220423_mirdb_mir_122_extract.csv", header = TRUE)

# miR-7 target genes (in both DB)
mir7_genes <- intersect(mir7target.df$V1, mir7mirdb.df$X4)

# miR_122 target genes (in both DB)
mir122_genes <- intersect(mir122target.df$V1, mir122mirdb.df$X4)


# calculate traits
mir7.df <- cor_mat_sub[rownames(cor_mat_sub) %in% mir7_genes, ] %>% colMeans
mir122.df <- cor_mat_sub[rownames(cor_mat_sub) %in% mir122_genes, ] %>% colMeans
mir122_specific.df <- cor_mat_sub[rownames(cor_mat_sub) %in% setdiff(mir122_genes, mir7_genes), ] %>% colMeans

KO_WT_UP = cor_mat_sub[rownames(cor_mat_sub) %in% res1.df[which(res1.df$log2FoldChange > 0.5 & res1.df$padj < 0.05),]$gene , ] %>% colMeans
KO_WT_DOWN = cor_mat_sub[rownames(cor_mat_sub) %in% res1.df[which(res1.df$log2FoldChange < 0.5 & res1.df$padj < 0.05),]$gene , ] %>% colMeans

mir7_WT_UP = cor_mat_sub[rownames(cor_mat_sub) %in% result2.df[which(result2.df$log2FoldChange > 0.5 & result2.df$padj < 0.05),]$gene , ] %>% colMeans
mir7_WT_DOWN = cor_mat_sub[rownames(cor_mat_sub) %in% result2.df[which(result2.df$log2FoldChange < 0.5 & result2.df$padj < 0.05),]$gene , ] %>% colMeans

mir7_KO_UP = cor_mat_sub[rownames(cor_mat_sub) %in% result3.df[which(result3.df$log2FoldChange > 0.5 & result3.df$padj < 0.05),]$gene , ] %>% colMeans
mir7_KO_DOWN = cor_mat_sub[rownames(cor_mat_sub) %in% result3.df[which(result3.df$log2FoldChange < 0.5 & result3.df$padj < 0.05),]$gene , ] %>% colMeans

mir7_KO_UP_WT = cor_mat_sub[rownames(cor_mat_sub) %in% res4.df[which(res4.df$log2FoldChange > 0.5 & res4.df$padj < 0.05),]$gene , ] %>% colMeans
mir7_KO_DOWN_WT = cor_mat_sub[rownames(cor_mat_sub) %in% res4.df[which(res4.df$log2FoldChange < 0.5 & res4.df$padj < 0.05),]$gene , ] %>% colMeans


nGenes = ncol(genes_t.df)
nSamples = nrow(genes_t.df)
moduleColors = netwk[[1]]

MEs0 = moduleEigengenes(genes_t.df, moduleColors)$eigengenes
MEs = orderMEs(MEs0)


datTraits <- data.frame("miR_7_5p" = mir7.df,
                        #miR_122_5p = mir122.df,
                        KO_WT_Up = KO_WT_UP,
                        KO__WT_Down = KO_WT_DOWN,
                        "WTm7oe_WT_Up" = mir7_WT_UP,
                        "WTm7oe_WT_Down" = mir7_WT_DOWN,
                        "KOm7oe_KO_Up" = mir7_KO_UP,
                        "KOm7oe_KO_Down" = mir7_KO_DOWN,
                        "KOm7oe_WT_Up" = mir7_KO_UP_WT,
                        "KOm7oe_WT_Down" = mir7_KO_DOWN_WT)

moduleTraitCor = cor(MEs, datTraits, use= "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

textMatrix= paste(signif(moduleTraitCor, 2), "\n(",
                  signif(moduleTraitPvalue, 1), ")", sep= "")
dim(textMatrix)= dim(moduleTraitCor)
par(mar= c(10, 10, 3, 3))

labeledHeatmap(Matrix= moduleTraitCor,
               xLabels= names(datTraits),
               yLabels= str_remove(names(MEs), "ME"),
               ySymbols= str_remove(names(MEs), "ME"),
               colorLabels= FALSE,
               colors= c(rev(paletteer_d("ggsci::light_blue_material")), paletteer_d("ggsci::orange_material")),
               textMatrix= textMatrix,
               setStdMargins= FALSE,
               cex.text = 1,
               cex.lab.y = 1,
               zlim= c(-1,1),
               main= paste("Module-trait correlation"))
# save images in WGCNA_Analysis_Final.Rmd file

bg_ids <- res4.df$gene  


Run_GO <- function(GO_targets1 = GO_targets2, exp = "KO_WT_Up"){
  
  topgo_BP <- topGOtable(GO_targets1, bg_ids,
                         ontology = "BP",
                         mapping = "org.Mm.eg.db",
                         geneID = "symbol")
  
  saveRDS(topgo_BP, file = paste0("/data/rajewsky/home/skim/Neuronal_Activity_Cledi/201911_m7_overexpression/publication_codes_figures/WGCNA_Analysis/GO_", exp, "_GO_BP.rda"))
  
  topgo_MF <- topGOtable(GO_targets1, bg_ids,
                         ontology = "MF",
                         mapping = "org.Mm.eg.db",
                         geneID = "symbol")
  
  saveRDS(topgo_MF, file = paste0("/data/rajewsky/home/skim/Neuronal_Activity_Cledi/201911_m7_overexpression/publication_codes_figures/WGCNA_Analysis/GO_", exp, "_GO_MF.rda"))
  
  topgo_CC <- topGOtable(GO_targets1, bg_ids,
                         ontology = "CC",
                         mapping = "org.Mm.eg.db",
                         geneID = "symbol")
  
  saveRDS(topgo_CC, file = paste0("/data/rajewsky/home/skim/Neuronal_Activity_Cledi/201911_m7_overexpression/publication_codes_figures/WGCNA_Analysis/GO_", exp, "_GO_CC.rda"))
  
  
}

GO_targets1 <- module_df2[module_df2$color_check %in% "turquoise", ]$gene_id
Run_GO(GO_targets1 = GO_targets1, exp = "turquoise")

GO_targets1 <- module_df2[module_df2$color_check %in% "green", ]$gene_id
Run_GO(GO_targets1 = GO_targets1, exp = "green")

GO_targets1 <- module_df2[module_df2$color_check %in% "yellow", ]$gene_id
Run_GO(GO_targets1 = GO_targets1, exp = "yellow")

GO_targets1 <- module_df2[module_df2$color_check %in% "blue", ]$gene_id
Run_GO(GO_targets1 = GO_targets1, exp = "blue")

GO_targets1 <- module_df2[module_df2$color_check %in% "brown", ]$gene_id
Run_GO(GO_targets1 = GO_targets1, exp = "brown")
