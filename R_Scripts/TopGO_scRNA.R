# Written by SJK 08.08 2019

# This R script is for DroNc-seq AD project
# This file aims for GO Analysis of scRNA-seq data
# Prerequisites are foreground genes which are interesting genes and background genes
# For background genes, I used FindMarkers function from Seurat with parameter logfc.threshold = 0.01 

# Ref1: https://github.com/karthikshekhar/CellTypeMIMB/blob/master/utilities.R
# Ref2: https://rpubs.com/kshekhar/349874
# Ref3: https://bioconductor.org/packages/release/bioc/vignettes/topGO/inst/doc/topGO.pdf

####################################################
# Modification 24. Jan. 2020
# previous using fg.genes and bg.genes before topGO
# I create fg.genes and bg.genes inside function
# by just provided DE results from Seurat MAST
####################################################

library(topGO)
library(org.Mm.eg.db)
library(ggplot2)

########################old codes################################

#scRNA_topGO <- function(fg.genes = NULL,
#                       bg.genes = NULL,
#                       organism = "Mouse", 
#                       stats.use = "fisher",
#                       algorithm.use = "elim",
#                      topnodes.print=20,
#                       num.char=100){

#if (is.null(fg.genes) | is.null(bg.genes)){
#    stop("Error : Both gene lists are empty")}
##################################################################

scRNA_topGO <- function(DE_results.df, 
                       logFC = 0.2,
                       InterestGenes = "all", 
                       significant = 0.01,
                       organism = "Mouse", 
                       stats.use = "fisher",
                       algorithm.use1 = "elim",
                       algorithm.use2 = "classic",
                       topnodes.print=50,
                       num.char=100){

if((InterestGenes == "up" &&  length(which(DE_results.df$avg_log2FC > logFC & DE_results.df$p_val_adj < 0.01)) > 0) | (InterestGenes == "down" &&  length(which(DE_results.df$avg_log2FC < -logFC & DE_results.df$p_val_adj < 0.01)) > 0) | (InterestGenes == "all" && length(which(abs(DE_results.df$avg_log2FC) > logFC & DE_results.df$p_val_adj < 0.01)) > 0)){

#########################################
# 1. Predefined list of interesting genes
#########################################

####################### old codes ###################################

#n <- length(bg.genes)
#geneList <- integer(n)
#names(geneList) <- bg.genes
#geneList[intersect(names(geneList), fg.genes)] <- 1
#geneList <- factor(geneList)
# geneList object is a named factor that indicates which genes are interesting and which not.
#####################################################################

# gene universe can be all the genes expressed in experiments 
  geneUniverse <- rownames(DE_results.df)

# extract gene of interest
  if(InterestGenes == "up"){
  genesOfInterest <- rownames(DE_results.df[which(DE_results.df$avg_log2FC > logFC & DE_results.df$p_val_adj < significant), ])
  }
  if(InterestGenes == "down"){
  genesOfInterest <- rownames(DE_results.df[which(DE_results.df$avg_log2FC < -logFC & DE_results.df$p_val_adj < significant), ])
  }
  if(InterestGenes == "all"){
  genesOfInterest <- c(rownames(DE_results.df[which(DE_results.df$avg_log2FC > logFC & DE_results.df$p_val_adj < significant), ]), rownames(DE_results.df[which(DE_results.df$avg_log2FC < -logFC & DE_results.df$p_val_adj < significant), ]))}

# export Genes of Interest
  genesOfInterest <<- genesOfInterest

# create gene list
  geneList <- factor(as.integer(geneUniverse %in% genesOfInterest))
  names(geneList) <- geneUniverse


##########################################################################
# 2. Performing test 
# ThetopGOdataobject contains the list of genes, the list of interesting genes and the gene scores
# gene scores (optional) can be p-value for differential expression or correlation with a phenotype
##########################################################################

onts = c( "MF", "BP", "CC" )
tab <- as.list(onts)
names(tab) <- onts

for(i in 1:3){

GOdata <- new("topGOdata",
                description = "GOanalysis",
                ontology = onts[i],
                allGenes = geneList,
                annot = annFUN.org,
                mapping = "org.Mm.eg.db",
                ID = "SYMBOL",
                nodeSize = 20) # nodeSize: prune the GO hierarchy less than nodeSize

# export GO ID from BP
if(i == 1){
  allGO_MF <<- genesInTerm(GOdata)
}

if(i == 2){
  allGO_BP <<- genesInTerm(GOdata)
}
if(i == 3){
  allGO_CC <<- genesInTerm(GOdata)
}



res.result1 <- runTest(GOdata, statistic = stats.use, algorithm = algorithm.use1)
res.result2 <- runTest(GOdata, statistic = stats.use, algorithm = algorithm.use2)


tab[[i]] <- cbind(data.frame(category = onts[i]), GenTable(GOdata, Fisher.elim = res.result1,
                        Fisher.classic = res.result2,
                        orderBy = "Fisher.elim" , topNodes = 30))


}

topGOResults <- plyr::rbind.fill(tab)
topGOResults.df <- as.data.frame(topGOResults)

# Visualization
# showSigOfNodes(GOdata, score(res.result), firstSigNodes = 5, useInfo ='all')

topGOResults.df$gene_ratio <- topGOResults.df$Significant / topGOResults.df$Annotated

# modification appropriate for plot
topGOResults.df$Fisher.elim <- as.numeric(topGOResults.df$Fisher.elim)
topGOResults.df$Fisher.classic <- as.numeric(topGOResults.df$Fisher.classic)
topGOResults.df$Term <- factor(topGOResults.df$Term, levels = rev(unique(topGOResults.df$Term)))

return(topGOResults.df)
}else{
 # if there are no candidate genes to do GO analysis, below message print
 message("no intersting GO")
 
}

}

##############################
# 1-2. edgeR GO function
##############################

scRNA_edgeR_topGO <- function(DE_results.df, 
                        logFC = 0.2,
                        InterestGenes = "all", 
                        significant = 0.01,
                        organism = "Mouse", 
                        stats.use = "fisher",
                        algorithm.use1 = "elim",
                        algorithm.use2 = "classic",
                        topnodes.print=100,
                        num.char=100){
  
  if((InterestGenes == "up" &&  length(which(DE_results.df$logFC > logFC & DE_results.df$FDR < 0.01)) > 0) | (InterestGenes == "down" &&  length(which(DE_results.df$logFC < -logFC & DE_results.df$FDR < 0.01)) > 0) | (InterestGenes == "all" && length(which(abs(DE_results.df$logFC) > logFC & DE_results.df$FDR < 0.01)) > 0)){
    
    #########################################
    # 1. Predefined list of interesting genes
    #########################################
    
    ####################### old codes ###################################
    
    #n <- length(bg.genes)
    #geneList <- integer(n)
    #names(geneList) <- bg.genes
    #geneList[intersect(names(geneList), fg.genes)] <- 1
    #geneList <- factor(geneList)
    # geneList object is a named factor that indicates which genes are interesting and which not.
    #####################################################################
    
    # gene universe can be all the genes expressed in experiments 
    geneUniverse <- rownames(DE_results.df)
    
    # extract gene of interest
    if(InterestGenes == "up"){
      genesOfInterest <- rownames(DE_results.df[which(DE_results.df$logFC > logFC & DE_results.df$FDR < significant), ])
    }
    if(InterestGenes == "down"){
      genesOfInterest <- rownames(DE_results.df[which(DE_results.df$logFC < -logFC & DE_results.df$FDR < significant), ])
    }
    if(InterestGenes == "all"){
      genesOfInterest <- c(rownames(DE_results.df[which(DE_results.df$logFC > logFC & DE_results.df$FDR < significant), ]), rownames(DE_results.df[which(DE_results.df$logFC < -logFC & DE_results.df$FDR < significant), ]))}
    
    # export Genes of Interest
    genesOfInterest <<- genesOfInterest
    
    # create gene list
    geneList <- factor(as.integer(geneUniverse %in% genesOfInterest))
    names(geneList) <- geneUniverse
    
    
    ##########################################################################
    # 2. Performing test 
    # ThetopGOdataobject contains the list of genes, the list of interesting genes and the gene scores
    # gene scores (optional) can be p-value for differential expression or correlation with a phenotype
    ##########################################################################
    
    onts = c( "MF", "BP", "CC" )
    tab <- as.list(onts)
    names(tab) <- onts
    
    for(i in 1:3){
      
      GOdata <- new("topGOdata",
                    description = "GOanalysis",
                    ontology = onts[i],
                    allGenes = geneList,
                    annot = annFUN.org,
                    mapping = "org.Mm.eg.db",
                    ID = "SYMBOL",
                    nodeSize = 20) # nodeSize: prune the GO hierarchy less than nodeSize
      
      # export GO ID from BP
      if(i == 1){
        allGO_MF <<- genesInTerm(GOdata)
      }
      
      if(i == 2){
        allGO_BP <<- genesInTerm(GOdata)
      }
      if(i == 3){
        allGO_CC <<- genesInTerm(GOdata)
      }
      
      
      
      res.result1 <- runTest(GOdata, statistic = stats.use, algorithm = algorithm.use1)
      res.result2 <- runTest(GOdata, statistic = stats.use, algorithm = algorithm.use2)
      
      
      tab[[i]] <- cbind(data.frame(category = onts[i]), GenTable(GOdata, Fisher.elim = res.result1,
                                                                 Fisher.classic = res.result2,
                                                                 orderBy = "Fisher.elim" , topNodes = 100))
      
      
    }
    
    topGOResults <- plyr::rbind.fill(tab)
    topGOResults.df <- as.data.frame(topGOResults)
    
    # Visualization
    # showSigOfNodes(GOdata, score(res.result), firstSigNodes = 5, useInfo ='all')
    
    topGOResults.df$gene_ratio <- topGOResults.df$Significant / topGOResults.df$Annotated
    
    # modification appropriate for plot
    topGOResults.df$Fisher.elim <- as.numeric(topGOResults.df$Fisher.elim)
    topGOResults.df$Fisher.classic <- as.numeric(topGOResults.df$Fisher.classic)
    topGOResults.df$Term <- factor(topGOResults.df$Term, levels = rev(unique(topGOResults.df$Term)))
    
    return(topGOResults.df)
  }else{
    # if there are no candidate genes to do GO analysis, below message print
    message("no intersting GO")
    
  }
  
}



##############################
# topGO results figure 1
##############################

scRNA_topGO_plot <- function(topGOResults.df, color1 = "chocolate2"){
# "dodgerblue3" for downregulated 

topGOResults.df$pval <- as.numeric(topGOResults.df$pval)
topGOResults.df$PValue <- -log10(topGOResults.df$pval)
topGOResults.df$Term <- factor(topGOResults.df$Term, levels = rev(unique(topGOResults.df$Term)))

ggplot(topGOResults.df, aes(x = Term, y = PValue)) + geom_bar(stat="identity", color = color1, fill = color1) + coord_flip() +
  ylab("-log10(p-value)") +
    theme(axis.title.y = element_blank(), axis.text = element_text(color = "black", size = 12), 
        axis.title.x = element_text(size = 12, color = "black"), title = element_text(size = 12),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), 
        axis.line = element_line(colour = "black"), 
        axis.ticks.y = element_blank())

}

##################################
# topGO results figure 2
##################################

scRNA_TopGO_plot2 <- function(topGOResults.df){

topGOResults.df$Fisher.elim <- as.numeric(topGOResults.df$Fisher.elim)
topGOResults.df$Fisher.elim <- -log10(topGOResults.df$Fisher.elim)
topGOResults.df$Term <- factor(topGOResults.df$Term, levels = rev(unique(topGOResults.df$Term)))
  
ggplot(topGOResults.df, aes(x=gene_ratio,
               y=Term,
               colour=Fisher.elim,
               size=Significant)) +
        geom_point() +
        expand_limits(x=0) +
        labs(x="gene ratio", y="GO term", colour="-log10(p-value)", size="Significant") +
  theme_minimal() + theme(axis.text = element_text(size = 10))

}

scRNA_TopGO_plot2 <- function(topGOResults.df, title = "title"){

  topGOResults.df$Fisher.elim <- as.numeric(topGOResults.df$Fisher.elim)
  topGOResults.df$Fisher.elim <- -log10(topGOResults.df$Fisher.elim)
  topGOResults.df$Term <- factor(topGOResults.df$Term, levels = rev(unique(topGOResults.df$Term)))

  ggplot(topGOResults.df, aes(x=gene_ratio,
                              y=Term,
                              colour=Fisher.elim,
                              size=Significant)) +
    geom_point() +
    expand_limits(x=0) +
    labs(x="gene ratio", y="GO term", colour="-log10(p-value)", size="Significant") +
    ggtitle(title) +
    theme_minimal() + theme(axis.text.x = element_text(size = 10, color = "black"), axis.text.y = element_text(size = 10, color = "black"), axis.title = element_text(size = 12, color = "black"), title = element_text(hjust = 0.5)) + guides(
      size = guide_legend(order = 1),
      fill = guide_legend(order = 0)
    )

}

##############################################
# convert genes from GO terms to proper format
##############################################

genes_from_GO_Term <- function(DE.df = DE_table, GO.df = GO_Results){
  final.df <- data.frame()

for(i in c(1:10)){
  genes <- GO.df[i, 11]
  genes2 <- unlist(strsplit(genes, ","))
  genes3 <- str_squish(genes2)
  DE.df1 <- DE.df[DE.df$gene %in% genes3, ]
  DE.df1$pathway <- GO.df[i, 3]
  final.df <- rbind(final.df, DE.df1)
}
  
  final.df$pathway <- factor(final.df$pathway, levels = rev(unique(GO.df$Term)))
  
  return(final.df)
}

##################################
# genes violin plot
##################################

Vln_plot <- function(df = DF){
  ggplot(df, aes(x=pathway, y=avg_log2FC)) +
    geom_violin(trim=FALSE, adjust = 3) + coord_flip() + geom_jitter(data = df[df$p_val_adj < 0.01, ], aes(color = avg_log2FC > 0),  position=position_jitter(0.1), alpha = 0.3) +
    theme(axis.title.y = element_blank()) +
    scale_color_manual(values = c("TRUE" = "red", "FALSE" = "blue")) + theme_minimal() +
    theme(axis.title.y = element_blank(),
          axis.title.x = element_text(size = 10, family = "helvetica", color = "black"),
          axis.text = element_text(size = 8, family = "helvetica", color = "black"),
          panel.grid.minor = element_blank(),
          legend.position = "none") + 
    geom_text_repel(data = df, aes(label = df$gene), size = 3.5)
  
}


