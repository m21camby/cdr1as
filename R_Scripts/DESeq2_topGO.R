# Created by SJK 22. Jan 2020
# This file is for analyzing GO analysis of DESeq2 results 


library(topGO)
library(org.Mm.eg.db)
library(clusterProfiler)

###############################################
# This function is for GO analysis using topGO 
###############################################

runTopGO <- function(DESeq2_results.df, InterestGenes = "all", significant = 0.01){
  
  # gene universe can be all the genes expressed in experiments 
  geneUniverse <- rownames(DESeq2_results.df) 

  # extract gene of interest
  if(InterestGenes == "up"){
  genesOfInterest <- rownames(DESeq2_results.df[which(DESeq2_results.df$log2FoldChange > 0.5 & DESeq2_results.df$padj < significant), ])  
  }
  if(InterestGenes == "down"){
  genesOfInterest <- rownames(DESeq2_results.df[which(DESeq2_results.df$log2FoldChange < -0.5 & DESeq2_results.df$padj < significant), ])
  }
  if(InterestGenes == "all"){
  genesOfInterest <- c(rownames(DESeq2_results.df[which(DESeq2_results.df$log2FoldChange > 0.5 & DESeq2_results.df$padj < significant), ]), rownames(DESeq2_results.df[which(DESeq2_results.df$log2FoldChange < -0.5 & DESeq2_results.df$padj < significant), ]))
  }

  # export Genes of Interest
  genesOfInterest <<- genesOfInterest
  
  # create gene list
  geneList <- factor(as.integer(geneUniverse %in% genesOfInterest))
  names(geneList) <- geneUniverse  

  # create topGOdata object run GO analysis
  onts = c( "MF", "BP", "CC" )
  tab <- as.list(onts)
  names(tab) <- onts

  for(i in 1:3){
  sampleGOdata <- new("topGOdata",
                    description = "GO_analysis", 
                    ontology = onts[i],
                    allGenes = geneList, 
                    nodeSize = 10,
                    annot=annFUN.org, mapping="org.Mm.eg.db", ID = "symbol")

  # export GO ID from BP
  if(i == 2){
  allGO <<- genesInTerm(sampleGOdata)
  }
  
  # run tests
  resultTopGO.elim <- runTest(sampleGOdata, algorithm = "elim", statistic = "Fisher" )
  resultTopGO.classic <- runTest(sampleGOdata, algorithm = "classic", statistic = "Fisher" )

  # save to table
   ## look at results
  tab[[i]] <- cbind(data.frame(category = onts[i]), GenTable(sampleGOdata, Fisher.elim = resultTopGO.elim, 
                        Fisher.classic = resultTopGO.classic,
                        orderBy = "Fisher.elim" , topNodes = 20))

}
  # save to dataframe
  topGOResults <- plyr::rbind.fill(tab)
  topGOResults.df <- as.data.frame(topGOResults)

  # calculate gene ratio
  # the ratio be-tween significantly differentially regulated genes inthe respective gene set measured by sequencingand the total number of genes belonging in that GOcategory. 
  # ref: https://www.cell.com/cell-reports/pdfExtended/S2211-1247(18)31242-7 
  topGOResults.df$gene_ratio <- topGOResults.df$Significant / topGOResults.df$Annotated
  topGOResults.df$gene_ratio <- round(topGOResults.df$gene_ratio, 4)
  # modification appropriate for plot
  topGOResults.df$Fisher.elim <- as.numeric(topGOResults.df$Fisher.elim)
  topGOResults.df$Fisher.classic <- as.numeric(topGOResults.df$Fisher.classic)
  topGOResults.df$Term <- factor(topGOResults.df$Term, levels = rev(unique(topGOResults.df$Term)))

  return(topGOResults.df)
}

###################################################
# This function is for figure of GO results
###################################################

figureTopGO <- function(res_GO.df){
  ggplot(res_GO.df, aes(x=gene_ratio, 
               y=Term, 
               colour=Fisher.elim, 
               size=Significant)) +
        geom_point() +
        expand_limits(x=0) +
        labs(x="gene ratio", y="GO term", colour="p value", size="Significant") + 
  theme_minimal() + theme(axis.text = element_text(size = 10))

}

#############################
# enrichGO function
#############################
runenrichGO <- function(DESeq2_results.df, InterestGenes = "all", significant = 0.01){
   # gene universe can be all the genes expressed in experiments 
  geneUniverse <- rownames(DESeq2_results.df)

  # extract gene of interest
  if(InterestGenes == "up"){
  genesOfInterest <- rownames(DESeq2_results.df[which(DESeq2_results.df$log2FoldChange > 0.5 & DESeq2_results.df$padj < significant), ])
  }
  if(InterestGenes == "down"){
  genesOfInterest <- rownames(DESeq2_results.df[which(DESeq2_results$log2FoldChange < -0.5 & DESeq2_results.df$padj < significant), ])
  }
  if(InterestGenes == "all"){
  genesOfInterest <- c(rownames(DESeq2_results.df[which(DESeq2_results.df$log2FoldChange > 0.5 & DESeq2_results.df$padj < significant), ]), rownames(DESeq2_results.df[which(DESeq2_results.df$log2FoldChange < -0.5 & DESeq2_results.df$padj < significant), ]))
  }

  # create gene list
  geneList <- factor(as.integer(geneUniverse %in% genesOfInterest))
  names(geneList) <- geneUniverse
 
  enrichGO.df <- enrichGO(gene = genesOfInterest, universe = names(gene_list),
                      OrgDb = org.Mm.eg.db, 
                      keyType = 'SYMBOL',
                      readable = T,
                      ont = "BP",
                      pvalueCutoff = 0.01, 
                      qvalueCutoff = 0.10)
}


############################
# references
############################

# ref: https://www.bioconductor.org/packages/devel/bioc/vignettes/topGO/inst/doc/topGO.pdf
# ref: http://avrilomics.blogspot.com/2015/07/using-topgo-to-test-for-go-term.html
# ref: https://www.huber.embl.de/users/klaus/Teaching/DESeq2Predoc2014.html#gene-ontology-enrichment-analysis
# ref: https://www.cell.com/cell-reports/pdfExtended/S2211-1247(18)31242-7


