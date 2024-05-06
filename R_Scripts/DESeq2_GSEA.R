

library(fgsea)

####################################
# run fgsea from DESeq object
####################################
# appy DESeq results shrinkage for fgsea

runfgsea <- function(DESeq_object, experiment){

  # shrinkage to reduce variance in low counts
  # e.g. coef = "condition_KO_Ctrl_vs_WT_Ctrl"
  res.shr <- lfcShrink(DESeq_object, coef = experiment, type = "apeglm")
  res.shr.df <- as.data.frame(res.shr)

  # create rnk file
  res.shr.df$Gene <- rownames(res.shr.df)
  res.shr.df$fcsign <- sign(res.shr.df$log2FoldChange)
  res.shr.df$logP = -log10(res.shr.df$pvalue)
  res.shr.df$metric <- res.shr.df$logP/res.shr.df$fcsign
  y <- res.shr.df[,c("Gene", "metric")]

  # Load the pathways into a named list
  pathways.hallmark <- gmtPathways("/data/rajewsky/home/skim/Neuronal_Activity_Cledi/ReactomePathways.gmt")

  # rank file modification 
  res_rank <- y[, c("Gene", "metric")]
  res_rank$Gene <- toupper(res_rank$Gene)
  ranks <- tibble::deframe(res_rank)

  fgseaRes <- fgsea(pathways.hallmark, ranks, minSize=15, maxSize=500, nperm=1000)

  return(fgseaRes)

}


