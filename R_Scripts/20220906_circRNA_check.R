library(ciRcus)
library(plyr)
library(dplyr)


annot.file = "/scratch/shared_cache/R/ciRcus/mm10.sqlite"
annot.list = loadAnnotation(annot.file)

# Used in Iter_circ which create uniquely circRNAs
Only_Uniq <- function(df){
  names <- c("gene_id", "circ_rna")
  colnames(df) <- names
  DF <- group_by(df, gene_id) # group by gene_id
  DF <- summarise_each(DF, funs(sum)) # sum by gene_id
  DF <- filter(DF, grepl("ENS", gene_id)) # remove intergenic and ambiguous
  DF.name <- deparse(substitute(df))
  assign(DF.name, DF)
  return(data.frame(DF))
}

# Return data frame everything about ciRcus
Iter_circ <- function(i){
  j <- paste("/data/rajewsky/home/skim/Neuronal_Activity_Cledi/201911_m7_overexpression/circRNA_analysis/", i , "/circ_splice_sites.bed", sep = "")
  k <- data.frame(sample=c(i),filename=c(j))
  l <- summarizeCircs(colData = k, wobble=1, keepCols=1:12)
  m <- annotateCircs(l, annot.list, assembly = c("mm10"), fixCoordIndexing = TRUE)
  n <- resTable(m)
  p <- n[ , c(1,2,3,4,5,6, 9)]
  #q <- Only_Uniq(p)
  q_name <- paste(i)
  #colnames(q) <- c("gene_id", q_name)
  #return(q)
  colnames(p) <- c("chr", "start", "end", "width", "strand", "gene_id", "counts")
  p$sample <- q_name
  
  return(p)
}



KOm7oe1 <- Iter_circ("KOm7oex1_S12")
KOm7oe2 <- Iter_circ("KOm7oex2_S14")
KOm7oe3 <- Iter_circ("KOm7oex3_S16")
KOm7oe4 <- Iter_circ("KOm7oex4_S18")
KOm7oe5 <- Iter_circ("KOm7oex5_S20")

KO1 <- Iter_circ("KOC1_S11")
KO2 <- Iter_circ("KOC2_S13")
KO3 <- Iter_circ("KOC3_S15")
KO4 <- Iter_circ("KOC4_S17")
KO5 <- Iter_circ("KOC5_S19")

WT1 <- Iter_circ("WTC1_S1")
WT2 <- Iter_circ("WTC2_S3")
WT3 <- Iter_circ("WTC3_S5")
WT4 <- Iter_circ("WTC4_S7")
WT5 <- Iter_circ("WTC5_S9")

WTm7oe1 <- Iter_circ("WTm7oex1_S2")
WTm7oe2 <- Iter_circ("WTm7oex2_S4")
WTm7oe3 <- Iter_circ("WTm7oex3_S6")
WTm7oe4 <- Iter_circ("WTm7oex4_S8")
WTm7oe5 <- Iter_circ("WTm7oex5_S10")


circ_final <- rbind(KOm7oe1, KOm7oe2, KOm7oe3, KOm7oe4,KOm7oe5,
                    KO1, KO2, KO3, KO4, KO5,
                    WT1, WT2, WT3, WT4, WT5,
                    WTm7oe1, WTm7oe2, WTm7oe3,WTm7oe4, WTm7oe5)


saveRDS(circ_final, file = "/data/rajewsky/home/skim/Neuronal_Activity_Cledi/201911_m7_overexpression/circRNA_analysis/circ_final.rda")
