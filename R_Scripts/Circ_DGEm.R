# ciRcus workflow

library(ciRcus)
library(rtracklayer)
library(readxl)
library(biomaRt)



# This function is for creating table of total number of BSJ and genes give rise to BSJ

Circ_Overview <- function(Sample_ID, PATH = PATH, ...){

  # Loading data
  PATH <- PATH
  File <- paste(PATH, Sample_ID , "/circ_splice_sites.bed", sep = "")
  File <- data.frame(sample=c(Sample_ID),filename=c(File))

  # Create empty DataFrame
  Result_DF <- data.frame(ID = character(), desc = character(), counts = integer(), qual = character(), stringsAsFactors = FALSE)

  # qualfilter 
  File_circ <- summarizeCircs(colData = File, wobble=1, keepCols=1:12, qualfilter = TRUE)
  File_circ <- annotateCircs(File_circ, annot.list, assembly = c("mm10"), fixCoordIndexing = TRUE)
  File_circ.df <- resTable(File_circ)

  Result_DF[1,1] <- Sample_ID 
  Result_DF[1,2] <- "Total_counts"
  Result_DF[1,3] <- sum(File_circ.df[, 9])
  Result_DF[1,4] <- "TRUE"

  Result_DF[2,1] <- Sample_ID
  Result_DF[2,2] <- "Total_genes"
  Result_DF[2,3] <- nrow(File_circ.df)
  Result_DF[2,4] <- "TRUE"

  # no qualfilter 
  File_circ <- summarizeCircs(colData = File, wobble=1, keepCols=1:12, qualfilter = FALSE)
  File_circ <- annotateCircs(File_circ, annot.list, assembly = c("mm10"), fixCoordIndexing = TRUE)
  File_circ.df <- resTable(File_circ)

  Result_DF[3,1] <- Sample_ID
  Result_DF[3,2] <- "Total_counts"
  Result_DF[3,3] <- sum(File_circ.df[, 9])
  Result_DF[3,4] <- "FALSE"

  Result_DF[4,1] <- Sample_ID
  Result_DF[4,2] <- "Total_genes"
  Result_DF[4,3] <- nrow(File_circ.df)
  Result_DF[4,4] <- "FALSE"


  return(Result_DF)
}

Circ_transcripts_Matrix_mm10 <- function(Sample_ID, PATH = PATH1, ...){

  # Loading data
  PATH <- PATH
  File <- paste(PATH, Sample_ID , "/circ_splice_sites.bed", sep = "")
  File <- data.frame(sample=c(Sample_ID),filename=c(File))

  # Follow ciRcus flow
  File_circ <- summarizeCircs(colData = File, wobble=1, keepCols=1:12, ...)
  File_circ <- annotateCircs(File_circ, annot.list, assembly = c("mm10"), fixCoordIndexing = TRUE)
  File_circ.df <- resTable(File_circ)

  # Cdr1 annotate
  File_circ.df$gene_id <- ifelse(grepl("intergenic", File_circ.df$gene_id) & grepl("^61", File_circ.df$start) & grepl("chrX", File_circ.df$chr) & grepl("+", File_circ.df$strand), "ENSMUSG00000090546", File_circ.df$gene_id)

  File_circ.df$start <- as.character(File_circ.df$start)
  File_circ.df$end <- as.character(File_circ.df$end)
  # Combine to create ID
  File_circ.df$transcripts <- apply(File_circ.df[ , c(6,1,2,3,5)] , 1 , paste , collapse = "_" )

  # Extract CDS circ
  File_circ.df <- File_circ.df[grepl("ENS*", File_circ.df$transcripts), ]
  # Extract two columns ID and counts
  File_circ.df <- File_circ.df[, c(14,9)]

  return(File_circ.df)
}

Circbase_mm10 <- function(){

  # From circbase excel download
  mmu_mm9 <- read_excel("/data/rajewsky/home/skim/organoids_Aga/mmu_mm9_Rybak2015.xlsx")
  mmu_mm9_2nd <- mmu_mm9[, c(1,2,3,4,5,11,12)]
  colnames(mmu_mm9_2nd) <- c("chr","start", "end","strand","circRNA_ID","transcript","gene_symbol")
  mmu_mm9.gr <- makeGRangesFromDataFrame(mmu_mm9_2nd, keep.extra.columns = TRUE)

  # for liftover from mm9 to mm10 
  ch = import.chain("/data/rajewsky/home/skim/bash_Script/mm9ToMm10.over.chain")

  seqlevelsStyle(mmu_mm9.gr) = "UCSC" 

  mm10.gr = liftOver(mmu_mm9.gr, ch)
  mm10.gr <- unlist(mm10.gr)

  mm10.df <- data.frame(seq = seqnames(mm10.gr), start = start(mm10.gr)+1, end = end(mm10.gr), strand = strand(mm10.gr), circbase = mm10.gr$circRNA_ID, transcripts = mm10.gr$transcript, gene_symbol = mm10.gr$gene_symbol)
  # Create for combining circbase and our data (circbase)
  mm10.df$comb <- paste0(mm10.df$seq, "_", mm10.df$start, "_", mm10.df$end,"_", mm10.df$strand)
  
  return(mm10.df)
}

Circ_trancripts_All_mm10_Samples <- function(exp_list, PATH = PATH){

  exp1 <- Circ_transcripts_Matrix_mm10(exp_list[1], PATH = PATH)
  
  exp_list2 <- exp_list[-1]

  for(i in exp_list2){
    my_name <- i
    file_name <- Circ_transcripts_Matrix_mm10(my_name, PATH = PATH)
    exp1 <- list(exp1, file_name) %>% purrr::reduce(full_join, by = "transcripts")
    
  }
  
  exp1[is.na(exp1)] <- 0
  return(exp1)
}

Comb_transcripts_circbase <- function(circ_transcripts_DGEm, mm10.df){
  # extract ID
  circ_transcripts_DGEm2 <- (circ_transcripts_DGEm %>% separate(transcripts, c("ensembl_gene_id", "chr", "start", "end", "strand"), "_"))[, c(1:5)]
  # combine ID
  circ_transcripts_DGEm2$comb <- paste0(circ_transcripts_DGEm2$chr, "_", circ_transcripts_DGEm2$start, "_", circ_transcripts_DGEm2$end,"_", circ_transcripts_DGEm2$strand)
  # extract information
  mm10.df_2nd <- mm10.df[,c(5,7,8)]
  # join together
  circ_transcripts_DGEm3 <- left_join(circ_transcripts_DGEm2, mm10.df_2nd, by = "comb")

  return(circ_transcripts_DGEm3)
}

# extract unidentified circ from circbase
Extract_Novel_circ_DF <- function(circ_transcripts_DGEm3){

  circ_transcripts_DGEm_NA <- circ_transcripts_DGEm3[is.na(circ_transcripts_DGEm3$circbase), ]
  colnames(circ_transcripts_DGEm_NA) <- c("ensembl_gene_id", "chr", "start", "end", "strand", "comb", "circbase", "gene_symbol")
  return(circ_transcripts_DGEm_NA)
}

# extract unidentified circ from circbase and convert to GenomicRanges
Extract_Novel_circ_GR <- function(circ_transcripts_DGEm_NA){
  circ_transcripts_DGEm_NA.gr <- makeGRangesFromDataFrame(circ_transcripts_DGEm_NA)
  return(circ_transcripts_DGEm_NA.gr)
}

# extract information of unidentified circ from circbase and convert to GenomicRanges 
Extract_Novel_circ_info <- function(circ_transcritps_DGEm_NA){
  
  # From biomart
  mart <- useDataset("mmusculus_gene_ensembl", useMart("ensembl"))
  #listAttributes(mart = mart)
  circ_Info.df <- getBM(attributes = c('ensembl_gene_id', 'chromosome_name', 'start_position',
                                     'end_position', 'strand','mgi_symbol'),
                    filters = 'ensembl_gene_id',
                    values = circ_transcripts_DGEm_NA$ensembl_gene_id,
                    mart = mart)

  # remove overlapping rows
  circ_Info.df <- circ_Info.df[!duplicated(circ_Info.df), ] #2560


  # Add chr to chromosome
  circ_Info.df$chromosome_name <- paste0("chr", circ_Info.df$chromosome_name)
  # Plus and Minus change
  circ_Info.df[circ_Info.df == 1] <- "+"
  circ_Info.df[circ_Info.df == -1] <- "-"
  colnames(circ_Info.df) <- c("ensembl_gene_id", "chr", "start", "end", "strand", 'mgi_symbol')

  # Change Cdr1as information (- to +)
  circ_Info.df[circ_Info.df$mgi_symbol == "Cdr1", ]$strand <- "+"
  
  # Remove duplicated ensembl_gene_id and mgi_symbol
  circ_Info.df <- circ_Info.df[!duplicated(circ_Info.df$ensembl_gene_id) & !duplicated(circ_Info.df$mgi_symbol), ]  
   
  #circ_Info.gr <- makeGRangesFromDataFrame(circ_Info.df)
  
  return(circ_Info.df)

}

# Identify novel circ information

Extract_Novel_circ_info_Final <- function(circ_transcripts_DGEm3){
   
  # Extract novel circ 
  circ_transcripts_DGEm_NA <- Extract_Novel_circ_DF(circ_transcripts_DGEm3)
  
  print(circ_transcripts_DGEm_NA[c(1:2), c(1:2)])
  # convert to GR
  circ_transcripts_DGEm_NA.gr <- Extract_Novel_circ_GR(circ_transcripts_DGEm_NA)
  print(circ_transcripts_DGEm_NA[c(1:2), c(1:2)])
  # novel circ information
  circ_Info.df <- Extract_Novel_circ_info(circ_transcripts_DGEm_NA)
  # Test
  print(nrow(circ_Info.df))

  circ_Info.gr <- makeGRangesFromDataFrame(circ_Info.df)

  overlap_results <- findOverlaps(circ_transcripts_DGEm_NA.gr, circ_Info.gr, type = "within")
  overlap_results.df <- as.data.frame(overlap_results)
  # Test
  #print(nrow(overlap_results.df))

  circ_transcripts_DGEm_NA$queryHits <- seq(1,nrow(circ_transcripts_DGEm_NA))
  overlap_results.df2 <- right_join(overlap_results.df, circ_transcripts_DGEm_NA, by = "queryHits")

  # Test
  #print(nrow(overlap_results.df2))
 
  circ_Info.df$subjectHits <- seq(1,nrow(circ_Info.df))
  overlap_results.df3 <- left_join(overlap_results.df2, circ_Info.df, by = "subjectHits")
  
  # Test
  #print(paste("3rd", nrow(overlap_results.df3)))

  overlap_results.df3 <- overlap_results.df3[, c(8,11,16)] 
  colnames(overlap_results.df3) <- c("comb", "ensembl_gend_id", "mgi_symbol")
  # Remove duplicated ensembl_gene_id and mgi_symbol
  overlap_results.df3 <- overlap_results.df3[!is.na(overlap_results.df3$mgi_symbol), ]
  # Remove No name
  overlap_results.df3 <- overlap_results.df3[!overlap_results.df3$mgi_symbol == "", ]

  return(overlap_results.df3)
}


# Novel circ identification (Naming none existence of circ in circbase)

Naming_none_circbase <- function(overlap_results.df3){

  overlap_results.df3_summary <- as.data.frame(table(overlap_results.df3$mgi_symbol))

  # 1st row circ
  c <- overlap_results.df3_summary[1 , ]$Var1
  DF <- filter(overlap_results.df3, mgi_symbol == c)
  DF$number <- seq(1, overlap_results.df3_summary[1 , ]$Freq)

  # rbind rest circ
  for(i in 2:nrow(overlap_results.df3_summary)){
  c <- overlap_results.df3_summary[i , ]$Var1
  DF2 <- filter(overlap_results.df3, mgi_symbol == c)
  DF2$number <- seq(1,overlap_results.df3_summary[i , ]$Freq)
  DF <- rbind(DF, DF2)
  }

  DF$circbase <- paste0(DF$mgi_symbol, "_", "novel", "_", DF$number)
  DF$number <- NULL

  return(DF)
}

Finalized_circ_table <- function(circ_transcripts_DGEm3, DF){
  
  circ_transcripts_DGEm3$circbase <- as.character(circ_transcripts_DGEm3$circbase)
  circ_transcripts_DGEm3$gene_symbol <- as.character(circ_transcripts_DGEm3$gene_symbol)

  for(i in 1:nrow(DF)){
    c <- which(grepl(DF[i, ]$comb, circ_transcripts_DGEm3$comb), arr.ind = TRUE) 
    circ_transcripts_DGEm3[c, ]$gene_symbol <- DF[i, ]$mgi_symbol
    circ_transcripts_DGEm3[c, ]$circbase <- DF[i, ]$circbase
  
  }

  circ_transcripts_DGEm3 <- circ_transcripts_DGEm3[!is.na(circ_transcripts_DGEm3$gene_symbol), ]
  return(circ_transcripts_DGEm3)
}

Finalized_circ_DGEm <- function(cDGE, circ_transcripts_DGEm4){
  
  circ_transcripts_DGEm4$transcripts <- paste0(circ_transcripts_DGEm4$ensembl_gene_id, "_", circ_transcripts_DGEm4$comb)

  circ_transcripts_DGEm5 <- circ_transcripts_DGEm4[, c("circbase", "transcripts")]
  cDGE2 <- inner_join(cDGE, circ_transcripts_DGEm5, by = "transcripts")
  rownames(cDGE2) <- cDGE2$circbase
  cDGE2$transcripts <- NULL
  cDGE2$circbase <- NULL
  
  return(cDGE2)
} 


#############################
# Extract linear counts 
#############################
# Below code is for linear counts matrix which is similar to circ count extraction with minor moidification

######################################
# linear count matrix for each sample
######################################
Circ_transcripts_linear_Matrix_mm10 <- function(Sample_ID, PATH = PATH1, ...){

  # Loading data
  PATH <- PATH
  File <- paste(PATH, Sample_ID , "/circ_splice_sites.bed", sep = "")
  File <- data.frame(sample=c(Sample_ID),filename=c(File))

  # Follow ciRcus flow
  File_circ <- summarizeCircs(colData = File, wobble=1, keepCols=1:12, ...)
  File_circ <- annotateCircs(File_circ, annot.list, assembly = c("mm10"), fixCoordIndexing = TRUE)
  File_circ.df <- resTable(File_circ)

  # Cdr1 annotate
  File_circ.df$gene_id <- ifelse(grepl("intergenic", File_circ.df$gene_id) & grepl("^61", File_circ.df$start) & grepl("chrX", File_circ.df$chr) & grepl("+", File_circ.df$strand), "ENSMUSG00000090546", File_circ.df$gene_id)

  File_circ.df$start <- as.character(File_circ.df$start)
  File_circ.df$end <- as.character(File_circ.df$end)
  # Combine to create ID
  File_circ.df$transcripts <- apply(File_circ.df[ , c(6,1,2,3,5)] , 1 , paste , collapse = "_" )

  # Extract CDS circ
  File_circ.df <- File_circ.df[grepl("ENS*", File_circ.df$transcripts), ]
  # Extract two columns ID and counts
 File_circ.df <- File_circ.df[, c(14,11,12)]

  File_circ.df2 <- data.frame(transcripts = as.character(), linear = as.numeric(), stringsAsFactors = FALSE)

  for (i in 1:nrow(File_circ.df)){
    File_circ.df2[i, 1] <- File_circ.df[i,1]
    File_circ.df2[i, 2] <- ifelse(File_circ.df[i, 2] > File_circ.df[i, 3], File_circ.df[i, 2] , File_circ.df[i, 3])
}
  return(File_circ.df2)

}

##################################################
# linear counts matrix for all sample by iteration
##################################################
Circ_trancripts_linear_All_mm10_Samples <- function(exp_list, PATH = PATH){

  exp1 <- Circ_transcripts_linear_Matrix_mm10(exp_list[1], PATH = PATH)

  exp_list2 <- exp_list[-1]

  for(i in exp_list2){
    my_name <- i
    file_name <- Circ_transcripts_linear_Matrix_mm10(my_name, PATH = PATH)
    exp1 <- list(exp1, file_name) %>% purrr::reduce(full_join, by = "transcripts")

  }

  exp1[is.na(exp1)] <- 0
  return(exp1)
}










