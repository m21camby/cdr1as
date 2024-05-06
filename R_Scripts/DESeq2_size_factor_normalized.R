
# 1. Calculate size factor
# dds2.ds <- estimateSizeFactors(dds2.ds)

# 2. extract size factor
# WT_size_factor <- sizeFactors(dds2.ds)

# 3. check size factor
# WT_size_factor
# WTC3      WTC4      WTC5      WTm7oex3  WTm7oex4  WTm7oex5 
# 0.9657168 1.1410779 0.6699582 1.2735552 0.6568746 1.6773274 

# 4. Prepare Digital gene expression matrix
# DGEm2 

Calculate_size_factor_normalized <- function(DGEm, size_factor){
  DGEm_normalized <- sweep(DGEm, MARGIN = 2, size_factor, FUN = "/")
  return(DGEm_normalized)
}

# 5. usage example
# DGEm2_normalized <- Calculate_size_factor_normalized(DGEm = DGEm2, size_factor = WT_size_factor)







