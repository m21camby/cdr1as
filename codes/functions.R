
#############################
# Normalized counts box plot
#############################

box_plot <- function(df = DGEm_normalized, gene = "Kdm5d"){
  sub.df <- df[rownames(df) %in% gene, ] %>% t %>% as.data.frame()
  colnames(sub.df) <- "gene"
  sub.df$sample <- c(rep("WT",3), rep("WTm7oe",3 ), rep("KO",4), rep("KOm7oe",4 ))
  sub.df$sample <- factor(sub.df$sample, levels = c("WT","WTm7oe","KO", "KOm7oe"))
  sub.df$rep <- str_sub(rownames(sub.df), -1, -1)
  g1 <- ggplot(sub.df, aes(x  = sample, y = gene, label = rep)) + geom_boxplot() + geom_jitter(shape=1, position=position_jitter(0.2), color = "darkorange", alpha = 0.7) + 
    scale_color_identity() + theme_cowplot() + geom_label() + ggtitle(gene)
  g1
}

#########################
# Run Random Forest
##########################

Run_RF <- function(interesting_gene = pheno_genes[1], 
                   outdir = "/data/rajewsky/home/skim/Neuronal_Activity_Cledi/201911_m7_overexpression/R_Scripts/20220420_Elisabeth_code_check/",
                   layer = "low"){
  
  
  dir.create(paste0(outdir, interesting_gene))
  
  evaluationAll <- matrix(NA, length(ugroup), 2)
  colnames(evaluationAll) <- c("RMSE", "Cosine similarity")
  rownames(evaluationAll) <- c(1:14)
  
  xRaw <- r[, !colnames(r) %in% interesting_gene ]
  if(layer == "low"){
    xRaw <- xRaw[, !colnames(xRaw) %in% c(mir7_genes, pheno_genes) ]
  }
  if(layer == "high"){
    xRaw <- xRaw[, !colnames(xRaw) %in% c(middle_genes, pheno_genes) ]
  }
  
  yRaw <- r[, interesting_gene]
  
  
  for (i in 1:length(ugroup)){
    #print(i)
    x <- xRaw
    y <- yRaw
    
    wtest <- which(group %in% c(ugroup[i]))
    
    
    print("start training")
    model <- randomForest(x[-wtest, , drop=F], y[-wtest], ntree=1000, nodesize=2, importance=T, keep.forest=T)
    
    saveRDS(model, paste0("/data/rajewsky/home/skim/Neuronal_Activity_Cledi/201911_m7_overexpression/R_Scripts/20220420_Elisabeth_code_check/",interesting_gene,"/", interesting_gene ,"_" , i,".rda"))
    pred <- predict(model, x[wtest, , drop=F])
    evaluation <- c(computeRMSE(pred, y[wtest]), computePairCorrelation(pred, y[wtest], y[-wtest]))
    evaluationAll[i,] <- evaluation
    
  }
  saveRDS(evaluationAll, paste0("/data/rajewsky/home/skim/Neuronal_Activity_Cledi/201911_m7_overexpression/R_Scripts/20220420_Elisabeth_code_check/",interesting_gene,"/", interesting_gene ,"_" , "evaluation",".rda"))
  
}

################################
# RMSE function in RF function
################################

computeRMSE <- function(pred, y){
  r <- sqrt(mean((pred-y)^2))
  return(r)
}

################################
# cosine similarity function in RF function
################################

computePairCorrelation <- function(pred, y, ytrain){
  ytm <- matrix(rep(ytrain, length(y)), length(y), length(ytrain), byrow=T) # rep not necessary, but make intention clearer
  gt <- c(matrix(rep(y, length(ytrain)), length(y), length(ytrain), byrow=F)-ytm) # ground truth # c() makes column-wise concatenation 
  pr <- c(matrix(rep(pred, length(ytrain)), length(pred), length(ytrain), byrow=F)-ytm)
  r <- (pr %*% gt)/(norm(pr, type="2")*norm(gt, type="2")) # cosine similarity # −1 meaning exactly opposite, to 1 meaning exactly the same, V and aV are maximally similar for any constant a
  return(r)
}

