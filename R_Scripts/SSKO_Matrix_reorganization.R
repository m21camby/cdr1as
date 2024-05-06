
library(DESeq2)
library(limma)

New_DGEm <- read.csv("/data/rajewsky/home/skim/Neuronal_Activity_Cledi/202207_m7_oex_5ssKO/final_DF_stranded.csv", sep = "\t", row.names = 1)

colnames(New_DGEm) <- c("SSKOm7oe4", "WTJ6", "WTJ4", "WTJ2", 
                        "SSKO4", "WTJm7oe2", "WTJ5", "WTJm7oe1",
                        "WTJm7oe6", "SSKO2", "WTJ3", "WTm7oe5",
                        "WTJm7oe5", "WTC5", "SSKOm7oe2", "SSKOm7oe1",
                        "WTJm7oe4", "SSKOm7oe3", "WTm7oe4", "WTJm7oe3",
                        "WTJ1","SSKO3","WTC4","SSKO1")


SS_DGEm <- New_DGEm[, c(21,11,3,
                        8,20,17,
                        10,22,5,
                        15,18,1)]


colnames(SS_DGEm) <- c("WTJ1","WTJ2","WTJ3",
                       "WTJm7oe1", "WTJm7oe2", "WTJm7oe3",
                       "SSKO1" ,"SSKO2" ,"SSKO3",
                       "SSKOm7oe1","SSKOm7oe2","SSKOm7oe3")

saveRDS(SS_DGEm, 
        file = "/data/rajewsky/home/skim/Neuronal_Activity_Cledi/202207_m7_oex_5ssKO/publication_code/SSKO_Matrix.rda")


