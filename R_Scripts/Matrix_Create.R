
New_DGEm <- read.csv("/data/rajewsky/home/skim/Neuronal_Activity_Cledi/202207_m7_oex_5ssKO/final_DF_stranded.csv", sep = "\t", row.names = 1)

colnames(New_DGEm) <- c("SSKOm7oe4", "WTJ6", "WTJ4", "WTJ2", 
                        "SSKO4", "WTJm7oe2", "WTJ5", "WTJm7oe1",
                        "WTJm7oe6", "SSKO2", "WTJ3", "WTm7oe5",
                        "WTJm7oe5", "WTC5", "SSKOm7oe2", "SSKOm7oe1",
                        "WTJm7oe4", "SSKOm7oe3", "WTm7oe4", "WTJm7oe3",
                        "WTJ1","SSKO3","WTC4","SSKO1")

DGEm <- readRDS(file = "/data/rajewsky/home/skim/Neuronal_Activity_Cledi/201911_m7_overexpression/publication_codes_figures/DGE_matrix.rda")


all_DGEm <- cbind(DGEm, New_DGEm[,c(12,14)])

all_DGEm <- all_DGEm[,c(1,2,3,16,
                        4,5,6,15,
                        7,8,9,10,
                        11,12,13,14)]

colnames(all_DGEm) <- c("WTC1", "WTC2", "WTC3", "WTC4",
                        "WTm7oe1", "WTm7oe2", "WTm7oe3", "WTm7oe4",
                        "KOC1",    "KOC2", "KOC3","KOC4",
                        "KOm7oe1", "KOm7oe2", "KOm7oe3", "KOm7oe4")

saveRDS(all_DGEm, file = "/data/rajewsky/projects/cdr1as_ko_snRNA/codes_github/cdr1as/DGEm/WTN_DGEm.rda")


