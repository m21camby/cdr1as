.libPaths(c("/data/rajewsky/home/skim/R/usr_lib_Seurat/", "/data/rajewsky/home/skim/R/usr_lib_v4/"))
#.libPaths()
unloadNamespace("mgcv")
unloadNamespace("Matrix")


library(ggplot2)
library(gridExtra)
library(ggrepel)
library(dplyr)
library(cowplot)
library(paletteer)
library(stringr)
library(scales)

cor_mat <- readRDS("/data/rajewsky/projects/cdr1as_ko_snRNA/codes_github/cdr1as/DGEm/WTN_batch_corr_nested_DGEm.rda")

res1_bp_up <- readRDS(file = paste0("/data/rajewsky/projects/cdr1as_ko_snRNA/codes_github/cdr1as/Results/", "GO_Analysis_","KO_WT_UP", "_BP.rda"))
res1_bp_up <- res1_bp_up[which(res1_bp_up$p.value_elim < 0.05 & res1_bp_up$Annotated < 1000 & res1_bp_up$Annotated > 30), ] 
res1_bp_down <- readRDS(file = paste0("/data/rajewsky/projects/cdr1as_ko_snRNA/codes_github/cdr1as/Results/", "GO_Analysis_","KO_WT_DOWN", "_BP.rda"))
res1_bp_down <- res1_bp_down[which(res1_bp_down$p.value_elim < 0.05 & res1_bp_down$Annotated < 1000 & res1_bp_down$Annotated > 30), ] 

res2_bp_up <- readRDS(file = paste0("/data/rajewsky/projects/cdr1as_ko_snRNA/codes_github/cdr1as/Results/", "GO_Analysis_","WTm7oe_WT_UP", "_BP.rda"))
res2_bp_up <- res2_bp_up[which(res2_bp_up$p.value_elim < 0.05 & res2_bp_up$Annotated < 1000 & res2_bp_up$Annotated > 30), ] 
res2_bp_down <- readRDS(file = paste0("/data/rajewsky/projects/cdr1as_ko_snRNA/codes_github/cdr1as/Results/", "GO_Analysis_","WTm7oe_WT_DOWN", "_BP.rda"))
res2_bp_down <- res2_bp_down[which(res2_bp_down$p.value_elim < 0.05 & res2_bp_down$Annotated < 1000 & res2_bp_down$Annotated > 30), ] 

res3_bp_up <- readRDS(file = paste0("/data/rajewsky/projects/cdr1as_ko_snRNA/codes_github/cdr1as/Results/", "GO_Analysis_","KOm7oe_KO_UP", "_BP.rda"))
res3_bp_up <- res3_bp_up[which(res3_bp_up$p.value_elim < 0.05 & res3_bp_up$Annotated < 1000 & res3_bp_up$Annotated > 30), ] 
res3_bp_down <- readRDS(file = paste0("/data/rajewsky/projects/cdr1as_ko_snRNA/codes_github/cdr1as/Results/", "GO_Analysis_","KOm7oe_KO_DOWN", "_BP.rda"))
res3_bp_down <- res3_bp_down[which(res3_bp_down$p.value_elim < 0.05 & res3_bp_down$Annotated < 1000 & res3_bp_down$Annotated > 30), ] 

res4_bp_up <- readRDS(file = paste0("/data/rajewsky/projects/cdr1as_ko_snRNA/codes_github/cdr1as/Results/", "GO_Analysis_","KOm7oe_WT_UP", "_BP.rda"))
res4_bp_up <- res4_bp_up[which(res4_bp_up$p.value_elim < 0.05 & res4_bp_up$Annotated < 1000 & res4_bp_up$Annotated > 30), ] 
res4_bp_down <- readRDS(file = paste0("/data/rajewsky/projects/cdr1as_ko_snRNA/codes_github/cdr1as/Results/", "GO_Analysis_","KOm7oe_WT_DOWN", "_BP.rda"))
res4_bp_down <- res4_bp_down[which(res4_bp_down$p.value_elim < 0.05 & res4_bp_down$Annotated < 1000 & res4_bp_down$Annotated > 30), ] 



Terms <- c("innate immune response", 'positive regulation of secretion by cell',
           "extracellular matrix organization", "regulation of axonogenesis",
           "G protein-coupled receptor signaling pathway", "positive regulation of lipid transport",
           "regulation of glucose transmembrane transport", "histone acetylation",
           "SMAD protein signal transduction",
           "regulated exocytosis", "regulation of calcium ion-dependent exocytosis",
           "behavioral fear response", "startle response", 
           "regulation of postsynaptic membrane potential","regulation of synaptic transmission, glutamatergic", "regulation of synaptic transmission, GABAergic")


res1_bp_down_sub <- res1_bp_down[res1_bp_down$Term %in% Terms, ]
res1_bp_down_sub$exp <- "KO_vs_WT_Down"

res2_bp_down_sub <- res2_bp_down[res2_bp_down$Term %in% Terms, ]
res2_bp_down_sub$exp <- "WTm7oe_vs_WT_Down"

res3_bp_down_sub <- res3_bp_down[res3_bp_down$Term %in% Terms, ]
res3_bp_down_sub$exp <- "KOm7oe_vs_KO_Down"

res4_bp_down_sub <- res4_bp_down[res4_bp_down$Term %in% Terms, ]
res4_bp_down_sub$exp <- "KOm7oe_vs_WT_Down"

res1_bp_up_sub <- res1_bp_up[res1_bp_up$Term %in% c(Terms), ]
res1_bp_up_sub$exp <- "KO_vs_WT_Up"

res2_bp_up_sub <- res2_bp_up[res2_bp_up$Term %in% c(Terms), ]
res2_bp_up_sub$exp <- "WTm7oe_vs_WT_Up"

res3_bp_up_sub <- res3_bp_up[res3_bp_up$Term %in% Terms, ]
res3_bp_up_sub$exp <- "KOm7oe_vs_KO_Up"

res4_bp_up_sub <- res4_bp_up[res4_bp_up$Term %in% Terms, ]
res4_bp_up_sub$exp <- "KOm7oe_vs_WT_Up"


GO.df <- rbind(res1_bp_down_sub,
               res2_bp_down_sub, 
               res3_bp_down_sub,
               res4_bp_down_sub, 
               res1_bp_up_sub,
               res2_bp_up_sub,
               res3_bp_up_sub,
               res4_bp_up_sub)


GO.df <- GO.df %>% dplyr::select(Term, Annotated, Significant, p.value_elim, exp)

GO.df$ratio <- GO.df$Significant / GO.df$Annotated

GO.df$exp <- factor(GO.df$exp, levels = c("KO_vs_WT_Up", "KO_vs_WT_Down",
                                          "WTm7oe_vs_WT_Up", "WTm7oe_vs_WT_Down",
                                          "KOm7oe_vs_KO_Up", "KOm7oe_vs_KO_Down",
                                          "KOm7oe_vs_WT_Up", "KOm7oe_vs_WT_Down"))

GO.df$p.value_elim <- as.numeric(GO.df$p.value_elim)
GO.df$p.value_elim <- -log10(GO.df$p.value_elim)
GO.df$Term <- factor(GO.df$Term, levels = rev(c(Terms[1],
                                            Terms[3],
                                            Terms[2],
                                            Terms[5],
                                            Terms[6],
                                            Terms[7], Terms[8],
                                            Terms[9],
                                            Terms[11],
                                            Terms[c(10,12,13,14,15,16)])))




g1 <- ggplot(GO.df, aes(x=exp,
                        y=Term,
                        colour=p.value_elim,
                        size=ratio)) +
  geom_point() +
  expand_limits(x=0) +
  labs(x="Expected", y="GO term", colour="-log10(p-value)", size="gene ratio") +
  theme_minimal() + theme(axis.title = element_blank(),
                          axis.text.y =  element_text(size = 10, color = "black"),
                          axis.text.x = element_text(angle = 90, size = 10, vjust = 0.5, color = "black")) +
  scale_color_gradientn("-log10(p-value)", colours = c(paletteer_d("ggsci::amber_material"))) + 
  guides(color = guide_colourbar(order = 1), size = guide_legend(order = 2))


g1 

ggsave(filename = "/data/rajewsky/projects/cdr1as_ko_snRNA/codes_github/cdr1as/Figures/GO_figures_dot.pdf",
       plot = g1,
       scale = 1, width = 8, height = 6, units = "in", device = cairo_pdf,
       dpi = 300)


# ------------------------- #
# heatmap
# ------------------------- #
Terms2 <- c(Terms[9], Terms[11],
           Terms[7], Terms[c(12,13,15,10)])


extract_genes <- function(term = "startle response"){
  
  final_genes <- c()
  GO <- c(list)
  GO[[1]] <- res1_bp_down
  GO[[2]] <- res2_bp_down
  GO[[3]] <- res3_bp_down
  GO[[4]] <- res4_bp_down
  GO[[5]] <- res1_bp_up
  GO[[6]] <- res2_bp_up
  GO[[7]] <- res3_bp_up
  GO[[8]] <- res4_bp_up
  
  for(i in c(1:8)){
    genes1 <- GO[[i]][GO[[i]]$Term %in% term, ]$genes 
    genes1 <- unlist(strsplit(genes1, ","))
    genes1 <- str_squish(genes1)
    final_genes <- c(final_genes, genes1)
    
  }
  
  return(final_genes)
}
gene1 <- extract_genes(term = Terms2[1])
gene2 <- extract_genes(term = Terms2[2])
gene3 <- extract_genes(term = Terms2[3])
gene4 <- extract_genes(term = Terms2[4])
gene5 <- extract_genes(term = Terms2[5])
gene6 <- extract_genes(term = Terms2[6])
gene7 <- extract_genes(term = Terms2[7])
#gene8 <- extract_genes(term = Terms2[8])

ZS.df <- apply(cor_mat, 1, function(x) (x - mean(x)) / sd(x))
ZS.df <- data.frame(t(ZS.df))
ZS.df$gene <- rownames(ZS.df)


ZS.df <- ZS.df %>% mutate(GO = case_when(gene %in% gene1 ~ "SMAD protein \nsignal \ntransduction" ,
                                         gene %in% gene2 ~ "regulation of \ncalcium \nion-dependent \nexocytosis",
                                         gene %in% gene3 ~ "regulation of \nglucose \ntransmembrane \ntransport" ,
                                         gene %in% gene4 ~ "behavioral \nfear response",
                                         gene %in% gene5 ~ "startle \nresponse",
                                         gene %in% gene6 ~ "regulation of \nsynaptic \ntransmission, \nglutamatergic",
                                         gene %in% gene7 ~ "regulated \nexocytosis"))

ZS.df <- ZS.df[!is.na(ZS.df$GO), ]

ZS.df2 <- tidyr::gather(ZS.df, sample, z_score, "WTC1":"KOm7oe4")

ZS.df2 <- ZS.df2 %>% mutate(exp = case_when(sample %in% c("WTC1", "WTC2", "WTC3","WTC4") ~ "WT",
                                            sample %in% c("WTm7oe1", "WTm7oe2", "WTm7oe3","WTm7oe4") ~ "WTm7oe",
                                            sample %in% c("KOC1", "KOC2", "KOC3", "KOC4") ~ "KO",
                                            sample %in% c("KOm7oe1", "KOm7oe2", "KOm7oe3", "KOm7oe4") ~ "KOm7oe"))

ZS.df2$exp <- factor(ZS.df2$exp, levels = c("WT", "WTm7oe", "KO", "KOm7oe"))


ZS.df2$GO <- factor(ZS.df2$GO, levels = c("SMAD protein \nsignal \ntransduction",
                                          "regulation of \ncalcium \nion-dependent \nexocytosis",
                                          "regulation of \nglucose \ntransmembrane \ntransport" ,
                                          "behavioral \nfear response",
                                          "startle \nresponse",
                                          "regulation of \nsynaptic \ntransmission, \nglutamatergic",
                                          "regulated \nexocytosis"))

# remove Mir7-1
ZS.df2 <- ZS.df2[!ZS.df2$gene %in% "Mir7-1", ]


g2 <- ggplot(data = ZS.df2, mapping = aes(x = sample, y = gene, fill = z_score)) +
  geom_tile() + 
  #scale_x_discrete(expand = c(0,0)) + 
  scale_y_discrete(position = "right",  expand = c(0,0))  + coord_flip() +  theme_classic() + theme(plot.title = element_text(size = 15, hjust = 0.5), 
                                                                                                    axis.title.y = element_blank(),
                                                                                                    axis.title.x = element_blank(),
                                                                                                    axis.text.y = element_blank(),
                                                                                                    axis.text.x = element_text(angle = 90, vjust = 11.5, size = 12, color = "black"),
                                                                                                    #axis.text.y = element_text(size = 12, color = "black"),
                                                                                                    axis.line = element_line(color = "white"),
                                                                                                    axis.ticks = element_line(color = "white"),
                                                                                                    legend.title = element_text(size = 11, color = "black"),
                                                                                                    legend.text = element_text(size = 11, color = "black"),
                                                                                                    panel.spacing.x = unit(0.05, "lines"), # change to 0.05 as Cledi asked for space
                                                                                                    panel.spacing.y = unit(0.05, "lines"), # change to 0.05 as Cledi asked for space
                                                                                                    panel.margin= unit(c(0), "lines"), 
                                                                                                    plot.margin = margin(0, 0, 0, 0, "cm"),
                                                                                                    strip.background = element_blank(),
                                                                                                    strip.text.x = element_text(size = 10)) + 
  #scale_fill_gradientn("expression", colours = c("#003366", "#FFCC66", "#990000"), limits = c(-1, 1), oob=squish) +
  scale_fill_gradientn("z-score", colours = c(rev(paletteer_d("ggsci::light_blue_material")), paletteer_d("ggsci::orange_material")), limits = c(-1.5, 1.5), oob=squish) +
  #scale_fill_gradientn(colours = cols, limits = c(-2, 2), oob=squish) + 
  facet_grid(exp ~ GO, scales = "free", space = "free", switch = 'x')
  
g2


ggsave(filename = "/data/rajewsky/projects/cdr1as_ko_snRNA/codes_github/cdr1as/Figures/GO_figures_heatmap.pdf",
       plot = g2,
       scale = 1, width = 11, height = 4.5, units = "in", device = cairo_pdf,
       dpi = 300)
