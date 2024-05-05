

# Snca, Cplx1, Prkcb, Pfn2, Zdhhc9, Cspa, Basp1, Phactr1, Wipf2

res1.df <- readRDS(file = "/data/rajewsky/projects/cdr1as_ko_snRNA/codes_github/cdr1as/Results/DESeq2_WTN_batch_corrected_res1.rda")
res2.df <- readRDS(file = "/data/rajewsky/projects/cdr1as_ko_snRNA/codes_github/cdr1as/Results/DESeq2_WTN_nested_corrected_res2.rda")
res3.df <- readRDS(file = "/data/rajewsky/projects/cdr1as_ko_snRNA/codes_github/cdr1as/Results/DESeq2_WTN_nested_corrected_res3.rda")
res4.df <- readRDS(file = "/data/rajewsky/projects/cdr1as_ko_snRNA/codes_github/cdr1as/Results/DESeq2_WTN_batch_corrected_res4.rda")


bar_plot <- function(df = res1.df, genes = miR7_target_genes){
  df_sub <- df[df$gene %in% genes, ]
  df_sub$sig <- ifelse(df_sub$padj < 0.05, TRUE, FALSE)
  tmp <- df_sub$gene[order(df_sub$log2FoldChange, decreasing = FALSE)]
  df_sub$gene <- factor(df_sub$gene, levels = tmp)
  return(df_sub)
}

pheno_genes <- c("Snca", "Cplx1", "Prkcb", "Pfn2", "Zdhhc9", "Cspa", "Basp1", "Phactr1", "Wipf2")

res1_sub <- bar_plot(df = res1.df, genes = pheno_genes)
res2_sub <- bar_plot(df = res2.df, genes = pheno_genes)
res3_sub <- bar_plot(df = res3.df, genes = pheno_genes)
res4_sub <- bar_plot(df = res4.df, genes = pheno_genes)


g1 <- ggplot(res1_sub) + 
  geom_bar(mapping = aes(x = gene, y = log2FoldChange, alpha= sig), fill = "#FF0066",stat = "identity", position=position_dodge()) +
  #  geom_point(mapping = aes(x = gene, y = (count) / 2000 , shape=exp, color=exp)) +
  scale_y_continuous(name = expression("log2FoldChange"), limits = c(-2, 2)) + theme_cowplot() +
  theme(plot.title = element_text(hjust = 0.5), axis.title.x = element_blank(), axis.text.x = element_text(angle = 90, hjust=0.95,vjust=0.3, size = 16), legend.position = c(0.8, 0.8))  +
  ggtitle("KO vs WT")


g2 <- ggplot(res2_sub) + 
  geom_bar(mapping = aes(x = gene, y = log2FoldChange, alpha= sig), fill = "#880E4F",stat = "identity", position=position_dodge()) +
  #  geom_point(mapping = aes(x = gene, y = (count) / 2000 , shape=exp, color=exp)) +
  scale_y_continuous(name = expression("log2FoldChange"), limits = c(-2, 2)) + theme_cowplot() +
  theme(plot.title = element_text(hjust = 0.5), axis.title.x = element_blank(),axis.text.x = element_text(angle = 90, hjust=0.95,vjust=0.3, size = 16), legend.position = c(0.8, 0.8))   +
  ggtitle("WTm7oe vs WT")

g3 <- ggplot(res3_sub) + 
  geom_bar(mapping = aes(x = gene, y = log2FoldChange, alpha= sig), fill = "#0099FF",stat = "identity", position=position_dodge()) +
  #  geom_point(mapping = aes(x = gene, y = (count) / 2000 , shape=exp, color=exp)) +
  scale_y_continuous(name = expression("log2FoldChange"), limits = c(-2, 2)) + theme_cowplot() +
  theme(plot.title = element_text(hjust = 0.5), axis.title.x = element_blank(),axis.text.x = element_text(angle = 90, hjust=0.95,vjust=0.3, size = 16), legend.position = c(0.8, 0.8)) +
  #scale_alpha_discrete(range = c(1, 0))   +
  ggtitle("KOm7oe vs KO")


g4 <- ggplot(res4_sub) + 
  geom_bar(mapping = aes(x = gene, y = log2FoldChange, alpha= sig), fill = "#0D47A1",stat = "identity", position=position_dodge()) +
  #  geom_point(mapping = aes(x = gene, y = (count) / 2000 , shape=exp, color=exp)) +
  scale_y_continuous(name = expression("log2FoldChange"), limits = c(-2, 2)) + theme_cowplot() +
  ggtitle("KOm7oe vs WT") + 
  theme(plot.title = element_text(hjust = 0.5), axis.title.x = element_blank(), axis.text.x = element_text(angle = 90, hjust=0.95,vjust=0.3, size = 16), legend.position = c(0.8, 0.8))

grid.arrange(g1, g2, g3, g4, ncol = 2)

c1 <- arrangeGrob(g1, g2, g3, g4, ncol = 2)

ggsave(filename = "/data/rajewsky/projects/cdr1as_ko_snRNA/codes_github/cdr1as/Figures/GRN_RF_SR_paper_genes_bar.pdf",
       plot = c1,
       scale = 1, width = 10, height = 7.5, units = "in", device = cairo_pdf,
       dpi = 300)
