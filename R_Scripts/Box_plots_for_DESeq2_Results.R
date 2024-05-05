
library(stringr)
library(ggplot2)
library(cowplot)

# data load FullKO exp
r.df <- readRDS(file = "/data/rajewsky/projects/cdr1as_ko_snRNA/codes_github/cdr1as/DGEm/WTN_batch_corr_nested_DGEm.rda")
DGEm_normalized <- readRDS(file = "/data/rajewsky/projects/cdr1as_ko_snRNA/codes_github/cdr1as/DGEm/WTN_batch_corr_nested_normalized_DGEm.rda")

# data load 5SSKO exp
r2.df <- readRDS(r.df, file = "/data/rajewsky/projects/cdr1as_ko_snRNA/codes_github/cdr1as/DGEm/WTJ_nested_DGEm.rda")
DGEm_normalized2 <- readRDS(file = "/data/rajewsky/projects/cdr1as_ko_snRNA/codes_github/cdr1as/DGEm/WTJ_nested_normalized_DGEm.rda")



box_plot <- function(df = DGEm_normalized, gene = "Kdm5d"){
  sub.df <- df[rownames(df) %in% gene, ] %>% t %>% as.data.frame()
  colnames(sub.df) <- "gene"
  sub.df$sample <- c(rep("WT",4), rep("WTm7oe",4 ), rep("KO",4), rep("KOm7oe",4 ))
  sub.df$sample <- factor(sub.df$sample, levels = c("WT","WTm7oe","KO", "KOm7oe"))
  sub.df$rep <- str_sub(rownames(sub.df), -1, -1)
  g1 <- ggplot(sub.df, aes(x  = sample, y = gene)) + geom_boxplot() + 
    geom_jitter(position=position_jitter(0.1), color = "darkorange", alpha = 0.6, size = 4) + 
    scale_color_identity() + theme_cowplot()  + ggtitle(gene) + ylab("normalized counts") + 
    theme(axis.title.x = element_blank(), plot.title = element_text(hjust = 0.5))
  
  g1
}

####################
# FullKO data
####################


# 3 genes that co-upregulated by KO vs WT and KOm7oe vs WT
## [1] "Crocc"    "Pisd-ps1" "Cdr1os"

# Cdr1os (C230004F18Rik)
b1 <- box_plot(df=DGEm_normalized, gene = "C230004F18Rik")
b2 <- box_plot(df=r.df, gene = "C230004F18Rik") + ylab("batch corrected normalized counts")
b3 <- arrangeGrob(b1, b2, ncol = 2)

ggsave(filename = "/data/rajewsky/projects/cdr1as_ko_snRNA/codes_github/cdr1as/Figures/Box_plot_FullKO_Cdr1os.pdf",
       plot = b3,
       scale = 1, width = 10, height = 4.5, units = "in", device = cairo_pdf,
       dpi = 300)

# Crocc
b1 <- box_plot(df=DGEm_normalized, gene = "Crocc")
b2 <- box_plot(df=r.df, gene = "Crocc") + ylab("batch corrected normalized counts")
b3 <- arrangeGrob(b1, b2, ncol = 2)

ggsave(filename = "/data/rajewsky/projects/cdr1as_ko_snRNA/codes_github/cdr1as/Figures/Box_plot_FullKO_Crocc.pdf",
       plot = b3,
       scale = 1, width = 10, height = 4.5, units = "in", device = cairo_pdf,
       dpi = 300)

# Pisd-ps1
b1 <- box_plot(df=DGEm_normalized, gene = "Pisd-ps1")
b2 <- box_plot(df=r.df, gene = "Pisd-ps1") + ylab("batch corrected normalized counts")
b3 <- arrangeGrob(b1, b2, ncol = 2)

ggsave(filename = "/data/rajewsky/projects/cdr1as_ko_snRNA/codes_github/cdr1as/Figures/Box_plot_FullKO_Pisd_ps1.pdf",
       plot = b3,
       scale = 1, width = 10, height = 4.5, units = "in", device = cairo_pdf,
       dpi = 300)

# 7 genes that co-downregulated by KO vs WT, WTm7oe vs WT and KOm7oe vs WT
## [1] "Ctss"  "C1qb"  "C1qc"  "C1qa"  "Lyz2"  "Trem2" "Mpeg1"

# Ctss
b1 <- box_plot(df=DGEm_normalized, gene = "Ctss")
b2 <- box_plot(df=r.df, gene = "Ctss") + ylab("batch corrected normalized counts")
b3 <- arrangeGrob(b1, b2, ncol = 2)

ggsave(filename = "/data/rajewsky/projects/cdr1as_ko_snRNA/codes_github/cdr1as/Figures/Box_plot_FullKO_Ctss.pdf",
       plot = b3,
       scale = 1, width = 10, height = 4.5, units = "in", device = cairo_pdf,
       dpi = 300)

# C1qb
b1 <- box_plot(df=DGEm_normalized, gene = "C1qb")
b2 <- box_plot(df=r.df, gene = "C1qb") + ylab("batch corrected normalized counts")
b3 <- arrangeGrob(b1, b2, ncol = 2)

ggsave(filename = "/data/rajewsky/projects/cdr1as_ko_snRNA/codes_github/cdr1as/Figures/Box_plot_FullKO_C1qb.pdf",
       plot = b3,
       scale = 1, width = 10, height = 4.5, units = "in", device = cairo_pdf,
       dpi = 300)

# C1qc
b1 <- box_plot(df=DGEm_normalized, gene = "C1qc")
b2 <- box_plot(df=r.df, gene = "C1qc") + ylab("batch corrected normalized counts")
b3 <- arrangeGrob(b1, b2, ncol = 2)

ggsave(filename = "/data/rajewsky/projects/cdr1as_ko_snRNA/codes_github/cdr1as/Figures/Box_plot_FullKO_C1qc.pdf",
       plot = b3,
       scale = 1, width = 10, height = 4.5, units = "in", device = cairo_pdf,
       dpi = 300)

# C1qa
b1 <- box_plot(df=DGEm_normalized, gene = "C1qa")
b2 <- box_plot(df=r.df, gene = "C1qa") + ylab("batch corrected normalized counts")
b3 <- arrangeGrob(b1, b2, ncol = 2)

ggsave(filename = "/data/rajewsky/projects/cdr1as_ko_snRNA/codes_github/cdr1as/Figures/Box_plot_FullKO_C1qa.pdf",
       plot = b3,
       scale = 1, width = 10, height = 4.5, units = "in", device = cairo_pdf,
       dpi = 300)

# Lyz2
b1 <- box_plot(df=DGEm_normalized, gene = "Lyz2")
b2 <- box_plot(df=r.df, gene = "Lyz2") + ylab("batch corrected normalized counts")
b3 <- arrangeGrob(b1, b2, ncol = 2)

ggsave(filename = "/data/rajewsky/projects/cdr1as_ko_snRNA/codes_github/cdr1as/Figures/Box_plot_FullKO_Lyz2.pdf",
       plot = b3,
       scale = 1, width = 10, height = 4.5, units = "in", device = cairo_pdf,
       dpi = 300)

# Trem2
b1 <- box_plot(df=DGEm_normalized, gene = "Trem2")
b2 <- box_plot(df=r.df, gene = "Trem2") + ylab("batch corrected normalized counts")
b3 <- arrangeGrob(b1, b2, ncol = 2)

ggsave(filename = "/data/rajewsky/projects/cdr1as_ko_snRNA/codes_github/cdr1as/Figures/Box_plot_FullKO_Trem2.pdf",
       plot = b3,
       scale = 1, width = 10, height = 4.5, units = "in", device = cairo_pdf,
       dpi = 300)

# Mpeg1
b1 <- box_plot(df=DGEm_normalized, gene = "Mpeg1")
b2 <- box_plot(df=r.df, gene = "Mpeg1") + ylab("batch corrected normalized counts")
b3 <- arrangeGrob(b1, b2, ncol = 2)

ggsave(filename = "/data/rajewsky/projects/cdr1as_ko_snRNA/codes_github/cdr1as/Figures/Box_plot_FullKO_Mpeg1.pdf",
       plot = b3,
       scale = 1, width = 10, height = 4.5, units = "in", device = cairo_pdf,
       dpi = 300)


# additional check
# 1700020I14Rik
b1 <- box_plot(df=DGEm_normalized, gene = "1700020I14Rik")
b2 <- box_plot(df=r.df, gene = "1700020I14Rik") + ylab("batch corrected normalized counts")
b3 <- arrangeGrob(b1, b2, ncol = 2)

ggsave(filename = "/data/rajewsky/projects/cdr1as_ko_snRNA/codes_github/cdr1as/Figures/Box_plot_FullKO_Cyrano.pdf",
       plot = b3,
       scale = 1, width = 10, height = 4.5, units = "in", device = cairo_pdf,
       dpi = 300)

####################
# 5SS KO data
####################
box_plot2 <- function(df = DGEm_normalized, gene = "Kdm5d"){
  sub.df <- df[rownames(df) %in% gene, ] %>% t %>% as.data.frame()
  colnames(sub.df) <- "gene"
  sub.df$sample <- c(rep("WT",3), rep("WTm7oe",3), rep("KO",3), rep("KOm7oe",3))
  sub.df$sample <- factor(sub.df$sample, levels = c("WT","WTm7oe","KO", "KOm7oe"))
  sub.df$rep <- str_sub(rownames(sub.df), -1, -1)
  g1 <- ggplot(sub.df, aes(x  = sample, y = gene)) + geom_boxplot() + 
    geom_jitter(position=position_jitter(0.1), color = "darkorange", alpha = 0.6, size = 4) + 
    scale_color_identity() + theme_cowplot()  + ggtitle(gene) + ylab("normalized counts") + 
    theme(axis.title.x = element_blank(), plot.title = element_text(hjust = 0.5))
  
  g1
}

# C230004F18Rik
b1 <- box_plot2(df=DGEm_normalized2, gene = "C230004F18Rik")
b2 <- box_plot2(df=r2.df, gene = "C230004F18Rik") + ylab("nested design normalized counts")
b3 <- arrangeGrob(b1, b2, ncol = 2)

ggsave(filename = "/data/rajewsky/projects/cdr1as_ko_snRNA/codes_github/cdr1as/Figures/Box_plot_5SSKO_Cdr1os.pdf",
       plot = b3,
       scale = 1, width = 10, height = 4.5, units = "in", device = cairo_pdf,
       dpi = 300)

# Pisd-ps1
b1 <- box_plot2(df=DGEm_normalized2, gene = "Pisd-ps1")
b2 <- box_plot2(df=r2.df, gene = "Pisd-ps1") + ylab("nested design normalized counts")
b3 <- arrangeGrob(b1, b2, ncol = 2)

ggsave(filename = "/data/rajewsky/projects/cdr1as_ko_snRNA/codes_github/cdr1as/Figures/Box_plot_5SSKO_Pisd_ps1.pdf",
       plot = b3,
       scale = 1, width = 10, height = 4.5, units = "in", device = cairo_pdf,
       dpi = 300)

# Crocc
b1 <- box_plot2(df=DGEm_normalized2, gene = "Crocc")
b2 <- box_plot2(df=r2.df, gene = "Crocc") + ylab("nested design normalized counts")
b3 <- arrangeGrob(b1, b2, ncol = 2)

ggsave(filename = "/data/rajewsky/projects/cdr1as_ko_snRNA/codes_github/cdr1as/Figures/Box_plot_5SSKO_Crocc.pdf",
       plot = b3,
       scale = 1, width = 10, height = 4.5, units = "in", device = cairo_pdf,
       dpi = 300)

## [1] "Ctss"  "C1qb"  "C1qc"  "C1qa"  "Lyz2"  "Trem2" "Mpeg1"
# Ctss
b1 <- box_plot2(df=DGEm_normalized2, gene = "Ctss")
b2 <- box_plot2(df=r2.df, gene = "Ctss") + ylab("nested design normalized counts")
b3 <- arrangeGrob(b1, b2, ncol = 2)

ggsave(filename = "/data/rajewsky/projects/cdr1as_ko_snRNA/codes_github/cdr1as/Figures/Box_plot_5SSKO_Ctss.pdf",
       plot = b3,
       scale = 1, width = 10, height = 4.5, units = "in", device = cairo_pdf,
       dpi = 300)

# C1qb
b1 <- box_plot2(df=DGEm_normalized2, gene = "C1qb")
b2 <- box_plot2(df=r2.df, gene = "C1qb") + ylab("nested design normalized counts")
b3 <- arrangeGrob(b1, b2, ncol = 2)

ggsave(filename = "/data/rajewsky/projects/cdr1as_ko_snRNA/codes_github/cdr1as/Figures/Box_plot_5SSKO_C1qb.pdf",
       plot = b3,
       scale = 1, width = 10, height = 4.5, units = "in", device = cairo_pdf,
       dpi = 300)

# C1qc
b1 <- box_plot2(df=DGEm_normalized2, gene = "C1qc")
b2 <- box_plot2(df=r2.df, gene = "C1qc") + ylab("nested design normalized counts")
b3 <- arrangeGrob(b1, b2, ncol = 2)

ggsave(filename = "/data/rajewsky/projects/cdr1as_ko_snRNA/codes_github/cdr1as/Figures/Box_plot_5SSKO_C1qc.pdf",
       plot = b3,
       scale = 1, width = 10, height = 4.5, units = "in", device = cairo_pdf,
       dpi = 300)

# C1qa
b1 <- box_plot2(df=DGEm_normalized2, gene = "C1qa")
b2 <- box_plot2(df=r2.df, gene = "C1qa") + ylab("nested design normalized counts")
b3 <- arrangeGrob(b1, b2, ncol = 2)

ggsave(filename = "/data/rajewsky/projects/cdr1as_ko_snRNA/codes_github/cdr1as/Figures/Box_plot_5SSKO_C1qa.pdf",
       plot = b3,
       scale = 1, width = 10, height = 4.5, units = "in", device = cairo_pdf,
       dpi = 300)

# Lyz2
b1 <- box_plot2(df=DGEm_normalized2, gene = "Lyz2")
b2 <- box_plot2(df=r2.df, gene = "Lyz2") + ylab("nested design normalized counts")
b3 <- arrangeGrob(b1, b2, ncol = 2)

ggsave(filename = "/data/rajewsky/projects/cdr1as_ko_snRNA/codes_github/cdr1as/Figures/Box_plot_5SSKO_Lyz2.pdf",
       plot = b3,
       scale = 1, width = 10, height = 4.5, units = "in", device = cairo_pdf,
       dpi = 300)

# Trem2
b1 <- box_plot2(df=DGEm_normalized2, gene = "Trem2")
b2 <- box_plot2(df=r2.df, gene = "Trem2") + ylab("nested design normalized counts")
b3 <- arrangeGrob(b1, b2, ncol = 2)

ggsave(filename = "/data/rajewsky/projects/cdr1as_ko_snRNA/codes_github/cdr1as/Figures/Box_plot_5SSKO_Trem2.pdf",
       plot = b3,
       scale = 1, width = 10, height = 4.5, units = "in", device = cairo_pdf,
       dpi = 300)

# Mpeg1
b1 <- box_plot2(df=DGEm_normalized2, gene = "Mpeg1")
b2 <- box_plot2(df=r2.df, gene = "Mpeg1") + ylab("nested design normalized counts")
b3 <- arrangeGrob(b1, b2, ncol = 2)

ggsave(filename = "/data/rajewsky/projects/cdr1as_ko_snRNA/codes_github/cdr1as/Figures/Box_plot_5SSKO_Mpeg1.pdf",
       plot = b3,
       scale = 1, width = 10, height = 4.5, units = "in", device = cairo_pdf,
       dpi = 300)

# additional check
# 1700020I14Rik
b1 <- box_plot2(df=DGEm_normalized2, gene = "1700020I14Rik")
b2 <- box_plot2(df=r2.df, gene = "1700020I14Rik") + ylab("nested design normalized counts")
b3 <- arrangeGrob(b1, b2, ncol = 2)

ggsave(filename = "/data/rajewsky/projects/cdr1as_ko_snRNA/codes_github/cdr1as/Figures/Box_plot_5SSKO_Cyranos.pdf",
       plot = b3,
       scale = 1, width = 10, height = 4.5, units = "in", device = cairo_pdf,
       dpi = 300)


