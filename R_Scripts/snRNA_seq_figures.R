.libPaths(c("/data/rajewsky/home/skim/R/usr_lib_Seurat/", "/data/rajewsky/home/skim/R/usr_lib_v4/"))
library(Seurat)
library(ggplot2)
library(dplyr)
library(ggrepel)
library(cowplot)
library(stringr)
library(gridExtra)

all.SO <- readRDS(file = "/data/rajewsky/projects/cdr1as_ko_snRNA/codes/snRNA_seq/merged_final.rda")
rownames(all.SO@meta.data) <- all.SO@assays$RNA@counts@Dimnames[[2]]
all.meta <- all.SO@meta.data

# doublets
scrublet.df <- all.meta[, c("seurat_clusters", "predicted_doublets")]
# zero counts (singlet)
scrublet.df1 <- scrublet.df %>% group_by(seurat_clusters) %>% summarise(zero_count = sum(predicted_doublets == "False"))
# non-zero counts (doublets)
scrublet.df2 <-  scrublet.df %>% group_by(seurat_clusters) %>% summarise(non_zero_count = sum(predicted_doublets == "True"))

# merge
scrublet.df3 <- cbind(scrublet.df1, non_zero_count = scrublet.df2$non_zero_count)
scrublet.df3$all <- scrublet.df3$zero_count + scrublet.df3$non_zero_count
scrublet.df3$doublets_percent <- scrublet.df3$non_zero_count/scrublet.df3$all

# 30, 39, 44 doublets
ggplot(scrublet.df3, aes(x = seurat_clusters, y = doublets_percent)) +
  geom_bar(stat="identity") + theme_cowplot()



# Cdr1as
Cdr1as_033_uniq <- read.csv("/data/rajewsky/projects/cdr1as_ko_snRNA/2nd_sequencing_run/circRNAs_analysis/20210424_extract_Cdr1as_CB_and_UB_SP081_033_unique.csv")
Cdr1as_034_uniq <- read.csv("/data/rajewsky/projects/cdr1as_ko_snRNA/2nd_sequencing_run/circRNAs_analysis/20210424_extract_Cdr1as_CB_and_UB_SP081_034_unique.csv")
Cdr1as_039_uniq <- read.csv("/data/rajewsky/projects/cdr1as_ko_snRNA/2nd_sequencing_run/circRNAs_analysis/20210424_extract_Cdr1as_CB_and_UB_SP081_039_unique.csv")

Cdr1as_033_uniq <- table(Cdr1as_033_uniq$CB) %>% as.data.frame
Cdr1as_033_uniq$Var1 <-  str_replace(Cdr1as_033_uniq$Var1, "CB:Z:", "")
Cdr1as_033_uniq$Var1 <-  str_replace(Cdr1as_033_uniq$Var1, "-1", "-1")

Cdr1as_034_uniq <- table(Cdr1as_034_uniq$CB) %>% as.data.frame
Cdr1as_034_uniq$Var1 <-  str_replace(Cdr1as_034_uniq$Var1, "CB:Z:", "")
Cdr1as_034_uniq$Var1 <-  str_replace(Cdr1as_034_uniq$Var1, "-1", "-2")

Cdr1as_039_uniq <- table(Cdr1as_039_uniq$CB) %>% as.data.frame
Cdr1as_039_uniq$Var1 <-  str_replace(Cdr1as_039_uniq$Var1, "CB:Z:", "")
Cdr1as_039_uniq$Var1 <-  str_replace(Cdr1as_039_uniq$Var1, "-1", "-5")

Cdr1as_all <- rbind(Cdr1as_033_uniq, Cdr1as_034_uniq, Cdr1as_039_uniq)
colnames(Cdr1as_all) <- c("cell_ID", "Cdr1as")


all.meta$cell_ID <- rownames(all.meta)
all.meta <- left_join(all.meta, Cdr1as_all, by = "cell_ID")
all.meta$Cdr1as[is.na(all.meta$Cdr1as)] <- 0
rownames(all.meta) <- all.meta$cell_ID
all.SO@meta.data <- all.meta

all.meta <- cbind(all.meta, all.SO@reductions$umap@cell.embeddings)



g1 <- ggplot(all.meta[all.meta$exp %in% "WT",], aes(x = UMAP_1, y = UMAP_2)) + 
  geom_point(data = all.meta[which(all.meta$exp %in% "WT" & all.meta$Cdr1as %in% 0), ], aes(x = UMAP_1, y = UMAP_2), size = 0.1, color = "gray") +  
  geom_point(data = all.meta[which(all.meta$exp %in% "WT" & all.meta$Cdr1as %in% 1), ], aes(x = UMAP_1, y = UMAP_2), size = 0.7, color = "red") + 
  theme_cowplot() + ggtitle("WT (Cdr1as)") +
  theme(axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(), 
        axis.title = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 50))

ggsave(filename = "/data/rajewsky/projects/cdr1as_ko_snRNA/codes_github/cdr1as/Figures_snRNA/WT_cdr1as.pdf",
       plot = g1,
       scale = 1, width = 5, height = 5, units = "in", device = cairo_pdf,
       dpi = 300)

g2 <- ggplot(all.meta[all.meta$exp %in% "FullKO",], aes(x = UMAP_1, y = UMAP_2)) + 
  geom_point(data = all.meta[which(all.meta$exp %in% "FullKO" & all.meta$Cdr1as %in% 0), ], aes(x = UMAP_1, y = UMAP_2), size = 0.1, color = "gray") +  
  geom_point(data = all.meta[which(all.meta$exp %in% "FullKO" & all.meta$Cdr1as %in% 1), ], aes(x = UMAP_1, y = UMAP_2), size = 0.7, color = "blue") + 
  theme_cowplot() + ggtitle("KO (Cdr1as)") +
  theme(axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(), 
        axis.title = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 50))

ggsave(filename = "/data/rajewsky/projects/cdr1as_ko_snRNA/codes_github/cdr1as/Figures_snRNA/KO_cdr1as.pdf",
       plot = g2,
       scale = 1, width = 5, height = 5, units = "in", device = cairo_pdf,
       dpi = 300)






dropviz_featureplot <- function(markers = "Astrocytes_Gja11", titletext = "Astrocytes dropviz markers"){
  FeaturePlot(object = all.SO, features = markers, cols = c("lightgray", "darkred"), order = TRUE) + 
    ggtitle(titletext) + theme_void() + theme(plot.title = element_text(size = 12, family = "helvetica", color = "black", hjust = 0.5),
                                              legend.text = element_text(size = 12, family = "helvetica", color = "black", hjust = 0.5))
}

f1 <- FeaturePlot(object = all.SO, features = 'Astrocytes_Gja11', cols = c("lightgrey", "darkred"))
f2 <- FeaturePlot(object = all.SO, features = 'Choroid_Plexus_Ttr1', cols = c("lightgrey", "darkred"))
f3 <- FeaturePlot(object = all.SO, features = 'Endothelial_Flt11', cols = c("lightgrey", "darkred"))
f4 <- FeaturePlot(object = all.SO, features = 'Ependyma1', cols = c("lightgrey", "darkred"))
f5 <- FeaturePlot(object = all.SO, features = 'Fibroblast_Dcn1', cols = c("lightgrey", "darkred"))
f6 <- FeaturePlot(object = all.SO, features = 'Interneuron_Gad21', cols = c("lightgrey", "darkred"))
f7 <- FeaturePlot(object = all.SO, features = 'Microglia_Macrophage_C1qb1', cols = c("lightgrey", "darkred"))
f8 <- FeaturePlot(object = all.SO, features = 'Mural_Rgs5Acta21', cols = c("lightgrey", "darkred"))
f9 <- FeaturePlot(object = all.SO, features = 'Neurogenesis_Sox41', cols = c("lightgrey", "darkred"))
f10 <- FeaturePlot(object = all.SO, features = 'Neuron_CA1_Subiculum_Postsubiculum_Entorhinal1', cols = c("lightgrey", "darkred"))
f11 <- FeaturePlot(object = all.SO, features = 'Neuron_CA2CA3_Pvrl31', cols = c("lightgrey", "darkred"))
f12 <- FeaturePlot(object = all.SO, features = 'Neuron_CajalRetzius_Lhx11', cols = c("lightgrey", "darkred"))
f13 <- FeaturePlot(object = all.SO, features = 'Neuron_Dentate_C1ql21', cols = c("lightgrey", "darkred"))
f14 <- FeaturePlot(object = all.SO, features = 'Oligodendrocyte_Tfr1', cols = c("lightgrey", "darkred"))
f15 <- FeaturePlot(object = all.SO, features = 'Polydendrocyte_Tnr1', cols = c("lightgrey", "darkred"))
f16 <- FeaturePlot(object = all.SO, features = 'Subiculum_Entorhinal_Nxph31', cols = c("lightgrey", "darkred"))
f17 <- FeaturePlot(object = all.SO, features = 'Subiculum_Slc17a61', cols = c("lightgrey", "darkred"))

c1 <- arrangeGrob(f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11, f12, f13, f14, f15, f16, f17, ncol = 2)

ggsave(filename = "/data/rajewsky/projects/cdr1as_ko_snRNA/codes_github/cdr1as/Figures_snRNA/addmodulescore.pdf",
       plot = c1,
       scale = 1, width = 10, height = 30, units = "in", device = cairo_pdf,
       dpi = 300)



# Celltype define
new.cluster.ids <- c("CA1", "Dentate_Gyrus", "Dentate_Gyrus", "Oligodendrocytes", "Oligodendrocytes", 
                     "Neurons", "Dentate_Gyrus", "Oligodendrocytes", "CA2_3", "Astrocytes",
                     
                     "subiculum", "CA1", "Oligodendrocytes", "subiculum", "Microglia",
                     "Microglia", "CA1", "OPC", "Inhibitory", "subiculum",
                     
                     "Inhibitory", "subiculum", "Inhibitory", "CA2_3", "subiculum",
                     "Inhibitory", "Neurons", "subiculum", "subiculum", "Fibroblast",
                     
                     "doublets", "Cajal_Retzius", "Microglia", "Inhibitory", "Mural",
                     "subiculum", "Choroid_Plexus", "Neurons","Neurons","doublets",
                     
                  
                     "subiculum", "Endothelial", "subiculum","subiculum","doublets",
                     "Fibroblast", "Inhibitory", "COP")


names(new.cluster.ids) <- levels(all.SO)
all.SO <- RenameIdents(all.SO, new.cluster.ids)
all.meta <- all.SO@meta.data
all.meta$cell_type <- all.SO@active.ident

DimPlot(all.SO, reduction = "umap", label = TRUE, pt.size = 0.5)
all.SO$cell_type <- all.SO@active.ident
all.SO$X <- NULL

saveRDS(all.SO, file = "/data/rajewsky/projects/cdr1as_ko_snRNA/codes/snRNA_seq/merged_final.rda")


all.SO <- readRDS(file = "/data/rajewsky/projects/cdr1as_ko_snRNA/codes/snRNA_seq/merged_final.rda")

all.SO <- subset(x = all.SO, subset = cell_type != "doublets")

all.SO@active.ident <- factor(all.SO@active.ident, levels = c("Dentate_Gyrus", "CA1", "CA2_3","subiculum","Neurons",
                                                              "Inhibitory", 
                                                              "Astrocytes", "Microglia", 
                                                              "Oligodendrocytes", "COP", "OPC",
                                                              "Fibroblast", "Cajal_Retzius","Mural", "Choroid_Plexus", "Endothelial"))

d1 <- DimPlot(all.SO, reduction = "umap", label = FALSE, pt.size = 0.5, 
        cols = c('Dentate_Gyrus' = 'darkred',
                 'CA1' = 'red', 
                 'CA2_3' = 'orange',
                 'subiculum' = '#CC6600',
                 'Neurons' = 'grey',
                 'Inhibitory' = '#CC79A7',
                 'Astrocytes' = '#33CC33',
                 'Microglia' = 'darkgreen',
                 'Oligodendrocytes' = '#00CCFF',
                 'COP' = '#CCFF00', 
                 'OPC' = 'purple',
                 'Fibroblast' = '#CCFFFF',
                 'Cajal_Retzius' = '#FFFFCC',
                 'Mural' = 'navy',
                 'Choroid_Plexus' = '#99CC00',
                 'Endothelial' = '#FFFF00',
                 'doublets' = '#333333'
                 )) +    theme(axis.line = element_blank(),
                               axis.ticks = element_blank(),
                               axis.text = element_blank(), 
                               axis.title = element_blank()) + 
  guides(color = guide_legend(override.aes = list(size = 7), ncol=1) )


ggsave(filename = "/data/rajewsky/projects/cdr1as_ko_snRNA/codes_github/cdr1as/Figures_snRNA/UMAP_cell_type.pdf",
       plot = d1,
       scale = 1, width = 7, height = 5, units = "in", device = cairo_pdf,
       dpi = 300)

# broad cell types
all.meta <- all.meta %>% mutate(cell_type_precise = case_when(cell_type %in% "CA1" ~ "CA1",
                                                              cell_type %in% "Dentate_Gyrus" ~ "Dentate_Gyrus",
                                                              cell_type %in% "Oligodendrocytes" ~ "Oligodendrocytes",
                                                              cell_type %in% "Neurons" ~ "Neurons",
                                                              cell_type %in% "CA2_3" ~ "CA2_3",
                                                              cell_type %in% "Astrocytes" ~ "Astrocytes",
                                                              cell_type %in% "subiculum" ~ "subiculum",
                                                              cell_type %in% "Microglia" ~ "Microglia",
                                                              cell_type %in% "OPC" ~ "OPC",
                                                              cell_type %in% "COP" ~ "COP",
                                                              cell_type %in% "Inhibitory" ~ "Inhibitory",
                                                              cell_type %in% c("Fibroblast", "doublets", "Cajal_Retzius","Mural","Choroid_Plexus", "Endothelial") ~ "Rest"))

all.SO$broad_cell_type <- all.meta$cell_type_precise

d1 <- DimPlot(all.SO, reduction = "umap", label = FALSE, group.by =  "broad_cell_type", pt.size = 0.5,
              cols = c('Dentate_Gyrus' = 'darkred',
                       'CA1' = 'red', 
                       'CA2_3' = '#CC79A7',
                       'subiculum' = '#CC6600',
                       'Neurons' = '#FFFFCC',
                       'Inhibitory' = 'orange',
                       'Astrocytes' = '#CCFF00',
                       'Microglia' = 'darkgreen',
                       'Oligodendrocytes' = '#00CCFF',
                       'COP' = '#33CC33', 
                       'OPC' = 'purple',
                       'Rest' = '#CCCCCC'
              )) +    theme(axis.line = element_blank(),
                            axis.ticks = element_blank(),
                            axis.text = element_blank(), 
                            axis.title = element_blank(),
                            plot.title = element_blank(),
                            legend.text = element_text(size = 20)) + 
  guides(color = guide_legend(override.aes = list(size = 9), ncol=1) )

ggsave(filename = "/data/rajewsky/projects/cdr1as_ko_snRNA/codes_github/cdr1as/Figures_snRNA/UMAP_broad_cell_type.pdf",
       plot = d1,
       scale = 1, width = 8, height = 5, units = "in", device = cairo_pdf,
       dpi = 300)


# genes check featureplot
g1 <-  FeaturePlot(all.SO, features = "mir-7b", order = TRUE, pt.size = 1, cols = c("gray", "red")) + 
      theme_cowplot() + ggtitle("pri-miR-7b") +
  scale_color_gradientn( colours = c('lightgrey', 'red'),  limits = c(0, 1), oob = scales::squish) + 
      theme(axis.line = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_blank(), 
          axis.title = element_blank(),
          plot.title = element_text(hjust = 0.5, size = 50),
          legend.position = "none")
  
ggsave(filename = "/data/rajewsky/projects/cdr1as_ko_snRNA/codes_github/cdr1as/Figures_snRNA/UMAP_pri_mir7b.pdf",
       plot = g1,
       scale = 1, width = 5, height = 5, units = "in", device = cairo_pdf,
       dpi = 300)

g1 <-  FeaturePlot(all.SO, features = "mir-7a-2", order = TRUE, pt.size = 1, cols = c("gray", "red")) + 
  theme_cowplot() + ggtitle("pri-miR-7a-2") +
  scale_color_gradientn( colours = c('lightgrey', 'red'),  limits = c(0, 1), oob = scales::squish) + 
  theme(axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(), 
        axis.title = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 50),
        legend.position = "none")

ggsave(filename = "/data/rajewsky/projects/cdr1as_ko_snRNA/codes_github/cdr1as/Figures_snRNA/UMAP_pri_mir7a2.pdf",
       plot = g1,
       scale = 1, width = 5, height = 5, units = "in", device = cairo_pdf,
       dpi = 300)

g1 <-  FeaturePlot(all.SO, features = "mir-146a", order = TRUE, pt.size = 1, cols = c("gray", "red")) + 
  theme_cowplot() + ggtitle("pri-miR-146a") +
  scale_color_gradientn( colours = c('lightgrey', 'red'),  limits = c(0, 1), oob = scales::squish) + 
  theme(axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(), 
        axis.title = element_blank(),
        plot.title = element_text(hjust = 0.5,size = 50),
        legend.position = "none")

ggsave(filename = "/data/rajewsky/projects/cdr1as_ko_snRNA/codes_github/cdr1as/Figures_snRNA/UMAP_pri_mir146a.pdf",
       plot = g1,
       scale = 1, width = 5, height = 5, units = "in", device = cairo_pdf,
       dpi = 300)

g1 <-  FeaturePlot(all.SO, features = "mir-34a", order = TRUE, pt.size = 1, cols = c("gray", "red")) + 
  theme_cowplot() + ggtitle("pri-miR-34a") +
  scale_color_gradientn( colours = c('lightgrey', 'red'),  limits = c(0, 1), oob = scales::squish) + 
  theme(axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(), 
        axis.title = element_blank(),
        plot.title = element_text(hjust = 0.5,size =50),
        legend.position = "none")

ggsave(filename = "/data/rajewsky/projects/cdr1as_ko_snRNA/codes_github/cdr1as/Figures_snRNA/UMAP_pri_mir34a.pdf",
       plot = g1,
       scale = 1, width = 5, height = 5, units = "in", device = cairo_pdf,
       dpi = 300)
###################################################
# Snca, Cplx1, Prkcb, Pfn2, Zdhhc9, Cspa, Basp1, Phactr1, Wipf2
###################################################

# Cplx1
g1 <-  FeaturePlot(all.SO, features = "Cplx1", order = FALSE, pt.size = 1, cols = c("gray", "red")) + 
  theme_cowplot() + ggtitle("Cplx1") +
  scale_color_gradientn( colours = c('lightgrey', 'red'),  limits = c(0, 3), oob = scales::squish) + 
  theme(axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(), 
        axis.title = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 55),
        legend.position = "none")

ggsave(filename = "/data/rajewsky/projects/cdr1as_ko_snRNA/codes_github/cdr1as/Figures_snRNA/UMAP_Cplx1.pdf",
       plot = g1,
       scale = 1, width = 5, height = 5, units = "in", device = cairo_pdf,
       dpi = 300)

# Snca
g1 <-  FeaturePlot(all.SO, features = "Snca", order = FALSE, pt.size = 1, cols = c("gray", "red")) + 
  theme_cowplot() + ggtitle("Snca") +
  scale_color_gradientn( colours = c('lightgrey', 'red'),  limits = c(0, 3), oob = scales::squish) + 
  theme(axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(), 
        axis.title = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 55),
        legend.position = "none")

ggsave(filename = "/data/rajewsky/projects/cdr1as_ko_snRNA/codes_github/cdr1as/Figures_snRNA/UMAP_Snca.pdf",
       plot = g1,
       scale = 1, width = 5, height = 5, units = "in", device = cairo_pdf,
       dpi = 300)

# Prkcb
g1 <-  FeaturePlot(all.SO, features = "Prkcb", order = FALSE, pt.size = 1, cols = c("gray", "red")) + 
  theme_cowplot() + ggtitle("Prkcb") +
  scale_color_gradientn( colours = c('lightgrey', 'red'),  limits = c(0, 3), oob = scales::squish) + 
  theme(axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(), 
        axis.title = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 55),
        legend.position = "none")

ggsave(filename = "/data/rajewsky/projects/cdr1as_ko_snRNA/codes_github/cdr1as/Figures_snRNA/UMAP_Prkcb.pdf",
       plot = g1,
       scale = 1, width = 5, height = 5, units = "in", device = cairo_pdf,
       dpi = 300)

# Pfn2
g1 <-  FeaturePlot(all.SO, features = "Pfn2", order = FALSE, pt.size = 1, cols = c("gray", "red")) + 
  theme_cowplot() + ggtitle("Pfn2") +
  scale_color_gradientn( colours = c('lightgrey', 'red'),  limits = c(0, 3), oob = scales::squish) + 
  theme(axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(), 
        axis.title = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 55),
        legend.position = "none")

ggsave(filename = "/data/rajewsky/projects/cdr1as_ko_snRNA/codes_github/cdr1as/Figures_snRNA/UMAP_Pfn2.pdf",
       plot = g1,
       scale = 1, width = 5, height = 5, units = "in", device = cairo_pdf,
       dpi = 300)

# Zdhhc9
g1 <-  FeaturePlot(all.SO, features = "Zdhhc9", order = FALSE, pt.size = 1, cols = c("gray", "red")) + 
  theme_cowplot() + ggtitle("Zdhhc9") +
  scale_color_gradientn( colours = c('lightgrey', 'red'),  limits = c(0, 3), oob = scales::squish) + 
  theme(axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(), 
        axis.title = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 55),
        legend.position = "none")

ggsave(filename = "/data/rajewsky/projects/cdr1as_ko_snRNA/codes_github/cdr1as/Figures_snRNA/UMAP_Zdhhc9.pdf",
       plot = g1,
       scale = 1, width = 5, height = 5, units = "in", device = cairo_pdf,
       dpi = 300)

# Cspa not exist


# Basp1
g1 <-  FeaturePlot(all.SO, features = "Basp1", order = FALSE, pt.size = 1, cols = c("gray", "red")) + 
  theme_cowplot() + ggtitle("Basp1") +
  scale_color_gradientn( colours = c('lightgrey', 'red'),  limits = c(0, 3), oob = scales::squish) + 
  theme(axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(), 
        axis.title = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 55),
        legend.position = "none")

ggsave(filename = "/data/rajewsky/projects/cdr1as_ko_snRNA/codes_github/cdr1as/Figures_snRNA/UMAP_Basp1.pdf",
       plot = g1,
       scale = 1, width = 5, height = 5, units = "in", device = cairo_pdf,
       dpi = 300)

# Phactr1
g1 <-  FeaturePlot(all.SO, features = "Phactr1", order = FALSE, pt.size = 1, cols = c("gray", "red")) + 
  theme_cowplot() + ggtitle("Phactr1") +
  scale_color_gradientn( colours = c('lightgrey', 'red'),  limits = c(0, 3), oob = scales::squish) + 
  theme(axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(), 
        axis.title = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 55),
        legend.position = "none")

ggsave(filename = "/data/rajewsky/projects/cdr1as_ko_snRNA/codes_github/cdr1as/Figures_snRNA/UMAP_Phactr1.pdf",
       plot = g1,
       scale = 1, width = 5, height = 5, units = "in", device = cairo_pdf,
       dpi = 300)

# Wipf2
g1 <-  FeaturePlot(all.SO, features = "Wipf2", order = FALSE, pt.size = 1, cols = c("gray", "red")) + 
  theme_cowplot() + ggtitle("Wipf2") +
  scale_color_gradientn( colours = c('lightgrey', 'red'),  limits = c(0, 3), oob = scales::squish) + 
  theme(axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(), 
        axis.title = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 55),
        legend.position = "none")

ggsave(filename = "/data/rajewsky/projects/cdr1as_ko_snRNA/codes_github/cdr1as/Figures_snRNA/UMAP_Wipf2.pdf",
       plot = g1,
       scale = 1, width = 5, height = 5, units = "in", device = cairo_pdf,
       dpi = 300)

all_WT.SO <- subset(all.SO, subset =  exp %in% c("WT"))
all_KO.SO <- subset(all.SO, subset =  exp %in% c("FullKO"))
  

# Cdr1os
g1 <-  FeaturePlot(all_WT.SO, features = "Cdr1os", order = FALSE, pt.size = 1, cols = c("gray", "red")) + 
  theme_cowplot() + ggtitle("WT (Cdr1os)") +
  scale_color_gradientn( colours = c('lightgrey', 'red'),  limits = c(0, 4), oob = scales::squish) + 
  theme(axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(), 
        axis.title = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 55),
        legend.position = "none")

ggsave(filename = "/data/rajewsky/projects/cdr1as_ko_snRNA/codes_github/cdr1as/Figures_snRNA/UMAP_WT_Cdr1os.pdf",
       plot = g1,
       scale = 1, width = 5, height = 5, units = "in", device = cairo_pdf,
       dpi = 300)

g1 <-  FeaturePlot(all_KO.SO, features = "Cdr1os", order = FALSE, pt.size = 1, cols = c("gray", "red")) + 
  theme_cowplot() + ggtitle("KO (Cdr1os)") +
  scale_color_gradientn( colours = c('lightgrey', 'red'),  limits = c(0, 4), oob = scales::squish) + 
  theme(axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(), 
        axis.title = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 55),
        legend.position = "none")

ggsave(filename = "/data/rajewsky/projects/cdr1as_ko_snRNA/codes_github/cdr1as/Figures_snRNA/UMAP_KO_Cdr1os.pdf",
       plot = g1,
       scale = 1, width = 5, height = 5, units = "in", device = cairo_pdf,
       dpi = 300)

all.SO@active.ident <- factor(all.SO@active.ident, levels = rev(c("Dentate_Gyrus", "CA1", "CA2_3","subiculum","Neurons",
                                                              "Inhibitory", 
                                                              "Astrocytes", "Microglia", 
                                                              "Oligodendrocytes", "COP", "OPC",
                                                              "Fibroblast", "Cajal_Retzius","Mural", "Choroid_Plexus", "Endothelial")))


features <- c("Slc17a7", "Cdh9", "Hcn1", "Grik4", "Dpp10",
              "Gad2", "Grip1",
              "Gja1",  "Wdr17", 
              "C1qb", "Dock8",
              "Mobp", "Mog",
              "Fyn", "Enpp6",
              "Pdgfra", "Cdo1",
              "Dcn", "Cped1",
              "Ndnf", "Cdh4",
              
              "Slc38a11","Ebf1",
              "Ttr", "Enpp2",
              
              "Mecom", "Ly6c1")
d1 <- DotPlot(all.SO, features = features) + RotatedAxis() +
  scale_color_gradientn( colours = c('lightgrey', 'red'),  limits = c(0, 4), oob = scales::squish) +
  theme(axis.title = element_blank(),
        axis.text = element_text(size = 14))


ggsave(filename = "/data/rajewsky/projects/cdr1as_ko_snRNA/codes_github/cdr1as/Figures_snRNA/dotplot_marker.pdf",
       plot = d1,
       scale = 1, width = 10, height = 4.5, units = "in", device = cairo_pdf,
       dpi = 300)
