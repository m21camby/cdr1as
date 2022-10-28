
library(Seurat , lib = "/data/local/rajewsky/home/skim/R/usr_lib_Seurat/")


KO_data <- Read10X("/data/rajewsky/projects/cdr1as_ko_snRNA/3rd_sequencing_run/SP081_aggr/outs/count/filtered_feature_bc_matrix/")

pri_miRNA <- Read10X("/data/rajewsky/projects/cdr1as_ko_snRNA/3rd_sequencing_run/FullKO_primiRNA_aggr_final/outs/count/filtered_feature_bc_matrix/")

# keep only intergenic counts
host_biotype <- read.csv("/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/miRNA_lncRNA/miRNA_gtf_whole_pri_miRNA_pre_miRNA/Extract_miRNA_host_gene_biotypes.csv", row.names = 1)
pri_miRNA_intergenic <- pri_miRNA[!rownames(pri_miRNA) %in% host_biotype$miRNA_ID, ]

# merge two counts
merged <- rbind(KO_data, pri_miRNA_intergenic)

all.SO <- CreateSeuratObject(counts =  merged,  min.cells = 3, min.features = 200, project = "FullKO")
# calculate MT genes
all.SO[["percent.mt"]] <- PercentageFeatureSet(object = all.SO, pattern = "^mt-")
gemgroup <- sapply(strsplit(rownames(all.SO@meta.data), split="-"), "[[", 2)
all.SO <- AddMetaData(object=all.SO, metadata=data.frame(gemgroup=gemgroup, row.names=rownames(all.SO@meta.data)))
all.SO <- subset(x = all.SO, subset = nCount_RNA < 30000 & nCount_RNA > 500 & percent.mt < 5)

all.meta.data <- all.SO@meta.data
all.meta.data <- all.meta.data %>% mutate(exp = case_when(gemgroup %in% c(1, 2, 5) ~ "WT",
                                                          gemgroup %in% c(3, 4, 6) ~ "FullKO"))

all.SO$exp <- all.meta.data$exp

all.SO <- SCTransform(all.SO, vars.to.regress = c("percent.mt", "nCount_RNA"), verbose = TRUE)
all.SO <- RunPCA(all.SO, verbose = FALSE)
all.SO <- FindNeighbors(object = all.SO, reduction = "pca", dims = 1:40)
all.SO <- FindClusters(object = all.SO, resolution = 1, verbose = FALSE)
all.SO <- RunUMAP(all.SO, dims = 1:40, verbose = FALSE)

