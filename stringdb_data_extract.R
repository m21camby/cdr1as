
library(STRINGdb)
library(data.table)
protein_link <- fread("/data/rajewsky/home/skim/Neuronal_Activity_Cledi/201911_m7_overexpression/R_Scripts/external_files/stringdb/10090.protein.links.v11.5.txt.gz")


protein_info <- fread("/data/rajewsky/home/skim/Neuronal_Activity_Cledi/201911_m7_overexpression/R_Scripts/external_files/stringdb/10090.protein.info.v11.5.txt.gz")

protein1_info <- protein_info[,c(1,2)]



input <- readRDS(file = "/data/rajewsky/home/skim/Neuronal_Activity_Cledi/201911_m7_overexpression/publication_codes_figures/DESeq2_and_RF_code_graph_GRN_input.rda")

##############
# subset info
##############

protein_info_sub <- protein_info %>% filter(preferred_name %in% c(input$regulatoryGene, input$targetGene))

##############
# subset link
##############

protein_link_sub <- protein_link %>% filter(protein1 %in% c(protein_info_sub$`#string_protein_id`))

colnames(protein1_info) <- c("protein1", "gene1")
protein_link_sub_1st <- left_join(protein_link_sub, protein1_info, by = "protein1")

colnames(protein1_info) <- c("protein2", "gene2")
protein_link_sub_1st <- left_join(protein_link_sub_1st, protein1_info, by = "protein2")

saveRDS(protein_link_sub_1st, file = "/data/rajewsky/home/skim/Neuronal_Activity_Cledi/201911_m7_overexpression/R_Scripts/20220609_stringdb_data_extract.rda")

  
  