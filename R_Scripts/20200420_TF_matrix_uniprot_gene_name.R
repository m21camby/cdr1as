
library("biomaRt")
library(stringr)

# ref: https://m.ensembl.org/info/data/biomart/biomart_r_package.html


TF_matrix <- read.csv("/data/rajewsky/home/skim/Neuronal_Activity_Cledi/TF_database/TF_target_matrix.txt", sep = "\t", row.names = 2)
TF_matrix$X <- NULL
uniprots <- str_replace(colnames(TF_matrix), ".txt", "") 
colnames(TF_matrix) <- uniprots
#test <- TF_matrix[c(1:10),c(1:10)]


#################################### 
# ENSEMBL to get information 
####################################

mart <- useMart(biomart = "ensembl")
listDatasets(mart)

ensembl <- useEnsembl(biomart= "ensembl", dataset="mmusculus_gene_ensembl")
ensembl.df <- listAttributes(ensembl)
genes.df <- getBM(attributes=c("uniprot_gn_symbol", "uniprot_gn_id", "mgi_symbol"),
                                 filters = "uniprot_gn_id",
                    values = uniprots, mart = ensembl, useCache = FALSE)
genes.df$mgi_symbol[duplicated(genes.df$mgi_symbol)]
table(genes.df$uniprot_gn_id %in% str_replace(colnames(TF_matrix), ".txt", ""))
str_replace(colnames(TF_matrix), ".txt", "")[str_replace(colnames(TF_matrix), ".txt", "") %in% genes.df$uniprot_gn_id]
genes.df <- genes.df[,c("uniprot_gn_id", "mgi_symbol")]


# this is used for below DB
not_in_list <- setdiff(str_replace(colnames(TF_matrix), ".txt", "") , genes.df$uniprot_gn_id)

###################################
# DB download from site and check not_in_list
###################################
uniprot <- read.csv("/data/rajewsky/home/skim/Neuronal_Activity_Cledi/TF_database/uniprot-filtered-organism_Mus+musculus.tab", sep = "\t", stringsAsFactors = FALSE)

uniprot_sub <- subset(uniprot, Entry %in% not_in_list)

# change multiple elements 
uniprot_sub[uniprot_sub == "Usp21 Usp23"] <- "Usp21"
uniprot_sub[uniprot_sub == "Smarca2 Baf190b Brm Snf2a Snf2l2"] <- "Smarca2"
uniprot_sub[uniprot_sub == "Tlx1 Hox11 Tlx-1"] <- "Tlx1"
uniprot_sub[uniprot_sub == "Ubtf Tcfubf Ubf-1 Ubf1"] <- "Ubtf"
uniprot_sub[uniprot_sub == "Neurog2 Ath4a Atoh4 Ngn2"] <- "Neurog2"
uniprot_sub[uniprot_sub == "Npas3 Mop6"] <- "Npas3"
uniprot_sub[uniprot_sub == "Tcf7 Tcf-1 Tcf1"] <- "Tcf7"
uniprot_sub[uniprot_sub == "Taf1 Ccg1"] <- "Taf1"
uniprot_sub[uniprot_sub == "Rpa2 Rpa34"] <- "Rpa2"
uniprot_sub[uniprot_sub == "Tead1 Tcf13 Tef-1 Tef1"] <- "Tead1"
uniprot_sub[uniprot_sub == "Twist2 Dermo1"] <- "Twist2"
uniprot_sub[uniprot_sub == "Tfe3 Tcfe3"] <- "Tfe3"
uniprot_sub[uniprot_sub == "Sox5 Sox-5"] <- "Sox5"
uniprot_sub[uniprot_sub == "Tcf7l1 Tcf3"] <- "Tcf7l1"
uniprot_sub[uniprot_sub == "Runx1 Aml1 Cbfa2 Pebp2ab"] <- "Runx1"
uniprot_sub[uniprot_sub == "Pparg Nr1c3"] <- "Pparg"
uniprot_sub[uniprot_sub == "Supt16h Fact140 Factp140 Supt16"] <- "Supt16h"
uniprot_sub[uniprot_sub == "Pex2 Paf1 Pmp35 Pxmp3"] <- "Pex2"
uniprot_sub[uniprot_sub == "Tcf7l2 Tcf4"] <- "Tcf7l2"
uniprot_sub[uniprot_sub == "Usp7 Hausp"] <- "Usp7"
uniprot_sub[uniprot_sub == "Klf6 Copeb Cpbp"] <- "Klf6"
uniprot_sub[uniprot_sub == "Phc1 Edr Edr1 Rae28"] <- "Phc1"
uniprot_sub[uniprot_sub == "Mybl1 Amyb"] <- "Mybl1"
uniprot_sub[uniprot_sub == "Otx2 Otx-2"] <- "Otx2"
uniprot_sub[uniprot_sub == "Setdb1 Eset Kiaa0067"] <- "Setdb1"
uniprot_sub[uniprot_sub == "Sall1 Sal3"] <- "Sall1"
uniprot_sub[uniprot_sub == "Sox3 Sox-3"] <- "Sox3"
uniprot_sub[uniprot_sub == "Hey1 Herp2 Hesr1 Hrt1"] <- "Hey1"
uniprot_sub[uniprot_sub == "Sox2 Sox-2"] <- "Sox2"
uniprot_sub[uniprot_sub == "Thoc5 Fmip Kiaa0983"] <- "Thoc5"
uniprot_sub[uniprot_sub == "Neo1 Ngn"] <- "Neo1"
uniprot_sub[uniprot_sub == "Pou2f1 Oct-1 Otf-1 Otf1"] <- "Pou2f1"
uniprot_sub[uniprot_sub == "Parp1 Adprp Adprt Adprt1"] <- "Parp1"
uniprot_sub[uniprot_sub == "Nelfb Cobra1 MNCb-5210"] <- "Nelfb"
uniprot_sub[uniprot_sub == "Pou2f3 Epoc1 Oct11 Otf-11 Otf11"] <- "Pou2f3"
uniprot_sub[uniprot_sub == "Nfatc1 Nfat2 Nfatc"] <- "Nfatc1"
uniprot_sub[uniprot_sub == "Rorc Nr1f3 Rorg Thor"] <- "Rorc"
uniprot_sub[uniprot_sub == "Zfp42 Rex-1 Rex1"] <- "Zfp42"
uniprot_sub[uniprot_sub == "Wt1 Wt-1"] <- "Wt1"
uniprot_sub[uniprot_sub == "Trim33 Kiaa1113"] <- "Trim33"
uniprot_sub[uniprot_sub == "Spen Kiaa0929 Mint Sharp"] <- "Spen"
uniprot_sub[uniprot_sub == "Nr3c1 Grl Grl1"] <- "Nr3c1"
uniprot_sub[uniprot_sub == "Runx2 Aml3 Cbfa1 Osf2 Pebp2a"] <- "Runx2"
uniprot_sub[uniprot_sub == "Rarb Nr1b2"] <- "Rarb"
uniprot_sub[uniprot_sub == "Nsd3 Whsc1l1"] <- "Nsd3"
uniprot_sub[uniprot_sub == "Prmt5 Jbp1 Skb1"] <- "Prmt5"
uniprot_sub[uniprot_sub == "Kdm3b Jhdm2b Jmjd1b Kiaa1082"] <- "Kdm3b"
uniprot_sub[uniprot_sub == "Kdm2a Fbl11 Fbxl11 Jhdm1a Kiaa1004"] <- "Kdm2a"
uniprot_sub[uniprot_sub == "Klf1 Elkf"] <- "Klf1"
uniprot_sub[uniprot_sub == "Gfi1 Gfi-1"] <- "Gfi1"
uniprot_sub[uniprot_sub == "Ep400 Kiaa1498"] <- "Ep400"
uniprot_sub[uniprot_sub == "Htt Hd Hdh"] <- "Htt"
uniprot_sub[uniprot_sub == "Crebbp Cbp"] <- "Crebbp"
uniprot_sub[uniprot_sub == "Hdac2 Yy1bp"] <- "Hdac2"
uniprot_sub[uniprot_sub == "Hnf4g Nr2a2"] <- "Hnf4g"
uniprot_sub[uniprot_sub == "Arntl Bmal1"] <- "Arntl"
uniprot_sub[uniprot_sub == "Ncapd3 Kiaa0056"] <- "Ncapd3"
uniprot_sub[uniprot_sub == "Dzip3 Kiaa0675"] <- "Dzip3"
uniprot_sub[uniprot_sub == "Capg Mbh1"] <- "Capg"
uniprot_sub[uniprot_sub == "Tfap2a Ap2tf Tcfap2a"] <- "Tfap2a"
uniprot_sub[uniprot_sub == "Ddx5 Tnz2"] <- "Ddx5"
uniprot_sub[uniprot_sub == "Ikzf1 Ikaros Lyf1 Zfpn1a1 Znfn1a1"] <- "Ikzf1"
uniprot_sub[uniprot_sub == "Hoxd11 Hox-4.6 Hoxd-11"] <- "Hoxd11"
uniprot_sub[uniprot_sub == "Ncor1 Rxrip13"] <- "Ncor1"
uniprot_sub[uniprot_sub == "H1-6 H1f6 H1ft H1t Hist1h1t"] <- "H1-6"
uniprot_sub[uniprot_sub == "Kmt2c Mll3"] <- "Kmt2c"
uniprot_sub[uniprot_sub == "Med23 Crsp3 Kiaa1216 Sur2"] <- "Med23"
uniprot_sub[uniprot_sub == "Cbfa2t3 Cbfa2t3h Mtgr2"] <- "Cbfa2t3"
uniprot_sub[uniprot_sub == "Ncoa3 Aib1 Pcip Rac3 Tram1"] <- "Ncoa3"
uniprot_sub[uniprot_sub == "Nelfe D17h6s45 Rd Rdbp"] <- "Nelfe"
uniprot_sub[uniprot_sub == "Dux Dux4"] <- "Dux"
uniprot_sub[uniprot_sub == "Foxd3 Hfh2"] <- "Foxd3"

uniprot_sub <- uniprot_sub[, c("Entry", "Gene.names")]

###############################
# merge ensembl and uniprot
###############################

colnames(genes.df) <- c("uniprot", "gene")
colnames(uniprot_sub) <- c("uniprot", "gene")

merged.df <- rbind(genes.df, uniprot_sub)

#check duplicated
merged.df$uniprot[duplicated(merged.df$uniprot)]
# [1] "P09602" "P50247" "Q3TTC2"
# remove
merged.df <- merged.df[!merged.df$gene %in% c("Gm10282", "Gm4737", "Mbtps2"), ]

merged.df$gene[duplicated(merged.df$gene)]
# [1]  "Zbtb16"
# remove
merged.df <- merged.df[!merged.df$uniprot %in% c("A3KMN0"), ]


# replace column name

# identify unique id of uniprot
unique_uniprot <- intersect(colnames(TF_matrix), merged.df$uniprot)

# subset matrix
TF_matrix_2nd <- TF_matrix[ ,colnames(TF_matrix) %in% unique_uniprot]

# subset merged.df
merged.df <- merged.df[merged.df$uniprot %in% unique_uniprot, ]

colnames(TF_matrix_2nd) <- merged.df$gene

#################################### 
# ENSEMBL to get information of row
####################################

genes2.df <- getBM(attributes=c("ensembl_gene_id", "mgi_symbol"),
                  filters = "ensembl_gene_id",
                  values = rownames(TF_matrix_2nd), mart = ensembl, useCache = FALSE)

# remove empty gene
genes2.df <- genes2.df[!genes2.df$mgi_symbol %in% c(""),]

#check duplicated
genes2.df$ensembl_gene_id[duplicated(genes2.df$ensembl_gene_id)]

#check duplicated
genes2.df$mgi_symbol[duplicated(genes2.df$mgi_symbol)]
#  [1] "Pakap"         "A530058N18Rik" "Pakap"         "Snora43"       "1700030C10Rik" "Gm23925"       "Pcdha11"       "Gm25203"       "Gm26265"       "Gm4430"        "4930594M22Rik" "Gm16364"       "Gm36638"       "Snord3b-ps1"   "Gm16701"      
# [16] "Gm28724"       "Gm7270"        "Gm23377"       "Gm20690"       "Gm18433"       "Rprl1"      

# remove duplicated genes
genes2.df <- genes2.df[!genes2.df$ensembl_gene_id %in% c("ENSMUSG00000090053","ENSMUSG00000089945","ENSMUSG00000092201","ENSMUSG00000097052", "ENSMUSG00000099759","ENSMUSG00000100992", 
                                                    "ENSMUSG00000102206", "ENSMUSG00000105037", "ENSMUSG00000105690", "ENSMUSG00000103532", "ENSMUSG00000085423","ENSMUSG00000087014",
                                                    "ENSMUSG00000102635", "ENSMUSG00000104921", "ENSMUSG00000102548","ENSMUSG00000101909", "ENSMUSG00000103218","ENSMUSG00000105187",
                                                    "ENSMUSG00000103007", "ENSMUSG00000104144", "ENSMUSG00000106222"),]


# replace rownames to gene name
TF_matrix_2nd$ensembl_gene_id <- rownames(TF_matrix_2nd)

TF_matrix_final <- inner_join(TF_matrix_2nd, genes2.df, by = "ensembl_gene_id")
rownames(TF_matrix_final) <- TF_matrix_final$mgi_symbol

TF_matrix_final$mgi_symbol <- NULL
TF_matrix_final$ensembl_gene_id <- NULL

saveRDS(TF_matrix_final, file = "/data/rajewsky/home/skim/Neuronal_Activity_Cledi/TF_database/20200420_TF_matrix_uniprot_gene_name_Final_matrix.rda")

# TF_matrix_final <- readRDS(file = "/data/rajewsky/home/skim/Neuronal_Activity_Cledi/TF_database/20200420_TF_matrix_uniprot_gene_name_Final_matrix.rda")

