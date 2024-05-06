# Cdr1as project

* This project is to understand `Cdr1as` and `miR-7` in mouse neurons
* This project co-works with `Cledi A. Cerda Jara`
* 4 data sets: 1)**WT** control, 2)**WT m7** overexpression, 3)**Cdr1as KO**, 4)**Cdr1as KO m7** overexpression


## 1. Data access
* processed data is available from NCBI GEO under accession: GSE224184
* You can check a Paper: [here](https://www.biorxiv.org/content/10.1101/2023.01.26.525729v1)

## 2. Paper overview

<img src="https://www.biorxiv.org/content/biorxiv/early/2023/01/26/2023.01.26.525729/F7.large.jpg">


## 3. Repository structures

* `DGEm` = Gene expression matrix used in the paper 
* `Figures` = Figures used in the paper 
* `Figures_snRNA` = Figures used in the project (Not in the paper) 
* `R_Scripts` = R codes used in the paper
* `Scripts` = bash and python scripts used in the paper
* `Results` = Differential expression and Gene ontology and Gene regulatory network used in the paper
* `ext_files` = External files used in the paper
* `srna` = small RNA analysis codes used in the paper


## 4. Ppreprocess 
* The sequencing had been done by NextSeq500 from single-end PolyA selected libraries (1x150) 

* Extract miR-7 targets gene lists: `cut -f1 TargetScan7.1__miR-7-5p.predicted_targets.txt | sed '1d' > TargetScan7.1__miR-7-5p.predicted_targets_lists.txt` 

* Check read strandness: `infer_experiment.py -r mm10_RefSeq.bed -i WTC1_S1sorted_Aligned.out.bam`

* Check read distribution: `read_distribution.py -r mm10_RefSeq.bed -i WTC1_S1sorted_Aligned.out.bam` 

* Mapping: 
		```
                STAR --runThreadN 24
	     	 --genomeDir /mm10_GRCm38.p6_star_2.7.0a\
		 --readFilesIn /sequencing/mouse/191031_NB501326_0339_AHCM7WBGXC/fastq/WTC1_S1_R1_001.fastq.gz\
		 --outFileNamePrefix WTC1_S1_quant\
		 --sjdbGTFfile /annotation/GRCm38/M21/GRCm38.M21.gtf\
 		 --readFilesCommand zcat\
		 --quantMode GeneCounts &
		```


 

