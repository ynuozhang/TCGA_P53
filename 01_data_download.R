library(dplyr)

gzfile <- "TCGA-BRCA.mutect2_MC3_non_silent.txt.gz"
download.file("https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/mc3_gene_level%2FBRCA_mc3_gene_level.txt.gz", destfile = gzfile)
mutype_file <- read.table( gzfile,
                           header = T,
                           sep = '	',
                           quote = '', 
                           fill = TRUE)
TP53 <- mutype_file[mutype_file$sample == 'tp53' | mutype_file$sample == 'TP53',]
TP53 <- na.omit(TP53)
TP53 <- TP53[, colSums(TP53 != 0) > 0]

cancer_type = 'BRCA'
# Downloads mutation categorical data
mc3 <- paste0("TCGA_", cancer_type, "_mutect2_MC3_public.txt.gz")
download.file("https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/mc3%2FBRCA_mc3.txt.gz", destfile = mc3)
mutype_file_detail <- read.table( mc3,
                                  header = T,
                                  sep = '	',
                                  quote = '', 
                                  fill = TRUE)

TP53_all_mut <- mutype_file_detail[mutype_file_detail$gene == 'tp53' | mutype_file_detail$gene == 'TP53', ]
TP53_R273 <- TP53_all_mut[grep('R273', TP53_all_mut$Amino_Acid_Change), ];dim(TP53_R273);View(TP53_R273)
TP53_R175 <- TP53_all_mut[grep('R175', TP53_all_mut$Amino_Acid_Change), ];dim(TP53_R175);View(TP53_R175)
TP53_R248 <- TP53_all_mut[grep('R248', TP53_all_mut$Amino_Acid_Change), ];dim(TP53_R248);View(TP53_R248)
save(TP53_all_mut, file = paste0('./data/',cancer_type,'_TP53_all','_mutect2.Rdata'))
save(TP53_R273, file = paste0('./data/',cancer_type,'_TP53_R273','_mutect2.Rdata'))
save(TP53_R175, file = paste0('./data/',cancer_type,'_TP53_R175','_mutect2.Rdata'))
save(TP53_R248, file = paste0('./data/',cancer_type,'_TP53_R248','_mutect2.Rdata'))

#Downloads clinical data
#only select primary tumor
clinical_data <- paste0("TCGA_", cancer_type, "phenotype_data.txt.gz")
download.file("https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/TCGA.BRCA.sampleMap%2FBRCA_clinicalMatrix", destfile = clinical_data)
clinical <- read.table(clinical_data,
                       header = T,
                       sep = '	',
                       quote = '',
                       fill = TRUE)
primary.tumor <- clinical[grep('(?i)primary', clinical$sample_type),]
save(primary.tumor, file = paste0('./data/',cancer_type,'_pheno_primary_tumor','.Rdata'))

  

#Downloads TCGA dataset
#HiSeq_data <- paste0("TCGA_", cancer_type, "IlluminaHiSeq_pancan_normalized.txt.gz")
HiSeq_data <- paste0("TCGA_", cancer_type, "IlluminaHiSeq_pancan.txt.gz")
#download.file("https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/TCGA.BRCA.sampleMap%2FHiSeqV2_PANCAN.gz", destfile = HiSeq_data)
download.file("https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/TCGA.BRCA.sampleMap%2FHiSeqV2.gz", destfile = HiSeq_data)
HiSeq <- read.table(HiSeq_data,
                    header = T,
                    sep = '	',
                    quote = '',
                    fill = TRUE)
HiSeq <- rename_with(HiSeq, ~ gsub(".", "-", .x, fixed = TRUE), .cols = 2:last_col())
HiSeq_for_primary_tumor <- select(HiSeq, c(sample,which(colnames(HiSeq) %in% primary.tumor$sampleID)))
save(HiSeq_for_primary_tumor, file = paste0('./data/',cancer_type,'_HiSeqcounts_all_not_normalized','.Rdata'))

#select expression profiles of genes of interest &
#patients contains TP53 mutation
#used SFRS2 instead of SRSF2 for search
gene.of.interest <- c("BRCA1", "SRSF2", "ATF4", "SLC7A11",
                      "ZFAF1", "WT1", "PBRM1", "BAP1",
                      "CDH1", "CIC", "RB1", "STK11", "SMAD4", 
                      "NF1", "TP53",'SFRS2')
#gene.of.interest <- c('PR264', 'SC-35', 'SC35', 'SFRS2', 'SFRS2A','SRp30b')

TP53_sample <- unique(TP53_all_mut$sample)
HiSeq_selected_gene <- filter(HiSeq_for_primary_tumor, sample %in% gene.of.interest); View(HiSeq_selected_gene)

HiSeq_selected_patients_withTP53mut <- select(HiSeq_selected_gene, c(sample, which(colnames(HiSeq_selected_gene) %in% TP53_sample)))
HiSeq_selected_patients_woTP53mut <- select(HiSeq_selected_gene, c(sample, !which(colnames(HiSeq_selected_gene) %in% TP53_sample)))
save(HiSeq_selected_patients_withTP53mut, file = paste0('./data/',cancer_type,'_HiSeqcounts_w_TP53_GOI','.Rdata'))
save(HiSeq_selected_patients_woTP53mut, file = paste0('./data/',cancer_type,'_HiSeqcounts_wo_TP53_GOI','.Rdata'))

HiSeq_selected_patients_R175 <- select(HiSeq_selected_patients_withTP53mut, 
                                       c(sample, which(colnames(HiSeq_selected_patients_withTP53mut) %in% TP53_R175$sample)))
HiSeq_selected_patients_R248 <- select(HiSeq_selected_patients_withTP53mut, 
                                       c(sample, which(colnames(HiSeq_selected_patients_withTP53mut) %in% TP53_R248$sample)))
HiSeq_selected_patients_R273 <- select(HiSeq_selected_patients_withTP53mut, 
                                       c(sample, which(colnames(HiSeq_selected_patients_withTP53mut) %in% TP53_R273$sample)))
save(HiSeq_selected_patients_R175, file = paste0('./data/',cancer_type,'_HiSeqcounts_R175','.Rdata'))
save(HiSeq_selected_patients_R248, file = paste0('./data/',cancer_type,'_HiSeqcounts_R248','.Rdata'))
save(HiSeq_selected_patients_R273, file = paste0('./data/',cancer_type,'_HiSeqcounts_R273','.Rdata'))

