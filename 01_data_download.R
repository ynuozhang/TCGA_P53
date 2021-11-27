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
clinical_data <- paste0("TCGA_", cancer_type, "phenotype_curated_survival_data.txt.gz")
download.file("https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/survival%2FBRCA_survival.txt", destfile = clinical_data)
clinical <- read.table(clinical_data,
                       header = T,
                       sep = '	',
                       quote = '',
                       fill = TRUE)
save(clinical, file = paste0('./data/',cancer_type,'_clinical','.Rdata'))


#Downloads TCGA dataset
HiSeq_data <- paste0("TCGA_", cancer_type, "IlluminaHiSeq_pancan_normalized.txt.gz")
download.file("https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/TCGA.BRCA.sampleMap%2FHiSeqV2_PANCAN.gz", destfile = HiSeq_data)
HiSeq <- read.table(HiSeq_data,
                    header = T,
                    sep = '	',
                    quote = '',
                    fill = TRUE)
save(HiSeq, file = paste0('./data/',cancer_type,'_HiSeqcounts_all','.Rdata'))

#select expression profiles of genes of interest &
#patients contains TP53 mutation

gene.of.interest <- c("BRCA1", "SRSF2", "ATF4", "SLC7A11",
                      "ZFAF1", "WT1", "PBRM1", "BAP1",
                      "CDH1", "CIC", "RB1", "STK11", "SMAD4", 
                      "NF1", "TP53")

TP53_sample <- unique(TP53_all_mut$sample)
HiSeq_selected_gene <- filter(HiSeq, sample %in% gene.of.interest); View(HiSeq_selected_gene)
HiSeq_selected_gene <- rename_with(HiSeq_selected_gene, ~ gsub(".", "-", .x, fixed = TRUE))

HiSeq_selected_patients_withTP53mut <- select(HiSeq_selected_gene, which(colnames(HiSeq_selected_gene) %in% TP53_sample))
save(HiSeq_selected_patients_withTP53mut, file = paste0('./data/',cancer_type,'_HiSeqcounts_w_TP53_GOI','.Rdata'))
