library(dplyr)

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
TP53_R273 <- TP53_all_mut[grep('R273', TP53_all_mut$Amino_Acid_Change), ];dim(TP53_R273)
TP53_R175 <- TP53_all_mut[grep('R175', TP53_all_mut$Amino_Acid_Change), ];dim(TP53_R175)
TP53_R248 <- TP53_all_mut[grep('R248', TP53_all_mut$Amino_Acid_Change), ];dim(TP53_R248)

#check if any sample may contain over 1 TP53 mutation
n_occur_TP53 <- data.frame(table(TP53_all_mut$sample))
n_occur_TP53_sample <- n_occur_TP53[n_occur_TP53$Freq > 1,]

if (length(n_occur_TP53_sample) != 0) {
  sample.with.multiple.TP53.mut <- n_occur_TP53_sample$Var1
  TP53_all_mut <- mutate(TP53_all_mut, category = case_when(
    TP53_all_mut$sample %in% sample.with.multiple.TP53.mut ~ 'Multiple_Mut',
    TRUE ~ 'Single_Mut'))
  print(paste0(length(sample.with.multiple.TP53.mut), ' samples contains more than 1 TP53 mutation'))
  verbose.or.not <- readline(prompt="check these samples in detail (only accept y/n): ")
  if (verbose.or.not == 'y'){
    print(n_occur_TP53_sample)
  } else{
    invisible()
  }
} else {
  print('no sample contains more than 1 TP53 mutation')
}
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
rownames(HiSeq) <- HiSeq[,1]
HiSeq <- HiSeq[,-1]

HiSeq <- round(((2**HiSeq) - 1), 0)

HiSeq <- rename_with(HiSeq, ~ gsub(".", "-", .x, fixed = TRUE))
HiSeq_for_primary_tumor <- select(HiSeq, which(colnames(HiSeq) %in% primary.tumor$sampleID))
save(HiSeq_for_primary_tumor, file = paste0('./data/',cancer_type,'_HiSeqcounts_all_not_normalized','.Rdata'))

#select expression profiles of genes of interest &
#patients contains TP53 mutation
#used SFRS2 instead of SRSF2 for search
gene.of.interest <- c("BRCA1", "SRSF2", "ATF4", "SLC7A11",
                      "ZFAF1", "WT1", "PBRM1", "BAP1",
                      "CDH1", "CIC", "RB1", "STK11", "SMAD4", 
                      "NF1", "TP53",'SFRS2')
#gene.of.interest <- c('PR264', 'SC-35', 'SC35', 'SFRS2', 'SFRS2A','SRp30b')

HiSeq_selected_gene <- filter(HiSeq_for_primary_tumor, rownames(HiSeq_for_primary_tumor) %in% gene.of.interest); View(HiSeq_selected_gene)

TP53_sample <- filter(TP53_all_mut, category == 'Single_Mut')$sample
TP53_multi_sample <- filter(TP53_all_mut, category == 'Multiple_Mut')$sample

HiSeq_selected_patients_withTP53mut <- select(HiSeq_selected_gene, which(colnames(HiSeq_selected_gene) %in% TP53_sample))
HiSeq_selected_patients_withMultipleTP53mut <- select(HiSeq_selected_gene, 
                                             which(colnames(HiSeq_selected_gene) %in% TP53_multi_sample))
HiSeq_selected_patients_woTP53mut <- select(HiSeq_selected_gene, which(!(colnames(HiSeq_selected_gene) %in% TP53_sample)
                                                                  & 
                                                                    !(colnames(HiSeq_selected_gene) %in% TP53_multi_sample)))
HiSeq_selected_patients_R175 <- select(HiSeq_selected_gene, 
                              which(colnames(HiSeq_selected_gene) %in% TP53_R175$sample))
HiSeq_selected_patients_R248 <- select(HiSeq_selected_gene, 
                              which(colnames(HiSeq_selected_gene) %in% TP53_R248$sample))
HiSeq_selected_patients_R273 <- select(HiSeq_selected_gene, 
                              which(colnames(HiSeq_selected_gene) %in% TP53_R273$sample))
HiSeq_selected_patients_R <- select(HiSeq_selected_gene, 
                           which((colnames(HiSeq_selected_gene) %in% TP53_R175$sample) | 
                                   (colnames(HiSeq_selected_gene) %in% TP53_R248$sample) |
                                   (colnames(HiSeq_selected_gene) %in% TP53_R273$sample)))
save(HiSeq_selected_patients_withTP53mut, HiSeq_selected_patients_withMultipleTP53mut, HiSeq_selected_patients_woTP53mut,
     HiSeq_selected_patients_R175,HiSeq_selected_patients_R248,HiSeq_selected_patients_R273,HiSeq_selected_patients_R,
     file = paste0('./data/','Step1_selectedgene','_HiSeqcounts','.Rdata'))


