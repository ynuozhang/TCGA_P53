library(edgeR)
library(dplyr)
library(stringr)
library(tidyr)

# FPKM trace back
#load("D:/TCGA_P53/data/BRCA_HiSeqcounts_w_TP53.Rdata")
#load("D:/TCGA_P53/data/BRCA_HiSeqcounts_wo_TP53.Rdata")
load("D:/Github/TCGA_P53/data/BRCA_HiSeqcounts.Rdata")
cancer_type = 'BRCA'
#Group set 1 and 2
#combine all groups

Group.all <- Reduce(function(a,b){
  temp <- merge(a,b, all=TRUE, by = 0)
  row.names(temp) <- temp[,'Row.names']
  temp[,!names(temp) %in% 'Row.names']
}, list(HiSeq_patients_woTP53mut, HiSeq_patients_withTP53mut, HiSeq_patients_withMultipleTP53mut,
        HiSeq_patients_R175, HiSeq_patients_R248, HiSeq_patients_R273))

Group.all.label <- c(rep("TP53 WT",times=dim(HiSeq_patients_woTP53mut)[2]),
                     rep("TP53 mutant",times=dim(HiSeq_patients_withTP53mut)[2]), 
                     rep("TP53 multiple mutant",times=dim(HiSeq_patients_withMultipleTP53mut)[2]),
                     rep("TP53 R175",times=dim(HiSeq_patients_R175)[2]),
                     rep("TP53 R248",times=dim(HiSeq_patients_R248)[2]),
                     rep("TP53 R273",times=dim(HiSeq_patients_R273)[2]))

Group.all.label = factor(Group.all.label,
                         levels = c("TP53 WT","TP53 mutant","TP53 multiple mutant",
                                    "TP53 R175", "TP53 R248","TP53 R273"))
Group.all.input <- DGEList(counts=Group.all, group=Group.all.label)

#Group1and2.label <- c(rep("TP53 WT",times=dim(HiSeq_patients_woTP53mut)[2]),
                      #rep("TP53 mutant",times=dim(HiSeq_patients_withTP53mut)[2]))
#Group1and2.label = factor(Group1and2.label,
                          #levels = c("TP53 WT","TP53 mutant"))
#Group1and2.input <- DGEList(counts=Group1and2, group=Group1and2.label)

#Group1and2.input.copy <- Group1and2.input
#filtering for less representative genes
#min(apply(Group1and2.input.copy$counts, 2, sum))
#max(apply(Group1and2.input.copy$counts, 2, sum))



keep <- rowSums(cpm(Group1and2.input.copy)>100) >= 266
Group1and2.input.copy <- Group1and2.input.copy[keep,]
Group1and2.input.copy$samples$lib.size <- colSums(Group1and2.input.copy$counts)
Group1and2.input.copy$samples

#start normalization by total count
#"The calcNormFactors() function normalizes for RNA composition 
#by finding a set of scaling factors for the library sizes that 
#minimize the log-fold changes between the samples for most genes." 
Group1and2.input.copy <- calcNormFactors(Group1and2.input.copy)

#Estimating the Dispersion
#Naive model
d1 <- estimateCommonDisp(Group1and2.input.copy, verbose=T)
#output Disp = 0.32147 , BCV = 0.567
#Bayes tagwise dispersion
d1 <- estimateTagwiseDisp(d1)
plotBCV(d1)
#GLM estimation of dispersion
design.mat <- model.matrix(~ 0 + Group1and2.input.copy$samples$group)
rownames(design.mat) <- colnames(Group1and2.input.copy)
colnames(design.mat) <- levels(Group1and2.input.copy$samples$group)
d2 <- estimateGLMCommonDisp(Group1and2.input.copy,design.mat)
d2 <- estimateGLMTrendedDisp(d2,design.mat)
# You can change method to "auto", "bin.spline", "power", "spline", "bin.loess".
# The default is "auto" which chooses "bin.spline" when > 200 tags and "power" otherwise.
d2 <- estimateGLMTagwiseDisp(d2,design.mat)
plotBCV(d2)

fit <- glmFit(d2, design.mat)
results <- glmLRT(fit, contrast = c(-1, 1)) 
nredgeR <- topTags(results, n = nrow(d2))
nredgeR <- as.data.frame(nredgeR)
colnames(nredgeR)[4] <- c("P.Value")

#data.dir <- file.path('./TCGA_P53')
source('00_drawing_function.R')

draw_heatmap(V.method = nredgeR, gene.dataset = Group1and2, 
             annotation.class = factor(Group1and2.label), name = 'BRCA_edgeR')

#check if gene of interest is within top100
gene.of.interest <- c("BRCA1", "SRSF2", "ATF4", "SLC7A11",
                      "ZFAF1", "WT1", "PBRM1", "BAP1",
                      "CDH1", "CIC", "RB1", "STK11", "SMAD4", 
                      "NF1", "TP53",'SFRS2')
look.for.gene.of.interest <- filter(HiSeq_for_primary_tumor, rownames(choose_matrix) %in% gene.of.interest)
if (look.for.gene.of.interest) {
  print('Gene of interest is within top100 regulated gene.')
} else {
  print('Gene of interest is within top100 regulated gene.
        Check if G.O.I is within other high expression gene')
  another.check <- nredgeR[rownames(nredgeR) %in% gene.of.interest]
  if (another.check) {
    print(another.check)
  } else {
    print('G.O.I is not within high epression gene. They are filtered before performing edgeR.')
  }
}
  



###################
library(ComplexHeatmap)
nrDEG_Z = nredgeR[order(nredgeR$logFC), ]
nrDEG_F = nredgeR[order(-nredgeR$logFC), ]
choose_gene = c(rownames(nrDEG_Z)[1:50], rownames(nrDEG_F)[1:50])
choose_matrix = Group1and2[choose_gene,]
choose_matrix = t(scale(t(choose_matrix)))

choose_matrix[choose_matrix > 2] = 2
choose_matrix[choose_matrix < -2] = -2

#annotation_col = data.frame(CellType = factor(Group1and2.label))
rownames(annotation_col) = colnames(Group1and2)
library(circlize)
col_fun = colorRamp2(c(-2, 0, 2), c("#2fa1dd", "white", "#f87669"))
top_annotation = HeatmapAnnotation(
  cluster = anno_block(gp = gpar(fill = c("#2fa1dd", "#f87669")),
                       labels = c('without.TP53','with.TP53'),
                       labels_gp = gpar(col = "white", fontsize = 12)))

Heatmap(choose_matrix, col = col_fun, top_annotation = top_annotation,
        show_heatmap_legend = T, column_split = factor(Group1and2.label),
        show_row_names = T, show_column_names = F,column_title = NULL, 
        cluster_columns = F, cluster_column_slices = F)

#################
gene.of.interest <- c("BRCA1", "SRSF2", "ATF4", "SLC7A11",
                      "ZFAF1", "WT1", "PBRM1", "BAP1",
                      "CDH1", "CIC", "RB1", "STK11", "SMAD4", 
                      "NF1", "TP53",'SFRS2')
HiSeq_selected_gene <- data.frame(GeneOfInterest = )
  filter(HiSeq_for_primary_tumor, rownames(choose_matrix) %in% gene.of.interest)



draw_heatmap(nrDEG = nrDEG_edgeR, type = 'edgeR')
draw_volcano(nrDEG = nrDEG_edgeR, type = 'edgeR')

