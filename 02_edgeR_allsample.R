library(edgeR)
library(dplyr)
library(stringr)
library(tidyr)

# FPKM trace back
load("D:/Github/TCGA_P53/data/BRCA_HiSeqcounts.Rdata")
cancer_type = 'BRCA'
#combine all groups

Group.all <- Reduce(function(a,b){
  temp <- merge(a,b, all=TRUE, by = 0)
  row.names(temp) <- temp[,'Row.names']
  temp[,!names(temp) %in% 'Row.names']
}, list(HiSeq_patients_woTP53mut, HiSeq_patients_withTP53mut, HiSeq_patients_withMultipleTP53mut,
        HiSeq_patients_R175, HiSeq_patients_R248, HiSeq_patients_R273))

Group.all.label <- c(rep("TP53_WT",times=dim(HiSeq_patients_woTP53mut)[2]),
                     rep("TP53_mutant",times=dim(HiSeq_patients_withTP53mut)[2]), 
                     rep("TP53_multiple_mutant",times=dim(HiSeq_patients_withMultipleTP53mut)[2]),
                     rep("TP53_R175",times=dim(HiSeq_patients_R175)[2]),
                     rep("TP53_R248",times=dim(HiSeq_patients_R248)[2]),
                     rep("TP53_R273",times=dim(HiSeq_patients_R273)[2]))

Group.all.label = factor(Group.all.label,
                         levels = c("TP53_WT","TP53_mutant","TP53_multiple_mutant",
                                    "TP53_R175", "TP53_R248","TP53_R273"))
Group.all.input <- DGEList(counts=Group.all, group=Group.all.label)

#filtering for less representative genes
Group.all.input.copy <- Group.all.input
keep <- filterByExpr(Group.all.input.copy)
Group.all.input.copy <- Group.all.input.copy[keep,]
Group.all.input.copy$samples$lib.size <- colSums(Group.all.input.copy$counts)

#start normalization by total count
#"The calcNormFactors() function normalizes for RNA composition 
#by finding a set of scaling factors for the library sizes that 
#minimize the log-fold changes between the samples for most genes." 
Group.all.input.copy <- calcNormFactors(Group.all.input.copy)

#Estimating the Dispersion
#GLM estimation of dispersion
design.mat <- model.matrix(~ 0 + Group.all.input.copy$samples$group)
rownames(design.mat) <- colnames(Group.all.input.copy)
colnames(design.mat) <- levels(Group.all.input.copy$samples$group)
d2 <- estimateGLMCommonDisp(Group.all.input.copy,design.mat)
d2 <- estimateGLMTrendedDisp(d2,design.mat)
d2 <- estimateGLMTagwiseDisp(d2,design.mat)

#specify pair-wise comparison
CONTRASTS <- makeContrasts( Group1vs2 = TP53_WT-TP53_mutant,
                            Group1vs3 = TP53_WT-TP53_multiple_mutant,
                            Group1vs4 = TP53_WT-TP53_R175,
                            Group1vs5 = TP53_WT-TP53_R248,
                            Group1vs6 = TP53_WT-TP53_R273,
                            levels = design.mat )

fit <- glmFit(d2, design.mat)
glm.Group1vs2 <- glmLRT(fit, contrast = CONTRASTS[,"Group1vs2"])
glm.Group1vs3 <- glmLRT(fit, contrast = CONTRASTS[,"Group1vs3"])
glm.Group1vs4 <- glmLRT(fit, contrast = CONTRASTS[,"Group1vs4"])
glm.Group1vs5 <- glmLRT(fit, contrast = CONTRASTS[,"Group1vs5"])
glm.Group1vs6 <- glmLRT(fit, contrast = CONTRASTS[,"Group1vs6"])

output.Group1vs2 <- as.data.frame(topTags(glm.Group1vs2, n = nrow(d2)))
output.Group1vs3 <- as.data.frame(topTags(glm.Group1vs3, n = nrow(d2)))
output.Group1vs4 <- as.data.frame(topTags(glm.Group1vs4, n = nrow(d2)))
output.Group1vs5 <- as.data.frame(topTags(glm.Group1vs5, n = nrow(d2)))
output.Group1vs6 <- as.data.frame(topTags(glm.Group1vs6, n = nrow(d2)))

colnames(output.Group1vs2)[4] <- c("P.Value")
colnames(output.Group1vs3)[4] <- c("P.Value")
colnames(output.Group1vs4)[4] <- c("P.Value")
colnames(output.Group1vs5)[4] <- c("P.Value")
colnames(output.Group1vs6)[4] <- c("P.Value")


#data.dir <- file.path('./TCGA_P53')
source('00_drawing_function.R')

Group1and2.label <- c(rep("TP53_WT",times=dim(HiSeq_patients_woTP53mut)[2]),
                      rep("TP53_mutant",times=dim(HiSeq_patients_withTP53mut)[2]))

Group1and2.label <- factor(Group1and2.label,
                         levels = c("TP53_WT","TP53_mutant"))

Group1and2 <- transform(merge(HiSeq_patients_woTP53mut, HiSeq_patients_withTP53mut, by=0, all=T),
                                     row.names = Row.names, Row.names = NULL)

draw_heatmap(V.method = output.Group1vs2, gene.dataset = Group1and2, 
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
nrDEG_Z = output.Group1vs2[order(output.Group1vs2$P.Value), ]

nrDEG_F = output.Group1vs2[order(-output.Group1vs2$P.Value), ]
choose_gene = c(rownames(nrDEG_Z)[1:50], rownames(nrDEG_F)[1:50])
choose_matrix = Group1and2[choose_gene,]
choose_matrix = log2(choose_matrix + 1)
choose_matrix = t(scale(t(choose_matrix)))

#choose_matrix[choose_matrix > 2] = 2
#choose_matrix[choose_matrix < -2] = -2

#annotation_col = data.frame(CellType = factor(Group1and2.label))
#rownames(annotation_col) = colnames(Group1and2)
library(circlize)
col_fun = colorRamp2(c(-2, 0, 2), c("#2fa1dd", "white", "#f87669"))
top_annotation = HeatmapAnnotation(
  cluster = anno_block(gp = gpar(fill = c("#2fa1dd", "#f87669")),
                       labels = c('without.TP53','with.TP53'),
                       labels_gp = gpar(col = "white", fontsize = 12)))

Heatmap(choose_matrix, col = col_fun, top_annotation = top_annotation,
        show_heatmap_legend = T, column_split = factor(Group1and2.label),
        show_row_names = T, show_column_names = F,column_title = NULL, 
        cluster_columns = F, cluster_column_slices = F, row_names_gp = gpar(fontsize = 6))

#################
gene.of.interest <- c("BRCA1", "SRSF2", "ATF4", "SLC7A11",
                      "ZFAF1", "WT1", "PBRM1", "BAP1",
                      "CDH1", "CIC", "RB1", "STK11", "SMAD4", 
                      "NF1", "TP53",'SFRS2')
HiSeq_selected_gene <- data.frame(GeneOfInterest = )
  filter(HiSeq_for_primary_tumor, rownames(choose_matrix) %in% gene.of.interest)



draw_heatmap(nrDEG = nrDEG_edgeR, type = 'edgeR')
draw_volcano(nrDEG = nrDEG_edgeR, type = 'edgeR')

