library(edgeR)
library(dplyr)
library(stringr)
library(tidyr)

# FPKM trace back
load("D:/TCGA_P53/data/BRCA_HiSeqcounts_w_TP53.Rdata")
load("D:/TCGA_P53/data/BRCA_HiSeqcounts_wo_TP53.Rdata")
cancer_type = 'BRCA'
#Group set 1 and 2
#combine groups: patients with TP53 and patients wo TP53
Group1and2 <- transform(merge(HiSeq_patients_woTP53mut, HiSeq_patients_withTP53mut, by=0, all=T),
                        row.names = Row.names, Row.names = NULL)
Group1and2.label <- c(rep("without.TP53",times=dim(HiSeq_patients_woTP53mut)[2]),
                      rep("with.TP53",times=dim(HiSeq_patients_withTP53mut)[2]))
Group1and2.label = factor(Group1and2.label,
                          levels = c("without.TP53","with.TP53"))
Group1and2.input <- DGEList(counts=Group1and2, group=Group1and2.label)

Group1and2.input.copy <- Group1and2.input
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

