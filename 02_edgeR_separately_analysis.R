library(edgeR)
library(dplyr)
library(stringr)
library(tidyr)
library(circlize)
library(ComplexHeatmap)

# FPKM trace back
#load("D:/TCGA_P53/data/BRCA_HiSeqcounts_w_TP53.Rdata")
#load("D:/TCGA_P53/data/BRCA_HiSeqcounts_wo_TP53.Rdata")
load("D:/Github/TCGA_P53/data/BRCA_HiSeqcounts.Rdata")
cancer_type = 'BRCA'
#Group set 1 and 2

Group1and2 <- transform(merge(HiSeq_patients_woTP53mut, HiSeq_patients_withTP53mut, by=0, all=T),
                        row.names = Row.names, Row.names = NULL)
Group1and2.label <- c(rep("TP53_WT",times=dim(HiSeq_patients_woTP53mut)[2]),
                      rep("TP53_mutant",times=dim(HiSeq_patients_withTP53mut)[2]))
Group1and2.label = factor(Group1and2.label,
                          levels = c("TP53_WT","TP53_mutant"))
Group1and2.input <- DGEList(counts=Group1and2, group=Group1and2.label)

Group1and2.input.copy <- Group1and2.input
#filtering for less representative genes
keep <- filterByExpr(Group1and2.input.copy)
Group1and2.input.copy <- Group1and2.input.copy[keep,]
Group1and2.input.copy$samples$lib.size <- colSums(Group1and2.input.copy$counts)

Group1and2.input.copy <- calcNormFactors(Group1and2.input.copy)

#Estimating the Dispersion
#Naive model
d1 <- estimateCommonDisp(Group1and2.input.copy, verbose=T)
#output Disp = 0.32147 , BCV = 0.567
#Bayes tagwise dispersion
d1 <- estimateTagwiseDisp(d1)
plotBCV(d1)
#GLM estimation of dispersion
Group1and2.design.mat <- model.matrix(~ 0 + Group1and2.input.copy$samples$group)
rownames(Group1and2.design.mat) <- colnames(Group1and2.input.copy)
colnames(Group1and2.design.mat) <- levels(Group1and2.input.copy$samples$group)
Group1and2.d2 <- estimateGLMCommonDisp(Group1and2.input.copy,Group1and2.design.mat)
Group1and2.d2 <- estimateGLMTrendedDisp(Group1and2.d2, Group1and2.design.mat)
# You can change method to "auto", "bin.spline", "power", "spline", "bin.loess".
# The default is "auto" which chooses "bin.spline" when > 200 tags and "power" otherwise.
Group1and2.d2 <- estimateGLMTagwiseDisp(Group1and2.d2,Group1and2.design.mat)
plotBCV(Group1and2.d2)

Group1and2.fit <- glmFit(Group1and2.d2, Group1and2.design.mat)
Group1and2.results <- glmLRT(Group1and2.fit, contrast = c(-1, 1)) 
Group1and2.edgeR <- topTags(Group1and2.results, n = nrow(Group1and2.d2))
Group1and2.edgeR <- as.data.frame(Group1and2.edgeR)
colnames(Group1and2.edgeR)[4] <- c("P.Value")

#Group1and3=======================================================================
Group1and3 <- transform(merge(HiSeq_patients_woTP53mut, HiSeq_patients_withMultipleTP53mut, by=0, all=T),
                        row.names = Row.names, Row.names = NULL)
Group1and3.label <- c(rep("TP53_WT",times=dim(HiSeq_patients_woTP53mut)[2]),
                      rep("TP53_multiple_mutant",times=dim(HiSeq_patients_withMultipleTP53mut)[2]))
Group1and3.label = factor(Group1and3.label,
                          levels = c("TP53_WT","TP53_multiple_mutant"))
Group1and3.input <- DGEList(counts=Group1and3, group=Group1and3.label)

Group1and3.input.copy <- Group1and3.input
#filtering for less representative genes
keep <- filterByExpr(Group1and3.input.copy)
Group1and3.input.copy <- Group1and3.input.copy[keep,]
Group1and3.input.copy$samples$lib.size <- colSums(Group1and3.input.copy$counts)

Group1and3.input.copy <- calcNormFactors(Group1and3.input.copy)

#Estimating the Dispersion
#GLM estimation of dispersion
Group1and3.design.mat <- model.matrix(~ 0 + Group1and3.input.copy$samples$group)
rownames(Group1and3.design.mat) <- colnames(Group1and3.input.copy)
colnames(Group1and3.design.mat) <- levels(Group1and3.input.copy$samples$group)
Group1and3.d2 <- estimateGLMCommonDisp(Group1and3.input.copy,Group1and3.design.mat)
Group1and3.d2 <- estimateGLMTrendedDisp(Group1and3.d2, Group1and3.design.mat)
# You can change method to "auto", "bin.spline", "power", "spline", "bin.loess".
# The default is "auto" which chooses "bin.spline" when > 200 tags and "power" otherwise.
Group1and3.d2 <- estimateGLMTagwiseDisp(Group1and3.d2,Group1and3.design.mat)
plotBCV(Group1and3.d2)

Group1and3.fit <- glmFit(Group1and3.d2, Group1and3.design.mat)
Group1and3.results <- glmLRT(Group1and3.fit, contrast = c(-1, 1)) 
Group1and3.edgeR <- topTags(Group1and3.results, n = nrow(Group1and3.d2))
Group1and3.edgeR <- as.data.frame(Group1and3.edgeR)
colnames(Group1and3.edgeR)[4] <- c("P.Value")

#Group1and4============================================================
Group1and4 <- transform(merge(HiSeq_patients_woTP53mut, HiSeq_patients_R175, by=0, all=T),
                        row.names = Row.names, Row.names = NULL)
Group1and4.label <- c(rep("TP53_WT",times=dim(HiSeq_patients_woTP53mut)[2]),
                      rep("TP53_R175",times=dim(HiSeq_patients_R175)[2]))
Group1and4.label = factor(Group1and4.label,
                          levels = c("TP53_WT","TP53_R175"))
Group1and4.input <- DGEList(counts=Group1and4, group=Group1and4.label)

Group1and4.input.copy <- Group1and4.input
#filtering for less representative genes
keep <- filterByExpr(Group1and4.input.copy)
Group1and4.input.copy <- Group1and4.input.copy[keep,]
Group1and4.input.copy$samples$lib.size <- colSums(Group1and4.input.copy$counts)

Group1and4.input.copy <- calcNormFactors(Group1and4.input.copy)

#Estimating the Dispersion
#GLM estimation of dispersion
Group1and4.design.mat <- model.matrix(~ 0 + Group1and4.input.copy$samples$group)
rownames(Group1and4.design.mat) <- colnames(Group1and4.input.copy)
colnames(Group1and4.design.mat) <- levels(Group1and4.input.copy$samples$group)
Group1and4.d2 <- estimateGLMCommonDisp(Group1and4.input.copy,Group1and4.design.mat)
Group1and4.d2 <- estimateGLMTrendedDisp(Group1and4.d2, Group1and4.design.mat)
# You can change method to "auto", "bin.spline", "power", "spline", "bin.loess".
# The default is "auto" which chooses "bin.spline" when > 200 tags and "power" otherwise.
Group1and4.d2 <- estimateGLMTagwiseDisp(Group1and4.d2,Group1and4.design.mat)
plotBCV(Group1and4.d2)

Group1and4.fit <- glmFit(Group1and4.d2, Group1and4.design.mat)
Group1and4.results <- glmLRT(Group1and4.fit, contrast = c(-1, 1)) 
Group1and4.edgeR <- topTags(Group1and4.results, n = nrow(Group1and4.d2))
Group1and4.edgeR <- as.data.frame(Group1and4.edgeR)
colnames(Group1and4.edgeR)[4] <- c("P.Value")

#Group1and5=====================================================================================
Group1and5 <- transform(merge(HiSeq_patients_woTP53mut, HiSeq_patients_R248, by=0, all=T),
                        row.names = Row.names, Row.names = NULL)
Group1and5.label <- c(rep("TP53_WT",times=dim(HiSeq_patients_woTP53mut)[2]),
                      rep("TP53_R248",times=dim(HiSeq_patients_R248)[2]))
Group1and5.label = factor(Group1and5.label,
                          levels = c("TP53_WT","TP53_R248"))
Group1and5.input <- DGEList(counts=Group1and5, group=Group1and5.label)

Group1and5.input.copy <- Group1and5.input
#filtering for less representative genes
keep <- filterByExpr(Group1and5.input.copy)
Group1and5.input.copy <- Group1and5.input.copy[keep,]
Group1and5.input.copy$samples$lib.size <- colSums(Group1and5.input.copy$counts)

Group1and5.input.copy <- calcNormFactors(Group1and5.input.copy)

#Estimating the Dispersion
#GLM estimation of dispersion
Group1and5.design.mat <- model.matrix(~ 0 + Group1and5.input.copy$samples$group)
rownames(Group1and5.design.mat) <- colnames(Group1and5.input.copy)
colnames(Group1and5.design.mat) <- levels(Group1and5.input.copy$samples$group)
Group1and5.d2 <- estimateGLMCommonDisp(Group1and5.input.copy,Group1and5.design.mat)
Group1and5.d2 <- estimateGLMTrendedDisp(Group1and5.d2, Group1and5.design.mat)
# You can change method to "auto", "bin.spline", "power", "spline", "bin.loess".
# The default is "auto" which chooses "bin.spline" when > 200 tags and "power" otherwise.
Group1and5.d2 <- estimateGLMTagwiseDisp(Group1and5.d2,Group1and5.design.mat)
plotBCV(Group1and5.d2)

Group1and5.fit <- glmFit(Group1and5.d2, Group1and5.design.mat)
Group1and5.results <- glmLRT(Group1and5.fit, contrast = c(-1, 1)) 
Group1and5.edgeR <- topTags(Group1and5.results, n = nrow(Group1and5.d2))
Group1and5.edgeR <- as.data.frame(Group1and5.edgeR)
colnames(Group1and5.edgeR)[4] <- c("P.Value")

#Group1and6================================================
Group1and6 <- transform(merge(HiSeq_patients_woTP53mut, HiSeq_patients_R273, by=0, all=T),
                        row.names = Row.names, Row.names = NULL)
Group1and6.label <- c(rep("TP53_WT",times=dim(HiSeq_patients_woTP53mut)[2]),
                      rep("TP53_R273",times=dim(HiSeq_patients_R273)[2]))
Group1and6.label = factor(Group1and6.label,
                          levels = c("TP53_WT","TP53_R273"))
Group1and6.input <- DGEList(counts=Group1and6, group=Group1and6.label)

Group1and6.input.copy <- Group1and6.input
#filtering for less representative genes
keep <- filterByExpr(Group1and6.input.copy)
Group1and6.input.copy <- Group1and6.input.copy[keep,]
Group1and6.input.copy$samples$lib.size <- colSums(Group1and6.input.copy$counts)

Group1and6.input.copy <- calcNormFactors(Group1and6.input.copy)

#Estimating the Dispersion
#GLM estimation of dispersion
Group1and6.design.mat <- model.matrix(~ 0 + Group1and6.input.copy$samples$group)
rownames(Group1and6.design.mat) <- colnames(Group1and6.input.copy)
colnames(Group1and6.design.mat) <- levels(Group1and6.input.copy$samples$group)
Group1and6.d2 <- estimateGLMCommonDisp(Group1and6.input.copy,Group1and6.design.mat)
Group1and6.d2 <- estimateGLMTrendedDisp(Group1and6.d2, Group1and6.design.mat)
# You can change method to "auto", "bin.spline", "power", "spline", "bin.loess".
# The default is "auto" which chooses "bin.spline" when > 200 tags and "power" otherwise.
Group1and6.d2 <- estimateGLMTagwiseDisp(Group1and6.d2,Group1and6.design.mat)
plotBCV(Group1and6.d2)

Group1and6.fit <- glmFit(Group1and6.d2, Group1and6.design.mat)
Group1and6.results <- glmLRT(Group1and6.fit, contrast = c(-1, 1)) 
Group1and6.edgeR <- topTags(Group1and6.results, n = nrow(Group1and6.d2))
Group1and6.edgeR <- as.data.frame(Group1and6.edgeR)
colnames(Group1and6.edgeR)[4] <- c("P.Value")

#Group1and7==================================================================
Group1and7 <- transform(merge(HiSeq_patients_woTP53mut, HiSeq_patients_R, by=0, all=T),
                        row.names = Row.names, Row.names = NULL)
Group1and7.label <- c(rep("TP53_WT",times=dim(HiSeq_patients_woTP53mut)[2]),
                      rep("TP53_Group3",times=dim(HiSeq_patients_R)[2]))
Group1and7.label = factor(Group1and7.label,
                          levels = c("TP53_WT","TP53_Group3"))
Group1and7.input <- DGEList(counts=Group1and7, group=Group1and7.label)

Group1and7.input.copy <- Group1and7.input
#filtering for less representative genes
keep <- filterByExpr(Group1and7.input.copy)
Group1and7.input.copy <- Group1and7.input.copy[keep,]
Group1and7.input.copy$samples$lib.size <- colSums(Group1and7.input.copy$counts)

Group1and7.input.copy <- calcNormFactors(Group1and7.input.copy)

#Estimating the Dispersion
#GLM estimation of dispersion
Group1and7.design.mat <- model.matrix(~ 0 + Group1and7.input.copy$samples$group)
rownames(Group1and7.design.mat) <- colnames(Group1and7.input.copy)
colnames(Group1and7.design.mat) <- levels(Group1and7.input.copy$samples$group)
Group1and7.d2 <- estimateGLMCommonDisp(Group1and7.input.copy,Group1and7.design.mat)
Group1and7.d2 <- estimateGLMTrendedDisp(Group1and7.d2, Group1and7.design.mat)
# You can change method to "auto", "bin.spline", "power", "spline", "bin.loess".
# The default is "auto" which chooses "bin.spline" when > 200 tags and "power" otherwise.
Group1and7.d2 <- estimateGLMTagwiseDisp(Group1and7.d2,Group1and7.design.mat)
#plotBCV(Group1and7.d2)

Group1and7.fit <- glmFit(Group1and7.d2, Group1and7.design.mat)
Group1and7.results <- glmLRT(Group1and7.fit, contrast = c(-1, 1)) 
Group1and7.edgeR <- topTags(Group1and7.results, n = nrow(Group1and7.d2))
Group1and7.edgeR <- as.data.frame(Group1and7.edgeR)
colnames(Group1and7.edgeR)[4] <- c("P.Value")

#data.dir <- file.path('./TCGA_P53')
source('00_drawing_function.R')
Group1and2.top100 <- top100_heatmap(V.method = Group1and2.edgeR, gene.dataset = Group1and2, 
                                        annotation.class = factor(Group1and2.label))
Group1and3.top100 <- top100_heatmap(V.method = Group1and3.edgeR, gene.dataset = Group1and3, 
                                    annotation.class = factor(Group1and3.label))
Group1and4.top100 <- top100_heatmap(V.method = Group1and4.edgeR, gene.dataset = Group1and4, 
                                    annotation.class = factor(Group1and4.label))
Group1and5.top100 <- top100_heatmap(V.method = Group1and5.edgeR, gene.dataset = Group1and5, 
                                    annotation.class = factor(Group1and5.label))
Group1and6.top100 <- top100_heatmap(V.method = Group1and6.edgeR, gene.dataset = Group1and6, 
                                    annotation.class = factor(Group1and6.label))
Group1and7.top100 <- top100_heatmap(V.method = Group1and7.edgeR, gene.dataset = Group1and7, 
                                    annotation.class = factor(Group1and7.label))

col_fun = colorRamp2(c(-2, 0, 2), c("#2fa1dd", "white", "#f87669"))
lgd = Legend(col_fun = col_fun, title = "scaled log2Expression", direction = 'horizontal')
ht_list_top100 <- Group1and2.top100  + Group1and3.top100  
ht_list_top1002 <- Group1and4.top100   + Group1and5.top100  + Group1and6.top100 
ht_list_top1003 <- Group1and7.top100

draw(ht_list_top100, auto_adjust = FALSE, 
     column_title = "BRCA WT TP53 vs Different Groups", column_title_gp = gpar(fontsize = 16))
draw(ht_list_top1002, auto_adjust = FALSE, 
     column_title = "BRCA WT TP53 vs Different Groups", column_title_gp = gpar(fontsize = 16))
draw(ht_list_top1003, auto_adjust = FALSE, 
     column_title = "BRCA WT TP53 vs Different Groups", column_title_gp = gpar(fontsize = 16))

draw(lgd, x = unit(1, "npc"), y = unit(1, "npc"), just = c("right", "top"))

#Handle data Group1and2 ================================================================
Group1and2.edgeR.filter <- filter(Group1and2.edgeR, abs(Group1and2.edgeR$logFC) >=1
                       & Group1and2.edgeR$P.Value < 0.05)
print(length(rownames(Group1and2.edgeR.filter)))
k1 = (Group1and2.edgeR.filter$P.Value < 0.05)&(Group1and2.edgeR.filter$logFC < -1)
k2 = (Group1and2.edgeR.filter$P.Value < 0.05)&(Group1and2.edgeR.filter$logFC > 1)
Group1and2.edgeR.filter$change = ifelse(k1,"DOWN", 'UP')

Group1and2.choose_gene = rownames(Group1and2.edgeR.filter)
Group1and2.matrix = Group1and2[Group1and2.choose_gene,]
Group1and2.matrix = log2(Group1and2.matrix + 1)
Group1and2.matrix = t(scale(t(Group1and2.matrix)))
Group1and2.matrix[Group1and2.matrix > 2] = 2
Group1and2.matrix[Group1and2.matrix < -2] = -2

Group1and2.all.heatmap <- draw_heatmap(choose_matrix = Group1and2.matrix, 
                                         annotation.class = factor(Group1and2.label))
draw(Group1and2.all.heatmap, auto_adjust = FALSE, 
     column_title = "BRCA WT TP53 vs BRCA TP53 Mutant", column_title_gp = gpar(fontsize = 16))
draw(lgd, x = unit(1, "npc"), y = unit(1, "npc"), just = c("right", "top"))


#Handle data Group1and7===============================================================
Group1and7.edgeR.filter <- filter(Group1and7.edgeR, abs(Group1and7.edgeR$logFC) >=1
                                  & Group1and7.edgeR$P.Value < 0.05)
print(length(rownames(Group1and7.edgeR.filter)))
k1 = (Group1and7.edgeR.filter$P.Value < 0.05)&(Group1and7.edgeR.filter$logFC < -1)
k2 = (Group1and7.edgeR.filter$P.Value < 0.05)&(Group1and7.edgeR.filter$logFC > 1)
Group1and7.edgeR.filter$change = ifelse(k1,"DOWN", 'UP')

Group1and7.choose_gene = rownames(Group1and7.edgeR.filter)
Group1and7.matrix = Group1and7[Group1and7.choose_gene,]
Group1and7.matrix = log2(Group1and7.matrix + 1)
Group1and7.matrix = t(scale(t(Group1and7.matrix)))
Group1and7.matrix[Group1and7.matrix > 2] = 2
Group1and7.matrix[Group1and7.matrix < -2] = -2

Group1and7.all.heatmap <- draw_heatmap(choose_matrix = Group1and7.matrix, 
                                       annotation.class = factor(Group1and7.label))
draw(Group1and7.all.heatmap, auto_adjust = FALSE, 
     column_title = "BRCA WT TP53 vs BRCA TP53 Rmutant", column_title_gp = gpar(fontsize = 16))
draw(lgd, x = unit(1, "npc"), y = unit(1, "npc"), just = c("right", "top"))

#Handle data Group1and3===============================================================
Group1and3.edgeR.filter <- filter(Group1and3.edgeR, abs(Group1and3.edgeR$logFC) >=1
                                  & Group1and3.edgeR$P.Value < 0.05)
print(length(rownames(Group1and3.edgeR.filter)))
k1 = (Group1and3.edgeR.filter$P.Value < 0.05)&(Group1and3.edgeR.filter$logFC < -1)
k2 = (Group1and3.edgeR.filter$P.Value < 0.05)&(Group1and3.edgeR.filter$logFC > 1)
Group1and3.edgeR.filter$change = ifelse(k1,"DOWN", 'UP')

Group1and3.choose_gene = rownames(Group1and3.edgeR.filter)
Group1and3.matrix = Group1and3[Group1and3.choose_gene,]
Group1and3.matrix = log2(Group1and3.matrix + 1)
Group1and3.matrix = t(scale(t(Group1and3.matrix)))
Group1and3.matrix[Group1and3.matrix > 2] = 2
Group1and3.matrix[Group1and3.matrix < -2] = -2

Group1and3.all.heatmap <- draw_heatmap(choose_matrix = Group1and3.matrix, 
                                       annotation.class = factor(Group1and3.label))
draw(Group1and3.all.heatmap, auto_adjust = FALSE, 
     column_title = "BRCA WT TP53 vs BRCA TP53 multiple mutations", column_title_gp = gpar(fontsize = 16))
draw(lgd, x = unit(1, "npc"), y = unit(1, "npc"), just = c("right", "top"))

#Handle data Group1and4===============================================================
Group1and4.edgeR.filter <- filter(Group1and4.edgeR, abs(Group1and4.edgeR$logFC) >=1
                                  & Group1and4.edgeR$P.Value < 0.05)
print(length(rownames(Group1and4.edgeR.filter)))
k1 = (Group1and4.edgeR.filter$P.Value < 0.05)&(Group1and4.edgeR.filter$logFC < -1)
k2 = (Group1and4.edgeR.filter$P.Value < 0.05)&(Group1and4.edgeR.filter$logFC > 1)
Group1and4.edgeR.filter$change = ifelse(k1,"DOWN", 'UP')

Group1and4.choose_gene = rownames(Group1and4.edgeR.filter)
Group1and4.matrix = Group1and4[Group1and4.choose_gene,]
Group1and4.matrix = log2(Group1and4.matrix + 1)
Group1and4.matrix = t(scale(t(Group1and4.matrix)))
Group1and4.matrix[Group1and4.matrix > 2] = 2
Group1and4.matrix[Group1and4.matrix < -2] = -2

#Handle data Group1and5===============================================================
Group1and5.edgeR.filter <- filter(Group1and5.edgeR, abs(Group1and5.edgeR$logFC) >=1
                                  & Group1and5.edgeR$P.Value < 0.05)
print(length(rownames(Group1and5.edgeR.filter)))
k1 = (Group1and5.edgeR.filter$P.Value < 0.05)&(Group1and5.edgeR.filter$logFC < -1)
k2 = (Group1and5.edgeR.filter$P.Value < 0.05)&(Group1and5.edgeR.filter$logFC > 1)
Group1and5.edgeR.filter$change = ifelse(k1,"DOWN", 'UP')

Group1and5.choose_gene = rownames(Group1and5.edgeR.filter)
Group1and5.matrix = Group1and5[Group1and5.choose_gene,]
Group1and5.matrix = log2(Group1and5.matrix + 1)
Group1and5.matrix = t(scale(t(Group1and5.matrix)))
Group1and5.matrix[Group1and5.matrix > 2] = 2
Group1and5.matrix[Group1and5.matrix < -2] = -2

#Handle data Group1and6===============================================================
Group1and6.edgeR.filter <- filter(Group1and6.edgeR, abs(Group1and6.edgeR$logFC) >=1
                                  & Group1and6.edgeR$P.Value < 0.05)
print(length(rownames(Group1and6.edgeR.filter)))
k1 = (Group1and6.edgeR.filter$P.Value < 0.05)&(Group1and6.edgeR.filter$logFC < -1)
k2 = (Group1and6.edgeR.filter$P.Value < 0.05)&(Group1and6.edgeR.filter$logFC > 1)
Group1and6.edgeR.filter$change = ifelse(k1,"DOWN", 'UP')

Group1and6.choose_gene = rownames(Group1and6.edgeR.filter)
Group1and6.matrix = Group1and6[Group1and6.choose_gene,]
Group1and6.matrix = log2(Group1and6.matrix + 1)
Group1and6.matrix = t(scale(t(Group1and6.matrix)))
Group1and6.matrix[Group1and6.matrix > 2] = 2
Group1and6.matrix[Group1and6.matrix < -2] = -2


save(Group1and2.edgeR.filter,Group1and7.edgeR.filter,Group1and3.edgeR.filter,
     Group1and4.edgeR.filter,Group1and5.edgeR.filter,Group1and6.edgeR.filter,
     Group1and2.matrix,Group1and7.matrix,Group1and3.matrix,
     Group1and4.matrix,Group1and5.matrix,Group1and6.matrix,
     file = paste0('./data/','step2','.Rdata'))

#check if gene of interest is within top100
gene.of.interest <- c("BRCA1", "SRSF2", "ATF4", "SLC7A11", 
                      "ZFAF1", "WT1", "PBRM1", "BAP1",
                      "CDH1", "CIC", "RB1", "STK11", "SMAD4", 
                      "NF1", "TP53",'SFRS2')
look.for.gene.of.interest <- filter(HiSeq_for_primary_tumor, rownames(Group1and2.choose_matrix) %in% gene.of.interest)
if (look.for.gene.of.interest) {
  print('Gene of interest is within top100 regulated gene.')
} else {
  print('Gene of interest is within top100 regulated gene.
        Check if G.O.I is within other high expression gene')
  another.check <- Group1and2.edgeR[rownames(Group1and2.edgeR) %in% gene.of.interest]
  if (another.check) {
    print(another.check)
  } else {
    print('G.O.I is not within high epression gene. They are filtered before performing edgeR.')
  }
}
  



###################

edgeR.filter <- filter(Group1and2.edgeR.filter, abs(Group1and2.edgeR.filter$logFC) >=1
                       & Group1and2.edgeR.filter$P.Value < 0.05)
nrDEG_Z = edgeR.filter[order(edgeR.filter$logFC), ]
nrDEG_F = edgeR.filter[order(-edgeR.filter$logFC), ]
choose_gene_top100 = c(rownames(nrDEG_Z)[1:50], rownames(nrDEG_F)[1:50])
choose_matrix_top100 = Group1and2[choose_gene_top100,]
choose_matrix_top100 = log2(choose_matrix_top100 + 1)
choose_matrix_top100 = t(scale(t(choose_matrix_top100)))
choose_matrix_top100[choose_matrix_top100 > 2] = 2
choose_matrix_top100[choose_matrix_top100 < -2] = -2

choose_gene = rownames(edgeR.filter)
choose_matrix = Group1and2[choose_gene,]
choose_matrix = log2(choose_matrix + 1)
choose_matrix = t(scale(t(choose_matrix)))

choose_matrix[choose_matrix > 2] = 2
choose_matrix[choose_matrix < -2] = -2

col_fun = colorRamp2(c(-2, 0, 2), c("#2fa1dd", "white", "#f87669"))
top_annotation = HeatmapAnnotation(
  cluster = anno_block(gp = gpar(fill = c("#2fa1dd", "#f87669")),
                       labels = unique(factor(Group1and2.label)),
                       labels_gp = gpar(col = "white", fontsize = 12)))

fig_top100 <-  Heatmap(choose_matrix_top100, col = col_fun, top_annotation = top_annotation,
                       show_heatmap_legend = T, column_split = factor(Group1and2.label), 
                       show_row_names = T, show_column_names = F,column_title = 'top100', 
                       cluster_columns = F, cluster_column_slices = T, row_names_gp = gpar(fontsize = 6))
fig_all <-  Heatmap(choose_matrix, col = col_fun, top_annotation = top_annotation,
                    show_heatmap_legend = F, column_split = factor(Group1and2.label),
                    show_row_names = T, show_column_names = F,column_title = 'All_significant',
                    cluster_columns = F, cluster_column_slices = T, row_names_gp = gpar(fontsize = 6))
ht_list <- fig_top100 + fig_all
draw(ht_list, auto_adjust = FALSE,
     column_title = "Group WT TP53 VS Group TP53 Mutant", column_title_gp = gpar(fontsize = 16))


Group1and2.edgeR.filter <- filter(Group1and2.edgeR, abs(Group1and2.edgeR$logFC) >=1
                                     & Group1and2.edgeR$P.Value < 0.05)
nrDEG_Z = Group1and2.edgeR.filter[order(Group1and2.edgeR.filter$logFC), ]
nrDEG_F = Group1and2.edgeR.filter[order(-Group1and2.edgeR.filter$logFC), ]
choose_gene = c(rownames(nrDEG_Z)[1:50], rownames(nrDEG_F)[1:50])
choose_matrix = Group1and2[choose_gene,]
choose_matrix = log2(choose_matrix + 1)
choose_matrix = t(scale(t(choose_matrix)))
choose_matrix[choose_matrix > 2] = 2
choose_matrix[choose_matrix < -2] = -2


```
nrDEG_Z = Group1and2.edgeR[order(Group1and2.edgeR$logFC), ]
nrDEG_F = Group1and2.edgeR[order(-Group1and2.edgeR$logFC), ]
choose_gene = c(rownames(nrDEG_Z)[1:50], rownames(nrDEG_F)[1:50])
choose_matrix = Group1and2[choose_gene,]
choose_matrix = log2(choose_matrix + 1)
choose_matrix = t(scale(t(choose_matrix)))

choose_matrix[choose_matrix > 2] = 2
choose_matrix[choose_matrix < -2] = -2
```

#annotation_col = data.frame(CellType = factor(Group1and2.label))
#rownames(annotation_col) = colnames(Group1and2)

col_fun = colorRamp2(c(-2, 0, 2), c("#2fa1dd", "white", "#f87669"))
top_annotation = HeatmapAnnotation(
  cluster = anno_block(gp = gpar(fill = c("#2fa1dd", "#f87669")),
                       labels = c('without.TP53','with.TP53'),
                       labels_gp = gpar(col = "white", fontsize = 12)))

a<-Heatmap(choose_matrix, col = col_fun, top_annotation = top_annotation,
        show_heatmap_legend = T, column_split = factor(Group1and2.label),row_title = 'group2',
        show_row_names = T, show_column_names = F,column_title = 'top100', 
        cluster_columns = F, cluster_column_slices = T, row_names_gp = gpar(fontsize = 6))
b<-Heatmap(choose_matrix, col = col_fun, top_annotation = top_annotation,
          show_heatmap_legend = F, column_split = factor(Group1and2.label), row_title = 'group1',
          show_row_names = T, show_column_names = F,column_title = 'top1000', 
          cluster_columns = F, cluster_column_slices = T, row_names_gp = gpar(fontsize = 6))

ht_list = a+b
draw(ht_list, row_title = "Three heatmaps, row title", row_title_gp = gpar(col = "red"),
     column_title = "Three heatmaps, column title", column_title_gp = gpar(fontsize = 16))

plot_grid(a,b,ncol = 1, nrow = 2)
somePDFPath = paste0('./figs/test.pdf')
pdf(file=somePDFPath)
par(mfrow=(c(1,3)))
plot(a,b)
dev.off()
for (i in seq(5,10))   
{   
  par(mfrow = c(2,1))
  VAR1=rnorm(i)  
  VAR2=rnorm(i)  
  plot(VAR1,VAR2)   
} 
dev.off() 
savePlot(a, filename = paste0('./figs/', name='BRCA', '/_heatmap_top100_logFC.png'))
#################
gene.of.interest <- c("BRCA1", "SRSF2", "ATF4", "SLC7A11",
                      "ZFAF1", "WT1", "PBRM1", "BAP1",
                      "CDH1", "CIC", "RB1", "STK11", "SMAD4", 
                      "NF1", "TP53",'SFRS2')
HiSeq_selected_gene <- data.frame(GeneOfInterest = )
  filter(HiSeq_for_primary_tumor, rownames(choose_matrix) %in% gene.of.interest)



draw_heatmap(nrDEG = nrDEG_edgeR, type = 'edgeR')
draw_volcano(nrDEG = nrDEG_edgeR, type = 'edgeR')

