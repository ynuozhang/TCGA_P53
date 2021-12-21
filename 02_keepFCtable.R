############select for only top 50
load("D:/TCGA_P53/data/step2.Rdata")
load("D:/TCGA_P53/data/step2-top100.Rdata")
library(dplyr)
library(circlize)
library(ggplot2)
library(ComplexHeatmap)

Group1and2 <- transform(merge(HiSeq_patients_woTP53mut, HiSeq_patients_withTP53mut, by=0, all=T),
                        row.names = Row.names, Row.names = NULL)
Group1and3 <- transform(merge(HiSeq_patients_woTP53mut, HiSeq_patients_withMultipleTP53mut, by=0, all=T),
                        row.names = Row.names, Row.names = NULL)
Group1and4 <- transform(merge(HiSeq_patients_woTP53mut, HiSeq_patients_R175, by=0, all=T),
                        row.names = Row.names, Row.names = NULL)
Group1and5 <- transform(merge(HiSeq_patients_woTP53mut, HiSeq_patients_R248, by=0, all=T),
                        row.names = Row.names, Row.names = NULL)
Group1and6 <- transform(merge(HiSeq_patients_woTP53mut, HiSeq_patients_R273, by=0, all=T),
                        row.names = Row.names, Row.names = NULL)
Group1and7 <- transform(merge(HiSeq_patients_woTP53mut, HiSeq_patients_R, by=0, all=T),
                        row.names = Row.names, Row.names = NULL)
save(Group1and2,Group1and3,Group1and4, Group1and5, Group1and6, Group1and7, file = 'Step2_gene_exp_count.Rdata')

choose_gene_1and2 = rownames(Group1and2.edgeR.filter)
choose_gene_1and3 = rownames(Group1and3.edgeR.filter)
choose_gene_1and4 = rownames(Group1and4.edgeR.filter)
choose_gene_1and5 = rownames(Group1and5.edgeR.filter)
choose_gene_1and6 = rownames(Group1and6.edgeR.filter)
choose_gene_1and7 = rownames(Group1and7.edgeR.filter)

gene.of.interest <- c("BRCA1", "SRSF2", "ATF4", "SLC7A11", 
                      "ZFAF1", "WT1", "PBRM1", "BAP1",
                      "CDH1", "CIC", "RB1", "STK11", "SMAD4", 
                      "NF1", "TP53",'SFRS2')
Group1and2.edgeR.filter[rownames(Group1and2.edgeR.filter) %in% gene.of.interest,]
Group1and3.edgeR.filter[rownames(Group1and3.edgeR.filter) %in% gene.of.interest,]
Group1and4.edgeR.filter[rownames(Group1and4.edgeR.filter) %in% gene.of.interest,]
Group1and5.edgeR.filter[rownames(Group1and5.edgeR.filter) %in% gene.of.interest,]
Group1and6.edgeR.filter[rownames(Group1and6.edgeR.filter) %in% gene.of.interest,]
Group1and7.edgeR.filter[rownames(Group1and7.edgeR.filter) %in% gene.of.interest,]

keep.FC.1and2 <- Group1and2.edgeR[rownames(Group1and2.edgeR) %in% gene.of.interest,][,1, drop=FALSE]
keep.FC.1and3 <- Group1and3.edgeR[rownames(Group1and3.edgeR) %in% gene.of.interest,][,1, drop=FALSE]
keep.FC.1and4 <- Group1and4.edgeR[rownames(Group1and4.edgeR) %in% gene.of.interest,][,1, drop=FALSE]
keep.FC.1and5 <- Group1and5.edgeR[rownames(Group1and5.edgeR) %in% gene.of.interest,][,1, drop=FALSE]
keep.FC.1and6 <- Group1and6.edgeR[rownames(Group1and6.edgeR) %in% gene.of.interest,][,1, drop=FALSE]
keep.FC.1and7 <- Group1and7.edgeR[rownames(Group1and7.edgeR) %in% gene.of.interest,][,1, drop=FALSE]

keep.FC <- Reduce(function(a,b){
  ans <- merge(a,b,by="row.names",all=T)
  row.names(ans) <- ans[,"Row.names"]
  ans[,!names(ans) %in% "Row.names"]
}, list(keep.FC.1and2,keep.FC.1and3,keep.FC.1and4,
        keep.FC.1and5,keep.FC.1and6,keep.FC.1and7))
#write.csv(keep.FC,"KEEPED-FC.csv", row.names = TRUE)

#################only select top 50
Group1and2.label <- c(rep("TP53_WT",times=dim(HiSeq_patients_woTP53mut)[2]),
                      rep("TP53_mutant",times=dim(HiSeq_patients_withTP53mut)[2]))
Group1and2.label = factor(Group1and2.label,
                          levels = c("TP53_WT","TP53_mutant"))
Group1and3.label <- c(rep("TP53_WT",times=dim(HiSeq_patients_woTP53mut)[2]),
                      rep("TP53_multiple_mutant",times=dim(HiSeq_patients_withMultipleTP53mut)[2]))
Group1and3.label = factor(Group1and3.label,
                          levels = c("TP53_WT","TP53_multiple_mutant"))
Group1and4.label <- c(rep("TP53_WT",times=dim(HiSeq_patients_woTP53mut)[2]),
                      rep("TP53_R175",times=dim(HiSeq_patients_R175)[2]))
Group1and4.label = factor(Group1and4.label,
                          levels = c("TP53_WT","TP53_R175"))
Group1and5.label <- c(rep("TP53_WT",times=dim(HiSeq_patients_woTP53mut)[2]),
                      rep("TP53_R248",times=dim(HiSeq_patients_R248)[2]))
Group1and5.label = factor(Group1and5.label,
                          levels = c("TP53_WT","TP53_R248"))
Group1and6.label <- c(rep("TP53_WT",times=dim(HiSeq_patients_woTP53mut)[2]),
                      rep("TP53_R273",times=dim(HiSeq_patients_R273)[2]))
Group1and6.label = factor(Group1and6.label,
                          levels = c("TP53_WT","TP53_R273"))
Group1and7.label <- c(rep("TP53_WT",times=dim(HiSeq_patients_woTP53mut)[2]),
                      rep("TP53_Group3",times=dim(HiSeq_patients_R)[2]))
Group1and7.label = factor(Group1and7.label,
                          levels = c("TP53_WT","TP53_Group3"))
save(Group1and2.label,Group1and3.label,Group1and4.label, 
     Group1and5.label, Group1and6.label, Group1and7.label, file = 'Step2_group_label.Rdata')

nrDEG_Z = Group1and2.edgeR.filter[order(Group1and2.edgeR.filter$logFC), ]
nrDEG_F = Group1and2.edgeR.filter[order(-Group1and2.edgeR.filter$logFC), ]
Group1and2.choose_gene_top50 = c(rownames(nrDEG_Z)[1:25], rownames(nrDEG_F)[1:25])

nrDEG_Z = Group1and3.edgeR.filter[order(Group1and3.edgeR.filter$logFC), ]
nrDEG_F = Group1and3.edgeR.filter[order(-Group1and3.edgeR.filter$logFC), ]
Group1and3.choose_gene_top50 = c(rownames(nrDEG_Z)[1:25], rownames(nrDEG_F)[1:25])
nrDEG_Z = Group1and4.edgeR.filter[order(Group1and4.edgeR.filter$logFC), ]
nrDEG_F = Group1and4.edgeR.filter[order(-Group1and4.edgeR.filter$logFC), ]
Group1and4.choose_gene_top50 = c(rownames(nrDEG_Z)[1:25], rownames(nrDEG_F)[1:25])
nrDEG_Z = Group1and5.edgeR.filter[order(Group1and5.edgeR.filter$logFC), ]
nrDEG_F = Group1and5.edgeR.filter[order(-Group1and5.edgeR.filter$logFC), ]
Group1and5.choose_gene_top50 = c(rownames(nrDEG_Z)[1:25], rownames(nrDEG_F)[1:25])
nrDEG_Z = Group1and6.edgeR.filter[order(Group1and6.edgeR.filter$logFC), ]
nrDEG_F = Group1and6.edgeR.filter[order(-Group1and6.edgeR.filter$logFC), ]
Group1and6.choose_gene_top50 = c(rownames(nrDEG_Z)[1:25], rownames(nrDEG_F)[1:25])
nrDEG_Z = Group1and7.edgeR.filter[order(Group1and7.edgeR.filter$logFC), ]
nrDEG_F = Group1and7.edgeR.filter[order(-Group1and7.edgeR.filter$logFC), ]
Group1and7.choose_gene_top50 = c(rownames(nrDEG_Z)[1:25], rownames(nrDEG_F)[1:25])

source('00_drawing_function.R')
Group1and2.top50 <- draw_top50_bar(filtered.edgeR = Group1and2.edgeR.filter, top50genes = Group1and2.choose_gene_top50,
                             title = 'BRCA WT TP53 vs TP53 Mutant top50 regulated genes')
Group1and2.top50
Group1and3.top50 <- draw_top50_bar(filtered.edgeR = Group1and3.edgeR.filter, top50genes = Group1and3.choose_gene_top50,
                                   title = 'BRCA WT TP53 vs TP53 Multiple Mutant top50 regulated genes')
Group1and3.top50
Group1and4.top50 <- draw_top50_bar(filtered.edgeR = Group1and4.edgeR.filter, top50genes = Group1and4.choose_gene_top50,
                                   title = 'BRCA WT TP53 vs TP53 R175 top50 regulated genes')
Group1and5.top50 <- draw_top50_bar(filtered.edgeR = Group1and5.edgeR.filter, top50genes = Group1and5.choose_gene_top50,
                                   title = 'BRCA WT TP53 vs TP53 R248 top50 regulated genes')
Group1and6.top50 <- draw_top50_bar(filtered.edgeR = Group1and6.edgeR.filter, top50genes = Group1and6.choose_gene_top50,
                                   title = 'BRCA WT TP53 vs TP53 R273 top50 regulated genes')
Group1and7.top50 <- draw_top50_bar(filtered.edgeR = Group1and7.edgeR.filter, top50genes = Group1and7.choose_gene_top50,
                                   title = 'BRCA WT TP53 vs TP53 R175+248+273 top50 regulated genes')
Group1and2.edgeR.top50 <- Group1and2.edgeR.filter[Group1and2.choose_gene_top50, ,drop=FALSE]
Group1and3.edgeR.top50 <- Group1and3.edgeR.filter[Group1and3.choose_gene_top50, ,drop=FALSE]
Group1and4.edgeR.top50 <- Group1and4.edgeR.filter[Group1and4.choose_gene_top50, ,drop=FALSE]
Group1and5.edgeR.top50 <- Group1and5.edgeR.filter[Group1and5.choose_gene_top50, ,drop=FALSE]
Group1and6.edgeR.top50 <- Group1and6.edgeR.filter[Group1and6.choose_gene_top50, ,drop=FALSE]
Group1and7.edgeR.top50 <- Group1and7.edgeR.filter[Group1and7.choose_gene_top50, ,drop=FALSE]
save(Group1and2.edgeR.top50,Group1and3.edgeR.top50,Group1and4.edgeR.top50, 
     Group1and5.edgeR.top50, Group1and6.edgeR.top50, Group1and7.edgeR.top50, file = 'Step2_top50gene.Rdata')

######################################################################

table.to.draw.bar <- Group1and2.edgeR.filter[choose_gene_top50, c(1, 6), drop=FALSE]

table.to.draw.bar$genes <- rownames(table.to.draw.bar)

g_kegg <- ggplot(table.to.draw.bar,
                 aes(x = reorder(genes, order(logFC, decreasing=F)), 
                     y = logFC, fill = change)) + 
  geom_bar(stat = "identity")+
  scale_fill_manual(values = c("#2fa1dd","#f87669"),guide = 'none')+
  scale_y_continuous( name = "logFC" )+
  coord_flip() + theme_bw() +
  theme(axis.title.y=element_blank())
g_kegg

choose_gene_top50 = c(rownames(nrDEG_Z)[1:25], rownames(nrDEG_F)[1:25])
choose_matrix_top50 = Group1and2[choose_gene_top50,]
#choose_matrix_top50 = log2(choose_matrix_top50 + 1)
#choose_matrix_top50 = as.data.frame(t(scale(t(choose_matrix_top50))))
#col_fun = colorRamp2(c(-2, 0, 2), c("#2fa1dd", "white", "#f87669"))

choose_matrix_top50$log_expression_mean <- rowMeans(choose_matrix_top50)
table.to.draw.bar <- choose_matrix_top50[, ncol(choose_matrix_top50),drop = FALSE]
table.to.draw.bar = log2(table.to.draw.bar + 1)

table.to.draw.bar <- mutate(table.to.draw.bar, group =
                              case_when(rownames(table.to.draw.bar) %in% rownames(nrDEG_Z)[1:25] ~ 1,
                            rownames(table.to.draw.bar) %in% rownames(nrDEG_F)[1:25] ~ -1)) %>%
  mutate(log_expression_mean = ifelse(group == -1, -1*log_expression_mean, log_expression_mean))

table.to.draw.bar$genes <- rownames(table.to.draw.bar)
table.to.draw.bar$log_expression_mean <- as.numeric(as.character(table.to.draw.bar$log_expression_mean))
table.to.draw.bar<-table.to.draw.bar[order(table.to.draw.bar$log_expression_mean),]

g_kegg <- ggplot(table.to.draw.bar,
                 aes(x = reorder(genes, order(log_expression_mean, decreasing=F)), 
                     y = log_expression_mean, fill = group)) + 
  geom_bar(stat = "identity")+
  scale_fill_gradient( low = "#2fa1dd", high = "#f87669", guide = 'none' )+
  scale_y_continuous( name = "log_expression_mean" )+
  coord_flip() + theme_bw() +
  theme(axis.title.y=element_blank())
g_kegg



