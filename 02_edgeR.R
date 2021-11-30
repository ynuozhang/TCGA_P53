library(edgeR)
library(dplyr)

# FPKM trace back
load("D:/TCGA_P53/data/BRCA_HiSeqcounts_w_TP53_GOI.Rdata")
load("D:/TCGA_P53/data/BRCA_HiSeqcounts_wo_TP53_GOI.Rdata")

rownames(HiSeq_selected_patients_withTP53mut) <- HiSeq_selected_patients_withTP53mut[,1]
HiSeq_selected_patients_withTP53mut <- HiSeq_selected_patients_withTP53mut[,-1]
rownames(HiSeq_selected_patients_woTP53mut) <- HiSeq_selected_patients_woTP53mut[,1]
HiSeq_selected_patients_woTP53mut <- HiSeq_selected_patients_woTP53mut[,-1]

patients_withTP53mut_count <- round(((2**HiSeq_selected_patients_withTP53mut) - 1), 0)
patients_woTP53mut_count <- round(((2**HiSeq_selected_patients_woTP53mut) - 1), 0)

#Group set 1 and 2
#combine groups: patients with TP53 and patients wo TP53
Group1and2 <- transform(merge(patients_woTP53mut_count, patients_withTP53mut_count, by=0, all=T),
                        row.names = Row.names, Row.names = NULL)
Group1and2.label <- c(rep("without.TP53",times=dim(patients_woTP53mut_count)[2]),
                      rep("with.TP53",times=dim(patients_withTP53mut_count)[2]))
Group1and2.label = factor(Group1and2.label,
                          levels = c("without.TP53","with.TP53"))
Group1and2.input <- DGEList(counts=Group1and2, group=Group1and2.label)

#skip filtering for less representative genes
#start normalization by total count
#"The calcNormFactors() function normalizes for RNA composition 
#by finding a set of scaling factors for the library sizes that 
#minimize the log-fold changes between the samples for most genes." 
Group1and2.input <- calcNormFactors(Group1and2.input)

#Estimating the Dispersion
