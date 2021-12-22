library(dplyr)
library(enrichR)
load("D:/Github/TCGA_P53/data/step2.Rdata")

nrDEG_Z = Group1and2.edgeR.filter[order(Group1and2.edgeR.filter$logFC), ]
nrDEG_F = Group1and2.edgeR.filter[order(-Group1and2.edgeR.filter$logFC), ]
Group1and2.choose_gene_top100 = c(rownames(nrDEG_Z)[1:50], rownames(nrDEG_F)[1:50])

nrDEG_Z = Group1and3.edgeR.filter[order(Group1and3.edgeR.filter$logFC), ]
nrDEG_F = Group1and3.edgeR.filter[order(-Group1and3.edgeR.filter$logFC), ]
Group1and3.choose_gene_top100 = c(rownames(nrDEG_Z)[1:50], rownames(nrDEG_F)[1:50])
nrDEG_Z = Group1and4.edgeR.filter[order(Group1and4.edgeR.filter$logFC), ]
nrDEG_F = Group1and4.edgeR.filter[order(-Group1and4.edgeR.filter$logFC), ]
Group1and4.choose_gene_top100 = c(rownames(nrDEG_Z)[1:50], rownames(nrDEG_F)[1:50])
nrDEG_Z = Group1and5.edgeR.filter[order(Group1and5.edgeR.filter$logFC), ]
nrDEG_F = Group1and5.edgeR.filter[order(-Group1and5.edgeR.filter$logFC), ]
Group1and5.choose_gene_top100 = c(rownames(nrDEG_Z)[1:50], rownames(nrDEG_F)[1:50])
nrDEG_Z = Group1and6.edgeR.filter[order(Group1and6.edgeR.filter$logFC), ]
nrDEG_F = Group1and6.edgeR.filter[order(-Group1and6.edgeR.filter$logFC), ]
Group1and6.choose_gene_top100 = c(rownames(nrDEG_Z)[1:50], rownames(nrDEG_F)[1:50])
nrDEG_Z = Group1and7.edgeR.filter[order(Group1and7.edgeR.filter$logFC), ]
nrDEG_F = Group1and7.edgeR.filter[order(-Group1and7.edgeR.filter$logFC), ]
Group1and7.choose_gene_top100 = c(rownames(nrDEG_Z)[1:50], rownames(nrDEG_F)[1:50])

Group1and2.edgeR.top100 <- Group1and2.edgeR.filter[Group1and2.choose_gene_top100, ,drop=FALSE]
Group1and3.edgeR.top100 <- Group1and3.edgeR.filter[Group1and3.choose_gene_top100, ,drop=FALSE]
Group1and4.edgeR.top100 <- Group1and4.edgeR.filter[Group1and4.choose_gene_top100, ,drop=FALSE]
Group1and5.edgeR.top100 <- Group1and5.edgeR.filter[Group1and5.choose_gene_top100, ,drop=FALSE]
Group1and6.edgeR.top100 <- Group1and6.edgeR.filter[Group1and6.choose_gene_top100, ,drop=FALSE]
Group1and7.edgeR.top100 <- Group1and7.edgeR.filter[Group1and7.choose_gene_top100, ,drop=FALSE]
save(Group1and2.edgeR.top100,Group1and3.edgeR.top100,Group1and4.edgeR.top100, 
     Group1and5.edgeR.top100, Group1and6.edgeR.top100, Group1and7.edgeR.top100, file = 'Step2_top100gene.Rdata')

dbs <- "KEGG_2021_Human"
listEnrichrSites()
setEnrichrSite("Enrichr")
websiteLive <- TRUE
gene_up <- rownames(Group1and2.edgeR.top100[Group1and2.edgeR.top100$change == 'UP', ])
gene_down <- rownames(Group1and2.edgeR.top100[Group1and2.edgeR.top100$change == 'DOWN', ])
if (websiteLive) {
  enriched_UP <- enrichr(rownames(Group1and2.edgeR.top100), databases =dbs)
}
enriched_UP[["KEGG_2021_Human"]]


gene.df <- bitr(rownames(Group1and2.edgeR.top100),fromType="SYMBOL",toType=c("ENTREZID"),
                OrgDb = org.Hs.eg.db)



rownames(gene.df) <- gene.df[,1]
gene.df <- gene.df[,-1, drop = FALSE]

Group1and2.edgeR.top100 = transform(merge(Group1and2.edgeR.top100, gene.df, by=0, all.x=F, all.y = T),
                                   row.names = Row.names, Row.names = NULL)

gene_up <- Group1and2.edgeR.top100[Group1and2.edgeR.top100$change == 'UP', 'ENTREZID']
gene_down <- Group1and2.edgeR.top100[Group1and2.edgeR.top100$change == 'DOWN','ENTREZID']
gene_all <- as.character(Group1and2.edgeR.top100$ENTREZID)

Group1and2.kegg <- draw_kegg(gene_up = gene_up, gene_down = gene_down, gene_all = gene_all,
                             title = 'test')
Group1and2.kegg



#Group1and3==========================================================================
Group1and3.gene.df <- bitr(rownames(Group1and3.edgeR.top100),fromType="SYMBOL",toType=c("ENTREZID"),
                           OrgDb = org.Hs.eg.db)
rownames(Group1and3.gene.df) <- Group1and3.gene.df[,1]
Group1and3.gene.df <- Group1and3.gene.df[,-1, drop = FALSE]

Group1and3.edgeR.top100 = transform(merge(Group1and3.edgeR.top100, Group1and3.gene.df, by=0, all.x=F, all.y = T),
                                   row.names = Row.names, Row.names = NULL)

Group1and3.gene_up <- Group1and3.edgeR.top100[Group1and3.edgeR.top100$change == 'UP', 'ENTREZID']
Group1and3.gene_down <- Group1and3.edgeR.top100[Group1and3.edgeR.top100$change == 'DOWN','ENTREZID']
Group1and3.gene_all <- as.character(Group1and3.edgeR.top100$ENTREZID)

Group1and3.kegg <- draw_kegg(gene_up = Group1and3.gene_up, gene_down = Group1and3.gene_down, gene_all = Group1and3.gene_all,
                             title = 'BRCA WT TP53 vs TP53 multiple mut \n Pathway Enrichment')
Group1and3.kegg
#Group1and4==========================================================================
Group1and4.gene.df <- bitr(rownames(Group1and4.edgeR.top100),fromType="SYMBOL",toType=c("ENTREZID"),
                           OrgDb = org.Hs.eg.db)
rownames(Group1and4.gene.df) <- Group1and4.gene.df[,1]
Group1and4.gene.df <- Group1and4.gene.df[,-1, drop = FALSE]

Group1and4.edgeR.top100 = transform(merge(Group1and4.edgeR.top100, Group1and4.gene.df, by=0, all.x=F, all.y = T),
                                   row.names = Row.names, Row.names = NULL)

Group1and4.gene_up <- Group1and4.edgeR.top100[Group1and4.edgeR.top100$change == 'UP', 'ENTREZID']
Group1and4.gene_down <- Group1and4.edgeR.top100[Group1and4.edgeR.top100$change == 'DOWN','ENTREZID']
Group1and4.gene_all <- as.character(Group1and4.edgeR.top100$ENTREZID)

Group1and4.kegg <- draw_kegg(gene_up = Group1and4.gene_up, gene_down = Group1and4.gene_down, gene_all = Group1and4.gene_all,
                             title = 'BRCA WT TP53 vs TP53 R175 \n Pathway Enrichment')
Group1and4.kegg
#Group1and5==========================================================================
Group1and5.gene.df <- bitr(rownames(Group1and5.edgeR.top100),fromType="SYMBOL",toType=c("ENTREZID"),
                           OrgDb = org.Hs.eg.db)
rownames(Group1and5.gene.df) <- Group1and5.gene.df[,1]
Group1and5.gene.df <- Group1and5.gene.df[,-1, drop = FALSE]

Group1and5.edgeR.top100 = transform(merge(Group1and5.edgeR.top100, Group1and5.gene.df, by=0, all.x=F, all.y = T),
                                   row.names = Row.names, Row.names = NULL)

Group1and5.gene_up <- Group1and5.edgeR.top100[Group1and5.edgeR.top100$change == 'UP', 'ENTREZID']
Group1and5.gene_down <- Group1and5.edgeR.top100[Group1and5.edgeR.top100$change == 'DOWN','ENTREZID']
Group1and5.gene_all <- as.character(Group1and5.edgeR.top100$ENTREZID)

Group1and5.kegg <- draw_kegg(gene_up = Group1and5.gene_up, gene_down = Group1and5.gene_down, gene_all = Group1and5.gene_all,
                             title = 'BRCA WT TP53 vs TP53 R248 \n Pathway Enrichment')
Group1and5.kegg

#Group1and6==========================================================================
Group1and6.gene.df <- bitr(rownames(Group1and6.edgeR.top100),fromType="SYMBOL",toType=c("ENTREZID"),
                           OrgDb = org.Hs.eg.db)
rownames(Group1and6.gene.df) <- Group1and6.gene.df[,1]
Group1and6.gene.df <- Group1and6.gene.df[,-1, drop = FALSE]

Group1and6.edgeR.top100 = transform(merge(Group1and6.edgeR.top100, Group1and6.gene.df, by=0, all.x=F, all.y = T),
                                   row.names = Row.names, Row.names = NULL)

Group1and6.gene_up <- Group1and6.edgeR.top100[Group1and6.edgeR.top100$change == 'UP', 'ENTREZID']
Group1and6.gene_down <- Group1and6.edgeR.top100[Group1and6.edgeR.top100$change == 'DOWN','ENTREZID']
Group1and6.gene_all <- as.character(Group1and6.edgeR.top100$ENTREZID)

Group1and6.kegg <- draw_kegg(gene_up = Group1and6.gene_up, gene_down = Group1and6.gene_down, gene_all = Group1and6.gene_all,
                             title = 'BRCA WT TP53 vs TP53 R273 \n Pathway Enrichment')
Group1and6.kegg

#Group1and7==========================================================================
Group1and7.gene.df <- bitr(rownames(Group1and7.edgeR.top100),fromType="SYMBOL",toType=c("ENTREZID"),
                           OrgDb = org.Hs.eg.db)
rownames(Group1and7.gene.df) <- Group1and7.gene.df[,1]
Group1and7.gene.df <- Group1and7.gene.df[,-1, drop = FALSE]

Group1and7.edgeR.top100 = transform(merge(Group1and7.edgeR.top100, Group1and7.gene.df, by=0, all.x=F, all.y = T),
                                   row.names = Row.names, Row.names = NULL)

Group1and7.gene_up <- Group1and7.edgeR.top100[Group1and7.edgeR.top100$change == 'UP', 'ENTREZID']
Group1and7.gene_down <- Group1and7.edgeR.top100[Group1and7.edgeR.top100$change == 'DOWN','ENTREZID']
Group1and7.gene_all <- as.character(Group1and7.edgeR.top100$ENTREZID)
Group1and7.kegg <- draw_kegg(gene_up = Group1and7.gene_up, gene_down = Group1and7.gene_down, gene_all = Group1and7.gene_all,
                             title = 'BRCA WT TP53 vs TP53 Rmut \n Pathway Enrichment')
Group1and7.kegg