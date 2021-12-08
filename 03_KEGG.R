library( "clusterProfiler" )
library( "org.Hs.eg.db" )
#library( "ggplot2" )
source('00_drawing_function.R')
load('data/step2.Rdata')
gene.df <- bitr(rownames(Group1and2.edgeR.filter),fromType="SYMBOL",toType=c("ENTREZID"),
                OrgDb = org.Hs.eg.db)
rownames(gene.df) <- gene.df[,1]
gene.df <- gene.df[,-1, drop = FALSE]

Group1and2.edgeR.filter = transform(merge(Group1and2.edgeR.filter, gene.df, by=0, all.x=F, all.y = T),
                                    row.names = Row.names, Row.names = NULL)

gene_up <- Group1and2.edgeR.filter[Group1and2.edgeR.filter$change == 'UP', 'ENTREZID']
gene_down <- Group1and2.edgeR.filter[Group1and2.edgeR.filter$change == 'DOWN','ENTREZID']
gene_all <- as.character(Group1and2.edgeR.filter$ENTREZID)

Group1and2.kegg <- draw_kegg(gene_up = gene_up, gene_down = gene_down, gene_all = gene_all,
                             title = 'BRCA WT TP53 vs TP53 Mutant (single SNP) \n Pathway Enrichment')

kk.up <- enrichKEGG(gene =  gene_up,
                    organism = 'hsa',
                    pvalueCutoff  =  0.8,
                    qvalueCutoff  =  0.8,
                    universe = gene_all)
                       
kk.down <- enrichKEGG(gene =  gene_down,
                      organism = 'hsa',
                      pvalueCutoff  =  0.8,
                      qvalueCutoff  =  0.8,
                      universe = gene_all)

kegg_down_dt <- as.data.frame( kk.down )
kegg_up_dt <- as.data.frame( kk.up )
down_kegg <- kegg_down_dt[ kegg_down_dt$pvalue < 0.05, ]
up_kegg <- kegg_up_dt[ kegg_up_dt$pvalue < 0.05, ]


if (dim(up_kegg)[1] == 0){
  print('No up regulated gene was enriched.')
  down_kegg$group = -1
  dat <- down_kegg
} else if (dim(down_kegg)[1] == 0) {
  print('No down regulated gene was enriched.')
  up_kegg$group = 1
  dat = up_kegg
} else {
  up_kegg$group = 1
  down_kegg$group = -1
  dat <- rbind(up_kegg, down_kegg)
}

dat$pvalue = -log10( dat$pvalue )
dat$pvalue = dat$pvalue * dat$group
dat = dat[ order( dat$pvalue, decreasing = F ), ]

title = 'BRCA WT TP53 vs TP53 Mutant (single SNP) \n Pathway Enrichment'
g_kegg <- ggplot(dat,
                 aes(x = reorder(Description, order(pvalue, decreasing=F)), y = pvalue, fill = group)) + 
  geom_bar(stat = "identity") +
  scale_fill_gradient( low = "#2fa1dd", high = "#f87669", guide = 'none') +
  scale_x_discrete( name = "Pathway names" ) +
  scale_y_continuous( name = "log10P-value" ) +
  coord_flip() + theme_bw() + theme(plot.title = element_text(hjust = 1, size=10),
                                    text = element_text(size=7)) +
  ggtitle(title)
filename <- paste('./figs/kegg_.png', collapse = NULL)
ggsave(g_kegg, width = 1920/72, height = 1080/72, filename = filename)
print( g_kegg )


