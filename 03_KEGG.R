library( "clusterProfiler" )
library( "org.Hs.eg.db" )
k1 = (DEG$pvalue < 0.05)&(DEG$log2FoldChange < -logFC_cutoff)
k2 = (DEG$pvalue < 0.05)&(DEG$log2FoldChange > logFC_cutoff)
DEG$change = ifelse(k1,"DOWN",ifelse(k2,"UP","NOT"))

kegg_plot <- function(type) {
  kk.up <- enrichKEGG(   gene          =  gene_up    ,
                         organism      =  'hsa'      ,
                         universe      =  gene_all   ,
                         pvalueCutoff  =  0.8        ,
                         qvalueCutoff  =  0.8        )
  kk.down <- enrichKEGG( gene          =  gene_down  ,
                         organism      =  'hsa'      ,
                         universe      =  gene_all   ,
                         pvalueCutoff  =  0.8        ,
                         qvalueCutoff  =  0.8        )
  library( "ggplot2" )
  kegg_down_dt <- as.data.frame( kk.down )
  kegg_up_dt <- as.data.frame( kk.up )
  down_kegg <- kegg_down_dt[ kegg_down_dt$pvalue < 0.05, ]
  down_kegg$group <- -1
  up_kegg <- kegg_up_dt[ kegg_up_dt$pvalue < 0.05, ]
  up_kegg$group <- 1
  dat = rbind( up_kegg, down_kegg )
  dat$pvalue = -log10( dat$pvalue )
  dat$pvalue = dat$pvalue * dat$group
  dat = dat[ order( dat$pvalue, decreasing = F ), ]
  g_kegg <- ggplot( dat, 
                    aes(x = reorder( Description, order( pvalue, decreasing = F )), 
                        y = pvalue, fill = group)) + 
    geom_bar( stat = "identity" ) + 
    scale_fill_gradient( low = "blue", high = "red", guide = FALSE ) + 
    scale_x_discrete( name = "Pathway names" ) +
    scale_y_continuous( name = "log10P-value" ) +
    coord_flip() + theme_bw() + 
    theme( plot.title = element_text( hjust = 0.5 ), 
           axis.text.x = element_text(size = 10),
           axis.text.y = element_text(size = 7)) +
    ggtitle( "Pathway Enrichment" ) 
  print( g_kegg )
  filename <- paste('./fig/kegg_up_down_', type, '.png', sep = "", collapse = NULL)
  ggsave( g_kegg, filename = filename )
}