# function to draw heatmap and volcano plot

draw_heatmap <- function(V.method, type, gene.dataset){
  library("pheatmap")
  
  nrDEG_Z = nrDEG[order(nrDEG$logFC), ]
  nrDEG_F = nrDEG[order(-nrDEG$logFC), ]
  choose_gene = c(rownames(nrDEG_Z)[1:50], rownames(nrDEG_F)[1:50])
  choose_matrix = gene.dataset[choose_gene,]
  choose_matrix = t(scale(t(choose_matrix)))
  
  choose_matrix[choose_matrix > 2] = 2
  choose_matrix[choose_matrix < -2] = -2
  
  annotation_col = data.frame( CellType = factor( group_list ) )
  rownames( annotation_col ) = colnames( gene.dataset )
  filename <- paste('./fig/', type, '_heatmap_top100_logFC.png',
                    sep = "", collapse = NULL)
  pheatmap( fontsize = 6, choose_matrix, annotation_col = annotation_col, 
            show_rownames = T, show_colnames = F,
            annotation_legend = T, cluster_cols = T, 
            filename = filename)
}

draw_volcano <- function(nrDEG, type){
  library( "ggplot2" )
  logFC_cutoff <- with( nrDEG, mean( abs( logFC ) ) + 2 * sd( abs( logFC ) ) )
  
  nrDEG$change = as.factor( ifelse( 
    nrDEG$P.Value < 0.01 & abs(nrDEG$logFC) > logFC_cutoff,
    ifelse( nrDEG$logFC > logFC_cutoff, 'UP', 'DOWN' ), 'NOT' ) )
  nrDEGfile <- paste('./data/', type, '_nrDEG_by_logFC.Rdata',
                     sep = "", collapse = NULL)
  save( nrDEG, file = nrDEGfile )
  
  this_tile <- paste0( 
    'Cutoff for logFC is ', round( logFC_cutoff, 3 ),
    '
The number of up gene is ', nrow(nrDEG[ nrDEG$change == 'UP', ] ),
    '
The number of down gene is ', nrow(nrDEG[ nrDEG$change == 'DOWN', ] ) )
  
  volcano = ggplot(data = nrDEG, 
                   aes( x = logFC, y = -log10(P.Value), color = change)) +
    geom_point( alpha = 0.4, size = 1.75 ) +
    theme_set( theme_set( theme_bw( base_size = 15 ) ) ) +
    xlab( "log2 fold change" ) + ylab( "-log10 p-value" ) +
    ggtitle( this_tile ) + 
    theme( plot.title = element_text( size = 15, hjust = 0.5 )) +
    scale_colour_manual( values = c('blue','black','red') )
  print( volcano )
  filename <- paste('./fig/', type, '_volcano_logFC.png',
                    sep = "", collapse = NULL)
  ggsave( volcano, filename = filename )
  dev.off()
}