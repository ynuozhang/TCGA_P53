# function to draw heatmap and volcano plots

top100_heatmap <- function(V.method, gene.dataset, annotation.class){
  library(ComplexHeatmap)
  library(circlize)
  
  edgeR.filter <- filter(V.method, abs(V.method$logFC) >=1
                         & V.method$P.Value < 0.05)
  nrDEG_Z = edgeR.filter[order(edgeR.filter$logFC), ]
  nrDEG_F = edgeR.filter[order(-edgeR.filter$logFC), ]
  choose_gene_top100 = c(rownames(nrDEG_Z)[1:50], rownames(nrDEG_F)[1:50])
  choose_matrix_top100 = gene.dataset[choose_gene_top100,]
  choose_matrix_top100 = log2(choose_matrix_top100 + 1)
  choose_matrix_top100 = t(scale(t(choose_matrix_top100)))
  choose_matrix_top100[choose_matrix_top100 > 2] = 2
  choose_matrix_top100[choose_matrix_top100 < -2] = -2
  
  col_fun = colorRamp2(c(-2, 0, 2), c("#2fa1dd", "white", "#f87669"))
  top_annotation = HeatmapAnnotation(
    cluster = anno_block(gp = gpar(fill = c("#2fa1dd", "#f87669")),
                         labels = unique(annotation.class),
                         labels_gp = gpar(col = "black", fontsize = 12)))
  
  fig_top100 <-  Heatmap(choose_matrix_top100, col = col_fun, top_annotation = top_annotation,
                  column_split = annotation.class, show_heatmap_legend = F,
                  show_row_names = T, show_column_names = F,column_title = 'top100', 
                  cluster_columns = F, cluster_column_slices = T, row_names_gp = gpar(fontsize = 6))
  
  return(fig_top100)
}

draw_heatmap <- function(choose_matrix, annotation.class){
  library(ComplexHeatmap)
  library(circlize)
  
  col_fun = colorRamp2(c(-2, 0, 2), c("#2fa1dd", "white", "#f87669"))
  top_annotation = HeatmapAnnotation(
    cluster = anno_block(gp = gpar(fill = c("#2fa1dd", "#f87669")),
                         labels = unique(annotation.class),
                         labels_gp = gpar(col = "black", fontsize = 12)))
  
  fig_all <-  Heatmap(choose_matrix, col = col_fun, top_annotation = top_annotation,
                      show_heatmap_legend = F, column_split = annotation.class,
                      show_row_names = T, show_column_names = F,column_title = 'All_significant',
                      cluster_columns = F, cluster_column_slices = T, row_names_gp = gpar(fontsize = 1))
  
  return(fig_all)
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