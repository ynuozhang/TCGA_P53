# function to draw heatmap and kegg plots

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

draw_selected_heatmap <- function(choose_matrix, annotation.class){
  library(ComplexHeatmap)
  library(circlize)
  
  col_fun = colorRamp2(c(-2, 0, 2), c("#2fa1dd", "white", "#f87669"))
  top_annotation = HeatmapAnnotation(
    cluster = anno_block(gp = gpar(fill = c("#2fa1dd", "#f87669")),
                         labels = unique(annotation.class),
                         labels_gp = gpar(col = "black", fontsize = 12)))
  
  fig_all <-  Heatmap(choose_matrix, col = col_fun, top_annotation = top_annotation,
                      show_heatmap_legend = F, column_split = annotation.class,
                      show_row_names = T, show_column_names = F,
                      cluster_columns = T, cluster_column_slices = F, row_names_gp = gpar(fontsize = 12))
  
  return(fig_all)
}

draw_kegg <- function(gene_up, gene_down, gene_all, title) {
  library(ggplot2)
  library(clusterProfiler)
  library(org.Hs.eg.db)
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
  
  dat$pvalue = -log10(dat$pvalue)
  dat$pvalue = dat$pvalue * dat$group
  dat = dat[ order(dat$pvalue, decreasing = F), ]
  
  if (length(unique(dat$group)) == 1){
    if (unique(dat$group) == 1) {
      g_kegg <- ggplot(dat,
                       aes(x = reorder(Description, order(pvalue, decreasing=F)), 
                           y = pvalue, fill = group)) + 
        geom_bar(stat = "identity") +
        scale_fill_gradient( na.value = "#f87669", guide = 'none') +
        scale_x_discrete( name = "Pathway names" ) +
        scale_y_continuous( name = "log10P-value" ) +
        coord_flip() + theme_bw() + theme(plot.title = element_text(hjust = 1, size=10),
                                          text = element_text(size=7)) +
        ggtitle(title)
    } else {
      g_kegg <- ggplot(dat,
                       aes(x = reorder(Description, order(pvalue, decreasing=F)), 
                           y = pvalue, fill = group)) + 
        geom_bar(stat = "identity") +
        scale_fill_gradient( na.value = "#2fa1dd", guide = 'none') +
        scale_x_discrete( name = "Pathway names" ) +
        scale_y_continuous( name = "log10P-value" ) +
        coord_flip() + theme_bw() + theme(plot.title = element_text(hjust = 1, size=10),
                                          text = element_text(size=7)) +
        ggtitle(title)
      
    } }else {
      g_kegg <- ggplot(dat,
                       aes(x = reorder(Description, order(pvalue, decreasing=F)), 
                           y = pvalue, fill = group)) + 
        geom_bar(stat = "identity") +
        scale_fill_gradient( low = "#2fa1dd", high = "#f87669", guide = 'none' ) +
        scale_x_discrete( name = "Pathway names" ) +
        scale_y_continuous( name = "log10P-value" ) +
        coord_flip() + theme_bw() + theme(plot.title = element_text(hjust = 1, size=10),
                                          text = element_text(size=7)) +
        ggtitle(title)
    }
  
  
  
  filename <- paste('./figs/kegg_', gsub("[\n]", "", title), '.png', sep = "", collapse = NULL)
  ggsave(g_kegg, filename = filename, height = 1181 , width = 1206, units = "px")
  return(g_kegg)
}



