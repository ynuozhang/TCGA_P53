rm(list = objects( all = TRUE ))

bioPackages <- 
  c( 
    "dplyr", "stringi", "purrr", ## ERROR
    "R.utils", "data.table", ## unzip and read table
    "GEOquery", ## download
    "FactoMineR", "factoextra", "ggfortify", ## PCA
    "pheatmap", ## heatmap
    "ggplot2", ## Volcano plot
    "limma", "DESeq2", "edgeR", ## DEG
    "clusterProfiler", "org.Hs.eg.db", ## annotation
    "pathview" ## kegg
  )

lapply( bioPackages, 
        function(bioPackage) {
          if ( !require( bioPackage, character.only = T ) ) {
            CRANpackages <- available.packages()
            
            ## install packages by CRAN
            if ( bioPackage %in% rownames( CRANpackages ) ) {
              install.packages( bioPackage )
              
            }else{
              ## install packages by bioconductor
              ## R version >= 3.5 ===> BiocManager
              if (!requireNamespace("BiocManager", quietly = TRUE))
                install.packages("BiocManager")
              BiocManager::install(bioPackage, update = TRUE, ask = FALSE)
            }
          }
        }
)
dir.create("data")
