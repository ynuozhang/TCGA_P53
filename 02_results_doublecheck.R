#double check results with webdata
brca.webdata <- read.table('TCGA-BRCA-PRIMARY-FILTERED-WEB.tsv',
                           header = T,
                           sep = '\t');View(brca.webdata)

sum(is.na(brca.webdata$TP53))

