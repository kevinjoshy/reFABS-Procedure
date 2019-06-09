library(org.Hs.eg.db)

gene.ann = read.csv('geneAnnotationsDF.csv',as.is=T,row.names=1)
x = mget(gene.ann$gene_ensID,org.Hs.egENSEMBL2EG,ifnotfound=NA)

entid = integer(length = nrow(gene.ann))

for(i in 1:length(entid))
{
  entid[i] = x[[i]][1]
}


gene.anns.ext = data.frame(gene.ann,entrezid = entid)
gene.anns.selected = gene.anns.ext[!is.na(gene.anns.ext$entrezid),]

# Writing entrezids into CSV
# write.csv(gene.anns.ext, 'geneAnnotationsDF_entrezID.csv')

# Writing entrezids without NA into CSV
# write.csv(gene.anns.selected, 'geneAnnotationsDF_Selected_entrezID.csv')
