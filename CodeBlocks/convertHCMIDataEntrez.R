convertHCMIData=function(){
  library(biomaRt)
  load("HCMI_GDC_gx_plus_annotations_20241123.RData")
  ensembl.genes <-gsub("\\..*","",rownames(fpkm_uq_unstrand)) 
  mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))

  genes <- getBM(
    filters="ensembl_gene_id",
    attributes=c("ensembl_gene_id", "entrezgene_id"),
    values=ensembl.genes,
    mart=mart)

  entrezids=genes[match(ensembl.genes, unlist(genes$ensembl_gene_id)),"entrezgene_id"]
  rm_genes=which(is.na(entrezids))
  rm_genes2=which(duplicated(entrezids))
  all_rms=union(rm_genes, rm_genes2)
  fpkm_uq_unstrand=fpkm_uq_unstrand[-all_rms, ]
  rownames(fpkm_uq_unstrand)=entrezids[-all_rms]
  rownames(fpkm_uq_unstrand)=paste("T", rownames(fpkm_uq_unstrand), sep="")
  return(fpkm_uq_unstrand)

}