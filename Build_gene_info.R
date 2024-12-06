library(data.table)
library(biomaRt)

genes <- fread("/scratch/mahone1/PrediXcan_Resilience/Scripts/gene_names_new.txt",  header = F, stringsAsFactors=F)

ensembl <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl",  host="grch37.ensembl.org")

gene_info <- getBM(attributes = c("hgnc_symbol", "chromosome_name", "start_position", "end_position", "ensembl_gene_id"),
                 filters = c("ensembl_gene_id"),
                 values = genes$V1,
                 mart = ensembl)
#this gets ids for all but 51 genes

gene_info$hgnc_symbol[gene_info$hgnc_symbol==""] <- "NA"

#save
write.table(gene_info, "/scratch/mahone1/PrediXcan_Resilience/gene_info.txt", row.names = F, quote = F)