args <- commandArgs(TRUE)
cohorts <- args[1] #comma-separated, no spaces
genes <- args[2] #comma-separated ensembl ids, no spaces, no version number
gene_label_type <- args[3] #names or ids (no version number)
loc <- args[4] #where to save it
lab <- args[5] #label for the file to indicate what it is

# cohorts <- "VMAP"
# genes <- "ENSG00000087245,ENSG00000149968,ENSG00000100985,ENSG00000035862,ENSG00000102265"
# gene_label_type <- "ids"
# loc <- "/nfs/mahone1/MMP3_predicted_expression/"
# lab <- "MMP_TIMP"
# i <- 1

#parse
cohorts <- unlist(strsplit(cohorts, ","))
genes <- unlist(strsplit(genes, ","))

suppressMessages({
  library(data.table)
  library(stringr)
  library(dplyr)
})

print(paste0("Extracting predicted expression for ", length(genes), " genes in ", paste(cohorts, collapse = ","), " in GTEx tissues."))

gtex_gene_info <- fread("/scratch/mahone1/scripts/gtex_genes.txt", header = T)
commonmind_gene_info <- fread("/scratch/mahone1/scripts/commonmind_genes.txt", header = T)
monocyte_gene_info <- fread("/scratch/mahone1/scripts/monocyte_genes.txt", header = T)

if(gene_label_type == "ids"){
  #gtex
  gtex_genes <- gtex_gene_info$ensembl_gene_id_vnum[gtex_gene_info$ensembl_gene_id %in% genes]
  names(gtex_genes) <- gtex_gene_info$hgnc_symbol[gtex_gene_info$ensembl_gene_id %in% genes]
  
  #commonmind
  commonmind_genes <- commonmind_gene_info$ensembl_gene_id[commonmind_gene_info$ensembl_gene_id %in% genes] 
  names(commonmind_genes) <- commonmind_gene_info$hgnc_symbol[commonmind_gene_info$ensembl_gene_id %in% genes]
  
  #monocyte
  monocyte_genes <- monocyte_gene_info$ensembl_gene_id_vnum[monocyte_gene_info$ensembl_gene_id %in% genes]
  names(monocyte_genes) <- monocyte_gene_info$hgnc_symbol[monocyte_gene_info$ensembl_gene_id %in% genes]
  
} else if(gene_label_type == "names"){
  #gtex
  gtex_genes <- gtex_gene_info$ensembl_gene_id_vnum[gtex_gene_info$hgnc_symbol %in% genes]
  names(gtex_genes) <- gtex_gene_info$hgnc_symbol[gtex_gene_info$hgnc_symbol %in% genes]
  
  #commonmind
  commonmind_genes <- commonmind_gene_info$ensembl_gene_id[commonmind_gene_info$hgnc_symbol %in% genes] 
  names(commonmind_genes) <- commonmind_gene_info$hgnc_symbol[commonmind_gene_info$hgnc_symbol %in% genes]
  
  #monocyte
  monocyte_genes <- monocyte_gene_info$ensembl_gene_id_vnum[monocyte_gene_info$hgnc_symbol %in% genes]
  names(monocyte_genes) <- monocyte_gene_info$hgnc_symbol[monocyte_gene_info$hgnc_symbol %in% genes]
  
} else {
  stop("Please specify valid gene label type!")
}

print(paste0(length(gtex_genes), " out of ", length(genes), " requested genes present in GTEx predicted expression models: ", paste0(names(gtex_genes), collapse = ",")))
print(paste0(length(commonmind_genes), " out of ", length(genes), " requested genes present in CommonMind predicted expression models: ", paste0(names(commonmind_genes), collapse = ",")))
print(paste0(length(monocyte_genes), " out of ", length(genes), " requested genes present in Monocyte predicted expression models: ", paste0(names(monocyte_genes), collapse = ",")))

tissues <- readLines("/scratch/mahone1/PrediXcan_Cognition/Scripts/tissue_names.txt")

for(i in 1:length(cohorts)){
    
  cohort <- cohorts[i]
  
  #get index of tissues with predicted expression
  index <- sapply(tissues, function(x) grepl(paste0(gtex_genes, collapse="|"), readLines(paste0("/nfs/DATA/", cohort, "/PrediXcan/GTEx/", x, "_predicted_expression.txt"), n=1)))
  
  #subset to a vector of only those tissues
  tissues_wexp <- tissues[index]
  print(paste0(length(tissues_wexp), " tissues have predicted expression of one or more genes requested."))
  
  #read in the predicted expression for the genes of interest in the tissue with predicted expression; 
  #data will be a list with one data.table per tissue with expression
  ids <- fread(paste0("/nfs/DATA/", cohort, "/PrediXcan/GTEx/Adipose_Subcutaneous_predicted_expression.txt"), select=c("FID", "IID"))
  data <- lapply(tissues_wexp, function(x) fread(paste0("/nfs/DATA/", cohort, "/PrediXcan/GTEx/", x, "_predicted_expression.txt"), select=c(gtex_genes)))
  
  #set the names of the data.tables to be the tissue names
  names(data) <- tissues_wexp
  
  #rename columns with gene_name (if there are no genes "named" NA)
  if(sum(names(gtex_genes)=="NA")==0){
    data <- lapply(names(data), function(x) setNames(data[[x]], names(gtex_genes)[gtex_genes %in% names(data[[x]])]))
  }
  
  #add back tissue names
  names(data) <- tissues_wexp
  
  #rename the gene columns to keep them distinct ENSG0000123_TISSUE_NAME (have to keep ids separately)
  data <- lapply(names(data), function(x) setNames(data[[x]], paste0(names(data[[x]]), "_", x)))
  
  #make into a single dataframe
  data <- do.call(cbind, data)
  
  #add back ids
  data <- cbind(ids, data)
  
  ## attempt to add in commonmind and monocyte
  if(file.exists(paste0("/nfs/DATA/", cohort, "/PrediXcan/CommonMind/CommonMind_predicted_expression.txt")) & length(commonmind_genes)>0){
    cm <- fread(paste0("/nfs/DATA/", cohort, "/PrediXcan/CommonMind/CommonMind_predicted_expression.txt"), select=c("FID", "IID", genes))
    names(cm)[3:length(names(cm))] <- paste0(names(commonmind_genes), "_CommonMind")
    data <- merge(data, cm, by = c("FID", "IID"), all = T)
  } else {
    print(paste0("No CommonMind predicted expression for ", cohort, "."))
  }
  
  if(file.exists(paste0("/nfs/DATA/", cohort, "/PrediXcan/Monocyte/Monocyte_predicted_expression.txt"))){
    monocyte <- fread(paste0("/nfs/DATA/", cohort, "/PrediXcan/Monocyte/Monocyte_predicted_expression.txt"), select=c("FID", "IID", monocyte_genes))
    names(monocyte)[3:length(names(monocyte))] <- paste0(names(monocyte_genes), "_Monocyte")
    data <- merge(data, monocyte, by = c("FID", "IID"), all = T)
  } else {
    print(paste0("No monocyte predicted expression for ", cohort, "."))
  }
  
  #drop genes that are all zeros
  print(paste0(paste0(names(data)[!lengths(sapply(data, table))>1], collapse = ","), " has only one value and is being dropped."))
  data <- data[,lengths(sapply(data, table))>1, with = FALSE]
  
  #print out which genes are in which tissues
  tmp <- data.frame(do.call(rbind, str_split(names(data)[3:length(names(data))], "_", 2)), stringsAsFactors = F)
  names(tmp) <- c("gene", "tissue")
  tmp %>% group_by(gene) %>% summarise( tissues = paste(tissue, collapse = ", ")) %>% as.data.frame() %>% print()
  
  #write out the predicted expression
  write.csv(data, paste0(loc, "/", cohort, "_", lab, "_predicted_expression.csv"), row.names = F, quote = F)
  
  print(paste0(paste(cohorts, collapse = ","), " predicted expression ", " saved to ", loc, "."))
}

