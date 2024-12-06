args <- commandArgs(TRUE)
results_file <- args[1]
geno_file_stem <- args[2]
snps_for_annotation <- args[3]
inferno_input_loc <- args[4]

#results file should be meta-analysis results from GWAMA (path and all)
#geno_file_stem should be a genotype file set stem (path and all) from which to pull position and MAF
#snps_for_annotation should be a comma-separated (no space) list of variants for annotation

library(data.table)

#fix the snps_for_annotation into a vector
snps_for_annotation <- unlist(strsplit(snps_for_annotation, ","))

#get inferno input file name
  #strip the file path
inferno_input_file_name <- sapply(strsplit(results_file, "/"), "[", length(strsplit(results_file, "/")[[1]]))
  #replace extension with .forINFERNO.txt and add the file path
extension <- sapply(strsplit(inferno_input_file_name, "[.]"), "[", length(strsplit(inferno_input_file_name, "[.]")[[1]]))
inferno_input_file_name <- paste0(inferno_input_loc, "/", gsub(extension, "forINFERNO.txt", inferno_input_file_name))

#create corresponding SNP file name
inferno_snp_file_name <- gsub(".forINFERNO", ".topSNPs.forINFERNO", inferno_input_file_name)

### Now actually create the files ###

#read in summary stats
results <-  fread(results_file)

#get minor allele frequency from frq file from the larger set
frq <- fread(paste0(geno_file_stem, ".frq"))
frq <- frq[,c("SNP", "MAF")]
names(frq) <- c("rs_number", "MAF")
results <- merge(results, frq, by = "rs_number")

#get chromosome and position from the bim file
bim <- fread(paste0(geno_file_stem, ".bim"))
bim <- bim[,c("V1", "V2", "V4")]
names(bim) <- c("chr", "rs_number", "pos")
results <- merge(results, bim, by = "rs_number")
results <- results[,c("rs_number", "chr", "pos", "reference_allele",  "other_allele", "MAF", "beta", "se", "p-value")]

#replace colon with dash for variants with only chr:bp for name
results$rs_number <- gsub(":", "-", results$rs_number)

#save input file with results
write.table(results, inferno_input_file_name, row.names = F, col.names = F, quote  = F, sep = "\t")

#create top SNP list
snps <- results[results$rs_number %in% snps_for_annotation,c("chr", "rs_number", "pos")]
snps$label <- paste0(snps$chr, "-", snps$pos)

#make the chromosome column be "chr#" not just the number
snps$chr <- paste0("chr", snps$chr)

#organize columns
snps <- snps[,c("chr", "rs_number", "label", "pos")]

#write SNP file
write.table(snps, inferno_snp_file_name, row.names = F, col.names = F, quote  = F, sep = "\t")

