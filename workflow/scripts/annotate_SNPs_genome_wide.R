#############################################################################
#
#    Annotating SNPs with MAF and genomic features for ase/imprinting
#    cross-referencing - GENOME WIDE
#
#############################################################################


#####  Load packages  #####

library(readr)
library(ggplot2)
# library(snplist)
library(biomaRt)
library(dplyr)

#####  Use Biormart to get gene list and genes boundaries based on hgnc_symbol  #####
genemart = useMart(host="www.ensembl.org",
                   biomart = "ENSEMBL_MART_ENSEMBL",
                   dataset="hsapiens_gene_ensembl")

gene_coord <- getBM(attributes=c("hgnc_symbol", 
                                 "chromosome_name",
                                 "start_position", 
                                 "end_position",
                                 "gene_biotype"),
                    filters = "biotype", # Not gene_biotype for some reason!
                    values = "protein_coding" , # 67159 genes without this line
                    mart = genemart)

colnames(gene_coord) <- c("geneID", "chr", "start", "stop", "gene_biotype") 

listFilters(genemart)
listAttributes(genemart)

# Potential issues:
# A: Some entries have no hgnc_symbol 
no_gene_symbol <- gene_coord[which(gene_coord$geneID==""),] 

# B: Some entires from chr patches are included
full_chr_list <- unique(gene_coord$chr)

# Filter - Remove mitochondrial reads, genes on patches and genes with no hngc symbol
genes_on_patches <- filter(gene_coord, grepl("_|K|G|MT", chr))
final_gene_list <- filter(gene_coord, !grepl("_|K|G|MT", chr))
final_gene_list <- final_gene_list[!(final_gene_list$geneID==""),]
final_gene_list <- final_gene_list[order(final_gene_list$geneID),]
final_chr_list <- unique(final_gene_list$chr) # Only chr 1-22, X and Y

# Inititalise log file
sink("gw_snps_gt0.05_in_genes.log", append=FALSE, split=TRUE)

cat(paste0("-----------------------------------------------------------------\n\n",
           "BIOMART - Get protein coding genes for genome-wide analysis\n\n",
           "-----------------------------------------------------------------\n\n",
           "Protein coding genes reported in biomart: ", nrow(gene_coord),
           "\nGenes with no gene symbol: ", nrow(no_gene_symbol),
           "\nGenes on patches, MT, G or K chr: ", nrow(genes_on_patches),
           "\nGenes taken futher analysis - after filtering: ", nrow(final_gene_list),
           "\nUnique chrs before filtering: \n\n"))

print(full_chr_list)

cat("\nUnique chrs after filtering: \n\n")

print(final_chr_list)

cat("\n\n-------------------------------------------------\n\n")

#####  Use Biormart to filter SNPs in gene boundaries on MAF > 0.05  #####
snpmart <- useMart(host="www.ensembl.org", 
                   biomart="ENSEMBL_MART_SNP", 
                   dataset="hsapiens_snp")

# dim(gene_coord)
final_gene_list_subset <- final_gene_list[1:5,] # For testing

cat(paste0("-----------------------------------------------------------------\n\n",
           "BIOMART - Get SNPs in coding genes with MAF >= 0.05\n\n",
           "-----------------------------------------------------------------\n\n"))

snp_error_df <- data.frame(geneID = character(), error = character(), stringsAsFactors=FALSE) 

for (i in 1:nrow(final_gene_list_subset)) {
  
  tryCatch({
    
    gene_ID <- final_gene_list_subset$geneID[i]
    gene_chr <- final_gene_list_subset$chr[i]
    gene_start <- final_gene_list_subset$start[i]
    gene_stop <- final_gene_list_subset$stop[i]
    
    cat("\n----------------------------------------------------------\n")
    cat(paste0("Retrieving SNPs for: ", gene_ID, " - chr", gene_chr,":", gene_start, "-", gene_stop))
    cat("\n----------------------------------------------------------\n")
    
    snps <- getBM(attributes=c("refsnp_id", "allele", "chr_name", "chrom_start", "chrom_end", "minor_allele_freq"),
                  filters = c("chr_name", "start", "end"),
                  values = list(gene_chr, gene_start, gene_stop), mart = snpmart)
    
    snps_maf <- snps[ which(snps$minor_allele_freq >= 0.05), ]
    
    cat(paste0("Total loci in ", gene_ID, ": ", nrow(snps)))
    cat(paste0("\nTotal loci after MAF filter (>= 0.05): ", nrow(snps_maf)))
    cat("\n----------------------------------------------------------\n")
    
    write.table(snps_maf, paste0("snps/", gene_ID, "_snps"), quote = FALSE, sep = "\t",
                col.names = FALSE, row.names = FALSE )
    
  }, error = function(e) {
    
    cat("ERROR: For ", gene_ID, "\n\n", conditionMessage(e), "\n\n")
    snp_error <- matrix(c(gene_ID, conditionMessage(e)), nrow = 1)
    snp_error_df <<- rbind(snp_error_df, snp_error)
    
  })
  
}  

sink()

# Dealing with timeouts - https://www.biostars.org/p/363530/
# listFilters(snpmart)
# listAttributes(snpmart)


#############################################################################
#############################################################################
