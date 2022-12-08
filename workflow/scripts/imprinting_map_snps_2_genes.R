#--------------------------------------------------------------------------------------
#
#    Map snps to genes
#
#--------------------------------------------------------------------------------------

##  Resources  ------------------------------------------------------------------------
# Note that snakemake is tracking the log file not the individual gene_snp output files

## Initialise R library  --------------------------------------------------------------
system('export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/scratch/c.c1477909/.conda/envs/snakemake/lib64/')
.libPaths( c( "/scratch/c.c1477909/R/library", .libPaths() ) )

##  Load Packages  --------------------------------------------------------------------
library(tidyverse)
library(biomaRt)
library(readxl)

## Set variables  ---------------------------------------------------------------------
IMP_GENES <- snakemake@input$genes
LOG_FILE <- snakemake@output
OUT_DIR <- snakemake@params

##  Load Data  ------------------------------------------------------------------------
genes <- read_excel(IMP_GENES, sheet = 'Human_All_Known IG_NoStatus') %>%
  filter(!grepl("placenta", Gene)) %>%
  drop_na(Gene) %>%
  arrange(Gene) %>%
  pull(Gene) 
  
##  Use Biormart to extract gene boundaries based on hgnc_symbol  --------------------
genemart = useMart(host="www.ensembl.org",
                   biomart = "ENSEMBL_MART_ENSEMBL",
                   dataset="hsapiens_gene_ensembl")

gene_coord <- getBM(attributes=c("hgnc_symbol", 
                                 "chromosome_name",
                                 "start_position", 
                                 "end_position",
                                 "gene_biotype"),
                    filters = c("hgnc_symbol"),
                    values = list(genes), mart = genemart)

#searchAttributes(mart = genemart)

## !! ONLY PICKING UP 217 OF 228 GENES - 11 Not found in BiomaRt !!  ##

# Note - there are 39 missing genes, however biomaRt detects more than one
# locus for some genes

# Extract and print gene summary info
colnames(gene_coord) <- c("geneID", "chr", "start", "stop", "type") 
genes_with_coords <- gene_coord$geneID
num_genes_with_coords <- length(unique(as.vector(genes_with_coords)))
missing_genes <- setdiff(genes, as.vector(genes_with_coords)) 

cat(paste0("-------------------------------------------------\n",
           "\nBiomart annotate genes\n\n",
           "Number of genes in initial list: ", length(genes),
           "\nNumber of genes dectected in biomart: ", nrow(gene_coord),
           "\nNumber of unique genes detected in biomart: ", num_genes_with_coords,
           "\nNumber of missing genes - not detected in biomart: ", length(missing_genes),
           "\n\n-------------------------------------------------\n\n"))

# # Create tables for gene lists created in this section
# write.table(genes, paste0(SCRATCH, "genes_input"), 
#             quote = FALSE, sep = "\t", col.names = FALSE, row.names = FALSE )
# write.table(genes_with_coords, paste0(SCRATCH, "genes_with_biomaRt_coords"), 
#             quote = FALSE, sep = "\t", col.names = FALSE, row.names = FALSE )
# write.table(missing_genes,  paste0(SCRATCH, "genes_without_biomaRt_coords"), 
#             quote = FALSE, sep = "\t", col.names = FALSE, row.names = FALSE )
# write.table(gene_coord,  paste0(SCRATCH, "genes_without_biomaRt_coords_df"), 
#             quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE )

#####  Use Biormart to filter SNPs in gene boundaries on MAF > 0.05  #####
snpmart <- useMart(host="www.ensembl.org", 
                   biomart="ENSEMBL_MART_SNP", 
                   dataset="hsapiens_snp")

sink(paste0(LOG_FILE), append = FALSE, split = TRUE)

for (i in 1:nrow(gene_coord)) {
  
  tryCatch({
    
    gene_ID <- gene_coord$geneID[i]
    gene_chr <- gene_coord$chr[i]
    gene_start <- gene_coord$start[i]
    gene_stop <- gene_coord$stop[i]
    
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
    
    write.table(snps_maf, paste0(OUT_DIR, gene_ID, "_snps"), quote = FALSE, sep = "\t",
                col.names = FALSE, row.names = FALSE )
    
  }, error = function(e) {
    
    cat("ERROR: For ", gene_ID, "\n", conditionMessage(e), "\n")
    
  })
  
}  

sink()


#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
