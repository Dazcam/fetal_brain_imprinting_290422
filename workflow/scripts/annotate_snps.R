#############################################################################
#
#    Annotating SNPs with MAF and genomic features for ase/imprinting
#    cross-referencing - With predetermined gene list
#
#############################################################################

# Note this version is different from cluster version

#####  Load packages  #####

library(readr)
library(ggplot2)
library(snplist)
library(biomaRt)

#####  References  #####

# Biomart vingette
# https://www.bioconductor.org/packages/devel/bioc/vignettes/biomaRt/inst/doc/biomaRt.html

# Variant annotation vingette
# https://bioconductor.org/packages/release/bioc/vignettes/VariantAnnotation/inst/doc/VariantAnnotation.pdf

#####  Set variables  #####

folder <- "~/Desktop/fine_mapping/"

genes <- read_delim(paste0(folder, "imprinted_gene_list.txt"), 
                    "\t", escape_double = FALSE, trim_ws = TRUE)

genes <- as.vector(genes$Gene) # Required for getBM (see gene_coord chunk)
genes <- sort(genes)

#####  Use Biormart to extract gene boundaries based on hgnc_symbol  #####
genemart = useMart(host="www.ensembl.org",
                   biomart = "ENSEMBL_MART_ENSEMBL",
                   dataset="hsapiens_gene_ensembl")

gene_coord <- getBM(attributes=c("hgnc_symbol", 
                                 "chromosome_name",
                                 "start_position", 
                                 "end_position"),
                    filters = c("hgnc_symbol"),
                    values = list(genes), mart = genemart)


## !! ONLY PICKING UP 217 OF 228 GENES - 11 Not found in BiomaRt !!  ##

# Note - there are 39 missing genes, however biomaRt detects more than one
# locus for some genes

colnames(gene_coord) <- c("geneID", "chr", "start", "stop") 
genes_with_coords <- gene_coord$hgnc_symbol
num_genes_with_coords <- length(unique(as.vector(gene_coord$geneID)))
missing_genes <- setdiff(genes, as.vector(genes_with_coords)) 

cat(paste0("-------------------------------------------------\n",
           "\nBiomart annotate genes\n\n",
           "Number of genes in initial list: ", length(genes),
           "\nNumber of genes dectected in biomart: ", nrow(gene_coord),
           "\nNumber of unique genes detected in biomart: ", num_genes_with_coords,
           "\nNumber of missing genes - not detected in biomart: ", length(missing_genes),
           "\n\n-------------------------------------------------\n\n"))

# Create tables for gene lists created in this section
write.table(genes, paste0(folder, "genes_input"), 
            quote = FALSE, sep = "\t", col.names = FALSE, row.names = FALSE )
write.table(genes_with_coords, paste0(folder, "genes_with_biomaRt_coords"), 
            quote = FALSE, sep = "\t", col.names = FALSE, row.names = FALSE )
write.table(missing_genes,  paste0(folder, "genes_without_biomaRt_coords"), 
            quote = FALSE, sep = "\t", col.names = FALSE, row.names = FALSE )


#####  Use Biormart to filter SNPs in gene boundaries on MAF > 0.05  #####
snpmart <- useMart(host="www.ensembl.org", 
                   biomart="ENSEMBL_MART_SNP", 
                   dataset="hsapiens_snp")

sink("annotated_imprinted_genes.log", append=TRUE, split=TRUE)

gene_coord

sink()

# dim(gene_coord)
# gene_coord <- gene_coord[11:135,] # For testing

sink("snps_gt0.05_in_genes.log", append=FALSE, split=TRUE)

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
  
    write.table(snps_maf, paste0("snps/", gene_ID, "_snps"), quote = FALSE, sep = "\t",
              col.names = FALSE, row.names = FALSE )
    
  }, error = function(e) {
    
    cat("ERROR: For ", gene_ID, "\n", conditionMessage(e), "\n")
    
  })
  
}  

sink()

# Dealing with timeouts - https://www.biostars.org/p/363530/
# listFilters(snpmart)
# listAttributes(snpmart)


# ------ Leave the feature annotation until after ase/imprinted gene SNP cross-ref

vcf <- read_delim(paste0("Desktop/fine_mapping/", chr, ".vcf"), 
                  "\t", escape_double = FALSE, trim_ws = TRUE)

# Use Variant Annotate to annotate variant to genomic feature
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(VariantAnnotation)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

rng <- GRanges(seqnames="14", ranges=IRanges(
  start=100779410,
  end=100861031,
  names="MEG3"))

tab <- TabixFile(paste0("Desktop/fine_mapping/", chr, ".vcf.gz"))
vcf_rng <- readVcf(tab, "hg38", param=rng)

seqlevels(vcf_rng) <- paste0("chr", chr)

rd <- rowRanges(vcf)
start <- c(100779410, )
stop <- c(100861031, )
gene <- c('MEG3', )

14: 100779410-
  
  var_list <- list(
    var1 = 1:10,          # 1, ..., 10
    var2 = letters[17:26] # q, ..., z
  )

for (i in 1:nrow(snps_maf)) {
  
  vcf_subset <- GRanges(seqnames=chr, ranges=IRanges(
    start=c(start),
    end=c(stop),
    names=c(gene)))
  
  chr <- snps_maf$chr_name[i])
start <- snps_maf$chrom_start[i]
stop <- snps_maf$chrom_end[i]

}

# Process each of var1 and var2
print(paste(
  'var1:', var_list$var1[i],
  '; var2:', var_list$var2[i]
))

}




AS3MT_ase <- read_delim("Desktop/fine_mapping/Results/AS3MT_ase.txt", 
                        "\t", escape_double = FALSE, trim_ws = TRUE)


martCheck(ENSEMBL_MART_SNP)

genes_failing_snp_search <- list()

for (i in 1:nrow(gene_coord)) {
  
  gene_ID <- gene_coord$geneID[i]
  print(gene_ID)
  genes_failing_snp_search <- list(genes_failing_snp_search, list(gene_ID))
  
}
  
sessionInfo()
#############################################################################
#############################################################################
