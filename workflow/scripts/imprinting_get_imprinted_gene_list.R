#--------------------------------------------------------------------------------------
#
#    Get imprinted genes
#
#--------------------------------------------------------------------------------------

##  Resources  ------------------------------------------------------------------------
# Note that snakemake is tracking the log file not the individual gene_snp output files

## Initialise R library  --------------------------------------------------------------
system('export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/scratch/c.c1477909/.conda/envs/snakemake/lib64/')
.libPaths( c( "/scratch/c.c1477909/R/library", .libPaths() ) )

##  Load Packages  --------------------------------------------------------------------
library(tidyverse)
library(readxl)

## Set variables  ---------------------------------------------------------------------
TUCCI_GENES <- snakemake@input$tucci_genes
GW_GENES <- snakemake@input$genomewide_genes
LOG_FILE <- toString(snakemake@log)
OUT_DIR <- snakemake@params

##  Load Data  ------------------------------------------------------------------------
tucci_genes <- read_excel(TUCCI_GENES, sheet = 'Human_All_Known IG_NoStatus') %>%
  filter(!grepl("placenta", Gene)) %>%
  drop_na(Gene) %>%
  arrange(Gene) %>%
  select(Gene) %>%
  write_tsv(paste0(OUT_DIR, 'Tucci_2019_imprinted_genes.txt'), col_names = FALSE)
  
genes_gw <- read_tsv(GW_GENES, col_names = FALSE) %>%
  rename(genes = X1)

tucci_genes_in_gw <- tucci_genes %>%
  filter(Gene %in% genes_gw$genes) %>%
  write_tsv(paste0(OUT_DIR, 'Tucci_2019_genes_in_gw_list.txt'), col_names = FALSE)

tucci_genes_not_in_gw <- tucci_genes %>%
  filter(!Gene %in% genes_gw$genes) %>%
  write_tsv(paste0(OUT_DIR, 'Tucci_2019_genes_not_in_gw_list.txt'), col_names = FALSE)

## Logging  --------------------------------------------------------------------------
sink(LOG_FILE)

cat("\n-------------------------------------------------------------------------\n")
cat("Cross reference Tucci 2019 genes with genome wide genes", append = TRUE)
cat("\n-------------------------------------------------------------------------\n", append = TRUE)

cat("\nNumber of Tucci genes:", nrow(tucci_genes), append = TRUE)
cat("\nNumber of genome wide genes genes:", nrow(genes_gw), append = TRUE)
cat("\nNumber of Tucci genes in genomewide list:", nrow(tucci_genes_in_gw), append = TRUE)
cat("\nNumber of Tucci genes NOT in genomewide list:", nrow(tucci_genes_not_in_gw), append = TRUE)

cat("\n\nTucci genes in genome wide list:\n\n", tucci_genes_in_gw %>% pull(Gene), append = TRUE)
cat("\n\nTucci genes in genome wide list:\n\n", tucci_genes_not_in_gw %>% pull(Gene), append = TRUE)

cat("\n-------------------------------------------------------------------------\n", append = TRUE)

sink()

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
