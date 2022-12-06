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
IMP_GENES <- snakemake@input$genes
LOG_FILE <- snakemake@output
OUT_DIR <- snakemake@params

##  Load Data  ------------------------------------------------------------------------
genes <- read_excel(IMP_GENES, sheet = 'Human_All_Known IG_NoStatus') %>%
  filter(!grepl("placenta", Gene)) %>%
  drop_na(Gene) %>%
  arrange(Gene) %>%
  select(Gene) %>%
  write_tsv(paste0(OUT_DIR, 'Tucci_2019_imprinted_genes.txt'), col_names = FALSE)
  
#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
