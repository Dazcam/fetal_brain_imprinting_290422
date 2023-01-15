#--------------------------------------------------------------------------------------
#
#    Imprinting - Create summary file
#
#--------------------------------------------------------------------------------------

## Requirements  ----------------------------------------------------------------------
# module load libgit2/1.1.0

## Info  ------------------------------------------------------------------------------

#  Run the analysis up until cluster QC cell removal 

## Initialise R library  --------------------------------------------------------------
.libPaths( c( "/scratch/c.c1477909/R/library", .libPaths() ) )
library(tidyverse)
library(openxlsx)

## Set variables  ---------------------------------------------------------------------
SNPS2GENE_DIR <- 'results/00LOG/13SNP2GENES_GW/'
CROSSREF_DIR <- 'results/15CROSSREF_TUCCI/'
ISGENEIMPDIR <- 'results/16ISGENEIMPRINTED_TUCCI/'
OUTFILE <- 'imprinting_summary.xlsx'


## Map variants to genes summary  -----------------------------------------------------
for (CHR in seq(1,22)) {
  
  FILE <- paste0(SNPS2GENE_DIR, 'snps2genes_gw_chr', CHR,'.log')
  CHR_SUMMARY <- readLines(FILE) 
  SUMMARY <- as.data.frame(CHR_SUMMARY[grep('^Total*', CHR_SUMMARY)]) %>%
    rename(Test = names(.)) %>%
    separate(Test, c("Discription", paste0('chr', CHR)), sep = ": ") %>%
    head(5)
  
  assign(paste0('chr', CHR, '_df'), SUMMARY)
  
}

snp2gene_df <- plyr::join_all(list(chr1_df, chr2_df, chr3_df, chr4_df, chr5_df,
              chr6_df, chr7_df, chr8_df, chr9_df, chr10_df,
              chr11_df, chr12_df, chr13_df, chr14_df, chr15_df,
              chr16_df, chr17_df, chr18_df, chr19_df, chr20_df,
              chr21_df, chr22_df), by = 'Discription', type = 'left') %>% 
  mutate_at(vars(contains('chr')), as.numeric) %>%
  mutate(Total = rowSums(select(., contains("chr"))))

## Crossref summary  -------------------------------------------------------------------
cross_ref_df <- read_tsv(paste0(CROSSREF_DIR, 'summary_6_tucci_ase_imprinting_genes_summary.txt'),
                         col_names = FALSE) 
colnames(cross_ref_df) <- c('Gene', 'Vars MAF >= 0.05', 'Vars MAF >= 0.05 exp in samples', 
                            'Vars MAF >= 0.05, 20 reads In 10 samples')

## Is gene imprinted summary  ----------------------------------------------------------
imp_df <- read_tsv(paste0(ISGENEIMPDIR, 'imp_test.txt'),
                         col_names = FALSE) 
colnames(imp_df) <- c('Gene', 'Variant', 'No. Samples', 'No. Samples Imp', 'Prop. Samples Imp', 
                      'Sum Samples Imp Prop. Imp', 'Prop. Samples Part Imp', 'Imp status')

df_list <- list(snp2gene_df, cross_ref_df, imp_df)

write.xlsx(df_list, OUTFILE, sheetName = c("map_snp_2_gene", "cross_ref", 
                                           "is_gene_imprinted"))

#-------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------
