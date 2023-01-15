#--------------------------------------------------------------------------------------
#
#    Imprinting - Fishers Exact testing for gene list overlap
#
#--------------------------------------------------------------------------------------

## Info  ------------------------------------------------------------------------------

#  Run fishers exact tests 

## Set log  ---------------------------------------------------------------------------
#log <- file(snakemake@log[[1]], open="wt")
#sink(log)

## Initialise R library  --------------------------------------------------------------
.libPaths( c( "/scratch/c.c1477909/R/library", .libPaths() ) )

##  Load Packages  --------------------------------------------------------------------
library(tidyverse)
library(GeneOverlap)
library(readxl)
library(rmarkdown)

## Set variables  ---------------------------------------------------------------------
cat('\nSetting variables ... \n')
NOVEL_GENES <- toString(snakemake@input$NOVEL_GENES)
OVERLAP_GENES <- toString(snakemake@input$OVERLAP_GENES)
TUCCI_GENES <- toString(snakemake@input$TUCCI_GENES)
GW_GENES <- toString(snakemake@input$GW_GENES)
NDD_GENES <- toString(snakemake@input$NDD_GENES)
SCZ_CV_GENES <- toString(snakemake@input$SCZ_CV_GENES)
SCZ_RV_GENES <- toString(snakemake@input$SCZ_RV_GENES)
ASD_GENES <- toString(snakemake@input$ASD_GENES)
DDG2p_GENES <- toString(snakemake@input$DDG2p_GENES)
MARKDOWN_FILE <- toString(snakemake@input$MARKDOWN_FILE) 
REPORT_DIR <- toString(snakemake@params$REPORT_DIR)
REPORT_FILE <- toString(snakemake@params$REPORT_FILE)

##  Load gene sets from this study   --------------------------------------------------
cat('\nLoading gene sets from this study ... \n')
novel_genes <- read_tsv(NOVEL_GENES, col_names = FALSE) %>% pull()
overlap_genes <- read_tsv(OVERLAP_GENES, col_names = FALSE) %>% pull()
tucci_only_genes <- read_tsv(TUCCI_GENES, col_names = FALSE) %>% pull()
all_genes <- read_tsv(GW_GENES, col_names = FALSE) %>% pull()

##  Load gene sets from public studies   ----------------------------------------------
cat('\nLoading public gene sets ... \n')
ndd_genes <- read_tsv(NDD_GENES) %>%
  pull(symbol)
scz_CV_genes <- read_excel(SCZ_CV_GENES, sheet = 'Prioritised') %>% pull(Symbol.ID)
scz_RV_genes <- read_tsv(SCZ_RV_GENES) %>%
  pull(symbol)
asd_genes <- read_tsv(ASD_GENES) %>%
  pull(gene)
ddg2p_genes <- read_tsv(DDG2p_GENES) %>%
  pull(gene.symbol)
comb_genes <- c(ndd_genes, asd_genes, ddg2p_genes)

##  Run Fisher's   --------------------------------------------------------------------
cat('\nRun Fishers ... \n')
imp_gene_list <- list(novel_genes = novel_genes, 
                      overlap_genes = overlap_genes, 
                      tucci_only_genes = tucci_only_genes)
test_gene_list <- list(ndd_genes = ndd_genes, 
                       scz_CV_genes = scz_CV_genes, 
                       scz_RV_genes = scz_RV_genes,
                       asd_genes = asd_genes, 
                       ddg2p_genes = ddg2p_genes)

gom.obj <- newGOM(imp_gene_list, test_gene_list, genome.size = length(all_genes))


##  Plot   ----------------------------------------------------------------------------
cat('\nPlotting ... \n')
drawHeatmap(gom.obj)
getMatrix(gom.obj, name = "pval")
getMatrix(gom.obj, "odds.ratio")

## Pull out intersections  ------------------------------------------------------------  
cat('\nPull out intersections ...\n')
inter.nl <- getNestedList(gom.obj, name = "intersection")
str(inter.nl)

## Create markdown doc  ---------------------------------------------------------------
cat('\nCreating markdown report ... \n')
render(MARKDOWN_FILE, output_file = REPORT_FILE, output_dir = REPORT_DIR)

cat('\nDONE.\n')

#sink()

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
