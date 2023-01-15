#--------------------------------------------------------------------------------------
#
#    ASE/imprinting gene SNP cross-referencing
#
#--------------------------------------------------------------------------------------

## Method  ----------------------------------------------------------------------------

# 1. Cross ref ASE results with MAF >= 0.05 SNPs for each gene in gene list
# 2. Summary files:
#    - Count total genes processed
#    - Count SNPs with MAF >= 0.05 in each gene?
#    - Count SNPs with MAF >= 0.05 have gEx values in at least 1 of the 120 samples
#      and have > 20 reads over SNP in >= 10 samples

# Final summary output

#  Gene
#  Variants with MAF >= 0.05
#  Variants with MAF >= 0.05 have gEx values in at least 1 of the 120 samples
#  No. of variants with >= 20 reads in 10 samples


## Set variables ----------------------------------------------------------------------
GENE=$1
INDIR=$2
OUTDIR=$3


## Get ase values over imprinted gene SNPs  -------------------------------------------

printf "\n|------------------------------------------------|\n"
printf "  Obtaining ase values over imprinted gene SNPs ... \n"
printf "|------------------------------------------------|\n\n"

printf "ASE results for: ${GENE}\n\n"

mkdir -p ${OUTDIR}${GENE}

for rsID in `cut -f4 ${INDIR}${GENE}_snps`; do

    printf "${GENE} - ${rsID}\n"

    printf "\n|-------------------------------------------|\n\n" >> ${OUTDIR}${GENE}/${rsID}_ase
    echo "ASE results for: "${GENE} - ${rsID} >> ${OUTDIR}${GENE}/${rsID}_ase
    echo"" >> ${OUTDIR}${GENE}/${rsID}_ase
    grep -wE ${rsID} ../results/11ASE_MAPPINGS/* >> ${OUTDIR}${GENE}/${rsID}_ase
    printf "\n|-------------------------------------------|\n" >> ${OUTDIR}${GENE}/${rsID}_ase

done

touch ${OUTDIR}${GENE}/touch.txt

printf "\nASE results for: ${GENE}  - DONE\n\n"

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
