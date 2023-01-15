#--------------------------------------------------------------------------------------
#
#    ASE/imprinting gene SNP cross-referencing
#
#--------------------------------------------------------------------------------------

## Method  ----------------------------------------------------------------------------

# 1. Cross ref ASE results with MAF >= 0.05 SNPs for each gene in Tucci gene list
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
FILE=$1
INDIR=$2
OUTDIR=$3

## Create summary files  ---------------------------------------------------------------
printf "\n|--------------------------------------------|\n"
printf "  Creating summary files ... \n"
printf "|--------------------------------------------|\n\n"

# Count total genes processed
echo "Total Genes processed:" $'\t' `ls -l ${OUTDIR} | grep '^d' | wc -l` > ${OUTDIR}summary_1_genes_processed.txt

# Count SNPs with MAF >= 0.05 in each gene?
cat ${FILE} | while read GENE; do

  wc -l ${INDIR}${GENE}_snps | awk -v OFS=' ' '{print $2 "\t" $1}' | cut -d/ -f 5 | sed 's/_snps//g' >> ${OUTDIR}summary_2_SNPS_MAF_0.05.txt

done

# Count SNPs with MAF >= 0.05 have gEx values in at least 1 of the 120 samples
cat ${FILE} | while read GENE; do

  echo $GENE
  echo $GENE >> ${OUTDIR}summary_3_gene_list.txt
  find ${OUTDIR}${GENE} -type f | xargs wc -l | sed '$ d' | awk '$1 > 7' | wc -l >> ${OUTDIR}summary_4_MAF_SNPs_exp_in_samples.txt

done

# Which variants have > 20 reads in any sample?

printf "\n|-------------------------------------------------|\n"
printf "  Does each SNP have >= 20 reads in any sample?\n"
printf "|-------------------------------------------------|\n\n"

printf "The following varinats have >= 20 reads in at least 1 sample:\n\n"


cat ${FILE} | while read GENE; do

  for SNP in `find ${OUTDIR}${GENE} -type f | xargs wc -l | awk '$1 > 7' | sed '$ d' | awk -v OFS=' ' '{print $2}'`; do
  
  SNPnoExt=${SNP%_*}
  echo ${SNPnoExt} | cut -d/ -f4,7 | sed 's/\// - /g'
  
    awk -F "\t" '{ 
    
    if($2+$3 >= 20)
      { print $0 } 
    
    }' ${SNP} | cut -d/ -f 4 > ${SNPnoExt}_20plus_reads_in_sample
    
    
  done
 
done


# Of those that have SNP > 20 reads, do at least 10 samples have 20 reads over this SNP?

printf "\n|-----------------------------------------------------------------------------------------|\n"
printf "  Of those that have SNP >= 20 reads, do at least 10 samples have 20 reads over this SNP?\n"
printf "|-----------------------------------------------------------------------------------------|\n\n"

cat ${FILE} | while read GENE; do

 # Total SNPs with 20 reads in 10 samples
 find ${OUTDIR}${GENE}/*_20plus_reads_in_sample -type f | xargs wc -l | awk '$1 >= 10' |\
  sed '$ d' | wc -l >> ${OUTDIR}summary_5_SNPnum_20readsIn10samples.txt

 # Get list of rsIDs for SNPs with 20 reads in 10 samples
 find ${OUTDIR}${GENE}/*_20plus_reads_in_sample -type f | xargs wc -l | awk '$1 >= 10' |\
 sed '$ d' | awk -v OFS=' ' '{print $2}' | cut -d'/' -f7 | sed 's/_20plus_reads_in_sample//g' >> ${OUTDIR}${GENE}/rsIDs_20readsIn10samples

done


# Create general summary
paste ${OUTDIR}summary_2_SNPS_MAF_0.05.txt \
${OUTDIR}summary_4_MAF_SNPs_exp_in_samples.txt \
${OUTDIR}summary_5_SNPnum_20readsIn10samples.txt > ${OUTDIR}summary_6_ase_imprinting_genes_summary.txt

printf "Done."

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
