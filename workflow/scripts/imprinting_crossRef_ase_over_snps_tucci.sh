#--------------------------------------------------------------------------------------
#
#    ASE/imprinting gene SNP cross-referencing - for Tucci 2019 gene list
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


## Get ase values over imprinted gene SNPs  -------------------------------------------

printf "\n|------------------------------------------------|\n"
printf "Geting ase values over imprinted gene SNPs ... \n"
printf "|------------------------------------------------|\n\n"

cat ${FILE} | while read GENE; do

  printf "ASE results for: ${GENE}\n\n"

  mkdir -p ${OUTDIR}${GENE}

  for rsID in `cut -f4 ${INDIR}${GENE}_snps`

  do

   printf "${GENE} - ${rsID}\n"

   printf "\n|-------------------------------------------|\n\n" >> ${OUTDIR}${GENE}/${rsID}_ase
   echo "ASE results for: "${GENE} - ${rsID} >> ${OUTDIR}${GENE}/${rsID}_ase
   echo"" >> ${OUTDIR}${GENE}/${rsID}_ase
   grep -wE ${rsID} ../results/11ASE_MAPPINGS/* >> ${OUTDIR}${GENE}/${rsID}_ase
   printf "\n|-------------------------------------------|\n" >> ${OUTDIR}${GENE}/${rsID}_ase

  done

  printf "\nASE results for: ${GENE}  - DONE\n\n"

done

## Create summary files  ---------------------------------------------------------------
printf "\n|--------------------------------------------|\n"
printf "Creating summary files ... \n"
printf "|--------------------------------------------|\n\n"

# Count total genes processed
echo "Total Genes processed:" $'\t' `ls -l ${OUTDIR} | grep '^d' | wc -l` > ${OUTDIR}summary_1_tucci_genes_processed.txt

# Count SNPs with MAF >= 0.05 in each gene?
cat ${FILE} | while read GENE; do

  wc -l ${INDIR}${GENE}_snps | awk -v OFS=' ' '{print $2 "\t" $1}' | cut -d/ -f 5 | sed 's/_snps//g' >> ${OUTDIR}summary_2_tucci_SNPS_MAF_0.05.txt

done

# Count SNPs with MAF >= 0.05 have gEx values in at least 1 of the 120 samples
cat ${FILE} | while read GENE; do

  echo $GENE
  echo $GENE >> ${OUTDIR}summary_3_tucci_gene_list.txt
  find ${OUTDIR}${GENE} -type f | xargs wc -l | sed '$ d' | awk '$1 > 7' | wc -l >> ${OUTDIR}summary_4_tucci_MAF_SNPs_exp_in_samples.txt

done

# Does each SNP have > 20 reads in any sample?

printf "\n|-------------------------------------------------|\n"
printf "Does each SNP have > 20 reads in any sample?\n"
printf "|-------------------------------------------------|\n\n"


cat ${FILE} | while read GENE; do

  for SNP in `find ${OUTDIR}${GENE} -type f | xargs wc -l | awk '$1 > 7' | sed '$ d' | awk -v OFS=' ' '{print $2}'`; do
  
  SNPnoExt=${SNP%_*}
  echo ${SNPnoExt}
  
    awk -F "\t" '{ 
    
    if($2+$3 >= 20)
      { print $0 } 
    
    }' ${SNP} > ${SNPnoExt}_20plus_reads_in_sample
    
    
  done
 
done


# Of those that have SNP > 20 reads, do at least 10 samples have 20 reads over this SNP?

printf "\n|--------------------------------------------------------------------------------------|\n"
printf "Of those that have SNP > 20 reads, do at least 10 samples have 20 reads over this SNP?\n"
printf "|--------------------------------------------------------------------------------------|\n\n"

cat ${FILE} | while read GENE; do

 # Total SNPs with 20 reads in 10 samples
 find ${OUTDIR}${GENE}/*_20plus_reads_in_sample -type f | xargs wc -l | awk '$1 >= 10' |\
  sed '$ d' | wc -l >> ${OUTDIR}summary_5_tucci_SNPnum_20readsIn10samples.txt

 # Get list of rsIDs for SNPs with 20 reads in 10 samples
 find ${OUTDIR}${GENE}/*_20plus_reads_in_sample -type f | xargs wc -l | awk '$1 >= 10' |\
 sed '$ d' | awk -v OFS=' ' '{print $2}' | cut -d'/' -f5 | sed 's/_20plus_reads_in_sample//g' >> ${OUTDIR}${GENE}/rsIDs_20readsIn10samples

done


# Create general summary
paste ${OUTDIR}summary_2_tucci_SNPS_MAF_0.05.txt \
${OUTDIR}summary_4_tucci_MAF_SNPs_exp_in_samples.txt \
${OUTDIR}summary_5_tucci_SNPnum_20readsIn10samples.txt > ${OUTDIR}summary_6_tucci_ase_imprinting_genes_summary.txt

printf "Done."

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------

