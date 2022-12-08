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


## Set variables ----------------------------------------------------------------------
OUTDIR=$1

## Get ase values over imprinted gene SNPs  -------------------------------------------
for SNPs_for_GENE in `ls ../results/12SNP2GENES/*_snps`
do
  
  PREFIX="$(basename ${SNPs_for_GENE} _snps)" 
  echo "ASE results for: "${PREFIX}
  
  mkdir -p ${OUTDIR}${PREFIX}
  
  for rsID in `cut -f1 ${SNPs_for_GENE}`
  
  do
   
   echo "ASE results for: "${PREFIX} - ${rsID}
   
   printf "\n|-------------------------------------------|\n\n" >> ${OUTDIR}${PREFIX}/${rsID}_ase
   echo "ASE results for: "${PREFIX} - ${rsID} >> ${OUTDIR}${PREFIX}/${rsID}_ase
   echo"" >> ${OUTDIR}${PREFIX}/${rsID}_ase
   grep -wE ${rsID} ../results/11ASE_MAPPINGS/* >> ${OUTDIR}${PREFIX}/${rsID}_ase
   printf "\n|-------------------------------------------|\n" >> ${OUTDIR}${PREFIX}/${rsID}_ase

  done
  
  echo "ASE results for: "${PREFIX} " - DONE\n"

done 


## Create summary files  --------------------------------------------------------------- 
# Count total genes processed
echo "Total Genes processed:" $'\t' `ls ../results/12SNP2GENES/*_snps | wc -l` > ${OUTDIR}summary_1_genes_processed.txt

# Count SNPs with MAF >= 0.05 in each gene?
wc -l ../results/12SNP2GENES/*_snps |\
awk -v OFS=' ' '{print $2 "\t" $1}' |\
sed '$ d' |\
cut -d/ -f 4 |\
sed 's/_snps//g' >> ${OUTDIR}summary_2_SNPS_MAF_0.05.txt

# Count SNPs with MAF >= 0.05 have gEx values in at least 1 of the 120 samples
for GENE in `ls ../results/12SNP2GENES/*_snps`
do
  
  PREFIX="$(basename ${GENE} _snps)"
  echo $PREFIX
  echo $PREFIX >> ${OUTDIR}summary_3_gene_list.txt
  find ${OUTDIR}${PREFIX} -type f | xargs wc -l | sed '$ d' | awk '$1 > 7' | wc -l >> ${OUTDIR}summary_4_MAF_SNPs_exp_in_samples.txt
 
done

# Does each SNP have > 20 reads in any sample?
for GENE in `ls ../results/12SNP2GENES/*_snps`; do 

  PREFIX="$(basename ${GENE} _snps)"

  for SNP in `find ${OUTDIR}${PREFIX} -type f | xargs wc -l | awk '$1 > 7' | sed '$ d' | awk -v OFS=' ' '{print $2}'`; do
  
  SNPnoExt=${SNP%_*}
  echo ${SNPnoExt}
  
    awk -F "\t" '{ 
    
    if($2+$3 >= 20)
      { print $0 } 
    
    }' ${SNP} > ${SNPnoExt}_20plus_reads_in_sample
    
    
  done
 
done


# Of those that have SNP > 20 reads, do at least 10 samples have 20 reads over this SNP?
for GENE in `ls ../results/12SNP2GENES/*_snps`; do 

  PREFIX="$(basename ${GENE} _snps)"
  

      # Total SNPs with 20 reads in 10 samples
      find ${OUTDIR}${PREFIX}/*_20plus_reads_in_sample -type f | xargs wc -l | awk '$1 >= 10' | sed '$ d' | wc -l >> ${OUTDIR}summary_5_SNPnum_20readsIn10samples.txt
  
      # Get list of rsIDs for SNPs with 20 reads in 10 samples
      find ${OUTDIR}${PREFIX}/*_20plus_reads_in_sample -type f | xargs wc -l | awk '$1 >= 10' | sed '$ d' |\
      awk -v OFS=' ' '{print $2}' | cut -d'/' -f3 | sed 's/_gt20reads//g' >> ${OUTDIR}${PREFIX}/rsIDs_20readsIn10samples

done

paste ${OUTDIR}summary_2_SNPS_MAF_0.05.txt \
${OUTDIR}summary_4_MAF_SNPs_exp_in_samples.txt \
${OUTDIR}summary_5_SNPnum_20readsIn10samples.txt > ${OUTDIR}summary_6_ase_imprinting_genes_summary.txt

# cut -f2 total_snps_maf_gt.05 | sort |  uniq -c
# sort snps_expressed | sort |  uniq -c
# sort sSNPnum_20readsIn10samples | sort |  uniq -c



# Quick check of ase scores for SNPs of interest

# Manually 

#GENE=KCNQ1OT1

#while read SNP; do

#cat snps/${GENE}/${SNP}_ase

#done < ${OUTDIR}${GENE}/rsIDs_20readsIn10samples


# Systematically
#while read GENE; do

#  while read SNP; do

#    cat snps/${GENE}/${SNP}_ase
    
#  done < ${OUTDIR}${GENE}/rsIDs_20readsIn10samples

#  sleep 10
  
#done < ${OUTDIR}summary_3_gene_list
