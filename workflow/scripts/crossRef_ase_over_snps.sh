#############################################################################
#
#    ASE/imprinting gene SNP cross-referencing
#
#############################################################################


#####  Get ase values over imprinted gene SNPs  #####

for GENE in `ls snps/`
do
  
  PREFIX="$(basename ${GENE} _snps)" 
  echo "ASE results for: "${PREFIX}
  
  mkdir snps/${PREFIX}
  
  for rsID in `cut -f1 snps/${GENE}`
  
  do
   
   echo "ASE results for: "${PREFIX} - ${rsID}
   
   printf "\n|-------------------------------------------|\n\n" >> snps/${PREFIX}/${rsID}_ase
   echo "ASE results for: "${PREFIX} - ${rsID} >> snps/${PREFIX}/${rsID}_ase
   echo"" >> snps/${PREFIX}/${rsID}_ase
   grep -wE ${rsID} 09ASE_MAPPINGS/* >> snps/${PREFIX}/${rsID}_ase
   printf "\n|-------------------------------------------|\n" >> snps/${PREFIX}/${rsID}_ase

  done
  
  echo "ASE results for: "${PREFIX} " - DONE"

done 


#####  Create summary file  #####

# How many genes in total were processed?
echo "Total Genes processed" $'\t' `ls snps/*_snps | wc -l` > total_genes

# How many SNPs with MAF >= 0.05 are there within each gene?
wc -l snps/*_snps | awk -v OFS=' ' '{print $2 "\t" $1}' | sed '$ d' | sed -E 's/_snps//g' |\
sed -E 's/snps\///g' >> total_snps_maf_gt.05

# How many SNPs with MAF >= 0.05 have gene expression values (GEV) in at least 1 of the 120 samples?

for GENE in `ls snps/*_snps`
do
  
  PREFIX="$(basename ${GENE} _snps)"
  echo $PREFIX
  echo $PREFIX >> gene_list
  find snps/${PREFIX} -type f | xargs wc -l | sed '$ d' | awk '$1 > 7' | wc -l >> SNPs_expressed
 
done


##### How many SNPs with MAF >= 0.05 have GEV in at least 1 of the 120 samples,  #####
##### and have > 20 reads over SNP in 10 or more samples?                        #####

# Does each SNP have > 20 reads in any sample?

for GENE in `ls snps/*_snps`; do 

  PREFIX="$(basename ${GENE} _snps)"

  for SNP in `find snps/${PREFIX} -type f | xargs wc -l | awk '$1 > 7' | sed '$ d' | awk -v OFS=' ' '{print $2}'`; do
  
  SNPnoExt=${SNP%%_*}
  echo ${SNPnoExt}
  
    awk -F "\t" '{ 
    
    if($2+$3 >= 20)
      { print $0 } 
    
    }' ${SNP} > ${SNPnoExt}_gt20reads 
    
    
  done
 
done


# Of those that have SNP > 20 reads, do at least 10 samples have 20 reads over this SNP?

for GENE in `ls snps/*_snps`; do 

  PREFIX="$(basename ${GENE} _snps)"
  

      # Total SNPs with 20 reads in 10 samples
      find snps/${PREFIX}/*_gt20reads -type f | xargs wc -l | awk '$1 >= 10' | sed '$ d' | wc -l >> SNPnum_20readsIn10samples
  
      # Get list of rsIDs for SNPs with 20 reads in 10 samples
      find snps/${PREFIX}/*_gt20reads -type f | xargs wc -l | awk '$1 >= 10' | sed '$ d' |\
      awk -v OFS=' ' '{print $2}' | cut -d'/' -f3 | sed 's/_gt20reads//g' >> snps/${PREFIX}/rsIDs_20readsIn10samples

done

paste total_snps_maf_gt.05 snps_expressed SNPnum_20readsIn10samples > ase_imprinting_genes_summary

# cut -f2 total_snps_maf_gt.05 | sort |  uniq -c
# sort snps_expressed | sort |  uniq -c
# sort sSNPnum_20readsIn10samples | sort |  uniq -c



# Quick check of ase scores for SNPs of interest

# Manually 

GENE=KCNQ1OT1

while read SNP; do

cat snps/${GENE}/${SNP}_ase

done < snps/${GENE}/rsIDs_20readsIn10samples


# Systematically
while read GENE; do

  while read SNP; do

    cat snps/${GENE}/${SNP}_ase
    
  done < snps/${GENE}/rsIDs_20readsIn10samples

  sleep 10
  
done < gene_list
