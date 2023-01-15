#--------------------------------------------------------------------------------------
#
#    Determine if gene are consistent with genomic imprinting
#
#--------------------------------------------------------------------------------------

## Method  ----------------------------------------------------------------------------
#  Consistent with genomic imprinting is defined as:

#  At least 90% of reads map to one of the two alleles
#  in 80% of our heterozygotes for each SNP

## Set variables ----------------------------------------------------------------------
INDIR=$1
OUTDIR=$2

## Set variables ----------------------------------------------------------------------
for GENE in `awk -F"\t" '$4>0' ${INDIR}summary_6_ase_imprinting_genes_summary.txt | awk '{print $1}'`; do

  echo "======================================================"
  printf "${GENE}\n"
  echo "======================================================"
  
  mkdir ${OUTDIR}${GENE}
  while read SNP; do

  SAMPLES=`grep -c :rs ${INDIR}${GENE}/${SNP}_ase`

  echo "------------------------------------------------------"
  printf "\nChecking if ${SNP} in ${GENE} is Imprinted\n"
  printf "\nTotal individuals expressing ${SNP} is: ${SAMPLES}\n\n"


  # 1st check if 90% of reads map to REF or ALT - use col 4 - REF expression ratio
  # Need awk due to floating pts

  IMP_UPPER=0.90
  IMP_LOWER=0.10
  PARTIAL_UPPER=0.6700000
  PARTIAL_LOWER=0.3300000

  grep :rs ${INDIR}${GENE}/${SNP}_ase |\
  awk -v imp_upper="$IMP_UPPER" -v imp_lower="$IMP_LOWER" -v snp="$SNP" -v gene="$GENE" \
  -v part_upper="$PARTIAL_UPPER" -v part_lower="$PARTIAL_LOWER" '{

  if ($4 >= imp_upper || $4 <= imp_lower)
          print $1"\t"$4"\tImprinted" > "/localscratch/scw1641/results/16ISGENEIMPRINTED/"gene"/"snp"_imp_test";
  else if ($4 >= part_upper || $4 <= part_lower)
      print $1"\t"$4"\tPartially imprinted" > "/localscratch/scw1641/results/16ISGENEIMPRINTED/"gene"/"snp"_imp_test";
  else
          print $1"\t"$4"\tBiallelic" > "/localscratch/scw1641/results/16ISGENEIMPRINTED/"gene"/"snp"_imp_test";

  }'

  # Count number of individals SNP is imprinted in then
  # divide by total individuals to get prop of hets where SNP is imprinted
  IMPRINTED=`grep -c Imprinted ${OUTDIR}${GENE}/${SNP}_imp_test`
  PROP_IMPRINTED=$(echo "${IMPRINTED}/${SAMPLES}" | bc -l) # Feed calc to bc as floating pt

  PART_IMP=`grep -c Partially ${OUTDIR}${GENE}/${SNP}_imp_test`
  SUM_IMP_PART_IMP=$(echo "${PART_IMP}+${IMPRINTED}" | bc -l)
  PROP_PART_IMP=$(echo "${SUM_IMP_PART_IMP}/${SAMPLES}" | bc -l)

  # Report to stout SNP stats
  printf "Total individuals where ${SNP} has enough reads to be consistent with genomic imprinting is: ${IMPRINTED}\n"
  printf '\nThe proportion of hets where SNP is imprinted is: %.3f\n' $(echo "${PROP_IMPRINTED}")

  printf "\nTotal individuals where ${SNP} has enough reads to be consistent with partial genomic imprinting is: $(echo "${SUM_IMP_PART_IMP}")\n"
  printf '\nThe proportion of hets where SNP is partially imprinted is: %.3f\n' $(echo "${PROP_PART_IMP}") # Feed calc to bc as floating


  # Report to file whether SNPs are imprinted, partilaly imprinted, or not imprinted

  if (( $(echo "${PROP_IMPRINTED} >= 0.8" |bc -l) )); then

    printf ${GENE}"\t"${SNP}"\t"${SAMPLES}"\t"${IMPRINTED}"\t"${PROP_IMPRINTED}"\t"${SUM_IMP_PART_IMP}"\t"${PROP_PART_IMP}"\tImprinted\n" >> ${OUTDIR}imp_test.txt

  elif (( $(echo "${PROP_PART_IMP} >= 0.8" |bc -l) )); then

    printf "${GENE}\t${SNP}\t${SAMPLES}\t${IMPRINTED}\t${PROP_IMPRINTED}\t${SUM_IMP_PART_IMP}\t${PROP_PART_IMP}\tPartially Imprinted\n" >> ${OUTDIR}imp_test.txt

  else

    printf ${GENE}"\t"${SNP}"\t"${SAMPLES}"\tNA\tNA\tNA\tNA\tNot Imprinted\n" >> ${OUTDIR}imp_test.txt

  fi

  printf "\n------------------------------------------------------\n\n"

  done < ${INDIR}${GENE}/rsIDs_20readsIn10samples

  printf "======================================================\n\n\n\n"

done

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------

