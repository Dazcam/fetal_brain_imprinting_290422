#--------------------------------------------------------------------------------------
#
#    ASE/imprinting mao variants with MAF >= 0.05 to genes 
#
#--------------------------------------------------------------------------------------

## Method  ----------------------------------------------------------------------------

# 1. Extract all snps with MAF >= 0.05
# 2. Extract all gene coords from gff3 and convert to bed (using bed ops)     
#    - Extracting genes labelled gene | ncRNA_gene in column 3 
#    - Retaining if gene has HGNC gene name
#    - If not storing genes that only have Ensemble gene names in different file
# 3. Convert SNPs into bed format (using bed ops)
# 4. Loop through genes and pull out all SNPs in edited gff3 file with MAF >= 0.05
# 5. Write summary report per chromosome


## Set variables  ---------------------------------------------------------------------
CHR=$1
OUT_DIR=$2

REF_DIR=$3

mkdir -p ${OUT_DIR}genes/

printf '%s\n' '' '|  ----------------------------------------------  |'
printf '%s\n' '|  MAPPING SNPS (MAF >= 0.05 TO GENES GENOME WIDE  |' ''
printf '%s\n' '' '|  ----------------------------------------------  |'

# Extract all snps with MAF >= 0.05
printf '%s\n' '' 'Extracting all snps with MAF >= 0.05 (using bcftools) ...' 
printf '%s\n' "In file: ${REF_DIR}homo_sapiens-chr${CHR}.vcf.gz ..." 
printf '%s\n' "Out file: ${OUT_DIR}homo_sapiens-chr${CHR}_MAF.05.vcf.gz ..." ""

#these are loaded as modules by snakemake
#module load htslib/1.9
#module load bcftools/1.16.0

bcftools view -i 'MAF[0]>=0.05' ${REF_DIR}homo_sapiens-chr${CHR}.vcf.gz -o ${REF_DIR}homo_sapiens-chr${CHR}_MAF.05.vcf.gz

# Extract all gene coords from gff3 and convert to bed (using bed ops)
# Pull out gene and ncRNA_gene in col3
printf '%s\n' '' 'Extracting all gene coords from gff3 and convert to bed ...' ''
zcat ${REF_DIR}Homo_sapiens.GRCh38.108.chromosome.${CHR}.gff3.gz |\
awk '$3 == "gene" || $3 == "ncRNA_gene"' |\
grep Name |\
gff2bed | awk -F'[;]' '{print $1 "\t" $2}' |\
cut -f 1-3,11 |\
sed 's/Name=//g' |\
sed 's/^/chr/g' > ${OUT_DIR}chr${CHR}_genes.txt

# Collect gene coords with no HGNC ID from gff3 and convert to bed (using bed ops)
printf '%s\n' '' 'Collect gene coords with no HGNC ID from gff3 and convert to bed ...' ''
zcat ${REF_DIR}Homo_sapiens.GRCh38.108.chromosome.${CHR}.gff3.gz |\
awk '$3 == "gene" || $3 == "ncRNA_gene"' |\
grep -v Name |\
gff2bed | cut -f 1-4 |\
sed 's/Name=//g' > ${OUT_DIR}chr${CHR}_genes_noHGNC.txt

# Convert SNPs into bed format (using bed ops)
printf '%s\n' '' 'Convert SNPs into bed format ...' ''
zcat ${REF_DIR}homo_sapiens-chr${CHR}_MAF.05.vcf.gz |\
vcf2bed | cut -f 1-4 |\
sed 's/^/chr/g' > ${OUT_DIR}chr${CHR}_snps_MAF.05.txt

# Loop through genes and pull out all SNPs MAF >= 0.05
printf '%s\n' '' "Loop through chr${CHR} genes and pull out all SNPs MAF >= 0.05 ..." ''
while read -r line; 
do 
   echo $line > ${OUT_DIR}temp_${CHR}
   GENE=$(echo $line | cut -d' ' -f4)
   echo $GENE
   bedops --element-of ${OUT_DIR}chr${CHR}_snps_MAF.05.txt ${OUT_DIR}temp_${CHR} | uniq > ${OUT_DIR}genes/${GENE}_snps;
done < ${OUT_DIR}chr${CHR}_genes.txt

rm ${OUT_DIR}temp_${CHR}

# Reporting
printf '%s\n' '' 'Prepping report files and info  ...' ''
VARIANTS=$(zcat ${REF_DIR}homo_sapiens-chr${CHR}.vcf.gz | wc -l)
VARIANTS_05=$(zcat ${REF_DIR}homo_sapiens-chr${CHR}_MAF.05.vcf.gz | wc -l)
GENES=$(zcat ${REF_DIR}Homo_sapiens.GRCh38.108.chromosome.${CHR}.gff3.gz | awk '$3 == "gene" || $3 == "ncRNA_gene"' | wc -l)
GENES_HGNC=$(zcat ${REF_DIR}Homo_sapiens.GRCh38.108.chromosome.${CHR}.gff3.gz | awk '$3 == "gene" || $3 == "ncRNA_gene"' | grep Name | wc -l)
GENES_NO_HGNC=$(zcat ${REF_DIR}Homo_sapiens.GRCh38.108.chromosome.${CHR}.gff3.gz | awk '$3 == "gene" || $3 == "ncRNA_gene"' | grep -v Name | wc -l)
GENES_PROCESSED=$(ls ${OUT_DIR}genes/ | wc -l)
GENES_NO_SNPS=$(wc -l ${OUT_DIR}genes/* | grep -w 0 | wc -l)

wc -l ${OUT_DIR}genes/* > ${OUT_DIR}chr${CHR}_snp_cnts_per_gene.txt
wc -l ${OUT_DIR}genes/* | grep -w 0 | cut -d/ -f5 | cut -d_ -f1 > ${OUT_DIR}chr${CHR}_genes_no_snps.txt

printf '%s\n' '' "----------- Chr${CHR} Report -----------" ''
printf '%s\n' "Total variants in VCF file: $VARIANTS" 
printf '%s\n' "Total variants in VCF file, MAF >= 0.05: $VARIANTS_05" 
printf '%s\n' "Total genes in gff3 file: $GENES" 
printf '%s\n' "Total HGNC genes: $GENES_HGNC" 
printf '%s\n' "Total genes with no HGNC ID: $GENES_NO_HGNC" 
printf '%s\n' "Total genes in output directory: $GENES_PROCESSED" 
printf '%s\n' "Total genes in output directory containing no variants, MAF >= 0.05: $GENES_NO_SNPS" ''
printf '%s\n' "------------------------------------" ''


printf '%s\n' '' '|  ----------------------------------------------  |'
printf '%s\n' '|  ---------------  JOB FINISHED  ---------------  |'
printf '%s\n' '|  ----------------------------------------------  |'
