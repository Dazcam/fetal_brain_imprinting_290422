## Migration from XXXXXX - 

### Data processing information

[fastq_prep.smk](workflow/rules/fastq_prep.smk)

+ move_and_rename_fastqs - copy fastq files from neurocluster databank
+ zip_fastqs - standardise format of fastq files (some were zipped and some were not)
+ merge_fastqs - merge fastqs into single files for R1 and R2 (some samples were sequenced over multiple lanes)
+ fastqc_pretrim - run QC on fastq files
+ multiQC_pretrim - collate fastqc reports
+ trim_fastq - trim adatpers from fastqs and run fastqc on trimmed fastqs
+ hard_trim_fastq - Trim fastq reads that are longer that specified threshold
+ remove_short_reads - Trim fastq reads that are shorter than specified threshold
+ read_length_dist_post_QC_and_trimGalore - Assess read lengths of reads in fastq files before and after adapter trimming
+ get_read_length_dist_post_hard_and_SrtRead_trim - Assess read lengths of reads in fastq files after long and short read trimming

[allele_specific.smk](workflow/rules/allele_specific.smk) 

+ build_static_index - build static index file (required for ASElux)
+ extract_sample_genotypes - extract genotype infromation for each individual donor at all SNPs
+ ase_align - run allele specific expression analysis


NOTE: In former version of this workflow I generated json files to deal with the complex merge_fastq process
which had multiple different number of input files for each run of merge fastq. The generation of these files
occured on the output of the zip_fastqs rule. However, at the time I couldn't work out a way to integrate the 
json generation step into the snakemake workflow. For now the json scripts and data are stored in resources/json/.  


ISSUES: Discrepency between zip_fastqs and merge_fastqs rule caused by absolute path being in json files and 
relative path being in snake rules. Jobs appear to work now. 

***

### **BiomaRt method**

For the original list of 228 imprinted genes:

1. Get chromosome, start/stop coordinates for all genes (189 genes left - note some genes had multiple entries after this stage)
2. Get all rsIDs with MAF >= 0.05 for each gene - (169 genes left)
    + Only 180 genes processed as some genes were too large for software - I can get the SNPs for these manually
    + A further 11 genes had no SNPs within that range 
3. For SNPs MAF >= 0.05 filter those with ASE reads in >= 1 sample - (164 genes left - this step needed to make programming simpler)
4. For SNPs MAF >= 0.05 filter those with >= 20 ASE reads in >= 10 samples - (103 genes left - note thats 20 reads total across alleles)

Is there a more efficient way to do this? Partiularly for the genome wide analysis??

### Running the pipeline to deal with scratch quota limits

The pipeline has to be run in different blocks in oredr to balance the competing requiremnets
of scratch quota limit which are breeched early on and a huge number of small jobs later on.

During the `fastq_prep.smk`, multiple versions of the fastq files can be generated at each stage. This results in
the quota limit being breeched and pipeline choking at random points. To resolve this I used the following:

+ The snakemake `temp()` function to delete olde versions of fastqs as we moved through the pipeline
+ Lowered the max number of jobs that the pipeline could run on Hawk at any one time from 500 to 50
+ Gave later jobs in the process higher priority than earlier jobs


At the later `annotation.smk` step for ASE and SNP cross referencing step, for the genome wide anaysis we need to run 19K x 120 jobs,
these take ~30-60s and require negligable resources. It is important to reinstate the job limit to 500 jobs to churn through these
jobs quicker.     


### **Ensembl data dump method**

+ [GFF3 README](https://ftp.ensembl.org/pub/release-108/gff3/homo_sapiens/README)
+ [GFF3 description](http://gmod.org/wiki/GFF3)


1. **move_and_rename_fastqs**
- Standarise fastq names
- Move 526 fastq files sequenced across 3 sites Cardiff, Edinburgh, Exeter
- Due to differing naming conventions used json to populate sample names
2. **zip_fastqs**
- Some files were fastqs zipped and some were not. Zipped all.
3. **merge_fastqs**
- Merged fastqs that were sequenced over multiple lanes into one file (output: R1 and R2 file for each sample)
4. **fastqc_pretrim**
- Initial fastq QC - check for adapters
5. **multiQC_pretrim**
- Pool fastqc reports
6. **trim_fastq**
- Remove sequencing adapters from fastq files
7. **hard_trim_fastq**
- For ASElux reads need to be uniform length
- Trim RHS of reads to following length: Cardiff, 55; Edinburgh, 105; Exeter, 55
8. **remove_short_reads**
- Remove reads which hard a shorter length than the lengths specified above
9. **read_length_dist_post_QC_and_trimGalore**
- QC to confirm read lengths distribution pre-hard trim
10. **get_read_length_dist_post_hard_and_SrtRead_trim**
- QC to	confirm	read lengths distribution post-hard trim
11. **build_static_index**
- Build a static index reference file for ASElux
12. **extract_sample_genotypes**
- Extract genotypes for 120 samples
13. **ase_align**
- Run ASElux to measure allele specific expression in 120 samples
14. **ase_get_SNP_refs**
- Download Ensmebl SNP database file
- `https://ftp.ensembl.org/pub/release-108/variation/vcf/homo_sapiens/homo_sapiens-chr{CHR}.vcf.gz`
15. **ase_get_GENE_refs**
- Download Ensmebl GENE database file
- `https://ftp.ensembl.org/pub/release-108/gff3/homo_sapiens/Homo_sapiens.GRCh38.108.chromosome.{CHR}.gff3.gz`
16. **ase_map_snps_to_genes_genomewide**
- Map SNPs with MAF >= 0.05 to genes in genomewide 
- Note: gene labels used in column 3 of gff file were `gene` and `ncRNA_gene`
17. **ase_get_imprinted_gene_list**
- Generate gene list for is gene imprinted step
18. **ase_cross_ref**
- Bottleneck: Cross ref ASE variants for 120 samples with 30K MAF >= 0.05 varaints for ~30K genes
19. **ase_is_gene_imprinted**
- Check if gene is consistent with genomic imprinting:
- At least 90% of reads map to one of the two alleles in 80% of our heterozygotes for each SNP

