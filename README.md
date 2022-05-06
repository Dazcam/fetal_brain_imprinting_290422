## Migration from XXXXXX - 

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
