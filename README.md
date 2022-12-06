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
