# A genome-wide screen of ASE - V2 - on going

This project was carried out in the Division of Psychological Medicine and Clinical Neurosciences (DPMCN). The paper is [here](https://www.biologicalpsychiatryjournal.com/article/S0006-3223(22)01404-4/fulltext). The workflow follows the the snakemake [distribution and reproducibility](https://snakemake.readthedocs.io/en/stable/snakefiles/deployment.html) recommendations. 

***

A snakemake pipeline to for a geneome wide screen of allele specific expression using genotype and GeX from 120 human brain bulk tissue samples. Utilising the following packages:

+ [Snakemake 6.6.1](https://snakemake.readthedocs.io/en/stable/)
+ [Fastqc 0.11.8](https://github.com/s-andrews/FastQC)
+ [Multiqc 1.7](https://multiqc.info)
+ [Trim Galore 0.6.10](https://github.com/FelixKrueger/TrimGalore) 
+ [BBMap 38.84](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/bbmap-guide/) 
+ [Bcftools 1.16](https://samtools.github.io/bcftools/bcftools.html)
+ [ASElux 1.0.2](https://github.com/abl0719/ASElux)

***

**Data**

+ [FastQ files](https://ega-archive.org/search-results.php?query=EGAS00001003214)

***

Papers for public data and software

+ [O'Brien et al. (2018)](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-018-1567-1#Sec24)
+ [Miao eat al. (2018)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5905663/) - ASElux
 
***

### **Ensembl data dump method**


+ [GFF3 README](https://ftp.ensembl.org/pub/release-108/gff3/homo_sapiens/README)
+ [GFF3 description](http://gmod.org/wiki/GFF3)

****

**Scripts**

[fastq_prep.smk](workflow/rules/fastq_prep.smk)

+ move_and_rename_fastqs
    + Move 526 fastq files (sequenced across 3 sites Cardiff, Edinburgh, Exeter) from neurocluster
    + Standarise fastq names
    + Due to differing naming conventions used json to populate sample names
+ zip_fastqs 
    + Some files were fastqs zipped and some were not. Zipped all.
+ merge_fastqs
    + Merged fastqs that were sequenced over multiple lanes into one file (output: R1 and R2 file for each sample)
+ fastqc_pretrim 
    + Initial fastq QC - check for adapters
+ multiQC_pretrim 
    + Collate fastqc reports
+ trim_fastq 
    + Remove sequencing adapters from reads and run fastqc on trimmed reads
+ hard_trim_fastq
    + For ASElux reads need to be uniform length
    + Trim RHS of reads to following length: Cardiff, 55; Edinburgh, 105; Exeter, 55
+ remove_short_reads
    + Remove reads which hard a shorter length than the lengths specified above
+ read_length_dist_post_QC_and_trimGalore 
    + QC to confirm read length distribution pre-hard trim
+ get_read_length_dist_post_hard_and_SrtRead_trim 
    + QC to	confirm	read lengths distribution post-hard trim

[allele_specific.smk](workflow/rules/allele_specific.smk) 

+ build_static_index 
    + Build a static index reference file for ASElux
+ extract_sample_genotypes - 
    + Extract genotype infromation for each individual donor at all SNPs
+ ase_align - run allele specific expression analysis
    + Run ASElux to measure allele specific expression in 120 samples
    
[imprinting_annotation_local.smk](workflow/rules/imprinting_annotation_local.smk)

+ ase_get_SNP_ref
    + Download Ensembl SNP database file
+ ase_get_GENE_ref
    + Download Ensembl GENE database file
+ ase_map_snps_to_genes_genomewide
    + Map SNPs with MAF >= 0.05 to genes genomewide 
    + Note: gene labels used in column 3 of gff file were `gene` and `ncRNA_gene`
+ ase_get_imprinted_gene_list
    + Generate gene list for is gene imprinted step
+ ase_cross_ref
    + Bottleneck: Cross ref ASE variants for 120 samples with 30K MAF >= 0.05 varaints for ~30K genes
+ ase_is_gene_imprinted
    + Check if gene is consistent with genomic imprinting which is defined as 
    + At least 90% of reads map to one of the two alleles in 80% of our heterozygotes for each SNP
    
**Note**: To improve performance the last two steps were run on GPU. There is code to run it on normal node and GPU node.

***

### **BiomaRt method**

The original idea was to measure ASE in our 120 samples in a list of 228 imprinted genes. The initial attempt was to use BiomaRt, but there were some issues with time outs. This method would also not be suitable for a genome wide screen.

1. Get chromosome, start/stop coordinates for all genes (189 genes left - note some genes had multiple entries after this stage)
2. Get all rsIDs with MAF >= 0.05 for each gene - (169 genes left)
    + Only 180 genes processed as some genes were too large for software - I can get the SNPs for these manually
    + A further 11 genes had no SNPs within that range 
3. For SNPs MAF >= 0.05 filter those with ASE reads in >= 1 sample - (164 genes left - this step needed to make programming simpler)
4. For SNPs MAF >= 0.05 filter those with >= 20 ASE reads in >= 10 samples - (103 genes left - note thats 20 reads total across alleles)

***

### Running the pipeline to deal with scratch quota limits

The pipeline has to be run in different blocks in order to balance the competing requirements
of the scratch file / memory quota limits which are breeched early on and a huge number of small jobs later on.

During the `fastq_prep.smk`, multiple versions of the fastq files can be generated at each stage. This results in
the quota limit being breeched and pipeline choking at random points. To resolve this I used the following:

+ The snakemake `temp()` function to delete olde versions of fastqs as we moved through the pipeline
+ Lowered the max number of jobs that the pipeline could run on Hawk at any one time from 500 to 50
+ Gave later jobs in the process higher priority than earlier jobs


At the later `annotation.smk` step for ASE and SNP cross referencing step, for the genome wide anaysis we need to run 19K x 120 jobs,
these take ~30-60s and require negligable resources. It is important to reinstate the job limit to 500 jobs to churn through these
jobs quicker.     

***










