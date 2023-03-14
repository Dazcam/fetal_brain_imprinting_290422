# A genome-wide screen of ASE 


### **Ensembl data dump method**

+ [GFF3 README](https://ftp.ensembl.org/pub/release-108/gff3/homo_sapiens/README)
+ [GFF3 description](http://gmod.org/wiki/GFF3)

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

***










