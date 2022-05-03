## Migration from XXXXXX - 

fastq_prep.smk

1. move_and_rename_fastqs - copy fastq files from neurocluster databank
2. zip_fastqs - standardise format of fastq files (some were zipped and some were not)
3. merge_fastqs - merge fastqs into single files for R1 and R2 (some samples were sequenced over multiple lanes)
4. 



NOTE: In former version of this workflow I generated json files to deal with the complex merge_fastq process
which had multiple different number of input files for each run of merge fastq. The generation of these files
occured on the output of the zip_fastqs rule. However, at the time I couldn't work out a way to integrate the 
json generation step into the snakemake workflow. For now the json scripts and data are stored in resources/json/.  

