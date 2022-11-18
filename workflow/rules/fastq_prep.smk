import numpy
import pandas as pd
import os
import json

# -------  LOAD SAMPLE LIST  ---------

sample_file = config['SAMPLE_LIST']
sample_df = pd.read_table(sample_file, sep="\t+", header=0)
sampleID = sample_df.sampleID
fastq = sample_df.fastq


#FILES = json.load(open(config['SAMPLES_JSON']))
MERGE_FILES = json.load(open(config['MERGE_FILES_JSON']))
#print()
ALL_SAMPLES = sorted(MERGE_FILES.keys())
ALL_MERGED_SAMPLES = sorted(MERGE_FILES.keys())
#print(ALL_MERGED_SAMPLES)


rule move_and_rename_fastqs:
    input:  fastq = lambda w: sample_df[sample_df.sampleID == w.sample].fastq.tolist()
    output: temp("../results/01RAW_fqs/{sample}")
    log:    "../results/log/01RAW_fqs/{sample}.log"
    shell:
            "cp {input.fastq} {output} 2> {log}"

rule zip_fastqs:
    input:  "../results/01RAW_fqs/{sample}"
    output: temp("../results/02ZIPPED_fqs/{sample}")
    params: outdir = "../results/02ZIPPED_fqs/"
    log:    "../results/log/02ZIPPED_fqs/{sample}.log"
    run:

        if wildcards.sample.endswith('.fastq'):
            shell("""gzip {input} 2> {log}""")
            shell("""mv {input}.gz {params.outdir} 2> {log}""")
            shell("""touch {output}""")
        else:
            shell("""mv {input} {params.outdir} 2> {log}""")
            shell("""touch {output}""")

rule merge_fastqs:
    input:  r1 = lambda wildcards: MERGE_FILES[wildcards.sampleID]['R1'],
            r2 = lambda wildcards: MERGE_FILES[wildcards.sampleID]['R2']
    output: r1 = temp("../results/03MRGD_fqs/{sampleID}_R1.fastq.gz"),
            r2 = temp("../results/03MRGD_fqs/{sampleID}_R2.fastq.gz")
    log:    r1 = "../results/log/03MRGD_fqs/{sampleID}_R1.log",
            r2 = "../results/log/03MRGD_fqs/{sampleID}_R2.log"
    params: indir = "../results/02ZIPPED_fqs/", 
            outdir = "../results//03MRGD_fqs/"
    run:
            if len(input.r1) > 1:
                print(wildcards.sampleID, ": has > 1 fastq file per read.\n", {input.r1})
                shell("cat {params.indir}{wildcards.sampleID}*_R1.fastq.gz > {output.r1} 2> {log.r1}"),
                shell("cat {params.indir}{wildcards.sampleID}*_R2.fastq.gz > {output.r2} 2> {log.r2}")

            else:
                print(wildcards.sampleID, ": has only 1 fastq file per read.\n", {input.r1})
                shell("mv {params.indir}{wildcards.sampleID}*_R1.fastq.gz {params.outdir}{wildcards.sampleID}_R1.fastq.gz 2> {log.r1}"),
                shell("mv {params.indir}{wildcards.sampleID}*_R2.fastq.gz {params.outdir}{wildcards.sampleID}_R2.fastq.gz 2> {log.r2}")


rule fastqc_pretrim:
    # Note: Running both fqc runs for R1 and R2 in one rule, not as efficient as it could be but compromise for now
    input:   rules.merge_fastqs.output.r1 
    output:  fastqs = "../results/04FASTQC/pre-trimmed/{sampleID}_R2_fastqc.zip",
    log:     "../results/log/04FASTQC/pre-trimmed/{sampleID}_fastqc"
    params:  outdir = "../results/04FASTQC/pre-trimmed/"
    message: "fastqc {input}"
    shell:
             """
             # fastqc works fine on .gz file as well
             module load FastQC
             fastqc -o {params.outdir} -f fastq ../results/03MRGD_fqs/{wildcards.sampleID}_R1.fastq.gz  2> {log}
             fastqc -o {params.outdir} -f fastq ../results/03MRGD_fqs/{wildcards.sampleID}_R2.fastq.gz  2> {log}
             """

rule multiQC_pretrim:
    input:   expand(rules.fastqc_pretrim.output, sampleID = ALL_MERGED_SAMPLES)
    output:  "../results/05MULTIQC/pre-trimmed/pre-trimmed.html"
    log:     "../results/log/05MULTIQC/pre-trimmed/multiqc.log"
    params:  indir = "../results/04FASTQC/pre-trimmed",
             outdir = "../results/05MULTIQC/pre-trimmed/",
             outname = "pre-trimmed"
    message: "multiqc for pre-trimmed fastqc"
    shell:
             """
             module load multiqc
             multiqc {params.indir} -o {params.outdir} -f -v -n {params.outname} 2> {log}
             """

rule trim_fastq:
    input:   r1 = "../results/03MRGD_fqs/{sampleID}_R1.fastq.gz",
             r2 = "../results/03MRGD_fqs/{sampleID}_R2.fastq.gz"
    output:  "../results/06TRIM_FQs/{sampleID}_R1_val_1.fq.gz",
             "../results/06TRIM_FQs/{sampleID}_R2_val_2.fq.gz",
             "../results/04FASTQC/trimmed_mrg/{sampleID}_R1_val_1_fastqc.zip",
             "../results/04FASTQC/trimmed_mrg/{sampleID}_R2_val_2_fastqc.zip"
    log:     "../results/log/06TRIM_FQs/{sampleID}.log"
    threads: 4
    params:  outdir_trim = "../results/06TRIM_FQs/",
             outdir_trim_fqc = "../results/04FASTQC/trimmed_mrg/"
    message: "Trim Galore: {input}"
    shell:
             """
             module load trimgalore
             trim_galore -e 0.1 -q 20  --paired --basename {wildcards.sampleID} \
             --illumina --output_dir {params.outdir_trim} -j 4 \
             --fastqc_args "-o {params.outdir_trim_fqc} -f fastq" {input.r1} {input.r2} 2> {log}
             """

rule hard_trim_fastq:
    input:   r1 = "../results/06TRIM_FQs/{sampleID}_R1_val_1.fq.gz", 
             r2 = "../results/06TRIM_FQs/{sampleID}_R2_val_2.fq.gz" 
    output:  r1 = temp("../results/07HARDTRIM_FQs/{sampleID}_R1_val_1.fq.gz"),
             r2 = temp("../results/07HARDTRIM_FQs/{sampleID}_R2_val_2.fq.gz")
    log:     "../results/log/07HARDTRIM_fqs/{sampleID}.log"
    params:  crdf = 54, extr = 79, edin = 104 # Co-ords are 0 based
    message: "\nHard trimming {wildcards.sampleID}\n"
    run:
            if "Crdf" in wildcards.sampleID:
                
                print("\n", wildcards.sampleID, "will be trimmed to a length of 55 reads\n")
                shell("../resources/bbmap/bbduk.sh in1={input.r1} in2={input.r2} out1={output.r1} out2={output.r2} ftr={params.crdf} ordered=t 2> {log}")
            
                
            elif "Extr" in wildcards.sampleID:

                print("\n", wildcards.sampleID, "will be trimmed to a length of 80 reads\n")
                shell("../resources/bbmap/bbduk.sh in1={input.r1} in2={input.r2} out1={output.r1} out2={output.r2} ftr={params.extr} ordered=t 2> {log}")       
           
            else:
              
                print("\n", wildcards.sampleID, "will be trimmed to a length of 105 reads\n")
                shell("../resources/bbmap/bbduk.sh in1={input.r1} in2={input.r2} out1={output.r1} out2={output.r2} ftr={params.edin} ordered=t 2> {log}")

rule remove_short_reads:
    input:   r1 = "../results/07HARDTRIM_FQs/{sampleID}_R1_val_1.fq.gz",
             r2 = "../results/07HARDTRIM_FQs/{sampleID}_R2_val_2.fq.gz"
    output:  r1 = "../results/07HARDTRIM_FQs/{sampleID}_R1_val_1.noSR.fq.gz",
             r2 = "../results/07HARDTRIM_FQs/{sampleID}_R2_val_2.noSR.fq.gz"
    log:     "../results/log/07HARDTRIM_FQs/{sampleID}.noSR.log"
    params:  crdf = 55, extr = 80, edin = 105
    message: "Removing short reads for {wildcards.sampleID}"
    run:
             if "Crdf" in wildcards.sampleID:

                 print("\nReads shorter than 55bps will be removed from: ", wildcards.sampleID, "\n")
                 shell("../resources/bbmap/reformat.sh in1={input.r1} in2={input.r2} out1={output.r1} out2={output.r2} minlength={params.crdf} 2> {log}")
             
             elif "Extr" in wildcards.sampleID:
                 
                 print("\nReads shorter than 80bps will be removed from: ", wildcards.sampleID, "\n")
                 shell("../resources/bbmap/reformat.sh in1={input.r1} in2={input.r2} out1={output.r1} out2={output.r2} minlength={params.extr} 2> {log}")
            
             else:

                 print("\nReads shorter than 105bps will be removed from: ", wildcards.sampleID, "\n")
       	       	 shell("../resources/bbmap/reformat.sh in1={input.r1} in2={input.r2} out1={output.r1} out2={output.r2} minlength={params.edin} 2> {log}")

rule read_length_dist_post_QC_and_trimGalore:
    input:   trm_r1 = "../results/06TRIM_FQs/{sampleID}_R1_val_1.fq.gz",
             trm_r2 = "../results/06TRIM_FQs/{sampleID}_R2_val_2.fq.gz",
    output:  trm_r1 = "../results/08READ_DISTRIBUTION/trimmed/{sampleID}_R1_readDist.txt",
             trm_r2 = "../results/08READ_DISTRIBUTION/trimmed/{sampleID}_R2_readDist.txt",
    log:     trm = "../results/log/08READ_DISTRIBUTION/trimmed/{sampleID}.log",
    message: "Extracting read distribution for {wildcards.sampleID}"
    shell:
             """
             printf "\n--------------------------\nExtracting read dist for: {wildcards.sampleID}_R1\n--------------------------\n\n" 2> {log.trm}
             gunzip -c {input.trm_r1} | awk 'NR%4==2{{print length($0)}}' | sort | uniq -c | sort -nr > {output.trm_r1} 2> {log.trm}
             echo "\n--------------------------\nExtracting read dist for: {wildcards.sampleID}_R2\n--------------------------\n\n" 2> {log.trm}
             gunzip -c {input.trm_r2} | awk 'NR%4==2{{print length($0)}}' | sort | uniq -c | sort -nr > {output.trm_r2} 2> {log.trm}
             """

rule get_read_length_dist_post_hard_and_SrtRead_trim:
    input:   hrd_trm_r1 = "../results/07HARDTRIM_FQs/{sampleID}_R1_val_1.fq.gz",
             hrd_trm_r2 = "../results/07HARDTRIM_FQs/{sampleID}_R2_val_2.fq.gz",
             hrd_trm_noSR_r1 = "../results/07HARDTRIM_FQs/{sampleID}_R1_val_1.noSR.fq.gz",
             hrd_trm_noSR_r2 = "../results/07HARDTRIM_FQs/{sampleID}_R2_val_2.noSR.fq.gz"
    output:  hrd_trm_r1 = "../results/08READ_DISTRIBUTION/hrd_trm/{sampleID}_R1_readDist.txt",
             hrd_trm_r2 = "../results/08READ_DISTRIBUTION/hrd_trm/{sampleID}_R2_readDist.txt",
             hrd_trm_noSR_r1 = "../results/08READ_DISTRIBUTION/hrd_trm_noSR/{sampleID}_R1_readDist.txt",
             hrd_trm_noSR_r2 = "../results/08READ_DISTRIBUTION/hrd_trm_noSR/{sampleID}_R2_readDist.txt",
    log:     hrd_trm = "../results/log/08READ_DISTRIBUTION/hrd_trm/{sampleID}.log",
             hrd_trm_noSR = "../results/log/08READ_DISTRIBUTION/hrd_trm_noSR/{sampleID}.log"

    message: "Extracting read distribution for {wildcards.sampleID}"

    run:
             print("\n--------------------------\nExtracting read dist for: {wildcards.sampleID}_R1\n--------------------------\n\n")
             shell("gunzip -c {input.hrd_trm_r1} | awk 'NR%4==2{{print length($0)}}' | sort | uniq -c | sort -nr > {output.hrd_trm_r1} 2> {log.hrd_trm}")

             print("\n--------------------------\nExtracting read dist for: {wildcards.sampleID}_R2\n--------------------------\n\n")
             shell("gunzip -c {input.hrd_trm_r2} | awk 'NR%4==2{{print length($0)}}' | sort | uniq -c | sort -nr > {output.hrd_trm_r2} 2> {log.hrd_trm}")

             print("\n--------------------------\nExtracting read dist for: {wildcards.sampleID}_R1\n--------------------------\n\n")
             shell("gunzip -c {input.hrd_trm_noSR_r1} | awk 'NR%4==2{{print length($0)}}' | sort | uniq -c | sort -nr > {output.hrd_trm_noSR_r1} 2> {log.hrd_trm_noSR}")

             print("\n--------------------------\nExtracting read dist for: {wildcards.sampleID}_R2\n--------------------------\n\n")
             shell("gunzip -c {input.hrd_trm_noSR_r2} | awk 'NR%4==2{{print length($0)}}' | sort | uniq -c | sort -nr > {output.hrd_trm_noSR_r2} 2> {log.hrd_trm_noSR}")
