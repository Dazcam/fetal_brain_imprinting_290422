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
                print(wildcards.sample, ": has > 1 fastq file per read.\n", {input.r1})
                shell("cat {params.indir}{wildcards.sampleID}*_R1.fastq.gz > {output.r1} 2> {log.r1}"),
                shell("cat {params.indir}{wildcards.sampleID}*_R2.fastq.gz > {output.r2} 2> {log.r2}")

            else:
                print(wildcards.sample, ": has only 1 fastq file per read.\n", {input.r1})
                shell("mv {params.indir}{wildcards.sampleID}*_R1.fastq.gz {params.outdir}{wildcards.sampleID}_R1.fastq.gz 2> {log.r1}"),
                shell("mv {params.indir}{wildcards.sampleID}*_R2.fastq.gz {params.outdir}{wildcards.sampleID}_R2.fastq.gz 2> {log.r2}")


rule fastqc_pretrim:
    # Note: Running both fqc runs for R1 and R2 in one rule, not as efficient as it could be but compromise for now
    input:   rules.merge_fastqs.output.r1 
    output:  "../results/04FASTQC/pre-trimmed/{sampleID}_R2_fastqc.zip"
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
    input: expand(rules.fastqc_pretrim.output, sampleID = ALL_MERGED_SAMPLES)
    output:  "../results/05MULTIQC/pre-trimmed/multiQC_log.html"
    log:     "../results/log/05MULTIQC/pre-trimmed/multiqc.log"
    params:  indir = "../results/04FASTQC/pre-trimmed", outfile = "../results/05MULTIQC/pre-trimmed"
    message: "multiqc for pre-trimmed fastqc"
    shell:
             """
             multiqc {params.indir} -o {output} -d -f -v -n {params.outfile} 2> {log}

             """
