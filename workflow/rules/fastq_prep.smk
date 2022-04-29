rule move_and_rename_fastqs:
    input:  fastq = lambda w: sample_df[sample_df.sampleID == w.sample].fastq.tolist()
    output: temp("01RAW_fqs/{sample}")
    log:    "/00logs/01RAW_fqs/{sample}.log"
    shell:
            "cp {input.fastq} {output} 2> {log}"

rule zip_fastqs:
    input:  "01RAW_fqs/{sample}"
    output: temp("02ZIPPED_fqs/{sample}")
    params: outdir = "02ZIPPED_fqs/"
    log:    "00logs/02ZIPPED_fqs/{sample}.log"
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
    output: r1 = temp("03MRGD_fqs/{sampleID}_R1.fastq.gz"),
            r2 = temp("03MRGD_fqs/{sampleID}_R2.fastq.gz")
    log:    r1 = "00logs/03MRGD_fqs/{sampleID}_R1.log",
            r2 = "00logs/03MRGD_fqs/{sampleID}_R2.log"
    params: indir = "02ZIPPED_fqs/", 
            outdir = "/03MRGD_fqs/"
    run:
            if len(input.r1) > 1:
                print(wildcards.sample, ": has > 1 fastq file per read.\n", {input.r1})
                shell("cat {params.indir}{wildcards.sampleID}*_R1.fastq.gz > {output.r1} 2> {log.r1}"),
                shell("cat {params.indir}{wildcards.sampleID}*_R2.fastq.gz > {output.r2} 2> {log.r2}")

            else:
                print(wildcards.sample, ": has only 1 fastq file per read.\n", {input.r1})
                shell("mv {params.indir}{wildcards.sampleID}*_R1.fastq.gz {params.outdir}{wildcards.sampleID}_R1.fastq.gz 2> {log.r1}"),
                shell("mv {params.indir}{wildcards.sampleID}*_R2.fastq.gz {params.outdir}{wildcards.sampleID}_R2.fastq.gz 2> {log.r2}")
