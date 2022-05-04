rule build_static_index:
    input:   ref = "../resources/refs/genome_sequence/" + config['REF_GENOME'],
             gtf = "../resources/refs/annotation_files/" + config['GENE_ANNOT']
    output:  "../results/09STATIC_INDEX/ASElux_static_index.annotation"
    log:     "../results/logs/09STATIC_INDEX/build_static_index.log"
    params:  outdir = "../results/09STATIC_INDEX/ASElux_static_index"
    shell:       
             """
             module load compiler/gnu/5
             ASElux build --gtf {input.gtf} \
             --ref {input.ref} \
             --out {params.outdir} 2> {log} 
             """
 
rule extract_sample_genotypes:
    input:   "../resources/vcf/genotypes.chr.vcf-v4.2"
    output:  "../results/10VCFs/{sampleID}.chr.vcf"
    log:     "../results/log/10VCFs/{sampleID}.log"
    params : jobname = "{sampleID}"
    message: "Extracting {wildcards.sampleID} genotype info"
    shell:
             """
             module load bcftools
             bcftools annotate -x INFO,^FORMAT/GT {input} \
             | bcftools view -Ov -s {wildcards.sampleID} > {output} 2> {log}     
             """

rule ase_align:
    input:   genotype = "../results/10VCFs/{sampleID}.chr.vcf",
             index = "../results/09STATIC_INDEX/ASElux_static_index.annotation",
             r1 = "../results/07HARDTRIM_FQs/{sampleID}_val_1.noSR.fq.gz",
             r2 = "../results/07HARDTRIM_FQs/{sampleID}_val_2.noSR.fq.gz"
    output:  "../results/11ASE_MAPPINGS/{sampleID}.ase.txt"    
    log:     "../results/log/11ASE_MAPPINGS/{sampleID}.ase.log"
    
    params:  jobname = "{sampleID}", crdf = 55, extr = 80, edin = 105,
             index = "../results/09STATIC_INDEX/ASElux_static_index"
    message: "Calculating allele specific expression for sample {wildcards.sampleID}"
    run: 
             if "Crdf" in wildcards.sampleID:

                 print("\nRunning ASE for: ", wildcards.sampleID, "\n")
                           
                 shell("""
                 
                 module load compiler/gnu/5
                 PREFIX_r1=$(echo "{input.r1}" | cut -f 1-4 -d '.')
                 PREFIX_r2=$(echo "{input.r2}" | cut -f 1-4 -d '.')
       	     
       	         gunzip -c {input.r1} > ${{PREFIX_r1}}
       	         gunzip -c {input.r2} > ${{PREFIX_r2}}
                 ASElux align --fq --pe --readLen {params.crdf} --index {params.index} \
                 --vcf {input.genotype} --seqFiles ${{PREFIX_r1}} ${{PREFIX_r2}} \
                 --out {output} 2> {log}
                 rm ${{PREFIX_r1}}
                 rm ${{PREFIX_r2}} 
                 
                 """)

             elif "Extr" in wildcards.sampleID:

                 print("\nRunning ASE for: ", wildcards.sampleID, "\n")

                 shell("""
                 module load compiler/gnu/5
                 PREFIX_r1=$(echo '{input.r1}' | cut -f 1-4 -d '.')
                 PREFIX_r2=$(echo '{input.r2}' | cut -f 1-4 -d '.')
                 gunzip -c {input.r1} > ${{PREFIX_r1}}
                 gunzip -c {input.r2} > ${{PREFIX_r2}}
                 ASElux align --fq --pe --readLen {params.extr} --index {params.index} \
                 --vcf {input.genotype} --seqFiles ${{PREFIX_r1}} ${{PREFIX_r2}} \
                 --out {output} 2> {log}
                 rm ${{PREFIX_r1}} 
                 rm ${{PREFIX_r2}} 
       	         
                 """)

             else:

                 print("\nRunning ASE for: ", wildcards.sampleID, "\n")

                 shell("""
                 module load compiler/gnu/5
                 PREFIX_r1=$(echo '{input.r1}' | cut -f 1-4 -d '.')
                 PREFIX_r2=$(echo '{input.r2}' | cut -f 1-4 -d '.')
                 gunzip -c {input.r1} > ${{PREFIX_r1}}
                 gunzip -c {input.r2} > ${{PREFIX_r2}}
                 ASElux align --fq --pe --readLen {params.edin} --index {params.index} \
                 --vcf {input.genotype} --seqFiles ${{PREFIX_r1}} ${{PREFIX_r2}} \
                 --out {output} 2> {log}
                 rm ${{PREFIX_r1}}
                 rm ${{PREFIX_r2}}
                  
                 """)
