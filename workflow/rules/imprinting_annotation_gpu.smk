# Needs to be run from /localscratch/scw1641
# salloc -p gpu_v100_test -t 2:00:00 --gres=gpu:0 --account=scw1641
# ssh ccs2101 

import pandas as pd

GENE_DF = pd.read_table("../results/14GETGENELIST/genomewide_genelist.txt")
GENE_LIST = list(GENE_DF.iloc[:,0])


rule ase_cross_ref_summary:
    input:  gene_file = "../results/14GETGENELIST/genomewide_genelist.txt",
            genes = expand("/localscratch/scw1641/results/15CROSSREF/{GENE}/touch.txt", GENE = GENE_LIST)
    output: "/localscratch/scw1641/results/15CROSSREF/summary_6_ase_imprinting_genes_summary.txt"
    params: indir = "../results/13SNP2GENES_GW/genes/",
            outdir = "/localscratch/scw1641/results/15CROSSREF/"
    log:    "../results/00LOG/15CROSSREF/ase_cross_ref_summary.log"
    shell:
            "scripts/imprinting_crossRef_summary.sh {input.gene_file} {params.indir} {params.outdir}  &> {log}"

rule ase_is_gene_imprinted:
    input:  "/localscratch/scw1641/results/15CROSSREF/summary_6_ase_imprinting_genes_summary.txt"
    output: "/localscratch/scw1641/results/16ISGENEIMPRINTED/imp_test.txt"
    params: indir = "/localscratch/scw1641/results/15CROSSREF/",
            outdir = "/localscratch/scw1641/results/16ISGENEIMPRINTED/"
    log:    "../results/00LOG/16ISGENEIMPRINTED/ase_is_gene_imprinted.log"
    shell:
            "scripts/imprinting_is_gene_imprinted.sh {params.indir} {params.outdir}  &> {log}"

rule ase_mv_GW_outfiles_from_server:
    input:  "/localscratch/scw1641/results/15CROSSREF/summary_6_ase_imprinting_genes_summary.txt",
            "/localscratch/scw1641/results/16ISGENEIMPRINTED/imp_test.txt"
    output: "../results/17GW_OUTPUT/summary_6_ase_imprinting_genes_summary.txt",
            "../results/17GW_OUTPUT/imp_test.txt",
    params: outfile = "/localscratch/scw1641/results/15CROSSREF/summary_*",
            outdir = "../results/17GW_OUTPUT/",
    log:    "../results/00LOG/17GW_OUTPUT/ase_mv_GW_outfiles_from_server.log"
    shell:
            "cp {params.outfile} {params.outdir} &> {log}; "
            "cp {input[1]} {params.outdir} &> {log}" 

rule ase_extract_imp_genes:
    input:  "../results/17GW_OUTPUT/imp_test.txt"
    output: "../results/17GW_OUTPUT/imp_or_part.txt",
            "../results/17GW_OUTPUT/imp_or_part_uniq_genes.txt"
    log:    "../results/00LOG/17GW_OUTPUT/ase_extract_imp_genes.log"
    shell:
            "grep -v Not {input} > {output} 2> {log}; "
            "cut -f1 {output[0]} | uniq | cut -d: -f2 > {output[1]} 2> {log}"


