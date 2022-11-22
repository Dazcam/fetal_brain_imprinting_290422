rule ase_map_snps_to_genes:
    input:  genes = "../resources/sheets/Tucci_2019_mmc1.xlsx"
#           ase = "../results/11ASE_MAPPINGS/{sampleID}.ase.txt"
    output: "../results/12SNP2GENES/snps2genes_maf_0.05.log"
    params: "../results/12SNP2GENES/"
    envmodules: "libgit2/1.1.0"
    script:
            "../scripts/imprinting_map_snps_2_genes.R"

rule ase_cross_ref:
    input:  "../results/12SNP2GENES/snps2genes_maf_0.05.log"
    output: "../results/13CROSSREF/summary_6_ase_imprinting_genes_summary.txt"
    params: "../results/13CROSSREF/"
    log:    "../results/00LOG/13CROSSREF/ase_cross_ref.log"
    shell:
            "scripts/imprinting_crossRef_ase_over_snps.sh {params} 2> {log}"
