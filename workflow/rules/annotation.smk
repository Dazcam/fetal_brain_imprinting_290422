rule ase_get_SNP_refs:
    output: "../results/13GETREFS/homo_sapiens-chr{CHR}.vcf.gz",
    params: OUT_DIR = "../results/13GETREFS/",
            REPO_SNPS_DIR = "https://ftp.ensembl.org/pub/release-108/variation/vcf/homo_sapiens/",
            SNPS = "homo_sapiens-chr{CHR}.vcf.gz",
    log:    "../results/00LOG/13GETREFS/get_SNP_refs_{CHR}.log"
    shell:
            "wget -nc {params.REPO_SNPS_DIR}{params.SNPS} -P {params.OUT_DIR} &> {log}"

rule ase_get_GENE_refs:
    output: "../results/13GETREFS/Homo_sapiens.GRCh38.108.chromosome.{CHR}.gff3.gz",
    params: OUT_DIR = "../results/13GETREFS/",
            REPO_GENES_DIR = "https://ftp.ensembl.org/pub/release-108/gff3/homo_sapiens/",
            GENES = "Homo_sapiens.GRCh38.108.chromosome.{CHR}.gff3.gz"
    log:    "../results/00LOG/13GETREFS/get_GENE_refs_{CHR}.log"
    shell:
            "wget -nc {params.REPO_GENES_DIR}{params.GENES} -P {params.OUT_DIR} &> {log}"

rule ase_map_snps_to_genes_genomewide:
    input:  "../results/13GETREFS/homo_sapiens-chr{CHR}.vcf.gz",
            "../results/13GETREFS/Homo_sapiens.GRCh38.108.chromosome.{CHR}.gff3.gz"
    output: "../results/14SNP2GENES_GW/chr{CHR}_genes_no_snps.txt"
    params: OUT_DIR = "../results/14SNP2GENES_GW/",
            REF_DIR = "../results/13GETREFS/"
    envmodules: "bcftools"
    log:    "../results/00LOG/14SNP2GENES_GW/snps2genes_gw_chr{CHR}.log"
    shell:
            "scripts/imprinting_map_snps_2_genes.sh {wildcards.CHR} {params.OUT_DIR} {params.REF_DIR} &> {log}"

rule ase_get_imprinted_gene_list:
    input:  genes = "../resources/sheets/Tucci_2019_mmc1.xlsx"
    output: "../results/15GETGENELIST/Tucci_2019_imprinted_genes.txt"
    params: "../results/15GETGENELIST/"
    envmodules: "libgit2/1.1.0"
    script:
            "../scripts/imprinting_get_imprinted_gene_list.R"

rule ase_get_genomewide_gene_list:
    input:  expand("../results/14SNP2GENES_GW/chr{CHR}_genes_no_snps.txt", CHR = range(1,23))
    output: "../results/15GETGENELIST/genomewide_genelist.txt"
    shell:
            "ls ../results/14SNP2GENES_GW/genes/ | cut -d_ -f1 > {output}"

rule ase_cross_ref:
    input:  "../results/12SNP2GENES/snps2genes_maf_0.05.log"
    output: "../results/13CROSSREF/summary_6_ase_imprinting_genes_summary.txt"
    params: "../results/13CROSSREF/"
    log:    "../results/00LOG/13CROSSREF/ase_cross_ref.log"
    shell:
            "scripts/imprinting_crossRef_ase_over_snps.sh {params} 2> {log}"
