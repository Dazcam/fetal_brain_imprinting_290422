rule ase_get_SNP_refs:
    output: "../results/12GETREFS/homo_sapiens-chr{CHR}.vcf.gz",
    params: OUT_DIR = "../results/12GETREFS/",
            REPO_SNPS_DIR = "https://ftp.ensembl.org/pub/release-108/variation/vcf/homo_sapiens/",
            SNPS = "homo_sapiens-chr{CHR}.vcf.gz",
    log:    "../results/00LOG/12GETREFS/get_SNP_refs_{CHR}.log"
    shell:
            "wget -nc {params.REPO_SNPS_DIR}{params.SNPS} -P {params.OUT_DIR} &> {log}"

rule ase_get_GENE_refs:
    output: "../results/12GETREFS/Homo_sapiens.GRCh38.108.chromosome.{CHR}.gff3.gz",
    params: OUT_DIR = "../results/12GETREFS/",
            REPO_GENES_DIR = "https://ftp.ensembl.org/pub/release-108/gff3/homo_sapiens/",
            GENES = "Homo_sapiens.GRCh38.108.chromosome.{CHR}.gff3.gz"
    log:    "../results/00LOG/12GETREFS/get_GENE_refs_{CHR}.log"
    shell:
            "wget -nc {params.REPO_GENES_DIR}{params.GENES} -P {params.OUT_DIR} &> {log}"

rule ase_map_snps_to_genes_genomewide:
    input:  "../results/12GETREFS/homo_sapiens-chr{CHR}.vcf.gz",
            "../results/12GETREFS/Homo_sapiens.GRCh38.108.chromosome.{CHR}.gff3.gz"
    output: "../results/13SNP2GENES_GW/chr{CHR}_genes_no_snps.txt"
    params: OUT_DIR = "../results/13SNP2GENES_GW/",
            REF_DIR = "../results/12GETREFS/"
    envmodules: "bcftools"
    log:    "../results/00LOG/13SNP2GENES_GW/snps2genes_gw_chr{CHR}.log"
    shell:
            "scripts/imprinting_map_snps_2_genes.sh {wildcards.CHR} {params.OUT_DIR} {params.REF_DIR} &> {log}"

rule ase_get_genomewide_gene_list:
    input:  expand("../results/13SNP2GENES_GW/chr{CHR}_genes_no_snps.txt", CHR = range(1,23))
    output: "../results/14GETGENELIST/genomewide_genelist.txt"
    shell:
            "ls ../results/13SNP2GENES_GW/genes/ | cut -d_ -f1 > {output}"

rule ase_get_imprinted_gene_list:
    input:  tucci_genes = "../resources/sheets/Tucci_2019_mmc1.xlsx",
            genomewide_genes = "../results/14GETGENELIST/genomewide_genelist.txt"
    output: "../results/14GETGENELIST/Tucci_2019_genes_in_gw_list.txt",
            "../results/14GETGENELIST/Tucci_2019_imprinted_genes.txt"
    log:    "../results/00LOG/14GETGENELIST/get_imprinted_gene_list.log"
    params: "../results/14GETGENELIST/"
    envmodules: "libgit2/1.1.0"
    script:
            "../scripts/imprinting_get_imprinted_gene_list.R"

rule ase_rename_tucci_genes:
    input:  "../results/14GETGENELIST/Tucci_2019_imprinted_genes.txt",
            "../resources/sheets/tucci_genes_alias.txt"
    output: "../results/14GETGENELIST/Tucci_2019_imprinted_genes_alias.txt"
    log:    "../results/00LOG/14GETGENELIST/ase_mv_GW_outfiles_from_server.log"
    shell:
            "scripts/imprinting_rename_tucci_genes.sh {input[0]} {input[1]} {output} &> {log}"
