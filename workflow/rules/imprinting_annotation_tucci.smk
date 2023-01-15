rule ase_tucci_gene_overlap:
    input:  "../results/17GW_OUTPUT/imp_or_part_uniq_genes.txt",
            "../results/14GETGENELIST/Tucci_2019_imprinted_genes_alias.txt"
    output: "../results/18TUCCI_ANN/imprinting_novel_genes.txt",
            "../results/18TUCCI_ANN/imprinting_tucci_only_genes.txt",
            "../results/18TUCCI_ANN/imprinting_overlap_genes.txt",
            "../results/18TUCCI_ANN/imprinting_all_overlap_genes.txt"
    log:    "../results/00LOG/18TUCCI_ANN/tucci_gene_overlap.log"
    shell:
            """

            comm -23 <(sort {input[0]}) <(sort {input[1]}) > {output[0]} 2> {log}
            comm -13 <(sort {input[0]}) <(sort {input[1]}) > {output[1]} 2> {log}
            comm -12 <(sort {input[0]}) <(sort {input[1]}) > {output[2]} 2> {log}
            printf 'novel\toverlap\ttucci_only\n' > {output[3]} 2> {log}
            paste {output[0]} {output[1]} {output[2]} >> {output[3]} 2> {log} 

            """ 

rule ase_gene_list_fishers:
    input:  NOVEL_GENES = "../results/18TUCCI_ANN/imprinting_novel_genes.txt",
            OVERLAP_GENES = "../results/18TUCCI_ANN/imprinting_overlap_genes.txt",
            TUCCI_GENES = "../results/18TUCCI_ANN/imprinting_tucci_only_genes.txt",
            GW_GENES = "../results/14GETGENELIST/genomewide_genelist.txt",
            NDD_GENES = "../resources/sheets/NDD_denovoWest_Ensemble.tsv",
            SCZ_CV_GENES = "../resources/sheets/Supplementary_Table_12_Prioritized_Genes_UPDATED.xlsx",
            SCZ_RV_GENES = "../resources/sheets/schema_gene_conversion.txt", 
            ASD_GENES = "../resources/sheets/Fu_preprint_ASD_FDR001.tsv",
            DDG2p_GENES = "../resources/sheets/DDG2P_mono_Ensemble.tsv",
            MARKDOWN_FILE = "scripts/imprinting_gene_list_fishers.Rmd"
    output: "../results/18TUCCI_ANN/imprinting_gene_list_fishers.html"  
    params: REPORT_DIR = "../results/18TUCCI_ANN/",
            REPORT_FILE = "imprinting_gene_list_fishers.html"
    log:    "../results/00LOG/18TUCCI_ANN/ase_gene_list_fishers.log"
    envmodules: "libgit2/1.1.0", "R/4.2.0", "pandoc/2.7.3"    
    script: 
            "../scripts/imprinting_gene_list_fishers.R"
          
