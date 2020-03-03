# Snakemake workflow to reproduce this analysis

rule all:
    input:
        'results/DE_genes_H9.tsv',
        'results/DE_genes_102.tsv',
        'results/DE_genes_106.tsv',
        'results/DE_upset.png',
        'results/KEGG_2016_GSEA.tsv',
        'results/GO_Biological_Process_2018_GSEA.tsv',
        'results/KEGG_2016_upset.png',
        'results/GO_Biological_Process_2018_upset.png'

rule generate_tsv:
    input:
        'raw_data/raw_data.tsv'
    output:
        'data/data.tsv'
    run:
        shell('python generate_tsv.py {input} {output}')

rule run_ebseq:
    input:
        'data/data.tsv'
    output:
        'results/EBSeq_DE_H9.tsv',
        'results/EBSeq_DE_102.tsv',
        'results/EBSeq_DE_106.tsv'
    run:
        shell('python de_analysis.py {input} results')

rule process_results:
    input:
        'results/EBSeq_DE_H9.tsv',
        'results/EBSeq_DE_102.tsv',
        'results/EBSeq_DE_106.tsv'
    output:
        'results/DE_genes_H9.tsv',
        'results/DE_genes_102.tsv',
        'results/DE_genes_106.tsv', 
        'results/DE_upset.png',
        'results/KEGG_2016_GSEA.tsv',
        'results/GO_Biological_Process_2018_GSEA.tsv',
        'results/KEGG_2016_upset.png',
        'results/GO_Biological_Process_2018_upset.png'
    run:
        shell('python analyze_de_results.py results')
