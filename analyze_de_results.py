import matplotlib as mpl
mpl.use('Agg')
import seaborn as sns
from matplotlib import pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
import upsetplot
from upsetplot import from_contents
import gseapy as gp
from os.path import join
import sys

PPDE_THRESH = 0.95
GSEA_THRESH = 0.05

def main():
    dirr = sys.argv[1]

    df_102 = pd.read_csv(
        join(dirr, 'EBSeq_DE_102.tsv'),
        index_col=0,
        sep='\t'
    )	
    df_106 = pd.read_csv(
        join(dirr, 'EBSeq_DE_106.tsv'),
        index_col=0,
        sep='\t'
    )
    df_h9 = pd.read_csv(
        join(dirr, 'EBSeq_DE_H9.tsv'),
        index_col=0,
        sep='\t'
    )

    # Filter genes by PPDE
    df_102 = _filt(df_102)
    df_106 = _filt(df_106)
    df_h9 = _filt(df_h9)
    condition_to_de_genes = {
        'H9': set(df_h9.index),
        '102': set(df_102.index),
        '106': set(df_106.index)
    }

    # Save gene lists and fold-changes to files
    df_102['PostFC'].sort_values(ascending=False).to_csv(
            join(dirr, 'DE_genes_102.tsv'), 
            sep='\t'
    )
    df_106['PostFC'].sort_values(ascending=False).to_csv(
        join(dirr, 'DE_genes_106.tsv'), 
        sep='\t'
    )
    df_h9['PostFC'].sort_values(ascending=False).to_csv(
        join(dirr, 'DE_genes_H9.tsv'), 
        sep='\t'
    )

    # Create Upset plot for DE genes
    memb = from_contents(condition_to_de_genes, id_column='gene')
    upsetplot.plot(memb, subset_size='count')
    plt.savefig(join(dirr, 'DE_upset.png'), format='png')

    _gsea_analysis(df_h9, df_102, df_106, 'KEGG_2016', dirr)
    _gsea_analysis(df_h9, df_102, df_106, 'GO_Biological_Process_2018', dirr)


def _gsea_analysis(df_h9, df_102, df_106, gene_sets, dirr):
    # Perform GSEA
    terms_102 = _gsea(list(df_102.index), gene_sets)
    terms_106 = _gsea(list(df_106.index), gene_sets)
    terms_h9 = _gsea(list(df_h9.index), gene_sets)
    condition_to_terms = {
        'H9': terms_h9,
        '102': terms_102,
        '106': terms_106
    }

    # Create Upset plot for KEGG terms
    memb = from_contents(condition_to_terms, id_column='gene')
    memb.to_csv(
        join(dirr, '{}_GSEA.tsv'.format(gene_sets)), 
        sep='\t'
    )
    upsetplot.plot(memb, subset_size='count')
    plt.savefig(
        join(dirr, '{}_upset.png'.format(gene_sets)), 
        format='png'
    )


def _filt(df):
    df_filt = df.dropna(subset=['PPDE', 'PPEE'])
    df_filt = df.loc[df['PPDE'] > PPDE_THRESH]
    df_filt['log2_PostFC'] = np.log2(df['PostFC']+1)
    return df_filt


def _gsea(genes, gene_set):
    enr = gp.enrichr(
        gene_list=genes,
        gene_sets=[gene_set],
        no_plot=True,
        cutoff=0.05
    )
    enr.results = enr.results[enr.results["Adjusted P-value"] < GSEA_THRESH]
    print(enr.results)
    sig_terms = set(enr.results['Term'])
    return sig_terms

if __name__ == "__main__":
    main()
