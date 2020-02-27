import numpy as np
import pandas as pd
import rpy2.robjects as robjects
from rpy2.robjects.packages import importr
from rpy2.robjects import r, pandas2ri
from rpy2.robjects.conversion import localconverter
import rpy2.robjects as ro

def main():
    df_counts = pd.read_csv('./data/data.tsv', sep='\t')
    df_counts_sub = df_counts[['102-D20', '102-D285', '106-D20', '106-D100']]
    print(df_counts)
    df_meta = pd.DataFrame(
        data=['Before', 'After', 'Before', 'After'],
        index=['102-D20', '102-D285', '106-D20', '106-D100'],
        columns=['condition']
    )
    print(df_meta)
    res = runDESeq2(df_counts_sub, df_meta)
    print(res)
    pvalues = list(subset_RS4(res, 'padj'))
    logfolds = list(subset_RS4(res, 'log2FoldChange'))
    
    df_res = pd.DataFrame(
        data=[
            (pv, lf)
            for pv, lf in zip(pvalues, logfolds)
        ],
        index=df_counts['Geneid'],
        columns=['p-value', 'log2 fold-change']
    )
    df_res = df_res.dropna()
    print(df_res)

    df_res = df_res.loc[df_res['p-value'] < 0.05]
    #df_res = df_res.loc[df_res['log2 fold-change'] >= 1]
    print(df_res)

    #print(df_res.keys())

    #df_res = df_res.dropna(subset=['pvalue'])
    #print(df_res)


def runDESeq2(df_counts, df_meta):
    with localconverter(ro.default_converter + pandas2ri.converter):
        r_cts = ro.conversion.py2rpy(df_counts)
        r_coldata = ro.conversion.py2rpy(df_meta)
    rstring="""
        function(cts, coldata) {
            library(DESeq2)
            dds <- DESeqDataSetFromMatrix(
                countData = cts,
                colData = coldata,
                design = ~ condition
            )
            dds <- DESeq(dds)
            res <- results(dds)
            res
        }
    """
    re_deseq_func = robjects.r(rstring)
    r_deseq_res = re_deseq_func(r_cts, r_coldata)
    with localconverter(ro.default_converter + pandas2ri.converter):
        df_res = ro.conversion.rpy2py(r_deseq_res)
    return df_res


def subset_RS4(rs4, subset):
    subset_func = r("""function(o, s){
        o[[s]]
    }
    """)
    r_res = subset_func(rs4, subset)
    with localconverter(ro.default_converter + pandas2ri.converter):
        res = ro.conversion.rpy2py(r_res)
    return res


if __name__ == "__main__":
    main()
