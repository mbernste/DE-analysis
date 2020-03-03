#   Run EBSeq from Python
#

import numpy as np
import pandas as pd
from rpy2.robjects import r, pandas2ri
from rpy2.robjects.conversion import localconverter
import rpy2.robjects as ro

def runEBSeq(df_counts, conditions):
    """
    Args:
        df_counts: Mx2 dataframe where the rows are genes and 
            the columns correspond to the two conditions
        conditions: 2-length list of the condition names. For 
            example: ['treated', 'untreated']
    Returns:
        Mx3 dataframe where the rows are genes and the columns
            'PPEE', 'PPDE', and 'PostFC' provide the posterior
            probabilies of equal mean, of differential mean, and
            the posterior fold-change respectively.
    """
    with localconverter(ro.default_converter + pandas2ri.converter):
        r_cts = ro.conversion.py2rpy(df_counts)
        r_genes = ro.conversion.py2rpy(list(df_counts.index))
        r_conds = ro.conversion.py2rpy(conditions)
    rstring="""
        function(counts, genes, conds) {
            library(EBSeq)
            conditions <- as.factor(unlist(conds))
            cts <- as.matrix(counts)
            rownames(cts) <- unlist(genes)
            colnames(cts) <- colnames(counts)
            Sizes <- MedianNorm(cts)
            EBOut <- EBTest(
                Data=cts,
                Conditions=conditions,
                sizeFactors=Sizes, 
                maxround=5
            )
            EBDERes <- GetDEResults(EBOut, FDR=0.05)
            df <- data.frame(EBDERes$PPMat)
            GeneFC <- PostFC(EBOut)
            df['PostFC'] <- GeneFC$PostFC[rownames(df)]
            df
        }
    """
    r_func = ro.r(rstring)
    r_res = r_func(r_cts, r_genes, r_conds)
    with localconverter(ro.default_converter + pandas2ri.converter):
        df_res = ro.conversion.rpy2py(r_res)
    return df_res


def main():
    df_counts = pd.read_csv('./data/data.tsv', sep='\t')
    df_counts_sub = df_counts.set_index('Geneid')[['H9-D20', 'H9-D285']]
    res = runEBSeq(df_counts_sub, ['Day0', 'Day285'])
    print(res)
    res.to_csv('EBSeq_DE_H9.tsv', sep='\t')



if __name__ == "__main__":
    main()
