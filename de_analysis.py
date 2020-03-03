import pandas as pd
import sys

import run_ebseq as re

def main():
    in_f = sys.argv[1]
    out_dir = sys.argv[2]
    df_counts = pd.read_csv('./data/data.tsv', sep='\t')
    _run_de('H9-D20', 'H9-D285', df_counts, join(out_dir, 'EBSeq_DE_H9.tsv'))
    _run_de('102-D20','102-D285', df_counts, join(out_dir, 'EBSeq_DE_102.tsv')
    _run_de('106-D20', '106-D100', df_counts, join(out_dir, 'EBSeq_DE_106.tsv'))

def _run_de(sample1, sample2, df_counts, out_f):
    df_counts_sub = df_counts.set_index('Geneid')[[sample1, sample2]]
    res = re.runEBSeq(df_counts_sub, ['time_1', 'time_2'])
    res.to_csv(out_f, sep='\t')

if __name__ == "__main__":
    main()
