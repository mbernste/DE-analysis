import pandas as pd
import sys

def main():
    in_f = sys.argv[1]
    out_f = sys.argv[2]

    df = pd.read_csv(in_f, sep='\t')
    df = df[[
        'Geneid', 
        'H9-D20',   
        'H9-D285',
        '102-D20',
        '102-D285',  
        '106-D20',
        '106-D100'
    ]]
    df.set_index('Geneid')
    df.to_csv(out_f, sep='\t')

if __name__ == "__main__":
    main()
