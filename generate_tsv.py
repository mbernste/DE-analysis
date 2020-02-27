import pandas as pd

SAMPLE_ID_TO_CATEGORY = {
    'H9-D20': 'H9',
    'H9-D285': 'H9',
    '102-D20': 'Patient',	
    '102-D285': 'Patient',
    '106-D20': 'Patient',
    '106-D100': 'Patient'
}

def main():
    df = pd.read_csv('./raw_data/raw_data.tsv', sep='\t')
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
    df.to_csv('./data/data.tsv', sep='\t')

    
    

if __name__ == "__main__":
    main()
