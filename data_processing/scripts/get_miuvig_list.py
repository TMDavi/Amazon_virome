import pandas as pd
import argparse

def main():
    parser = argparse.ArgumentParser("Extract the miuvigs from the clustered file")
    parser.add_argument('-i', '--input', help="ANI clustered tsv table")
    parser.add_argument('-o','--output', help="The list of miuvigs")
    
    args = parser.parse_args()

    df = pd.read_csv(args.input, sep='\t', names=['miuvig_id', 'cluster'])
    miuvigs = df['miuvig_id'].to_frame()

    miuvigs.to_csv(args.output, sep='\t', index=False)

if __name__ == '__main__':  
    main()

