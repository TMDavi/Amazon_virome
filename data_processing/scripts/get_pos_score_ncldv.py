import argparse
import pandas as pd

def main():
    parser = argparse.ArgumentParser("Get the contigs with viralrecall positive scores")
    parser.add_argument('-i', '--input', required=True, help='Viral Recall sumary table')
    parser.add_argument('-o', '--output', required=True, help='Contig list with contigs with positive scores')
    args = parser.parse_args()

    df = pd.read_csv(args.input, sep='\t')
    df = df[df['score'] > 0]
    df = df.rename(columns={'replicon': 'contig_id'})
    ncldv = df['contig_id'].to_frame()
    ncldv.to_csv(args.output, sep='\t', index=False)

if __name__ == '__main__':
    main()