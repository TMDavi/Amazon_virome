import pandas as pd
import argparse
import re
import os
from crass_extract import parse_hmm, parse_domout, extract_genes

def filter_results(input_file, threshold_file, output_file):
    """
    Filters the results based on e-value and score.
    """

    thresholds = pd.read_csv(threshold_file, sep="\t")
    print("Filtering results...")
    df = pd.read_csv(input_file, sep="\t")
    df = df[(df['evalue'] < 1e-5) & (df['score'] > 70)]

    mask = []  
    for index, row in df.iterrows():
        gene = row['qname'].split('_')[0]
        if gene in thresholds['Gene'].values:
            cutoff = thresholds.loc[thresholds['Gene'] == gene, 'ali_cut'].values[0]

            mask.append(row['ali_len'] >= cutoff)
        else:
            mask.append(True)
    df = df[mask]


    df['contig_name'] = df['contig_id'].astype(str).str.replace(r'_\d+$', '', regex=True)
    df = df.sort_values(by='score', ascending=False)
    df = df.drop_duplicates(subset=['contig_name', 'qname'], keep='first')
    df = df.drop('contig_name',axis=1)
    df.to_csv(output_file, sep="\t", index=False)

def count_gene_hits(input_file, output_filtered_contigs_path=""):
    """
    Counts the number of gene hits per contig and filters those with more than 4 hits.
    """
    print("Counting gene hits...")
    df = pd.read_csv(input_file, sep="\t")
    df['contig_id'] = df['contig_id'].astype(str).apply(lambda x: re.sub(r'_\d+$', '', x))
    gene_counts = df.groupby('contig_id').size().reset_index(name='gene_hits')
    filtered = gene_counts[gene_counts['gene_hits'] > 1]
    filtered = filtered.sort_values(by='gene_hits', ascending=False)

    if output_filtered_contigs_path=="":
        pass
    else:
        filtered.to_csv(output_filtered_contigs_path, sep="\t", index=False)

    filtered = filtered['contig_id'].tolist()
    print(f"Filtered contigs with more than 2 gene hit: {len(filtered)}")
    return filtered

def main():
    parser = argparse.ArgumentParser(description="Extract putative Crassphage sequences from HMMER results.")
    parser.add_argument("-i", "--input", required=True, help="Input FASTA file with contigs.")
    parser.add_argument("-t", "--threads", type=int, default=1, help="Number of threads for HMMER search.")
    parser.add_argument("-o", "--output_dir", required=True, help="Output folder for filtered results.")
    args = parser.parse_args()

    threshold_file = "/MP_Data/Dados_David/viral_analysis/viral_HMMs/giantDB/cutoffs.tsv"
    if not os.path.exists(threshold_file):
        raise FileNotFoundError(f"Threshold file {threshold_file} does not exist.")

    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)

    ncldv_db = "/MP_Data/Dados_David/viral_analysis/viral_HMMs/giantDB/NCLDV_gvogs.hmm"

    parse_hmm(args.input, args.threads, ncldv_db, f"{args.output_dir}/hmmer_results.domout")
    parse_domout(f"{args.output_dir}/hmmer_results.domout", f"{args.output_dir}/parsed_results.txt")
    filter_results(f"{args.output_dir}/parsed_results.txt", threshold_file,f"{args.output_dir}/filtered_results.txt")
    extract_genes(f"{args.output_dir}/filtered_results.txt", 
                  args.input, 
                  count_gene_hits(f"{args.output_dir}/filtered_results.txt"),
                  f"{args.output_dir}/extracted_genes.faa")
    count_gene_hits(f"{args.output_dir}/filtered_results.txt", f"{args.output_dir}/putative_ncldv_contig_list.txt")

if __name__ == "__main__":
    main()