from Bio import SeqIO
import argparse
import os

def extract_contigs(fasta_file, file_list): 

    seen_ids = set()

    with open(fasta_file, "w") as mergedcontigs:
        for file in file_list:
            if not os.path.exists(file):
                print(f"File {file} does not exist.")
                continue

            for record in SeqIO.parse(file, "fasta"):
                ide = record.id.split("||")[0]
                if ide not in seen_ids:
                    SeqIO.write(record, mergedcontigs, "fasta")
                    seen_ids.add(ide)
                else:
                    print(f"Duplicate ID skipped: {ide}")

    print(f"\nTotal unique headers written: {len(seen_ids)}")

def main():
    parser = argparse.ArgumentParser(description="combine different fasta files with repeated ids")
    parser.add_argument("-o","--outdir", help="Output file with unique headers.")
    parser.add_argument("-f","--file_list", nargs='+', help="Virsorter2 file")
    
    args = parser.parse_args()

    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)

    extract_contigs(f'{args.outdir}/merged_viral_contigs.fasta', args.file_list)

if __name__ == "__main__":
    main()