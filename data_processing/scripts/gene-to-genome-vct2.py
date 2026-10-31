from Bio import SeqIO
import re
import argparse

def genetogenome(faafile, outfile):

    with open(outfile, "w") as out:
        out.write("contig_id,genome_id,keywords\n")
        for record in SeqIO.parse(faafile, "fasta"):
            gene_id = record.id
            genome_id = re.sub(r'_\d+$', '', gene_id)
            keywords = " "
            out.write(f"{gene_id},{genome_id},{keywords}\n")

def main():
    parser = argparse.ArgumentParser("produce the gene-to-genome table for VContact2")
    parser.add_argument("-i", "--input", help="Prodigal file")
    parser.add_argument("-o", "--output", help="Gene-to_genome table")
    args = parser.parse_args()

    genetogenome(args.input, args.output)

if __name__ == "__main__":
    main()