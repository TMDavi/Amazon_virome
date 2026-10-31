configfile: "config.yaml"
import os

rule all:
    input:
        expand("results/{sample}/iphop/Host_prediction_to_genome_m90.csv", sample=config["samples"])

rule iphop:
    input:
        #contigs = "results/{sample}/merged_viral_step2/merged_viral_contigs.fasta"
        contigs = "results/{sample}/final_filtered/filtered_genomad.fasta"
    params:
        outdir = "results/{sample}/iphop",
        db_dir = config["iphop_db_dir"]
    output:
        "results/{sample}/iphop/Host_prediction_to_genome_m90.csv"
    conda: 
        "iphop_env"
    threads: 50
    shell:
        """
        iphop predict --fa_file {input.contigs} --db_dir {params.db_dir} --out_dir {params.outdir} --num_threads {threads}
        """