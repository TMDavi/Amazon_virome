configfile: 'config.yaml'

rule all:
    input:
        #"results/clustering/all_viral_contigs.nhr",
        "results/clustering/all_viral_contigs_blastn.tsv",
        "results/miuvigs.fasta",
        expand("results/{sample}/miuvigs/miuvigs.fasta", sample=config['samples'].keys())

rule merge_contigs:
    input:
        phages="results/{sample}/final_filtered/filtered_genomad.fasta",
        crass="results/{sample}/crass_id/putative_crass_contigs.fasta",
        ncldv="results/{sample}/ncldv_id/putative_ncldv_contigs.fasta",
        virophage="results/{sample}/virophage_id/putative_virophage_contigs.fasta",
        plv="results/{sample}/virophage_id/putative_plv_contigs.fasta"
    output:
        "results/{sample}/merged_viral_step2/merged_viral_contigs.fasta"
    params:
        outdir = "results/{sample}/merged_viral_step2"
    shell:
        """
        python scripts/mergecontigs.py -o {params.outdir} -f {input.phages} {input.crass} {input.ncldv} {input.virophage} {input.plv}
        """

rule rename_contigs:
    input:
        assembly = "results/{sample}/merged_viral_step2/merged_viral_contigs.fasta"
    output:
        temp("results/{sample}/merged_viral_step2/merged_viral_contigs_renamed.fasta")
    params:
        sample_name = lambda wildcards: wildcards.sample
    shell:
        "python scripts/rename_contigs.py -i {input.assembly} -o {output} -p {params.sample_name}"

rule concatenate_contigs:
    input:
        expand("results/{sample}/merged_viral_step2/merged_viral_contigs_renamed.fasta", sample=list(config["samples"].keys()))
    output:
        "results/clustering/concatenated_merged_viral_contigs.fasta"
    shell:
        "cat {input} > {output}"

rule makeblastdb:
    input:
        "results/clustering/concatenated_merged_viral_contigs.fasta"
    output:
        "results/clustering/all_viral_contigs.nhr"
    params:
        db_name = "results/clustering/all_viral_contigs"
    conda:
        "viral-id-sop"
    shell:
        """
        makeblastdb -in {input} -dbtype nucl -out {params.db_name}
        """

rule blastn:
    input:
        query = "results/clustering/concatenated_merged_viral_contigs.fasta",
        nhr = "results/clustering/all_viral_contigs.nhr"
    output:
        "results/clustering/all_viral_contigs_blastn.tsv"
    params:
        db = "results/clustering/all_viral_contigs"
    conda:
        "/MP_Data/mambaforge/envs/viral-id-sop"
    threads: 10
    shell:
        """
        blastn -query {input.query} -db {params.db} -outfmt '6 std qlen slen' -max_target_seqs 10000 -out {output} -num_threads {threads}
        """

rule anicalc:
    input:
        "results/clustering/all_viral_contigs_blastn.tsv"
    output:
        "results/clustering/all_viral_contigs_ani.tsv"
    shell:
        """
        python programs/anicalc.py -i {input} -o {output}
        """

rule aniclust:
    input:
        fasta = "results/clustering/concatenated_merged_viral_contigs.fasta",
        ani = "results/clustering/all_viral_contigs_ani.tsv"
    output:
        "results/clustering/all_viral_contigs_clusters.tsv"
    shell:
        """
        python programs/aniclust.py --fna {input.fasta} --ani {input.ani} --out {output} --min_ani 95 --min_tcov 85 --min_qcov 0
        """

rule get_miuvig_list:
    input:
        "results/clustering/all_viral_contigs_clusters.tsv"
    output:
        "results/clustering/miuvig_list.tsv"
    shell:
        """
        python scripts/get_miuvig_list.py -i {input} -o {output}
        """

rule filter_miuvigs:
    input:
        header_file = "results/clustering/miuvig_list.tsv", 
        fasta_file = "results/clustering/concatenated_merged_viral_contigs.fasta"
    params:
        temp_dir_parent = "results/clustering"
    output:
        "results/miuvigs.fasta"
    shell:
        """
        python scripts/filter_contigs.py -i {input.fasta_file} -hd {input.header_file} -o {output} --temp_dir_parent {params.temp_dir_parent} --num_chunks 1 -t 1
        """
rule redistribute_miuvigs:
    input:
        "results/miuvigs.fasta"
    output:
        expand("results/{sample}/miuvigs/miuvigs.fasta", sample=config['samples'].keys())
    params:
        samplelist = list(config["samples"].keys())
    shell:
        """
        python scripts/distribute_miuvigs.py -m {input} -s {params.samplelist}  
        """
    