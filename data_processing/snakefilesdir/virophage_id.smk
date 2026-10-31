configfile: 'config.yaml'

rule all:
    input:
        expand("results/{sample}/virophage_plv_id/putative_virophage_contigs.fasta", sample=config["samples"].keys()),
        expand("results/{sample}/virophage_plv_id/putative_plv_contigs.fasta", sample=config["samples"].keys())
       

rule run_virophage_affi:
    input:
        assembly = lambda wildcards: os.path.abspath(config["samples"][wildcards.sample])
    output:
        "results/{sample}/virophage_plv_id/MP_Data_final_affiliation.tsv"
    params:
        outdir = "results/{sample}/virophage_plv_id/"
    conda:
        "virophage_classification"
    shell:
        """
        ./programs/ICTV_VirophageSG/run_virophage_affi.pl -i {input.assembly} -o {params.outdir}
        """

rule get_virophage_and_plv_lists:
    input:
        "results/{sample}/virophage_plv_id/MP_Data_final_affiliation.tsv"
    output:
        "results/{sample}/virophage_plv_id/putative_virophage_list.tsv",
        "results/{sample}/virophage_plv_id/putative_plv_list.tsv"
    params:
        outdir = "results/{sample}/virophage_plv_id/"
    conda:
        "virophage_classification"
    shell:
        """
        python scripts/virophage_process.py -i {input} -o {params.outdir} 
        """

rule filter_virophage_contigs:
    input:
        header_file = "results/{sample}/virophage_plv_id/putative_virophage_list.tsv", 
        fasta_file = lambda wildcards: os.path.abspath(config["samples"][wildcards.sample])
    params:
        temp_dir_parent = "results/{sample}/virophage_plv_id"
    output:
        "results/{sample}/virophage_plv_id/putative_virophage_contigs.fasta"
    threads: 20
    shell:
        """
        python scripts/filter_contigs.py -i {input.fasta_file} -hd {input.header_file} -o {output} --temp_dir_parent {params.temp_dir_parent} --num_chunks 20 -t {threads}
        """
rule filter_plv_contigs:
    input:
        header_file = "results/{sample}/virophage_plv_id/putative_plv_list.tsv", 
        fasta_file = lambda wildcards: os.path.abspath(config["samples"][wildcards.sample])
    params:
        temp_dir_parent = "results/{sample}/virophage_plv_id"
    output:
        "results/{sample}/virophage_plv_id/putative_plv_contigs.fasta"
    threads: 20
    shell:
        """
        python scripts/filter_contigs.py -i {input.fasta_file} -hd {input.header_file} -o {output} --temp_dir_parent {params.temp_dir_parent} --num_chunks 20 -t {threads}
        """
rule merge_contigs:
    input:
        viro = "results/{sample}/virophage_plv_id/putative_virophage_contigs.fasta",
        plv = "results/{sample}/virophage_plv_id/putative_plv_contigs.fasta"
    output:
        "results/{sample}/virophage_plv_id/putative_virophage_plv_contigs.fasta"
    shell:
        """
        cat {input.viro} {input.plv} >> {output}
        """

rule checkv:
    input:
        assembly = "results/{sample}/virophage_plv_id/putative_virophage_plv_contigs.fasta"
    output:
        "results/{sample}/virophage_plv_id/checkv/viruses.fna"
    params:
        db = "/MP_Data/database/checkv-db-v1.5",
        outdir = "results/{sample}/virophage_plv_id/checkv/"
    conda:
         "viral-id-sop"
    threads: 20
    shell:
        """
        checkv end_to_end {input} {params.outdir} -t {threads} -d {params.db}
        """

rule genomad:
    input:
        filtered = "results/{sample}/virophage_plv_id/checkv/viruses.fna"
    output:
        "results/{sample}/virophage_plv_id/genomad/viruses_annotate/viruses_genes.tsv"
    params:
        outdir = "results/{sample}/virophage_plv_id/genomad/",
        db = "/MP_Data/database/genomad_db" #adcionar no config 
    conda:
         "genomad_env"
    threads: 40
    shell:
        """
        genomad end-to-end --cleanup {input.filtered} {params.outdir} {params.db} --threads {threads}
        """