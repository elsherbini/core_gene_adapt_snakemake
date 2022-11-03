import pandas as pd

genome_info = pd.read_csv("input/kwon_lab_genomes_busco_sgb.csv")
target_taxa = ["Lactobacillus_crispatus", "Lactobacillus_iners", "Lactobacillus_gasseri", "Lactobacillus_jensenii", "Lactobacillus_mulieris"]
ALL_GENOMES = [g.replace("_genome.fna", "") for g in list(genome_info[(genome_info.species_name.isin(target_taxa)) & (genome_info.busco_single_copy > 85)].genome)]

lactobacillales_hmm_ids = glob_wildcards("input/lactobacillales_odb9/hmms/{hmm_id}.hmm").hmm_id

def get_orf_files_from_group(wildcards):
    if wildcards.group == "noniners":
        genome_group = list(genome_info[(genome_info.species_name.isin(["Lactobacillus_crispatus", "Lactobacillus_gasseri", "Lactobacillus_jensenii", "Lactobacillus_mulieris"])) & (genome_info.busco_single_copy > 85)].genome)
    elif wildcards.group == "iners":
        genome_group = list(genome_info[(genome_info.species_name.isin(["Lactobacillus_iners"])) & (genome_info.busco_single_copy > 85)].genome)
    if wildcards.type == "aa":
        return ["temp/rename_fasta/{}.faa".format(g.replace("_genome.fna", ""), g.replace("_genome.fna", "")) for g in genome_group]
    elif wildcards.type == "nt":
        return ["temp/rename_fasta/{}.ffn".format(g.replace("_genome.fna", ""), g.replace("_genome.fna", "")) for g in genome_group]

def get_hmm_files_from_taxon(wildcards):
    if wildcards.taxon == "firmicutes":
        return expand("input/firmicutes_odb9/hmms/{hmm_id}.hmm", hmm_id=firmicutes_hmm_ids)
    elif wildcards.taxon == "lactobacillales":
        return expand("input/lactobacillales_odb9/hmms/{hmm_id}.hmm", hmm_id=lactobacillales_hmm_ids)

def get_group_genomes(wildcards):
    if wildcards.group == "noniners":
        return ["/n/groups/kwon/joseph/projects/femmcat_assemblies_20220726/input/{}_genome.fna".format(g.replace("_genome.fna", "")) for g in list((genome_info[(genome_info.species_name.isin(["Lactobacillus_crispatus", "Lactobacillus_gasseri", "Lactobacillus_jensenii", "Lactobacillus_mulieris"])) & (genome_info.busco_single_copy > 85)]).groupby("species_name").head(10).reset_index(drop=True).genome)]
    elif wildcards.group == "iners":
        return ["/n/groups/kwon/joseph/projects/femmcat_assemblies_20220726/input/{}_genome.fna".format(g.replace("_genome.fna", "")) for g in list(genome_info[(genome_info.species_name.isin(["Lactobacillus_iners"])) & (genome_info.busco_single_copy > 85)].genome)]

localrules: target, make_concatenated_orf_file, cat_hmm_file, cat_group_genomes

rule target:
    input:
        expand("output/adapt/{group}/{gene}_{taxon}_guides.tsv", group = ["iners","noniners"], taxon = ["lactobacillales"], gene = lactobacillales_hmm_ids)

rule download_bakta:
    output:
        directory("bakta_db")
    conda:
        "envs/bakta.yaml"
    shell:
        "bakta_db download --output {output}"

rule annotate_bakta:
    input:
        genome_fasta = "/n/groups/kwon/joseph/projects/femmcat_assemblies_20220726/input/{genome}_genome.fna",
        db = "bakta_db"
    output:
        directory("output/annotations/{genome}")
    conda:
        "envs/bakta.yaml"
    threads:
        8
    shell:
        "bakta {input.genome_fasta} --db {input.db}/db --threads {threads} --output {output} --prefix {wildcards.genome} --keep-contig-headers"

rule align_genes:
    input:
        "temp/hmmer/{group, \w+}-{taxon}/{group}-{taxon}_dummy_file.txt"
    output:
        "temp/hmmer/{group, \w+}-{taxon, \w+}/aligned/{gene_name}_aligned.fasta"
    conda:
        "envs/mafft.yaml"
    threads:
        8
    shell:
        "mkdir -p temp/hmmer/{wildcards.group}-{wildcards.taxon}/aligned/; mafft-linsi --thread {threads} temp/hmmer/{wildcards.group}-{wildcards.taxon}/unaligned/{wildcards.gene_name}_unaligned.fasta > {output}"

rule make_unaligned_gene_files:
    input:
        singletons = "temp/hmmer/hmmer_out/{group, \w+}-{taxon}_singletons.txt",
        concatenated_orfs = "temp/hmmer/{group, \w+}_concatenated_orf.nt.fasta"
    output:
        "temp/hmmer/{group,\w+}-{taxon,\w+}/{group}-{taxon}_dummy_file.txt"
    conda:
        "envs/biopython.yaml"
    shell:
        "python scripts/make_phylo_tree/make_ribo_files.py -m {input[singletons]} -c {input[concatenated_orfs]} -o temp/hmmer/{wildcards.group}-{wildcards.taxon}/unaligned/ -d {output} "

rule filter_hmm_results:
    input:
        "temp/hmmer/hmmer_out/{group, \w+}-{taxon}_result.tbl"
    output:
        "temp/hmmer/hmmer_out/{group, \w+}-{taxon}_singletons.txt"
    conda:
        "envs/biopython.yaml"
    shell:
        "python scripts/make_phylo_tree/check_hmmsearch.py -i {input} -s {output}"

rule cat_hmm_file:
    input:
        get_hmm_files_from_taxon
    output:
        "temp/hmmer/cat_hmm_{taxon}.hmm"
    shell:
        "cat {input} > {output}"

rule hmmsearch:
    input:
        concatenated_orfs = "temp/hmmer/{group, \w+}_concatenated_orf.aa.fasta",
        hmmfile = "temp/hmmer/cat_hmm_{taxon}.hmm"
    output:
        text = "temp/hmmer/hmmer_out/{group, \w+}-{taxon}_result.txt",
        table = "temp/hmmer/hmmer_out/{group, \w+}-{taxon}_result.tbl",
        domains = "temp/hmmer/hmmer_out/{group, \w+}-{taxon}_result.dom.tbl"
    conda:
        "envs/hmmer.yaml"
    threads:
        8
    shell:
        "hmmsearch --cpu {threads} -o {output[text]} --tblout {output[table]} --domtblout {output[domains]} {input[hmmfile]} {input[concatenated_orfs]}"


rule rename_fastas_faa:
    input:
        fasta = "output/annotations/{isolate}/{isolate}.faa"
    output:
        fasta = "temp/rename_fasta/{isolate}.faa"
    conda:
        "envs/biopython.yaml"
    script:
        "scripts/rename_fastas.py"

rule rename_fastas_ffn:
    input:
        fasta = "output/annotations/{isolate}/{isolate}.ffn"
    output:
        fasta = "temp/rename_fasta/{isolate}.ffn"
    conda:
        "envs/biopython.yaml"
    script:
        "scripts/rename_fastas.py"

rule make_concatenated_orf_file:
    input:
        get_orf_files_from_group
    output:
        "temp/hmmer/{group, \w+}_concatenated_orf.{type}.fasta"
    shell:
        "cat {input} > {output}"

rule cat_group_genomes:
    input:
        get_group_genomes
    output:
        "temp/{group}_genomes_for_adapt.fna"
    shell:
        "cat {input} > {output}"


def get_adapt_inputs_noniners(wildcards):
    return {
        "alignment": "temp/hmmer/noniners-{}/aligned/{}_aligned.fasta".format(wildcards.taxon, wildcards.gene),
        "dont_hit": "temp/{}_genomes_for_adapt.fna".format("iners")
    }

def get_adapt_inputs_iners(wildcards):
    return {
        "alignment": "temp/hmmer/iners-{}/aligned/{}_aligned.fasta".format(wildcards.taxon, wildcards.gene),
        "dont_hit": "temp/{}_genomes_for_adapt.fna".format("noniners")
    }

rule run_adapt_noniners:
    input:
        unpack(get_adapt_inputs_noniners)
    output:
        "output/adapt/noniners/{gene}_{taxon, \w+}_guides.tsv"
    conda:
        "envs/adapt.yaml"
    shell:
        "design.py complete-targets fasta {input.alignment} -o output/adapt/noniners/{wildcards.gene}_{wildcards.taxon}_guides --obj maximize-activity --predict-cas13a-activity-model --best-n-targets 5 --verbose --specific-against-fastas {input.dont_hit}"


rule run_adapt_iners:
    input:
        unpack(get_adapt_inputs_iners)
    output:
        "output/adapt/iners/{gene}_{taxon, \w+}_guides.tsv"
    conda:
        "envs/adapt.yaml"
    shell:
        "design.py complete-targets fasta {input.alignment} -o output/adapt/iners/{wildcards.gene}_{wildcards.taxon}_guides --obj maximize-activity --predict-cas13a-activity-model --best-n-targets 5 --verbose --specific-against-fastas {input.dont_hit}"
