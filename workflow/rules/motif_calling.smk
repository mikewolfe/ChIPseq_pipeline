## HELPER FUNCTIONS inherited from parent SnakeFile:
# samples(pep)
# lookup_sample_metadata(sample, key, pep)

# GLOBAL CONSTANTS inherited from parent Snakefile:
# RES - resolution of coverage
# WITHIN - normalization to perform within samples
 
rule clean_motif_calling:
    shell:
        "rm -fr results/motif_calling"


def determine_motif_calling_files(config, pep):
    outfiles = []
    for model in lookup_in_config(config, ["motif_calling"], []):
        motif_caller = lookup_in_config(config, ["motif_calling", model, "motif_caller"], 
        err = "Need motif caller specified for motif_caller model %s in config file. I.e. \nmotif_calling:\n\t%s:\n\t\tmotif_caller: 'macs2'"%(model,model))
        these_samples = filter_samples(pep, \
        lookup_in_config(config, ["motif_calling", model, "filter"], "not input_sample.isnull()"))
        if motif_caller == "streme":
            for sample in these_samples:
                outfiles.append("results/motif_calling/%s/fimo/%s.tsv"%(model, sample))
        if motif_caller == "meme":
            for sample in these_samples:
                outfiles.append("results/motif_calling/%s/fimo/%s.tsv"%(model,sample))
    return outfiles

rule run_motif_calling:
    input:
        determine_motif_calling_files(config, pep)


def determine_filesig_input(sample, model, config, pep):
   file_sig = lookup_in_config(config, ["motif_calling", model, "filesignature"])
   return file_sig%(sample)

rule get_peak_seqs:
    input: 
        bed=lambda wildcards: determine_filesig_input(wildcards.sample, wildcards.model, config, pep),
        fasta = lambda wildcards: "results/alignment/combine_fasta/%s/%s.fa"%(\
        lookup_sample_metadata(wildcards.sample, "genome", pep),\
        lookup_sample_metadata(wildcards.sample, "genome", pep))

    output:
        "results/motif_calling/{model}/get_peak_seqs/{sample}_peak_seqs.fa"
    log:
        stdout="results/motif_calling/logs/{model}/get_peak_seqs/{sample}.log",
        stderr="results/motif_calling/logs/{model}/get_peak_seqs/{sample}.err"
    params:
        get_peak_seqs_param_string = lambda wildcards: lookup_in_config_persample(config,\
        pep, ["motif_calling", wildcards.model, "get_peak_seqs_param_string"], wildcards.sample,\
        "--min_size 30 ") 
    conda:
        "../envs/motif_calling.yaml"
    shell:
        "python3 workflow/scripts/bed_to_fasta.py {input.bed} {input.fasta} {output} "
        "--meme_header "
        "{params.get_peak_seqs_param_string} "
        "> {log.stdout} 2> {log.stderr}"


rule get_markov:
    input:
        infasta = "results/alignment/combine_fasta/{genome}/{genome}.fa"
    output:
        "results/motif_calling/{model}/get_markov/{genome}.txt"
    log:
        stdout="results/motif_calling/logs/{model}/get_markov/{genome}.log",
        stderr="results/motif_calling/logs/{model}/get_markov/{genome}.err"
    conda:
        "../envs/motif_calling.yaml"
    params: 
        meme_param_string = lambda wildcards: lookup_in_config(config,\
        ["motif_calling", wildcards.model, "get_markov_param_string"],\
        "-dna -m 0 ") 
    shell:
        "fasta-get-markov {input.infasta} {output} > {log.stdout} 2> {log.stderr}"


rule streme_call_motifs:
    input:
        inseqs ="results/motif_calling/{model}/get_peak_seqs/{sample}_peak_seqs.fa",
        bg = lambda wildcards: "results/motif_calling/{model}/get_markov/%s.txt"%(\
        lookup_sample_metadata(wildcards.sample, "genome", pep))
    output:
        "results/motif_calling/{model}/streme/{sample}/sequences.tsv",
        "results/motif_calling/{model}/streme/{sample}/streme.txt",
        "results/motif_calling/{model}/streme/{sample}/streme.html"
    log:
        stdout="results/motif_calling/logs/{model}/streme/{sample}_streme.log",
        stderr="results/motif_calling/logs/{model}/streme/{sample}_streme.err"
    params:
        streme_param_string = lambda wildcards: lookup_in_config_persample(config,\
        pep, ["motif_calling", wildcards.model, "streme_param_string"], wildcards.sample,\
        "--dna --minw 8 --maxw 30") 
    conda:
        "../envs/motif_calling.yaml"
    shell:
        "streme --p {input.inseqs} "
        "--oc results/motif_calling/{wildcards.model}/streme/{wildcards.sample}/ "
        "--bfile {input.bg} "
        "{params.streme_param_string} "
        "> {log.stdout} 2> {log.stderr}"


rule rename_streme_output:
    input:
        motifs="results/motif_calling/{model}/streme/{sample}/streme.txt",
        seqs= "results/motif_calling/{model}/streme/{sample}/sequences.tsv"
    output:
        outmotifs="results/motif_calling/{model}/streme/renamed_output/{sample}.txt",
        outseqs="results/motif_calling/{model}/streme/renamed_output/{sample}.tsv"

    threads: 1
    shell:
        "cp {input.motifs} {output.outmotifs} && cp {input.seqs} {output.outseqs}"

rule meme_call_motifs:
    input:
        inseqs = "results/motif_calling/{model}/get_peak_seqs/{sample}_peak_seqs.fa",
        bg = lambda wildcards: "results/motif_calling/{model}/get_markov/%s.txt"%(\
        lookup_sample_metadata(wildcards.sample, "genome", pep))
    output:
        "results/motif_calling/{model}/meme/{sample}/meme.html",
        "results/motif_calling/{model}/meme/{sample}/meme.txt"
    log:
        stdout="results/motif_calling/logs/{model}/meme/{sample}_meme.log",
        stderr="results/motif_calling/logs/{model}/meme/{sample}_meme.err"
    params:
        meme_param_string = lambda wildcards: lookup_in_config_persample(config,\
        pep, ["motif_calling", wildcards.model, "meme_param_string"], wildcards.sample,\
        "-dna -revcomp -minw 8 -maxw 32 -mod anr -nmotifs 2") 
    conda:
        "../envs/motif_calling.yaml"
    shell:
        "meme {input.inseqs} "
        "-oc results/motif_calling/{wildcards.model}/meme/{wildcards.sample}/ "
        "-bfile {input.bg} "
        "{params.meme_param_string} "
        "> {log.stdout} 2> {log.stderr}"

rule rename_meme_output:
    input:
        motifs="results/motif_calling/{model}/meme/{sample}/meme.txt"
    output:
        outmotifs="results/motif_calling/{model}/meme/renamed_output/{sample}.txt"
    threads: 1
    shell:
        "cp {input.motifs} {output.outmotifs} "

def find_motif_file(sample, model, motif_caller):
    if motif_caller == "meme":
        out = "results/motif_calling/%s/meme/renamed_output/%s.txt"%(model, sample)
    elif motif_caller == "streme":
        out = "results/motif_calling/%s/streme/renamed_output/%s.txt"%(model, sample)
    else:
        raise ValueError("Don't know where motif file is for sample %s from motif caller %s"%(sample, motif_caller))
    return out
    

rule fimo_find_motifs:
    input:
        motif_file = lambda wildcards: find_motif_file(wildcards.sample, wildcards.model,\
        lookup_in_config(config, ["motif_calling", wildcards.model, "motif_caller"])),
        fasta = lambda wildcards: "results/alignment/combine_fasta/%s/%s.fa"%(\
        lookup_sample_metadata(wildcards.sample, "genome", pep),\
        lookup_sample_metadata(wildcards.sample, "genome", pep)),
        bg = lambda wildcards: "results/motif_calling/{model}/get_markov/%s.txt"%(\
        lookup_sample_metadata(wildcards.sample, "genome", pep))
    output:
        "results/motif_calling/{model}/fimo/{sample}.tsv"
    log:
        stderr = "results/motif_calling/logs/{model}/fimo/{sample}.err"
    conda:
        "../envs/motif_calling.yaml"
    params:
        fimo_param_string = lambda wildcards: lookup_in_config_persample(config,\
        pep, ["motif_calling", wildcards.model, "fimo_param_string"], wildcards.sample,\
        "")  
    shell:
       "fimo --bfile {input.bg} "
       "--text "
       "{params.fimo_param_string} "
       "{input.motif_file} {input.fasta} " 
       "> {output} 2> {log.stderr}"
