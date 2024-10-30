

def get_sampling_genome(sample, pep):
    reference = lookup_sample_metadata(sample, "genome", pep)
    return "results/alignment/combine_fasta/%s/%s.fa"%(reference,reference)

rule sample_fastqs:
    message: "Sampling for test sample {wildcards.sample}"
    input:
        genome = lambda wildcards: get_sampling_genome(wildcards.sample, pep),
        bed = lambda wildcards: lookup_sample_metadata(wildcards.sample,
        "region_weights", pep)
    output:
        "results/test/input_fastqs/{sample}_R1.fastq.gz",
        "results/test/input_fastqs/{sample}_R2.fastq.gz"
    log:
        stdout="results/test/logs/sample_fastqs/{sample}.log",
        stderr="results/test/logs/sample_fastqs/{sample}.err"
    params:
        sample_param_string = lambda wildcards: lookup_in_config_persample(config,
        pep, ["test", "sample_fastqs", "sample_param_string"], wildcards.sample,
        default = "--rng_seed 42")
    conda:
        "../envs/test.yaml"
    shell:
        "python3 workflow/scripts/FastqSim.py Regions {input.genome} "
        "results/test/input_fastqs/{wildcards.sample} "
        "--scores {input.bed} "
        "{params.sample_param_string} "
        "> {log.stdout} 2> {log.stderr}"
