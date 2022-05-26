## HELPER FUNCTIONS inherited from parent SnakeFile:
# samples(pep)
# lookup_sample_metadata(sample, key, pep)

# GLOBAL CONSTANTS inherited from parent Snakefile:
# RES - resolution of coverage
# WITHIN - normalization to perform within samples
 
rule clean_peak_calling:
    shell:
        "rm -fr results/peak_calling"


def determine_peak_calling_files(config, pep):
    outfiles = []
    for model in lookup_in_config(config, ["peak_calling"], []):
        peak_caller = lookup_in_config(config, ["peak_calling", model, "peak_caller"], 
        err = "Need peak caller specified for peak_caller model %s in config file. I.e. \npeak_calling:\n\t%s:\n\t\tpeak_caller: 'macs2'"%(model,model))
        these_samples = filter_samples(pep, \
        lookup_in_config(config, ["peak_calling", model, "filter"], "not input_sample.isnull()"))
        if peak_caller == "macs2_broad":
            for sample in these_samples:
                outfiles.append("results/peak_calling/%s/macs2/%s_peaks.broadPeak"%(model, sample))
        if peak_caller == "macs2":
            for sample in these_samples:
                outfiles.append("results/peak_calling/%s/macs2/%s_peaks.narrowPeak"%(model, sample))
        if peak_caller == "cmarrt":
            for sample in these_samples:
                outfiles.append("results/peak_calling/%s/cmarrt/%s.narrowPeak"%(model,sample))
        if peak_caller == "macs3":
            for sample in these_samples:
                outfiles.append("results/peak_calling/%s/macs3/%s_peaks.xls"%(model, sample))
    return outfiles

rule run_peak_calling:
    input:
        determine_peak_calling_files(config, pep)

def determine_peak_coverage_files(config, pep):
    outfiles = []
    for model in lookup_in_config(config, ["peak_calling"], []):
        if lookup_in_config(config, ["peak_calling", model, "peak_coverage"], ""):
            these_samples = filter_samples(pep, \
            lookup_in_config(config, ["peak_calling", model, "filter"], "not input_sample.isnull()"))
            outfiles.extend(["results/peak_calling/%s/peak_coverage/%s_peak_coverage.tsv.gz"%(model, sample) for sample in these_samples])
    return outfiles

rule get_peak_coverage:
    input:
        determine_peak_coverage_files(config, pep)


def determine_cmarrt_input(sample, model, config, pep):
   file_sig = lookup_in_config(config, ["peak_calling", model, "filesignature"], 
   "results/coverage_and_norm/bwtools_compare/%s_BPM_subinp.bw")
   return file_sig%(sample)

rule cmarrt_call_peaks:
    input:
        lambda wildcards: determine_cmarrt_input(wildcards.sample, wildcards.model, config, pep)
    output:
        "results/peak_calling/{model}/cmarrt/{sample}.narrowPeak"
    log:
        stdout="results/peak_calling/logs/{model}/cmarrt/{sample}_cmarrt.log",
        stderr="results/peak_calling/logs/{model}/cmarrt/{sample}_cmarrt.err"
    params:
        cmarrt_param_string = lambda wildcards: lookup_in_config_persample(config,\
        pep, ["peak_calling", wildcards.model, "cmarrt_param_string"], "--resolution 25 --plots"),
        window_size = lambda wildcards: lookup_in_config_persample(config,\
        pep, ["peak_calling", wildcards.model, "cmarrt_window_size"], 25),
        outpre = lambda wildcards: "results/peak_calling/%s/cmarrt/%s"%(wildcards.model, wildcards.sample),
        resolution = RES,
    conda:
        "../envs/peak_calling.yaml"
    shell:
        "python3 "
        "workflow/CMARRT_python/run_cmarrt_bigwig.py "
        " {input} {params.window_size} "
        "-o {params.outpre} "
        "--resolution {params.resolution} "
        "{params.cmarrt_param_string} > {log.stdout} 2> {log.stderr}"


def get_macs2_matching_input(sample, pep):
   input_sample = lookup_sample_metadata(sample, "input_sample", pep)
   return "results/alignment/bowtie2/%s_sorted.bam"%(input_sample)

rule macs2_call_peaks:
    input:
        ext = "results/alignment/bowtie2/{sample}_sorted.bam",
        inp = lambda wildcards: get_macs2_matching_input(wildcards.sample, pep),
        genome_size= lambda wildcards: determine_effective_genome_size_file(wildcards.sample, config, pep)
    output:
        "results/peak_calling/{model}/macs2/{sample}_summits.bed",
        "results/peak_calling/{model}/macs2/{sample}_peaks.narrowPeak"
    log:
        stdout="results/peak_calling/logs/{model}/macs2/{sample}_macs2.log",
        stderr="results/peak_calling/logs/{model}/macs2/{sample}_macs2.err"
    params:
        macs2_param_string = lambda wildcards: lookup_in_config_persample(config,\
        pep, ["peak_calling", wildcards.model, "macs2_param_string"], wildcards.sample,\
        "-f BAMPE ") 
    conda:
        "../envs/peak_calling.yaml"
    shell:
        "macs2 callpeak -t {input.ext} -c {input.inp} -n {wildcards.sample} "
        "--outdir results/peak_calling/{wildcards.model}/macs2/ "
        "--call-summits "
        "{params.macs2_param_string} "
        "-g $(cat {input.genome_size}) > {log.stdout} 2> {log.stderr}"


rule macs2_broad_call_peaks:
    input:
        ext = "results/alignment/bowtie2/{sample}_sorted.bam",
        inp = lambda wildcards: get_macs2_matching_input(wildcards.sample, pep),
        genome_size= lambda wildcards: determine_effective_genome_size_file(wildcards.sample, config, pep)
    output:
        "results/peak_calling/{model}/macs2/{sample}_peaks.broadPeak"
    log:
        stdout="results/peak_calling/logs/{model}/macs2/{sample}_macs2.log",
        stderr="results/peak_calling/logs/{model}/macs2/{sample}_macs2.err"
    params:
        macs2_param_string = lambda wildcards: lookup_in_config_persample(config,\
        pep, ["peak_calling", wildcards.model, "macs2_param_string"], wildcards.sample,\
        "-f BAMPE ") 
    conda:
        "../envs/peak_calling.yaml"
    shell:
        "macs2 callpeak -t {input.ext} -c {input.inp} -n {wildcards.sample} "
        "--outdir results/peak_calling/{wildcards.model}/macs2/ "
        "--broad "
        "{params.macs2_param_string} "
        "-g $(cat {input.genome_size}) > {log.stdout} 2> {log.stderr}"

rule macs3_call_peaks:
    input:
        ext = "results/alignment/bowtie2/{sample}_sorted.bam",
        inp = lambda wildcards: get_macs2_matching_input(wildcards.sample, pep),
        genome_size= lambda wildcards: determine_effective_genome_size_file(wildcards.sample, config, pep)
    output:
        "results/peak_calling/{model}/macs3/{sample}_peaks.xls",
        "results/peak_calling/{model}/macs3/{sample}_summits.bed",
        "results/peak_calling/{model}/macs3/{sample}.narrowPeak"
    log:
        stdout="results/peak_calling/logs/{model}/macs3/{sample}_macs3.log",
        stderr="results/peak_calling/logs/{model}/macs3/{sample}_macs3.err"
    params:
        macs3_param_string = lambda wildcards: lookup_in_config_persample(config,\
        pep, ["peak_calling", wildcards.model, "macs3_param_string"], wildcards.sample,\
        "-f BAMPE --broad") 
    conda:
        "../envs/peak_calling.yaml"
    shell:
        "macs3 callpeak -t {input.ext} -c {input.inp} -n {wildcards.sample} "
        "--outdir results/peak_calling/{wildcards.model}/macs3/ "
        "--call-summits "
        "{params.macs3_param_string} "
        "-g $(cat {input.genome_size}) > {log.stdout} 2> {log.stderr}" 

rule peak_coverage:
    input: 
        coverage = lambda wildcards: lookup_in_config_persample(config, pep,\
        ["peak_calling", wildcards.model,"peak_coverage","cov_filesig"],\
        wildcards.sample, "results/coverage_and_norm/deeptools_coverage/%s_raw.bw")%(wildcards.sample),
        regions = lambda wildcards: lookup_in_config_persample(config, pep,\
        ["peak_calling", wildcards.model,"peak_coverage", "peaks_filesig"], wildcards.sample)%(wildcards.sample) 
    output:
        outtext="results/peak_calling/{model}/peak_coverage/{sample}_peak_coverage.tsv.gz"
    log:
        stdout="results/peak_calling/logs/{model}/peak_coverage/{sample}.log",
        stderr="results/peak_calling/logs/{model}/peak_coverage/{sample}.err"
    params:
        upstream = lambda wildcards: lookup_in_config(config, ["peak_calling", wildcards.model, "peak_coverage", "upstream"], 0),
        downstream = lambda wildcards: lookup_in_config(config, ["peak_calling", wildcards.model, "peak_coverage", "downstream"], 0),
        res = RES
    threads:
        5
    conda:
        "../envs/coverage_and_norm.yaml"
    shell:
        "python3 workflow/scripts/bwtools.py query "
        "{output.outtext} "
        "{input.coverage} "
        "--res {params.res} "
        "--regions {input.regions} "
        "--upstream {params.upstream} "
        "--downstream {params.downstream} "
        "--samp_names {wildcards.sample} "
        "--summarize identity "
        "--gzip "
        "> {log.stdout} 2> {log.stderr} "
