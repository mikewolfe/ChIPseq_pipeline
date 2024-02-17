## HELPER FUNCTIONS inherited from parent SnakeFile:
# samples(pep)
# lookup_sample_metadata(sample, key, pep)
# GLOBAL CONSTANTS inherited from parent Snakefile:
# RES - resolution of coverage
# WITHIN - normalization to perform within samples
# ENDING - log2ratio or log2ratioRZ


# OVERALL RULES

rule clean_coverage_and_norm:
    shell:
        "rm -fr results/coverage_and_norm/"

rule run_coverage_and_norm:
    input:       
        expand("results/coverage_and_norm/bwtools_compare/{sample}_{within}_{ending}.bw",\
        sample = determine_extracted_samples(pep),\
        within = WITHIN,\
        ending = ENDING)


def determine_group_norm_files(config, pep):
    outfiles = []
    models = lookup_in_config(config, ["coverage_and_norm", "group_norm"], "")
    for model in models: 
        these_samples = filter_samples(pep, \
            lookup_in_config(config, ["coverage_and_norm", "group_norm", model, "filter"], "not input_sample.isnull()"))
        for norm_type in lookup_in_config(config, ["coverage_and_norm", "group_norm", model, "methods"], []):
            outfiles.extend(["results/coverage_and_norm/group_norm/%s/%s_%s.bw"%(model, sample, norm_type) for sample in these_samples]) 
        for dwnsample in lookup_in_config(config, ["coverage_and_norm", "group_norm", model, "dwnsample"], []):
            outfiles.extend(["results/coverage_and_norm/group_norm/%s/dwnsample/%s_%s.bw"%(model, sample, dwnsample) for sample in these_samples]) 

    return outfiles

rule run_group_norm:
    input:       
        determine_group_norm_files(config, pep)

rule get_raw_coverage:
    input:
        expand("results/coverage_and_norm/deeptools_coverage/{sample}_raw.bw", sample = samples(pep))

def raw_or_smoothed(sample, pep, config):
    filtered = filter_samples(pep, \
    lookup_in_config(config, ["coverage_and_norm", "smooth_samples", "filter"], \
    "not input_sample.isnull() and input_sample.isnull()"))
    if sample in filtered:
        out = "results/coverage_and_norm/bwtools_smooth/%s_smooth.bw"%sample
    else:
        out = "results/coverage_and_norm/deeptools_coverage/%s_raw.bw"%sample
    return out
        


def masked_regions_file_for_deeptools(config, sample, pep):
    genome = lookup_sample_metadata(sample, "genome", pep) 
    outfile = determine_masked_regions_file(config, genome)
    if outfile is not None:
        out = "--blackListFileName %s"%(outfile)
    else:
        out = ""
    return out

rule deeptools_coverage_raw:
    input:
        inbam="results/alignment/bowtie2/{sample}_sorted.bam",
        ind="results/alignment/bowtie2/{sample}_sorted.bam.bai"
    output:
        "results/coverage_and_norm/deeptools_coverage/{sample}_raw.bw"
    log:
        stdout="results/coverage_and_norm/logs/deeptools_coverage/{sample}_raw.log",
        stderr="results/coverage_and_norm/logs/deeptools_coverage/{sample}_raw.err"
    params:
        resolution = RES,
        masked_regions = lambda wildcards: masked_regions_file_for_deeptools(config, wildcards.sample, pep),
        bamCoverage_param_string= lambda wildcards: lookup_in_config_persample(config,\
        pep, ["coverage_and_norm", "deeptools_coverage", "bamCoverage_param_string"], wildcards.sample,\
        "--samFlagInclude 67 --extendReads")

    threads:
        5
    conda:
        "../envs/coverage_and_norm.yaml"
    shell:
        "bamCoverage --bam {input.inbam} --outFileName {output} "
        "--outFileFormat 'bigwig' "
        "--numberOfProcessors {threads} --binSize {params.resolution} "
        "{params.masked_regions} "
        "{params.bamCoverage_param_string} "
        "> {log.stdout} 2> {log.stderr}"

rule bwtools_smooth:
    input:
        "results/coverage_and_norm/deeptools_coverage/{sample}_raw.bw"
    output:
        "results/coverage_and_norm/bwtools_smooth/{sample}_smooth.bw"
    log:
        stdout="results/coverage_and_norm/logs/bwtools_smooth/{sample}_smoothed.log",
        stderr="results/coverage_and_norm/logs/bwtools_smooth/{sample}_smoothed.err"
    params:
        resolution = RES,
        dropNaNsandInfs = determine_dropNaNsandInfs(config),
        smooth_params = lambda wildcards: lookup_in_config_persample(config,\
        pep, ["coverage_and_norm", "bwtools_smooth", "param_string"], wildcards.sample,\
        "--operation 'flat_smooth' --wsize 50 --edge 'wrap'")
    threads:
        1
    conda:
        "../envs/coverage_and_norm.yaml"
    shell:
        "python3 "
        "workflow/scripts/bwtools.py manipulate "
        "{input} {output} "
        "--res {params.resolution} "
        "{params.smooth_params} "
        "{params.dropNaNsandInfs} "
        "> {log.stdout} 2> {log.stderr}"


rule deeptools_coverage:
    input:
        inbam="results/alignment/bowtie2/{sample}_sorted.bam",
        ind="results/alignment/bowtie2/{sample}_sorted.bam.bai",
        genome_size= lambda wildcards: determine_effective_genome_size_file(wildcards.sample, config, pep)
    output:
        "results/coverage_and_norm/deeptools_coverage/{sample}_{norm}.bw"
    log:
        stdout="results/coverage_and_norm/logs/deeptools_coverage/{sample}_{norm}.log",
        stderr="results/coverage_and_norm/logs/deeptools_coverage/{sample}_{norm}.err"
    params:
        resolution = RES,
        masked_regions = lambda wildcards: masked_regions_file_for_deeptools(config, wildcards.sample, pep),
        bamCoverage_param_string= lambda wildcards: lookup_in_config_persample(config,\
        pep, ["coverage_and_norm", "deeptools_coverage", "bamCoverage_param_string"], wildcards.sample,\
        "--samFlagInclude 66 --extendReads --exactScaling --minMappingQuality 10")
    wildcard_constraints:
        norm="RPKM|CPM|BPM|RPGC"
    threads:
        5
    conda:
        "../envs/coverage_and_norm.yaml"
    shell:
        "bamCoverage --bam {input.inbam} --outFileName {output} --outFileFormat 'bigwig' "
        "--numberOfProcessors {threads} --binSize {params.resolution} "
        "--normalizeUsing {wildcards.norm} "
        "--effectiveGenomeSize $(cat {input.genome_size}) "
        "{params.masked_regions} "
        "{params.bamCoverage_param_string} "
        "> {log.stdout} 2> {log.stderr}"

# helper for finding correct input sample
def get_log2ratio_matching_input(sample, norm, pep):
   input_sample = lookup_sample_metadata(sample, "input_sample", pep)
   return "results/coverage_and_norm/deeptools_coverage/%s_%s.bw"%(input_sample, norm)


rule bwtools_log2ratio:
    input:
        ext= "results/coverage_and_norm/deeptools_coverage/{sample}_{norm}.bw",
        inp= lambda wildcards: get_log2ratio_matching_input(wildcards.sample,wildcards.norm, pep)
    output:
        "results/coverage_and_norm/bwtools_compare/{sample}_{norm}_log2ratio.bw"
    wildcard_constraints:
        norm="RPKM|CPM|BPM|RPGC|median"
    params: 
        resolution = RES,
        dropNaNsandInfs = determine_dropNaNsandInfs(config)
    log:
        stdout="results/coverage_and_norm/logs/bwtools_compare/{sample}_{norm}_log2ratio.log",
        stderr="results/coverage_and_norm/logs/bwtools_compare/{sample}_{norm}_log2ratio.err"
    threads:
        1
    conda:
        "../envs/coverage_and_norm.yaml"
    shell:
        "python3 "
        "workflow/scripts/bwtools.py compare {input.ext} {input.inp} {output} "
        "--operation 'log2ratio' "
        "--res {params.resolution} "
        "{params.dropNaNsandInfs} "
        "> {log.stdout} 2> {log.stderr}"


rule bwtools_recipratio:
    input:
        ext= "results/coverage_and_norm/deeptools_coverage/{sample}_{norm}.bw",
        inp= lambda wildcards: get_log2ratio_matching_input(wildcards.sample,wildcards.norm, pep)
    output:
        "results/coverage_and_norm/bwtools_compare/{sample}_{norm}_recipratio.bw"
    wildcard_constraints:
        norm="RPKM|CPM|BPM|RPGC|median"
    params: 
        resolution = RES,
        dropNaNsandInfs = determine_dropNaNsandInfs(config)
    log:
        stdout="results/coverage_and_norm/logs/bwtools_compare/{sample}_{norm}_recipratio.log",
        stderr="results/coverage_and_norm/logs/bwtools_compare/{sample}_{norm}_recipratio.err"
    threads:
        1
    conda:
        "../envs/coverage_and_norm.yaml"
    shell:
        "python3 "
        "workflow/scripts/bwtools.py compare {input.ext} {input.inp} {output} "
        "--operation 'recipratio' "
        "--res {params.resolution} "
        "{params.dropNaNsandInfs} "
        "> {log.stdout} 2> {log.stderr}"

def get_input_bam(sample, pep):
    input_sample = lookup_sample_metadata(sample, "input_sample", pep)
    return "results/alignment/bowtie2/%s_sorted.bam"%(input_sample),

def get_input_bai(sample, pep):
    input_sample = lookup_sample_metadata(sample, "input_sample", pep)
    return "results/alignment/bowtie2/%s_sorted.bam.bai"%(input_sample),


rule deeptools_SES_log2ratio:
    input:
        ext= "results/alignment/bowtie2/{sample}_sorted.bam",
        ind_ext= "results/alignment/bowtie2/{sample}_sorted.bam.bai",
        inp = lambda wildcards: get_input_bam(wildcards.sample, pep),
        inp_ext = lambda wildcards: get_input_bai(wildcards.sample, pep)

    output:
        "results/coverage_and_norm/bwtools_compare/{sample}_SES_log2ratio.bw"
    log:
        stdout="results/coverage_and_norm/logs/deeptools_compare/{sample}_SES_log2ratio.log",
        stderr="results/coverage_and_norm/logs/deeptools_compare/{sample}_SES_log2ratio.err"
    threads:
        5
    params: 
        resolution = RES,
        bamCoverage_param_string= lambda wildcards: lookup_in_config_persample(config,\
        pep, ["coverage_and_norm", "deeptools_coverage", "bamCoverage_param_string"], wildcards.sample,\
        "--samFlagInclude 66 --extendReads --exactScaling --minMappingQuality 10")
    conda:
        "../envs/coverage_and_norm.yaml"
    shell:
        "bamCompare --bamfile1 {input.ext} --bamfile2 {input.inp} --outFileName {output} "
        "--outFileFormat 'bigwig' --scaleFactorsMethod 'SES' --operation 'log2' "
        "--binSize {params.resolution} "
        "--numberOfProcessors {threads} "
        "{params.bamCoverage_param_string} "
        "> {log.stdout} 2> {log.stderr}"

rule bwtools_median:
    input:
        lambda wildcards: raw_or_smoothed(wildcards.sample, pep, config)
    output:
        "results/coverage_and_norm/deeptools_coverage/{sample}_median.bw"
    params:
        resolution = RES,
        dropNaNsandInfs = determine_dropNaNsandInfs(config),
        pseudocount = lambda wildcards: lookup_in_config_persample(config,\
        pep, ["coverage_and_norm", "bwtools_median", "pseudocount"], wildcards.sample, 0)
    log:
        stdout="results/coverage_and_norm/logs/bwtools_manipulate/{sample}_median.log",
        stderr="results/coverage_and_norm/logs/bwtools_manipulate/{sample}_median.err"
    conda:
        "../envs/coverage_and_norm.yaml"
    shell:
       "python3 "
       "workflow/scripts/bwtools.py manipulate "
       "{input} {output} "
       "--res {params.resolution} --operation Median_norm "
       "--pseudocount {params.pseudocount} "
       "{params.dropNaNsandInfs} "
       "> {log.stdout} 2> {log.stderr}"

rule bwtools_spike_scale:
    input:
        infile = lambda wildcards: raw_or_smoothed(wildcards.sample, pep, config),
        fixed_regions = lambda wildcards: lookup_in_config_persample(config,\
        pep, ["coverage_and_norm", "bwtools_spike_scale", "fixed_regions"],
        wildcards.sample)
    output:
        "results/coverage_and_norm/deeptools_coverage/{sample}_spike.bw"
    params:
        resolution = RES,
        dropNaNsandInfs = determine_dropNaNsandInfs(config),
        summary_func = lambda wildcards: lookup_in_config_persample(config,\
        pep, ["coverage_and_norm", "bwtools_spike_scale", "summary_func"],
        wildcards.sample, "mean")
    log:
        stdout="results/coverage_and_norm/logs/bwtools_manipulate/{sample}_spike.log",
        stderr="results/coverage_and_norm/logs/bwtools_manipulate/{sample}_spike.err"
    conda:
        "../envs/coverage_and_norm.yaml"
    shell:
       "python3 "
       "workflow/scripts/bwtools.py manipulate "
       "{input.infile} {output} "
       "--res {params.resolution} --operation spike_scale "
       "--fixed_regions {input.fixed_regions} "
       "{params.dropNaNsandInfs} "
       "--summary_func {params.summary_func} "
       "> {log.stdout} 2> {log.stderr}"

rule bwtools_fixed_subtract:
    input:
        infile = "results/coverage_and_norm/bwtools_compare/{sample}_{norm}_{norm_btwn}.bw",
        fixed_regions = lambda wildcards: lookup_in_config_persample(config,\
        pep, ["coverage_and_norm", "bwtools_fixed_subtract", "fixed_regions"],\
        wildcards.sample)
    output:
        "results/coverage_and_norm/bwtools_compare/{sample}_{norm}_{norm_btwn}_fixedsub.bw"
    params:
        resolution = RES,
        dropNaNsandInfs = determine_dropNaNsandInfs(config),
        summary_func = lambda wildcards: lookup_in_config_persample(config,\
        pep, ["coverage_and_norm", "bwtools_fixed_subtract", "summary_func"],
        wildcards.sample, "mean")
    wildcard_constraints:
        norm="RPKM|CPM|BPM|RPGC|count|SES|median|spike",
        norm_btwn="log2ratio|recipratio|subinp|ratio"
    log:
        stdout="results/coverage_and_norm/logs/bwtools_manipulate/{sample}_{norm}_{norm_btwn}_fixedsub.log",
        stderr="results/coverage_and_norm/logs/bwtools_manipulate/{sample}_{norm}_{norm_btwn}_fixedsub.err"
    conda:
        "../envs/coverage_and_norm.yaml"
    shell:
       "python3 "
       "workflow/scripts/bwtools.py manipulate "
       "{input.infile} {output} "
       "--res {params.resolution} --operation fixed_subtract "
       "--fixed_regions {input.fixed_regions} "
       "{params.dropNaNsandInfs} "
       "--summary_func {params.summary_func} "
       "> {log.stdout} 2> {log.stderr}"

rule bwtools_fixed_scale:
    input:
        infile = "results/coverage_and_norm/bwtools_compare/{sample}_{norm}_{norm_btwn}_{bckgrd_sub}.bw",
        fixed_regions = lambda wildcards: lookup_in_config_persample(config,\
        pep, ["coverage_and_norm", "bwtools_fixed_scale", "fixed_regions"],\
        wildcards.sample)
    output:
        "results/coverage_and_norm/bwtools_compare/{sample}_{norm}_{norm_btwn}_{bckgrd_sub}_fixedscale.bw"
    params:
        resolution = RES,
        dropNaNsandInfs = determine_dropNaNsandInfs(config),
        summary_func = lambda wildcards: lookup_in_config_persample(config,\
        pep, ["coverage_and_norm", "bwtools_fixed_scale", "summary_func"],
        wildcards.sample, "mean")
    wildcard_constraints:
        norm="RPKM|CPM|BPM|RPGC|count|SES|median|spike",
        norm_btwn="log2ratio|recipratio|subinp|ratio",
        bckgrd_sub = "fixedsub|querysub"
    log:
        stdout="results/coverage_and_norm/logs/bwtools_manipulate/{sample}_{norm}_{norm_btwn}_{bckgrd_sub}_fixedscale.log",
        stderr="results/coverage_and_norm/logs/bwtools_manipulate/{sample}_{norm}_{norm_btwn}_{bckgrd_sub}_fixedscale.err"
    conda:
        "../envs/coverage_and_norm.yaml"
    shell:
       "python3 "
       "workflow/scripts/bwtools.py manipulate "
       "{input.infile} {output} "
       "--res {params.resolution} --operation fixed_scale "
       "--fixed_regions {input.fixed_regions} "
       "{params.dropNaNsandInfs} "
       "--summary_func {params.summary_func} "
       "> {log.stdout} 2> {log.stderr}"

rule bwtools_scale_max:
    input:
        "results/coverage_and_norm/bwtools_compare/{sample}_{norm}_{norm_btwn}_minbg.bw"
    output:
        "results/coverage_and_norm/bwtools_compare/{sample}_{norm}_{norm_btwn}_minbg_scalemax.bw"
    params:
        resolution = RES,
        dropNaNsandInfs = determine_dropNaNsandInfs(config)
    log:
        stdout="results/coverage_and_norm/logs/bwtools_manipulate/{sample}_{norm}_{norm_btwn}_minbg_scalemax.log",
        stderr="results/coverage_and_norm/logs/bwtools_manipulate/{sample}_{norm}_{norm_btwn}_minbg_scalemax.err"
    wildcard_constraints:
        norm="RPKM|CPM|BPM|RPGC|count|SES|median|spike",
        norm_btwn="log2ratio|recipratio|subinp|ratio"
    conda:
        "../envs/coverage_and_norm.yaml"
    shell:
       "python3 "
       "workflow/scripts/bwtools.py manipulate "
       "{input} {output} "
       "--res {params.resolution} --operation scale_max "
       "{params.dropNaNsandInfs} "
       "> {log.stdout} 2> {log.stderr}"

rule bwtools_query_subtract:
    input:
        infile="results/coverage_and_norm/bwtools_compare/{sample}_{norm}_{norm_btwn}.bw",
        query_regions = lambda wildcards: lookup_in_config_persample(config,\
        pep, ["coverage_and_norm", "bwtools_query_subtract", "query_regions"],\
        wildcards.sample,\
        "results/alignment/process_genbank/{genome}/{genome}.bed".format(genome = lookup_sample_metadata(wildcards.sample, "genome", pep)))
    output:
        "results/coverage_and_norm/bwtools_compare/{sample}_{norm}_{norm_btwn}_querysub.bw"
    params:
        resolution = RES,
        dropNaNsandInfs = determine_dropNaNsandInfs(config),
        number_of_regions = lambda wildcards: lookup_in_config_persample(config,\
        pep, ["coverage_and_norm", "bwtools_query_subtract", "number_of_regions"],\
        wildcards.sample, 20),
        summary_func = lambda wildcards: lookup_in_config_persample(config,\
        pep, ["coverage_and_norm", "bwtools_query_subtract", "summary_func"],
        wildcards.sample, "mean")
    log:
        stdout="results/coverage_and_norm/logs/bwtools_manipulate/{sample}_{norm}_{norm_btwn}_querysub.log",
        stderr="results/coverage_and_norm/logs/bwtools_manipulate/{sample}_{norm}_{norm_btwn}_querysub.err"
    wildcard_constraints:
        norm="RPKM|CPM|BPM|RPGC|count|SES|median|spike",
        norm_btwn="log2ratio|recipratio|subinp|ratio"
    conda:
        "../envs/coverage_and_norm.yaml"
    shell:
       "python3 "
       "workflow/scripts/bwtools.py manipulate "
       "{input.infile} {output} "
       "--res {params.resolution} --operation query_subtract "
       "{params.dropNaNsandInfs} "
       "--number_of_regions {params.number_of_regions} "
       "--query_regions {input.query_regions} "
       "--summary_func {params.summary_func} "
       "> {log.stdout} 2> {log.stderr}"


rule bwtools_query_scale:
    input:
        infile="results/coverage_and_norm/bwtools_compare/{sample}_{norm}_{norm_btwn}_{bckgrd_sub}.bw",
        query_regions = lambda wildcards: lookup_in_config_persample(config,\
        pep, ["coverage_and_norm", "bwtools_query_scale", "query_regions"],\
        wildcards.sample,\
        "results/alignment/process_genbank/{genome}/{genome}.bed".format(genome = lookup_sample_metadata(wildcards.sample, "genome", pep)))
    output:
        "results/coverage_and_norm/bwtools_compare/{sample}_{norm}_{norm_btwn}_{bckgrd_sub}_queryscale.bw"
    params:
        resolution = RES,
        dropNaNsandInfs = determine_dropNaNsandInfs(config),
        number_of_regions = lambda wildcards: lookup_in_config_persample(config,\
        pep, ["coverage_and_norm", "bwtools_query_scale", "number_of_regions"],\
        wildcards.sample, 20),
        summary_func = lambda wildcards: lookup_in_config_persample(config,\
        pep, ["coverage_and_norm", "bwtools_query_scale", "summary_func"],
        wildcards.sample, "mean")
    wildcard_constraints:
        norm="RPKM|CPM|BPM|RPGC|count|SES|median|spike",
        norm_btwn="log2ratio|recipratio|subinp|ratio",
        bckgrd_sub="querysub|fixedsub"
    log:
        stdout="results/coverage_and_norm/logs/bwtools_manipulate/{sample}_{norm}_{norm_btwn}_{bckgrd_sub}_queryscale.log",
        stderr="results/coverage_and_norm/logs/bwtools_manipulate/{sample}_{norm}_{norm_btwn}_{bckgrd_sub}_queryscale.err"
    conda:
        "../envs/coverage_and_norm.yaml"
    shell:
       "python3 "
       "workflow/scripts/bwtools.py manipulate "
       "{input.infile} {output} "
       "--res {params.resolution} --operation query_scale "
       "{params.dropNaNsandInfs} "
       "--number_of_regions {params.number_of_regions} "
       "--query_regions {input.query_regions} "
       "--summary_func {params.summary_func} "
       "> {log.stdout} 2> {log.stderr}"

rule bwtools_RobustZ:
    input:
        "results/coverage_and_norm/bwtools_compare/{sample}_{norm}_{norm_btwn}.bw"
    output:
        "results/coverage_and_norm/bwtools_compare/{sample}_{norm}_{norm_btwn}RZ.bw"
    log:
        stdout="results/coverage_and_norm/logs/bwtools_manipulate/{sample}_{norm}_{norm_btwn}RZ.log",
        stderr="results/coverage_and_norm/logs/bwtools_manipulate/{sample}_{norm}_{norm_btwn}RZ.err"
    wildcard_constraints:
        norm="RPKM|CPM|BPM|RPGC|count|SES|median|spike",
        norm_btwn="log2ratio|recipratio|subinp"
    params:
        resolution = RES,
        dropNaNsandInfs = determine_dropNaNsandInfs(config)
    conda:
        "../envs/coverage_and_norm.yaml"
    shell:
       "python3 "
       "workflow/scripts/bwtools.py manipulate "
       "{input} {output} "
       "--res {params.resolution} "
       "{params.dropNaNsandInfs} "
       "--operation RobustZ > {log.stdout} 2> {log.stderr}"

rule bwtools_subinp:
    input:
        ext = "results/coverage_and_norm/deeptools_coverage/{sample}_{norm}.bw",
        inp = lambda wildcards: get_log2ratio_matching_input(wildcards.sample, wildcards.norm, pep)
    output:
        "results/coverage_and_norm/bwtools_compare/{sample}_{norm}_subinp.bw"
    wildcard_constraints:
        norm = "RPKM|CPM|BPM|RPGC|median|spike"
    params:
        resolution = RES,
        dropNaNsandInfs = determine_dropNaNsandInfs(config)
    log:
        stdout="results/coverage_and_norm/logs/bwtools_compare/{sample}_{norm}_subinp.log",
        stderr="results/coverage_and_norm/logs/bwtools_compare/{sample}_{norm}_subinp.err"
    threads:
        1
    conda:
        "../envs/coverage_and_norm.yaml"
    shell:
        "python3 "
        "workflow/scripts/bwtools.py compare {input.ext} {input.inp} {output} "
        "--operation 'subtract' "
        "--res {params.resolution} "
        "{params.dropNaNsandInfs} "
        "> {log.stdout} 2> {log.stderr}"


rule bwtools_ratio:
    input:
        ext = "results/coverage_and_norm/deeptools_coverage/{sample}_{norm}.bw",
        inp = lambda wildcards: get_log2ratio_matching_input(wildcards.sample, wildcards.norm, pep)
    output:
        "results/coverage_and_norm/bwtools_compare/{sample}_{norm}_ratio.bw"
    wildcard_constraints:
        norm = "RPKM|CPM|BPM|RPGC|median|spike"
    params:
        resolution = RES,
        dropNaNsandInfs = determine_dropNaNsandInfs(config)
    log:
        stdout="results/coverage_and_norm/logs/bwtools_compare/{sample}_{norm}_ratio.log",
        stderr="results/coverage_and_norm/logs/bwtools_compare/{sample}_{norm}_ratio.err"
    threads:
        1
    conda:
        "../envs/coverage_and_norm.yaml"
    shell:
        "python3 "
        "workflow/scripts/bwtools.py compare {input.ext} {input.inp} {output} "
        "--operation 'divide' "
        "--res {params.resolution} "
        "{params.dropNaNsandInfs} "
        "> {log.stdout} 2> {log.stderr}"

def pull_bws_for_bwtools_multicompare(modelname, config, pep):
    groupA_samples = filter_samples(pep, \
    lookup_in_config(config, ['coverage_and_norm', 'bwtools_multicompare', modelname, 'filter_groupA']))

    groupB_filter = lookup_in_config(config, ['coverage_and_norm', 'bwtools_multicompare', modelname, 'filter_groupB'], "")
    file_sig = lookup_in_config(config, ['coverage_and_norm', 'bwtools_multicompare', modelname, 'filesignature'],\
    "results/coverage_and_norm/deeptools_coverage/%s_raw.bw")
    groupA = [file_sig%(sample) for sample in groupA_samples]
    if groupB_filter is not "":
        for sample in filter_samples(pep, groupB_filter):
            groupA.append(file_sig%(sample))
    return groupA

def grouped_bws_for_bwtools_multicompare(modelname, config, pep, group_name):
    group_sample_filter = lookup_in_config(config, ['coverage_and_norm', 'bwtools_multicompare', modelname, 'filter_%s'%group_name], "")
    file_sig = lookup_in_config(config, ['coverage_and_norm', 'bwtools_multicompare', modelname, 'filesignature'],\
    "results/coverage_and_norm/deeptools_coverage/%s_raw.bw")
    if group_sample_filter is not "":
        group = " ".join([file_sig%(sample) for sample in filter_samples(pep, group_sample_filter)])
        group = "--%s "%(group_name) + group
    else:
        group = " "
    return group



rule bwtools_multicompare:
    input:
        lambda wildcards: pull_bws_for_bwtools_multicompare(wildcards.model, config, pep)
    output:
        "results/coverage_and_norm/bwtools_multicompare/{model}.bw"
    params:
        resolution = RES,
        groupA = lambda wildcards: grouped_bws_for_bwtools_multicompare(wildcards.model, config, pep, "groupA"),
        groupB = lambda wildcards: grouped_bws_for_bwtools_multicompare(wildcards.model, config, pep, "groupB"),
        within_op = lambda wildcards: lookup_in_config(config, ["coverage_and_norm", 'bwtools_multicompare', wildcards.model, "within_operation"], "mean"),
        btwn_op = lambda wildcards: lookup_in_config(config, ["coverage_and_norm", 'bwtools_multicompare', wildcards.model, "between_operation"], "subtract"),
        dropNaNsandInfs = determine_dropNaNsandInfs(config)
    log:
        stdout="results/coverage_and_norm/logs/bwtools_multicompare/{model}.log",
        stderr="results/coverage_and_norm/logs/bwtools_multicompare/{model}.err"
    threads:
        1
    conda:
        "../envs/coverage_and_norm.yaml"
    shell:
        "python3 "
        "workflow/scripts/bwtools.py multicompare {output} {params.groupA} {params.groupB} "
        "--within_operation {params.within_op} "
        "--between_operation {params.btwn_op} "
        "--res {params.resolution} "
        "{params.dropNaNsandInfs} "
        "> {log.stdout} 2> {log.stderr}"

def determine_multicompare_models(config):
    outfiles = []
    if lookup_in_config(config, ["coverage_and_norm", "bwtools_multicompare"], ""):
        outfiles.extend(["results/coverage_and_norm/bwtools_multicompare/%s.bw"%model\
        for model in config["coverage_and_norm"]["bwtools_multicompare"]])
    return outfiles

rule run_bwtools_multicompare:
    input:
        determine_multicompare_models(config)


def pull_bws_for_group_norm_models(modelname, config, pep, ext_or_inp = "ext"):
    these_samples = filter_samples(pep, \
    lookup_in_config(config, ["coverage_and_norm", "group_norm", modelname, "filter"], "not input_sample.isnull()"))
    file_sig = lookup_in_config(config, ["coverage_and_norm", "group_norm", modelname, "filesignature"],\
    "results/coverage_and_norm/deeptools_coverage/%s_raw.bw")

    if ext_or_inp == "ext":
        files = [file_sig%(sample) for sample in these_samples]
    else: 
        files = [file_sig%(lookup_sample_metadata(sample, "input_sample", pep)) for sample in these_samples]

    return files

def pull_labels_for_group_norm_models(modelname, config, pep):
    these_samples = filter_samples(pep, \
    lookup_in_config(config, ["coverage_and_norm", "group_norm", modelname, "filter"], "not input_sample.isnull()"))
    return " ".join(these_samples)

def pull_expected_regions_group_norm_models(modelname, config, pep):
    spikeregions = lookup_in_config(config, ["coverage_and_norm", "group_norm", modelname, "regions"], "")
    if spikeregions is not "":
        out = "--expected_regions %s"%(spikeregions)
    else:
        out = ""
    return out

rule group_norm_table:
    input:
        inextbws= lambda wildcards: pull_bws_for_group_norm_models(wildcards.model,config, pep, ext_or_inp = "ext"),
        ininpbws= lambda wildcards: pull_bws_for_group_norm_models(wildcards.model,config, pep, ext_or_inp = "inp"),
        inbed= lambda wildcards: lookup_in_config(config, ["coverage_and_norm", "group_norm", wildcards.model, "regions"], "config/config.yaml"),
        fragtable = "results/quality_control/frags_per_contig/all_samples.tsv",
        md= lambda wildcards: lookup_in_config(config, ["coverage_and_norm", "group_norm", wildcards.model, "metadata"], None),
    output:
        outtext="results/coverage_and_norm/group_norm/{model}/{model}_scale_factors.tsv"
    log:
        stdout="results/coverage_and_norm/logs/group_norm/{model}/scale_factors.log",
        stderr="results/coverage_and_norm/logs/group_norm/{model}/scale_factors.err"
    params:
        labels = lambda wildcards: pull_labels_for_group_norm_models(wildcards.model, config, pep),
        pseudocount = lambda wildcards: lookup_in_config(config, ["coverage_and_norm", "group_norm", wildcards.model, "pseudocount"], 0.1),
        spikecontigs = lambda wildcards: lookup_in_config(config, ["coverage_and_norm", "group_norm", wildcards.model, "spikecontigs"], None),
        spikeregions = lambda wildcards: pull_expected_regions_group_norm_models(wildcards.model, config, pep),
        res = RES,
        upstream = lambda wildcards: lookup_in_config(config, ["coverage_and_norm", "group_norm", wildcards.model, "upstream"], 0),
        downstream = lambda wildcards: lookup_in_config(config, ["coverage_and_norm", "group_norm", wildcards.model, "downstream"], 0)
    threads:
        5
    conda:
        "../envs/coverage_and_norm.yaml"
    shell:
        "python3 workflow/scripts/bwtools.py normfactor "
        "{output.outtext} "
        "{input.fragtable} "
        "{input.md} "
        "--ext_bws {input.inextbws} "
        "--inp_bws {input.ininpbws} "
        "--pseudocount {params.pseudocount} "
        "{params.spikeregions} "
        "--spikecontigs {params.spikecontigs} "
        "--res {params.res} "
        "--samples {params.labels} "
        "--upstream {params.upstream} "
        "--downstream {params.downstream} "
        "> {log.stdout} 2> {log.stderr} "


def get_sample_for_scale_byfactor(modelname, sample, config, pep): 
    file_sig = lookup_in_config(config, ["coverage_and_norm", "group_norm", modelname, "filesignature"],\
    "results/coverage_and_norm/deeptools_coverage/%s_raw.bw")
    return file_sig%(sample)

rule bwtools_scale_byfactor:
    input:
        infile = lambda wildcards: get_sample_for_scale_byfactor(wildcards.model, wildcards.sample, config, pep),
        sf_tab="results/coverage_and_norm/group_norm/{model}/{model}_scale_factors.tsv"
    output:
        "results/coverage_and_norm/group_norm/{model}/{sample}_{norm}.bw"
    params:
        resolution = RES,
        dropNaNsandInfs = determine_dropNaNsandInfs(config),
        pseudocount = lambda wildcards: lookup_in_config(config, ["coverage_and_norm", "group_norm", wildcards.model, "pseudocount"], 0.1)
    wildcard_constraints:
        norm="total_frag_sfs|spike_frag_sfs|nonspike_frag_sfs|deseq2_sfs|deseq2_spike_sfs|regress_rpm_sfs|regress_total_rpm_sfs|regress_spike_rpm_sfs|tmm_rpm_sfs|tmm_spike_rpm_sfs"
    log:
        stdout="results/coverage_and_norm/logs/group_norm/{model}/{sample}_{norm}.log",
        stderr="results/coverage_and_norm/logs/group_norm/{model}/{sample}_{norm}.err"
    conda:
        "../envs/coverage_and_norm.yaml"
    shell:
       "python3 "
       "workflow/scripts/bwtools.py manipulate "
       "{input.infile} {output} "
       "--res {params.resolution} --operation scale_byfactor "
       "--pseudocount {params.pseudocount} "
       "--scalefactor_table {input.sf_tab} "
       "--scalefactor_id {wildcards.sample} {wildcards.norm} "
       "{params.dropNaNsandInfs} "
       "> {log.stdout} 2> {log.stderr}"

rule bwtools_dwnsample:
    input:
        inplus = lambda wildcards: get_sample_for_scale_byfactor(wildcards.model, wildcards.sample, config, pep),
        sf_tab="results/coverage_and_norm/group_norm/{model}/{model}_scale_factors.tsv"
    output:
        outplus="results/coverage_and_norm/group_norm/{model}/dwnsample/{sample}_plus_{norm}.bw",
        outminus="results/coverage_and_norm/group_norm/{model}/dwnsample/{sample}_minus_{norm}.bw"
    params:
        resolution = RES,
        dropNaNsandInfs = determine_dropNaNsandInfs(config)
    wildcard_constraints:
        norm="total_frag|spike_frag|nonspike_frag"
    log:
        stdout="results/coverage_and_norm/logs/group_norm/{model}/dwnsample/{sample}_{norm}.log",
        stderr="results/coverage_and_norm/logs/group_norm/{model}/dwnsample/{sample}_{norm}.err"
    conda:
        "../envs/coverage_and_norm.yaml"
    shell:
       "python3 "
       "workflow/scripts/bwtools.py manipulate "
       "{input.inplus} {output.outplus} "
       "--minus_strand {input.inminus} --minus_strand_out {output.outminus} "
       "--res {params.resolution} --operation downsample "
       "--scalefactor_table {input.sf_tab} "
       "--scalefactor_id {wildcards.sample} {wildcards.norm} "
       "{params.dropNaNsandInfs} "
       "> {log.stdout} 2> {log.stderr}"
