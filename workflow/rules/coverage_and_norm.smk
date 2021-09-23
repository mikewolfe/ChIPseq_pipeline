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

rule get_raw_coverage:
    input:
        expand("results/coverage_and_norm/deeptools_coverage/{sample}_raw.bw", sample = samples(pep))


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
        "--samFlagInclude 66 --extendReads --exactScaling --minMappingQuality 10")

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

# Old way to get the log2ratio directly using deeptools
#rule deeptools_log2ratio:
#    input:
#        ext= "results/coverage_and_norm/deeptools_coverage/{sample}_{norm}.bw",
#        inp= lambda wildcards: get_log2ratio_matching_input(wildcards.sample,wildcards.norm, pep)
#    output:
#        "results/coverage_and_norm/deeptools_log2ratio/{sample}_{norm}_log2ratio.bw"
#    wildcard_constraints:
#        norm="RPKM|CPM|BPM|RPGC|median"
#    params: 
#        resolution = RES,
#        pseudocount = determine_pseudocount(config)
#    log:
#        stdout="results/coverage_and_norm/logs/deeptools_log2ratio/{sample}_{norm}_log2ratio.log",
#        stderr="results/coverage_and_norm/logs/deeptools_log2ratio/{sample}_{norm}_log2ratio.err"
#    threads:
#        5
#
#    conda:
#        "../envs/coverage_and_norm.yaml"
#    shell:
#        "bigwigCompare -b1 {input.ext} -b2 {input.inp} --outFileName {output} "
#        "--operation 'log2' "
#        "--binSize {params.resolution} "
#        "--pseudocount {params.pseudocount} "
#        "--skipZeroOverZero --skipNonCoveredRegions "
#        "--numberOfProcessors {threads} > {log.stdout} 2> {log.stderr}"


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
        "results/coverage_and_norm/deeptools_coverage/{sample}_raw.bw"
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
        infile = "results/coverage_and_norm/deeptools_coverage/{sample}_raw.bw",
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
        infile = "results/coverage_and_norm/bwtools_compare/{sample}_{norm}_{norm_btwn}_fixedsub.bw",
        fixed_regions = lambda wildcards: lookup_in_config_persample(config,\
        pep, ["coverage_and_norm", "bwtools_fixed_scale", "fixed_regions"],\
        wildcards.sample)
    output:
        "results/coverage_and_norm/bwtools_compare/{sample}_{norm}_{norm_btwn}_fixedsub_fixedscale.bw"
    params:
        resolution = RES,
        dropNaNsandInfs = determine_dropNaNsandInfs(config),
        summary_func = lambda wildcards: lookup_in_config_persample(config,\
        pep, ["coverage_and_norm", "bwtools_fixed_scale", "summary_func"],
        wildcards.sample, "mean")
    wildcard_constraints:
        norm="RPKM|CPM|BPM|RPGC|count|SES|median|spike",
        norm_btwn="log2ratio|recipratio|subinp|ratio"
    log:
        stdout="results/coverage_and_norm/logs/bwtools_manipulate/{sample}_{norm}_{norm_btwn}_fixedsub_fixedscale.log",
        stderr="results/coverage_and_norm/logs/bwtools_manipulate/{sample}_{norm}_{norm_btwn}_fixedsub_fixedscale.err"
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
        infile="results/coverage_and_norm/bwtools_compare/{sample}_{norm}_{norm_btwn}_querysub.bw",
        query_regions = lambda wildcards: lookup_in_config_persample(config,\
        pep, ["coverage_and_norm", "bwtools_query_scale", "query_regions"],\
        wildcards.sample,\
        "results/alignment/process_genbank/{genome}/{genome}.bed".format(genome = lookup_sample_metadata(wildcards.sample, "genome", pep)))
    output:
        "results/coverage_and_norm/bwtools_compare/{sample}_{norm}_{norm_btwn}_querysub_queryscale.bw"
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
        norm_btwn="log2ratio|recipratio|subinp|ratio"
    log:
        stdout="results/coverage_and_norm/logs/bwtools_manipulate/{sample}_{norm}_{norm_btwn}_querysub_queryscale.log",
        stderr="results/coverage_and_norm/logs/bwtools_manipulate/{sample}_{norm}_{norm_btwn}_querysub_queryscale.err"
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
