## HELPER FUNCTIONS inherited from parent SnakeFile:
# samples(pep)
# lookup_sample_metadata(sample, key, pep)

ruleorder: bowtie2_map > bowtie2_map_se

rule clean_alignment:
    shell:
        "rm -fr results/alignment/"

rule run_alignment:
    input:
        expand("results/alignment/bowtie2/{sample}_sorted.bam.bai", sample = samples(pep))

rule run_annotations:
    input:
        expand("results/alignment/combine_bed/{genome}/{genome}.bed", genome = config["reference"].keys()),
        expand("results/alignment/combine_bed/{genome}/{genome}_annotations.tsv", genome = config["reference"].keys())

def get_genome_fastas(config, genome):
    return config["reference"][genome]["fastas"]

def get_genome_annotations(config, genome, ext = "bed"):
    out = []
    for key in lookup_in_config(config, ["reference", genome, "genbanks"]).keys():
       out.append("results/alignment/process_genbank/%s/%s.%s"%(genome, key, ext))
    return out

def get_bt2_index(sample, pep):
    reference = lookup_sample_metadata(sample, "genome", pep)
    return "results/alignment/bowtie2_index/%s/%s"%(reference,reference)

def get_bt2_index_file(sample, pep):
    reference = lookup_sample_metadata(sample, "genome", pep)
    return "results/alignment/bowtie2_index/%s/%s.1.bt2"%(reference,reference)
    

rule pull_genbank:
    message: "Download genbank for {wildcards.accession}"
    output:
        "resources/genbanks/{accession}.gbk"
    log:
        stdout="results/alignment/logs/pull_genbank/{accession}.log",
        stderr="results/alignment/logs/pull_genbank/{accession}.err"
    threads: 1
    conda:
        "../envs/alignment.yaml"
    shell:
        "ncbi-acc-download {wildcards.accession} --out {output} > {log.stdout} "
        "2> {log.stderr}"

def determine_input_genbank(genome, contig, config):
    outfile = config["reference"][genome]["genbanks"][contig]
    return outfile

rule process_genbank:
    message: "Processing genbank for {wildcards.genome} {wildcards.contig}"
    input:
        lambda wildcards: determine_input_genbank(wildcards.genome, wildcards.contig, config)
    output:
        "results/alignment/process_genbank/{genome}/{contig}.{outfmt}"
    log:
        stdout="results/alignment/logs/process_genbank/{genome}/{contig}_{outfmt}.log",
        stderr="results/alignment/logs/process_genbank/{genome}/{contig}_{outfmt}.err"
    params:
        features =  lookup_in_config(config, ["alignment", "process_genbank", "features", "value"], "CDS tRNA rRNA ncRNA"),
        qual_name = lookup_in_config(config, ["alignment", "process_genbank", "qual_name", "value"], "gene")
    threads: 1
    conda:
        "../envs/alignment.yaml"
    shell:
         "python3 workflow/scripts/parse_genbank.py {input} "
         "--outfmt {wildcards.outfmt} "
         "--features {params.features} "
         "--qual_name {params.qual_name} "
         "--chrm '{wildcards.contig}'  "
         " > {output} 2> {log.stderr}"

def masked_regions_file_for_combine_fastas(config, genome):
    outfile = determine_masked_regions_file(config, genome)
    if outfile is not None:
        out = "--masked_regions %s"%(outfile)
    else:
        out = ""
    return out

rule combine_fastas:
    message: "Generating fasta for genome {wildcards.genome}"
    input:
        lambda wildcards: get_genome_fastas(config, wildcards.genome)
    output:
        "results/alignment/combine_fasta/{genome}/{genome}.fa",
        "results/alignment/combine_fasta/{genome}/{genome}_contig_sizes.tsv",
        "results/alignment/combine_fasta/{genome}/{genome}_mappable_size.txt"
    log:
        stdout="results/alignment/logs/combine_fasta/{genome}/{genome}.log",
        stderr="results/alignment/logs/combine_fasta/{genome}/{genome}.err"
    threads: 1
    params:
        masked_regions = lambda wildcards: masked_regions_file_for_combine_fastas(config, wildcards.genome)
    conda:
        "../envs/alignment.yaml"
    shell:
         "python3 workflow/scripts/combine_fasta.py "
         "results/alignment/combine_fasta/{wildcards.genome}/{wildcards.genome} "
         "{params.masked_regions} "
         "{input} > {log.stdout} 2> {log.stderr}"

rule combine_beds:
    message: "Generating bed for genome {wildcards.genome}"
    input:
        lambda wildcards: get_genome_annotations(config, wildcards.genome, ext = "bed")
    output:
        "results/alignment/combine_bed/{genome}/{genome}.bed",
    log:
        stdout="results/alignment/logs/combine_bed/{genome}/{genome}.log",
        stderr="results/alignment/logs/combine_bed/{genome}/{genome}.err"
    threads: 1
    conda:
        "../envs/alignment.yaml"
    shell:
         "python3 workflow/scripts/combine_bed.py "
         "results/alignment/combine_bed/{wildcards.genome}/{wildcards.genome} "
         "{input} > {log.stdout} 2> {log.stderr}"

rule get_annotation_table:
    input:
        lambda wildcards: get_genome_annotations(config, wildcards.genome, ext = "tsv")
    output:
        "results/alignment/combine_bed/{genome}/{genome}_annotations.tsv",
    log:
        stdout="results/alignment/logs/get_annotation_table/{genome}/{genome}.log",
        stderr="results/alignment/logs/get_annotation_table/{genome}/{genome}.err"
    threads: 1
    run:
        shell("cat %s > {output}"%(input[0]))
        if len(input) > 1:
            for inf in input[1:]:
                shell("cat %s | tail -n +2 >> {output}"%(inf)) 
    
rule bowtie2_index:
    input:
        "results/alignment/combine_fasta/{genome}/{genome}.fa"
    output:
        "results/alignment/bowtie2_index/{genome}/{genome}.1.bt2"
    threads:
        5
    log:
        stdout="results/alignment/logs/bowtie2_index/{genome}.log",
        stderr="results/alignment/logs/bowtie2_index/{genome}.err" 
    conda:
        "../envs/alignment.yaml"
    shell:
        "bowtie2-build --threads {threads} "
        "{input} "
        "results/alignment/bowtie2_index/{wildcards.genome}/{wildcards.genome} "
        "> {log.stdout} 2> {log.stderr}"

rule bowtie2_map:
    input:
        in1="results/preprocessing/trimmomatic/{sample}_trim_paired_R1.fastq.gz",
        in2="results/preprocessing/trimmomatic/{sample}_trim_paired_R2.fastq.gz",
        bt2_index_file= lambda wildcards: get_bt2_index_file(wildcards.sample,pep)
    output:
        temp("results/alignment/bowtie2/{sample}.bam")
    log:
        stderr="results/alignment/logs/bowtie2/{sample}_bt2.log" 
    params:
        bt2_index= lambda wildcards: get_bt2_index(wildcards.sample,pep),
        bowtie2_param_string= lambda wildcards: lookup_in_config_persample(config,\
        pep, ["alignment", "bowtie2_map", "bowtie2_param_string"], wildcards.sample,\
        "--end-to-end --very-sensitive --phred33"),
        samtools_view_param_string = lambda wildcards: lookup_in_config_persample(config,\
        pep, ["alignment", "bowtie2_map", "samtools_view_param_string"],\
        wildcards.sample, "-b")
    threads: 
        5
    conda:
        "../envs/alignment.yaml"
    shell:
        "bowtie2 -x {params.bt2_index} -p {threads} "
        "-1 {input.in1} -2 {input.in2}  "
        "{params.bowtie2_param_string} 2> {log.stderr} "
        "| samtools view {params.samtools_view_param_string} > {output}"

rule bowtie2_map_se:
    input:
        in1="results/preprocessing/trimmomatic/{sample}_trim_R0.fastq.gz",
        bt2_index_file= lambda wildcards: get_bt2_index_file(wildcards.sample,pep)
    output:
        temp("results/alignment/bowtie2/{sample}.bam")
    log:
        stderr="results/alignment/logs/bowtie2/{sample}_bt2.log" 
    params:
        bt2_index= lambda wildcards: get_bt2_index(wildcards.sample,pep),
        bowtie2_param_string= lambda wildcards: lookup_in_config_persample(config,\
        pep, ["alignment", "bowtie2_map_se", "bowtie2_param_string"], wildcards.sample,\
        "--end-to-end --very-sensitive --phred33"),
        samtools_view_param_string = lambda wildcards: lookup_in_config_persample(config,\
        pep, ["alignment", "bowtie2_map_se", "samtools_view_param_string"],\
        wildcards.sample, "-b")
    threads: 
        5
    conda:
        "../envs/alignment.yaml"
    shell:
        "bowtie2 -x {params.bt2_index} -p {threads} "
        "-U {input.in1} "
        "{params.bowtie2_param_string} 2> {log.stderr} "
        "| samtools view {params.samtools_view_param_string} > {output}"

rule bam_sort:
    input:
        "results/alignment/bowtie2/{sample}.bam"
    output:
        "results/alignment/bowtie2/{sample}_sorted.bam"
    log:
        stderr="results/alignment/logs/bowtie2/{sample}_bt2_sort.log"
    conda:
        "../envs/alignment.yaml"
    shell:
        "samtools sort {input} > {output} 2> {log.stderr}"

rule bam_index:
    input:
        "results/alignment/bowtie2/{sample}_sorted.bam"
    output:
        "results/alignment/bowtie2/{sample}_sorted.bam.bai"
    log:
        stdout="results/alignment/logs/bowtie2/{sample}_bt2_index.log",
        stderr="results/alignment/logs/bowtie2/{sample}_bt2_index.log"
    conda:
        "../envs/alignment.yaml"
    shell:
        "samtools index {input} {output} > {log.stdout} 2> {log.stderr}"
